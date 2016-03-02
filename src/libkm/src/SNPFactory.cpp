/* 
 * File:   SNPFactory.cpp
 * Author: veraalva
 * 
 * Created on February 25, 2016, 9:22 AM
 */

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <string>
#include <memory>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <algorithm>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "bmath.h"
#include "svm.h"
#include "Global.h"
#include "KmersFactory.h"
#include "FastaFactory.h"
#include "SNPFactory.h"
#include "SVMPredict.h"

using namespace std;
using namespace fasta;
using namespace kmers;
using namespace snp;
using namespace svm;

SNP::SNP() {
    this->scoreDrop = 0.0;
    this->neighborKmer = 0.0;
    this->overlapKmer = 0.0;
}

SNP::~SNP() {
    if (this->seq) free(seq);
}

void SNP::CalculateKmerDescriptors(kmers::KmersFactory& kmersFactory) {
    unsigned long int i, startPos, endPos;
    double overlapMutated = 0;
    char t;
    unsigned long int len = static_cast<unsigned long int> (strlen(this->seq));

    if (static_cast<int> (this->pos - Global::instance()->GetOrder() - 1) <= 0) {
        startPos = 0;
    } else {
        startPos = this->pos - Global::instance()->GetOrder() + 1;
    }
    if (this->pos + Global::instance()->GetOrder() >= len) {
        endPos = len - Global::instance()->GetOrder();
    } else {
        endPos = this->pos;
    }

    for (i = 0; i <= len - Global::instance()->GetOrder(); i++) {
        t = this->seq[i + Global::instance()->GetOrder()];
        this->seq[i + Global::instance()->GetOrder()] = '\0';
        if (i >= startPos && i <= endPos) {
            this->overlapKmer += kmersFactory.GetKmerSig(this->seq + i);
            this->seq[this->pos] = this->alt;
            overlapMutated += kmersFactory.GetKmerSig(this->seq + i);
            this->seq[this->pos] = this->ref;
        } else {
            this->neighborKmer += kmersFactory.GetKmerSig(this->seq + i);
            if (isinf(this->neighborKmer)) {
                cout << this->GetId() << endl;
                cout << (this->seq + i) << endl;
                cout << "Sig: " << kmersFactory.GetKmerSig(this->seq + i) << endl;
                cout << endl;
                exit(0);
            }
        }
        this->seq[i + Global::instance()->GetOrder()] = t;
    }

    this->scoreDrop = this->overlapKmer - overlapMutated;
}

SNPFactory::SNPFactory() {
}

SNPFactory::SNPFactory(const SNPFactory& orig) {
}

SNPFactory::~SNPFactory() {
    for (auto it = this->snps.begin(); it != this->snps.end(); ++it) {
        delete (*it);
    }
}

int SNPFactory::ProcessSNPFromFiles(char* snpFileName, unsigned long int neighbors, FastaFactory &chrFactory, kmers::KmersFactory& kmersFactory, svm::SVMPredict& svmPredict) {
    FILE *snpFile = (FILE *) checkPointerError(fopen(snpFileName, "r"), "Can't open bed file", __FILE__, __LINE__, -1);

    int i, count = 0;
    size_t bufferSize, read, backupLineSize;
    char *buffer, *newLine, *str, *backupLine, *completeLine, *seq;
    char **fields = NULL;
    size_t fieldsSize = 0;
    Fasta *f = NULL;
    SNP *snp;
    int snpPos, startPos, endPos;

    int featNumber = 3;
    double target_label = 0.0;
    struct svm_node *x = (struct svm_node *) allocate(sizeof (struct svm_node) * (featNumber + 1), __FILE__, __LINE__);
    x[featNumber].index = -1;

    double *mean = (double *) allocate(sizeof (double *) * featNumber, __FILE__, __LINE__);
    double *sd = (double *) allocate(sizeof (double *) * featNumber, __FILE__, __LINE__);
    for (i = 0; i < featNumber; i++) {
        mean[i] = 0.0;
        sd[i] = 0.0;
    }

    backupLineSize = bufferSize = 100000000;
    buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
    backupLine = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);

    *backupLine = '\0';
    while (!feof(snpFile)) {
        read = fread(buffer, sizeof (char), bufferSize, snpFile);
        buffer[read] = '\0';
        str = buffer;
        while (*str != '\0') {
            newLine = strchr(str, '\n');
            if (newLine) *newLine = '\0';
            if (*backupLine != '\0') {
                if (strlen(backupLine) + strlen(str) + 1 > backupLineSize) {
                    backupLineSize += backupLineSize;
                    backupLine = (char *) reallocate(backupLine, sizeof (char) * (backupLineSize + 1), __FILE__, __LINE__);
                }
                strcat(backupLine, str);
                completeLine = backupLine;
            } else {
                completeLine = str;
            }

            if (*completeLine != '#') {
                fieldsSize = splitString(&fields, completeLine, "\t");
                if (fieldsSize != 5) {
                    printLog(stderr, "Input SNP file with a wrong format ", __FILE__, __LINE__, -1);
                }
                if (!f || f->GetId().compare(fields[0]) != 0) {
                    f = chrFactory.GetFastaFromID(fields[0]);
                }
                if (f) {
                    if (Global::instance()->isDebug3()) {
                        cout << completeLine << endl;
                    }
                    snp = new SNP();
                    snp->SetId(fields[2]);
                    snp->SetChr(fields[0]);

                    snp->SetRef(fields[3][0]);
                    snp->SetAlt(fields[4][0]);
                    snpPos = atoi(fields[1]) - 1;
                    if (snpPos - 1 >= 0 && snpPos - 1 < f->GetLength()) {
                        snp->SetChrPos(snpPos);
                    } else {
                        printLog(stderr, "SNP position out of range", __FILE__, __LINE__, -1);
                    }
                    if (static_cast<int> (snpPos - neighbors - 1) < 0) {
                        startPos = 0;
                        snp->SetPos(snpPos);
                    } else {
                        startPos = snpPos - neighbors;
                        snp->SetPos(neighbors);
                    }
                    if (snpPos + neighbors >= f->GetLength()) {
                        endPos = f->GetLength() - 1;
                    } else {
                        endPos = snpPos + neighbors;
                    }

                    seq = f->GetSubStr(startPos, endPos - startPos + 1);
                    toUpper(&seq);
                    snp->SetSeq(&seq);

                    if (Global::instance()->isDebug3()) {
                        cout << snp->GetRef() << " " << snp->GetSeq()[snp->GetPos()]
                                << " " << snp->GetPos()
                                << " " << snp->GetLength()
                                << " " << snp->GetSeq()
                                << endl;
                    }

                    if (snp->GetSeq()[snp->GetPos()] != snp->GetRef()) {
                        fprintf(stderr, "\nERROR1:\n\n%c != %c\n\n", snp->GetSeq()[snp->GetPos()], fields[3][0]);
                        fprintf(stderr, "%c\t%c\t%lu\n\n", snp->GetRef(), snp->GetSeq()[snp->GetPos()], snp->GetPos());
                        fprintf(stderr, "%s\n\n", snp->GetSeq());
                        if (Global::instance()->isDebug3()) {
                            printLog(stderr, "SNP is not in the chromosome position provided", __FILE__, __LINE__, -1);
                        }
                        delete snp;
                    } else {
                        count++;
                        snp->CalculateKmerDescriptors(kmersFactory);

                        /*
                         * Calculating overall sum for mean and sd
                         */
                        mean[0] += snp->GetScoreDrop();
                        mean[1] += snp->GetOverlapKmer();
                        mean[2] += snp->GetNeighborKmer();


                        this->snps.push_back(move(snp));
                    }
                }
                freeArrayofPointers((void **) fields, fieldsSize);
            }
            *backupLine = '\0';
            if (newLine) str = newLine + 1;
            else str += strlen(str);
        }

        if (strlen(str) > 0) {
            if (strlen(backupLine) + strlen(str) + 1 > backupLineSize) {
                backupLineSize += backupLineSize;
                backupLine = (char *) reallocate(backupLine, sizeof (char) * (backupLineSize + 1), __FILE__, __LINE__);
            }
            strcat(backupLine, str);
        }
    }

    /*
     * Calculating final mean
     */
    for (i = 0; i < featNumber; i++) {
        mean[i] = mean[i] / static_cast<double> (this->snps.size());
    }

    /*
     * Summing standard deviation 
     */
    for (auto it = this->snps.begin(); it != this->snps.end(); ++it) {
        SNP *s = *it;
        sd[0] += (s->GetScoreDrop() - mean[0])*(s->GetScoreDrop() - mean[0]);
        sd[1] += (s->GetOverlapKmer() - mean[1])*(s->GetOverlapKmer() - mean[1]);
        sd[2] += (s->GetNeighborKmer() - mean[2])*(s->GetNeighborKmer() - mean[2]);
    }

    /*
     * Calculating final standard deviation
     */
    for (i = 0; i < featNumber; i++) {
        sd[i] = sqrt(sd[i] / static_cast<double> (this->snps.size()));
    }

    /*
     * Calculating ZScore terms
     */
    for (auto it = this->snps.begin(); it != this->snps.end(); ++it) {
        SNP *s = *it;
        x[0].index = 1;
        x[0].value = (s->GetScoreDrop() - mean[0]) / sd[0];
        x[1].index = 2;
        x[1].value = (s->GetOverlapKmer() - mean[1]) / sd[1];
        x[2].index = 3;
        x[2].value = (s->GetNeighborKmer() - mean[2]) / sd[2];

        svmPredict.SVMPredictCalulation(x, target_label);
        s->SetProbPos(svmPredict.GetProb_estimates()[0]);
    }

    free(x);
    free(mean);
    free(sd);
    free(buffer);
    free(backupLine);
    fclose(snpFile);
    return count;
}


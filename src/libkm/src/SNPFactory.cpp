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
#include <cmath>
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
#include "FimoFactory.h"
#include "TFBSFactory.h"
#include "SNPFactory.h"
#include "SVMPredict.h"

using namespace std;
using namespace fasta;
using namespace kmers;
using namespace snp;
using namespace svm;
using namespace fimo;
using namespace tfbs;

SNP::SNP() {

}

SNP::~SNP() {
    if (seq) free(seq);
}

void SNP::CalculateKmerDescriptors(kmers::KmersFactory& kmersFactory, unsigned long int featNumber) {
    unsigned long int i, startPos, endPos;
    double overlapMutated = 0.0;
    char t;
    unsigned long int len = static_cast<unsigned long int> (strlen(seq));
    FILE *overlapingKmersFile;
    FILE *neighborsKmersFile;
    FILE *overlapingKmersMutatedFile;

    for (i = 0; i < featNumber; i++) {
        descriptors.push_back(0.0000);
    }

    if (static_cast<int> (pos - Global::instance()->GetOrder() - 1) <= 0) {
        startPos = 0;
    } else {
        startPos = pos - Global::instance()->GetOrder() + 1;
    }
    if (pos + Global::instance()->GetOrder() >= len) {
        endPos = len - Global::instance()->GetOrder();
    } else {
        endPos = pos;
    }

    if (Global::instance()->isDebug3()) {
        overlapingKmersFile = fopen("CAPE_debug_3_overlaping_kmers.txt", "a+");
        neighborsKmersFile = fopen("CAPE_debug_3_neighbors_kmers.txt", "a+");
        overlapingKmersMutatedFile = fopen("CAPE_debug_3_overlaping_kmers_mutated.txt", "a+");
    }

    for (i = 0; i <= len - Global::instance()->GetOrder(); i++) {
        t = seq[i + Global::instance()->GetOrder()];
        seq[i + Global::instance()->GetOrder()] = 0;

        if (i >= startPos && i <= endPos) {

            descriptors[1] += kmersFactory.GetKmerSig(seq + i);
            if (Global::instance()->isDebug3()) {
                fprintf(overlapingKmersFile, "%lu\t%s\t%.15f\t%.15f\n", i, seq + i, kmersFactory.GetKmerSig(seq + i), descriptors[1]);
            }

            seq[pos] = alt;
            overlapMutated += kmersFactory.GetKmerSig(seq + i);

            if (Global::instance()->isDebug3()) {
                fprintf(overlapingKmersMutatedFile, "%lu\t%s\t%.15f\t%.15f\n", i, seq + i, kmersFactory.GetKmerSig(seq + i), overlapMutated);
            }

            seq[pos] = ref;
        } else {
            if (Global::instance()->isDebug3()) {
                fprintf(neighborsKmersFile, "%lu\t%s\t%.5f\n", i, seq + i, kmersFactory.GetKmerSig(seq + i));
            }
            descriptors[2] += kmersFactory.GetKmerSig(seq + i);
            if (std::isinf(descriptors[2])) {
                cout << GetId() << endl;
                cout << (seq + i) << endl;
                cout << "Sig: " << kmersFactory.GetKmerSig(seq + i) << endl;
                cout << endl;
                exit(0);
            }
        }
        seq[i + Global::instance()->GetOrder()] = t;
    }
    if (Global::instance()->isDebug3()) {
        fprintf(overlapingKmersFile, "91-100\tsum\t%.15f\n", descriptors[1]);
        fprintf(overlapingKmersMutatedFile, "91-100\tsum\t%.15f\n", overlapMutated);
        fclose(overlapingKmersFile);
        fclose(neighborsKmersFile);
        fclose(overlapingKmersMutatedFile);
    }
    descriptors[0] = std::fabs(descriptors[1] - overlapMutated);
}

SNPFactory::SNPFactory() {
}

SNPFactory::SNPFactory(const SNPFactory& orig) {
}

SNPFactory::~SNPFactory() {
    for (auto it = snps.begin(); it != snps.end(); ++it) {
        delete (*it);
    }
}

void SNPFactory::ReadSNPFromFile(char* snpFileName, unsigned long int neighbors, FastaFactory &chrFactory) {
    FILE *snpFile = (FILE *) checkPointerError(fopen(snpFileName, "r"), "Can't open SNP file", __FILE__, __LINE__, -1);

    size_t bufferSize, read, backupLineSize;
    char *buffer, *newLine, *str, *backupLine, *completeLine, *seq;
    char **fields = NULL;
    size_t fieldsSize = 0;
    Fasta *f = NULL;
    SNP *snp;
    int snpPos, startPos, endPos;

    backupLineSize = bufferSize = 100000000;
    buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
    backupLine = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);

    *backupLine = 0;
    while (!feof(snpFile)) {
        read = fread(buffer, sizeof (char), bufferSize, snpFile);
        buffer[read] = 0;
        if (feof(snpFile)) {
            if (buffer[read - 1] != '\n') {
                buffer[read] = '\n';
                buffer[read + 1] = 0;
            }
        }
        str = buffer;
        while ((newLine = strchr(str, '\n')) != NULL) {
            *newLine = 0;
            if (*backupLine != 0) {
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
                    if (snpPos - 1 >= 0 && snpPos - 1 < static_cast<int> (f->GetLength())) {
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
                        fprintf(stderr, "ID: %s\n", snp->GetId().c_str());
                        fprintf(stderr, "%s\n\n", snp->GetSeq());
                        if (Global::instance()->isDebug3()) {
                            printLog(stderr, "SNP is not in the chromosome position provided", __FILE__, __LINE__, -1);
                        }
                        delete snp;
                    } else {
                        snps.push_back(move(snp));
                    }
                }
                freeArrayofPointers((void **) fields, fieldsSize);
            }
            *backupLine = 0;
            str = newLine + 1;
        }

        if (strlen(str) > 0) {
            if (strlen(backupLine) + strlen(str) + 1 > backupLineSize) {
                backupLineSize += backupLineSize;
                backupLine = (char *) reallocate(backupLine, sizeof (char) * (backupLineSize + 1), __FILE__, __LINE__);
            }
            strcat(backupLine, str);
        }
    }

    free(buffer);
    free(backupLine);
    fclose(snpFile);
}

void SNPFactory::WriteEnhansersFastaFile(char* fastaFile, bool binary) {
    FastaFactory fastaFactory;
    Fasta *f = NULL;
    SNP *snp = NULL;
    char *seq = NULL;
    string id;

    for (auto it = snps.begin(); it != snps.end(); ++it) {
        snp = *it;
        f = new Fasta();
        id = snp->GetId();
        f->SetId(id);
        seq = strdup(snp->GetSeq());
        f->SetLength(strlen(seq));
        f->SetSeq(&seq);
        fastaFactory.GetFastaMap().insert(pair<string, Fasta*>(id, move(f)));
    }

    fastaFactory.WriteSequencesToFile(fastaFile, binary);
}

int SNPFactory::ProcessSNPFromFile(char* snpFileName, unsigned long int neighbors, FastaFactory &chrFactory, KmersFactory& kmersFactory, SVMPredict& svmPredict, FimoFactory & fimoFactory, TFBSFactory & tFBSFactory) {
    FILE *snpFile = (FILE *) checkPointerError(fopen(snpFileName, "r"), "Can't open bed file", __FILE__, __LINE__, -1);

    double overlapValue, neighborSum;
    unsigned long int i, count = 0;
    size_t bufferSize, read, backupLineSize;
    char *buffer, *newLine, *str, *backupLine, *completeLine, *seq;
    char **fields = NULL;
    size_t fieldsSize = 0;
    Fasta *f = NULL;
    SNP *snp;
    int snpPos, startPos, endPos;
    FILE *featuresFile;
    FILE *zscoreFile;

    std::pair<std::string, double> cPair;

    unsigned long int featNumber = 3;
    double target_label = 0.0;

    if (Global::instance()->isDebug3()) {
        featuresFile = fopen("CAPE_debug_3_features.txt", "w");
        zscoreFile = fopen("CAPE_debug_3_zscores.txt", "w");
    }

    if (!fimoFactory.GetSnpIDMap().empty() || tFBSFactory.isReady()) {
        featNumber = 5;
    }

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

    cout.precision(4);
    *backupLine = 0;
    while (!feof(snpFile)) {
        read = fread(buffer, sizeof (char), bufferSize, snpFile);
        buffer[read] = 0;
        if (feof(snpFile)) {
            if (buffer[read - 1] != '\n') {
                buffer[read] = '\n';
                buffer[read + 1] = 0;
            }
        }
        str = buffer;
        while ((newLine = strchr(str, '\n')) != NULL) {
            *newLine = 0;
            if (*backupLine != 0) {
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
                    snp = new SNP();
                    snp->SetId(fields[2]);
                    snp->SetChr(fields[0]);

                    snp->SetRef(fields[3][0]);
                    snp->SetAlt(fields[4][0]);
                    snpPos = atoi(fields[1]) - 1;
                    if (snpPos - 1 >= 0 && snpPos - 1 < static_cast<int> (f->GetLength())) {
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
                        snp->CalculateKmerDescriptors(kmersFactory, featNumber);

                        /*
                         * Calculating overall sum for mean and sd
                         */
                        for (i = 0; i < 3; i++) {
                            mean[i] += snp->GetDescriptors()[i];
                        }

                        if (featNumber == 5) {
                            if (!fimoFactory.GetSnpIDMap().empty()) {
                                auto fimoMapIt = fimoFactory.GetSnpIDMap().find(snp->GetId());
                                if (fimoMapIt != fimoFactory.GetSnpIDMap().end()) {
                                    snp->GetDescriptors()[3] = fimoMapIt->second[0];
                                    snp->GetDescriptors()[4] = fimoMapIt->second[1];
                                }
                            } else {
                                overlapValue = neighborSum = 0;
                                tFBSFactory.ExtractTFBSFromFile(startPos, endPos, f);
                                for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
                                    TFBS *t = *it;
                                    cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), expressionCode);
                                    if (fabs(cPair.second) >= 1.0) {
                                        if (t->GetStart() <= snpPos && t->GetEnd() >= snpPos) {
                                            if (overlapValue < cPair.second) {
                                                overlapValue = cPair.second;
                                            }
                                        } else {
                                            neighborSum += cPair.second;
                                        }
                                    }
                                }
                                snp->GetDescriptors()[3] = overlapValue;
                                snp->GetDescriptors()[4] = neighborSum;
                            }
                            mean[3] += snp->GetDescriptors()[3];
                            mean[4] += snp->GetDescriptors()[4];
                        }
                        snps.push_back(move(snp));
                    }
                }
                freeArrayofPointers((void **) fields, fieldsSize);
            }
            *backupLine = 0;
            str = newLine + 1;
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
        mean[i] = mean[i] / static_cast<double> (snps.size());
    }

    /*
     * Summing standard deviation 
     */
    for (auto it = snps.begin(); it != snps.end(); ++it) {
        SNP *s = *it;
        for (i = 0; i < featNumber; i++) {
            sd[i] += (s->GetDescriptors()[i] - mean[i])*(s->GetDescriptors()[i] - mean[i]);
        }
    }

    /*
     * Calculating final standard deviation
     */
    for (i = 0; i < featNumber; i++) {
        sd[i] = sqrt(sd[i] / static_cast<double> (snps.size()));
    }

    /*
     * Calculating ZScore terms
     */
    for (auto it = snps.begin(); it != snps.end(); ++it) {
        SNP *s = *it;

        if (Global::instance()->isDebug3()) {
            fprintf(featuresFile, "%s\t", s->GetId().c_str());
            fprintf(zscoreFile, "%s\t", s->GetId().c_str());
        }
        for (i = 0; i < featNumber; i++) {
            x[i].index = i + 1;
            x[i].value = (s->GetDescriptors()[i] - mean[i]) / sd[i];
            if (Global::instance()->isDebug3()) {
                fprintf(featuresFile, "%.4f\t", s->GetDescriptors()[i]);
                fprintf(zscoreFile, "%.4f\t", x[i].value);
            }
        }
        if (Global::instance()->isDebug3()) {
            fprintf(featuresFile, "\n");
            fprintf(zscoreFile, "\n");
        }
        svmPredict.SVMPredictCalulation(x, target_label);
        s->SetProbPos(svmPredict.GetProb_estimates()[0]);
    }

    if (Global::instance()->isDebug3()) {
        fclose(featuresFile);
        fclose(zscoreFile);
    }
    free(x);
    free(mean);
    free(sd);
    free(buffer);
    free(backupLine);
    fclose(snpFile);

    cout.precision(2);

    return count;
}


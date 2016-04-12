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
#include "Exceptions.h"
#include "FileParserFactory.h"
#include "KmersFactory.h"
#include "FastaFactory.h"
#include "FimoFactory.h"
#include "TFBSFactory.h"
#include "SNPFactory.h"
#include "SVMPredict.h"

using namespace std;
using namespace parsers;
using namespace sequence;
using namespace kmers;
using namespace snp;
using namespace svm;
using namespace fimo;
using namespace tfbs;

SNP::SNP() {
    this->length = 0;
    this->seq = NULL;
    this->pos = 0;
    this->chrPos = 0;
    this->ref = 0;
    this->alt = 0;
    this->probPos = 0.0;
}

SNP::~SNP() {
    if (seq) free(seq);
}

void SNP::calculateKmerDescriptors(kmers::KmersFactory& kmersFactory, unsigned long int featNumber) {
    unsigned long int i, startPos, endPos;
    double overlapMutated = 0.0;
    unsigned long int len = static_cast<unsigned long int> (strlen(seq));

    for (i = 0; i < featNumber; i++) {
        descriptors.push_back(0.0000);
    }

    if (static_cast<int> (pos - Global::instance()->getOrder() - 1) <= 0) {
        startPos = 0;
    } else {
        startPos = pos - Global::instance()->getOrder() + 1;
    }
    if (pos + Global::instance()->getOrder() >= len) {
        endPos = len - Global::instance()->getOrder();
    } else {
        endPos = pos;
    }

    for (i = 0; i <= len - Global::instance()->getOrder(); i++) {
        char t = seq[i + Global::instance()->getOrder()];
        seq[i + Global::instance()->getOrder()] = 0;

        if (i >= startPos && i <= endPos) {

            descriptors[1] += kmersFactory.getKmerSig(seq + i);
            seq[pos] = alt;
            overlapMutated += kmersFactory.getKmerSig(seq + i);
            seq[pos] = ref;
        } else {
            descriptors[2] += kmersFactory.getKmerSig(seq + i);
        }
        seq[i + Global::instance()->getOrder()] = t;
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

void SNPFactory::parseSNPFile(std::string snpFileName, unsigned long int neighbors, FastaFactory &chrFactory) {
    FileParserFactory fParser(snpFileName);
    Seq *f = NULL;
    SNP *snp;
    int snpPos, startPos, endPos;

    while (fParser.iterate('#', "\t")) {
        if (fParser.getNWords() != 5) {
            printLog(stderr, "Input SNP file with a wrong format ", __FILE__, __LINE__, -1);
        }
        if (!f || f->getId().compare(fParser.getWords()[0]) != 0) {
            f = chrFactory.getSequenceFromID(fParser.getWords()[0]);
        }
        if (f) {
            snp = new SNP();
            snp->setId(fParser.getWords()[2]);
            snp->setChr(fParser.getWords()[0]);

            snp->setRef(fParser.getWords()[3][0]);
            snp->setAlt(fParser.getWords()[4][0]);
            snpPos = atoi(fParser.getWords()[1]) - 1;
            if (snpPos - 1 >= 0 && snpPos - 1 < static_cast<int> (f->getLength())) {
                snp->setChrPos(snpPos);
            } else {
                printLog(stderr, "SNP position out of range", __FILE__, __LINE__, -1);
            }
            if (static_cast<int> (snpPos - neighbors - 1) < 0) {
                startPos = 0;
                snp->setPos(snpPos);
            } else {
                startPos = snpPos - neighbors;
                snp->setPos(neighbors);
            }
            if (snpPos + neighbors >= f->getLength()) {
                endPos = f->getLength() - 1;
            } else {
                endPos = snpPos + neighbors;
            }

            char *seq = f->getSubStr(startPos, endPos - startPos + 1);
            toUpper(&seq);
            snp->setSeq(&seq);

            if (Global::instance()->isDebug3()) {
                cout << snp->getRef() << " " << snp->getSeq()[snp->getPos()]
                        << " " << snp->getPos()
                        << " " << snp->getLength()
                        << " " << snp->getSeq()
                        << endl;
            }

            if (snp->getSeq()[snp->getPos()] != snp->getRef()) {
                fprintf(stderr, "\nERROR1:\n\n%c != %c\n\n", snp->getSeq()[snp->getPos()], fParser.getWords()[3][0]);
                fprintf(stderr, "%c\t%c\t%lu\n\n", snp->getRef(), snp->getSeq()[snp->getPos()], snp->getPos());
                fprintf(stderr, "ID: %s\n", snp->getId().c_str());
                fprintf(stderr, "%s\n\n", snp->getSeq());
                if (Global::instance()->isDebug3()) {
                    printLog(stderr, "SNP is not in the chromosome position provided", __FILE__, __LINE__, -1);
                }
                delete snp;
            } else {
                snps.push_back(snp);
            }
        }
    }
}

void SNPFactory::writeEnhansersFastaFile(std::string fastaFile, bool binary) {
    FastaFactory fastaFactory;
    SNP *snp = NULL;
    char *seq = NULL;
    string id;

    for (auto it = snps.begin(); it != snps.end(); ++it) {
        snp = *it;
        Seq *f = new Seq();
        id = snp->getId();
        f->setId(id);
        seq = strdup(snp->getSeq());
        f->setLength(strlen(seq));
        f->setSeq(&seq);
        fastaFactory.getSequenceContainter().insert(pair<string, Seq*>(id, move(f)));
    }

    fastaFactory.writeSequencesToFile(fastaFile, binary);
}

int SNPFactory::processSNPFromFile(std::string snpFileName, unsigned long int neighbors, FastaFactory &chrFactory, KmersFactory& kmersFactory, SVMPredict& svmPredict, FimoFactory & fimoFactory, TFBSFactory & tFBSFactory) {
    FileParserFactory fParser(snpFileName);
    double overlapValue, neighborSum;
    unsigned long int i, count = 0;

    Seq *f = NULL;
    SNP *snp;
    int snpPos, startPos, endPos;
    FILE *featuresFile = stderr;
    FILE *zscoreFile = stderr;

    std::pair<std::string, double> cPair;

    unsigned long int featNumber = 3;
    double target_label = 0.0;

    if (Global::instance()->isDebug3()) {
        featuresFile = fopen("CAPE_debug_3_features.txt", "w");
        zscoreFile = fopen("CAPE_debug_3_zscores.txt", "w");
    }

    if (!fimoFactory.getSnpIDContainer().empty() || tFBSFactory.isReady()) {
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

    cout.precision(4);
    while (fParser.iterate('#', "\t")) {
        if (fParser.getNWords() != 5) {
            printLog(stderr, "Input SNP file with a wrong format ", __FILE__, __LINE__, -1);
        }
        if (!f || f->getId().compare(fParser.getWords()[0]) != 0) {
            f = chrFactory.getSequenceFromID(fParser.getWords()[0]);
        }
        if (f) {
            snp = new SNP();
            snp->setId(fParser.getWords()[2]);
            snp->setChr(fParser.getWords()[0]);

            snp->setRef(fParser.getWords()[3][0]);
            snp->setAlt(fParser.getWords()[4][0]);
            snpPos = atoi(fParser.getWords()[1]) - 1;
            if (snpPos - 1 >= 0 && snpPos - 1 < static_cast<int> (f->getLength())) {
                snp->setChrPos(snpPos);
            } else {
                printLog(stderr, "SNP position out of range", __FILE__, __LINE__, -1);
            }
            if (static_cast<int> (snpPos - neighbors - 1) < 0) {
                startPos = 0;
                snp->setPos(snpPos);
            } else {
                startPos = snpPos - neighbors;
                snp->setPos(neighbors);
            }
            if (snpPos + neighbors >= f->getLength()) {
                endPos = f->getLength() - 1;
            } else {
                endPos = snpPos + neighbors;
            }

            char *seq = f->getSubStr(startPos, endPos - startPos + 1);
            toUpper(&seq);
            snp->setSeq(&seq);

            if (snp->getSeq()[snp->getPos()] != snp->getRef()) {
                fprintf(stderr, "\nERROR1:\n\n%c != %c\n\n", snp->getSeq()[snp->getPos()], fParser.getWords()[3][0]);
                fprintf(stderr, "%c\t%c\t%lu\n\n", snp->getRef(), snp->getSeq()[snp->getPos()], snp->getPos());
                fprintf(stderr, "%s\n\n", snp->getSeq());
                if (Global::instance()->isDebug3()) {
                    printLog(stderr, "SNP is not in the chromosome position provided", __FILE__, __LINE__, -1);
                }
                delete snp;
            } else {
                count++;
                snp->calculateKmerDescriptors(kmersFactory, featNumber);

                /*
                 * Calculating overall sum for mean and sd
                 */
                for (i = 0; i < 3; i++) {
                    mean[i] += snp->getDescriptors()[i];
                }

                if (featNumber == 5) {
                    if (!fimoFactory.getSnpIDContainer().empty()) {
                        auto fimoMapIt = fimoFactory.getSnpIDContainer().find(snp->getId());
                        if (fimoMapIt != fimoFactory.getSnpIDContainer().end()) {
                            snp->getDescriptors()[3] = fimoMapIt->second[0];
                            snp->getDescriptors()[4] = fimoMapIt->second[1];
                        }
                    } else {
                        overlapValue = neighborSum = 0;
                        tFBSFactory.extractTFBSFromFile(startPos, endPos, f);
                        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
                            TFBS *t = *it;
                            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), expressionCode);
                            if (fabs(cPair.second) >= 1.0) {
                                if (t->getStart() <= snpPos && t->getEnd() >= snpPos) {
                                    if (overlapValue < cPair.second) {
                                        overlapValue = cPair.second;
                                    }
                                } else {
                                    neighborSum += cPair.second;
                                }
                            }
                        }
                        snp->getDescriptors()[3] = overlapValue;
                        snp->getDescriptors()[4] = neighborSum;
                    }
                    mean[3] += snp->getDescriptors()[3];
                    mean[4] += snp->getDescriptors()[4];
                }
                snps.push_back(move(snp));
            }
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
            sd[i] += (s->getDescriptors()[i] - mean[i])*(s->getDescriptors()[i] - mean[i]);
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
            fprintf(featuresFile, "%s\t", s->getId().c_str());
            fprintf(zscoreFile, "%s\t", s->getId().c_str());
        }
        for (i = 0; i < featNumber; i++) {
            x[i].index = i + 1;
            x[i].value = (s->getDescriptors()[i] - mean[i]) / sd[i];
            if (Global::instance()->isDebug3()) {
                fprintf(featuresFile, "%.4f\t", s->getDescriptors()[i]);
                fprintf(zscoreFile, "%.4f\t", x[i].value);
            }
        }
        if (Global::instance()->isDebug3()) {
            fprintf(featuresFile, "\n");
            fprintf(zscoreFile, "\n");
        }
        svmPredict.svmPredictCalulation(x, target_label);
        s->setProbPos(svmPredict.getProbEstimates()[0]);
    }

    if (Global::instance()->isDebug3()) {
        fclose(featuresFile);
        fclose(zscoreFile);
    }
    free(x);
    free(mean);
    free(sd);
    cout.precision(2);

    return count;
}


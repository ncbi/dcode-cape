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
#include <fstream>

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

    for (i = 0; i < featNumber; i++) {
        descriptors.push_back(0.0000);
    }

    if (static_cast<int> (pos - Global::instance()->getOrder() - 1) <= 0) {
        startPos = 0;
    } else {
        startPos = pos - Global::instance()->getOrder() + 1;
    }
    if (pos + Global::instance()->getOrder() >= length) {
        endPos = length - Global::instance()->getOrder();
    } else {
        endPos = pos;
    }

    for (i = 0; i <= length - Global::instance()->getOrder(); i++) {
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

    try {
        while (fParser.iterate('#', "\t")) {
            if (fParser.getNWords() != 5) {
                cerr << "Input SNP file with a wrong format " << endl;
                exit(-1);
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
                    cerr << "SNP position out of range" << endl;
                    exit(-1);
                }
                if (static_cast<int> (snpPos - neighbors - 1) < 0) {
                    startPos = 0;
                    snp->setPos(snpPos);
                } else {
                    startPos = snpPos - static_cast<int> (neighbors);
                    snp->setPos(neighbors);
                }
                if (snpPos + neighbors >= f->getLength()) {
                    endPos = f->getLength() - 1;
                } else {
                    endPos = snpPos + static_cast<int> (neighbors);
                }

                char *seq = f->getSubStr(startPos, endPos - startPos + 1);
                toUpper(&seq);
                snp->setSeq(&seq);
                snp->setLength(endPos - startPos + 1);

                if (Global::instance()->isDebug3()) {
                    cout << snp->getRef() << " " << snp->getSeq()[snp->getPos()]
                            << " " << snp->getPos()
                            << " " << snp->getLength()
                            << " " << snp->getSeq()
                            << endl;
                }

                if (snp->getSeq()[snp->getPos()] != snp->getRef()) {
                    cerr << "\nERROR1:\n\n" << snp->getSeq()[snp->getPos()] << " != " << fParser.getWords()[3][0] << "\n\n";
                    cerr << snp->getRef() << "\t" << snp->getSeq()[snp->getPos()] << "\t" << snp->getPos() << "\n\n";
                    cerr << snp->getSeq() << "\n\n";
                    if (Global::instance()->isDebug3()) {
                        cerr << "SNP is not in the chromosome position provided" << endl;
                        exit(-1);
                    }
                    delete snp;
                } else {
                    snps.push_back(snp);
                }
            }
        }
    } catch (exceptions::FileNotFoundException ex) {
        cerr << ex.what() << endl;
        cerr << "Error parsing file" << endl;
        exit(-1);
    } catch (exceptions::ErrorReadingFromFileException ex) {
        cerr << ex.what() << endl;
        cerr << "Error parsing file" << endl;
        exit(-1);
    } catch (exceptions::NotFoundException ex) {
        cerr << ex.what() << endl;
        cerr << "Error retrieving sequences" << endl;
        exit(-1);
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
        seq = strndup(snp->getSeq(), snp->getLength());
        f->setLength(snp->getLength());
        f->setSeq(&seq);
        fastaFactory.getSequenceContainter().insert(pair<string, Seq*>(id, f));
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

    ofstream featuresFile;
    ofstream zscoreFile;

    std::pair<std::string, double> cPair;

    unsigned long int featNumber = 3;
    double target_label = 0.0;

    if (Global::instance()->isDebug3()) {
        featuresFile.open("CAPE_debug_3_features.txt");
        if (!featuresFile.is_open()) {
            cerr << "Can't open features debug file named: CAPE_debug_3_features.txt" << endl;
            exit(-1);
        }
        featuresFile.precision(8);
        zscoreFile.open("CAPE_debug_3_zscores.txt");
        if (!zscoreFile.is_open()) {
            cerr << "Can't open zscores debug file named: CAPE_debug_3_zscores.txt" << endl;
            exit(-1);
        }
        zscoreFile.precision(8);
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
    try {
        while (fParser.iterate('#', "\t")) {
            if (fParser.getNWords() != 5) {
                cerr << "Input SNP file with a wrong format" << endl;
                exit(-1);
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
                    cerr << "SNP position out of range" << endl;
                    exit(-1);
                }
                if (static_cast<int> (snpPos - neighbors - 1) < 0) {
                    startPos = 0;
                    snp->setPos(snpPos);
                } else {
                    startPos = snpPos - static_cast<int> (neighbors);
                    snp->setPos(neighbors);
                }
                if (snpPos + neighbors >= f->getLength()) {
                    endPos = f->getLength() - 1;
                } else {
                    endPos = snpPos + static_cast<int> (neighbors);
                }

                try {
                    char *seq = f->getSubStr(startPos, endPos - startPos + 1);
                    toUpper(&seq);
                    snp->setSeq(&seq);
                    snp->setLength(endPos - startPos + 1);
                } catch (exceptions::OutOfRangeException ex) {
                    cerr << "Sequence length to retrieve is out of range" << endl;
                    exit(-1);
                }
                if (snp->getSeq()[snp->getPos()] != snp->getRef()) {
                    cerr << "\nERROR1:\n\n" << snp->getSeq()[snp->getPos()] << " != " << fParser.getWords()[3][0] << "\n\n";
                    cerr << snp->getRef() << "\t" << snp->getSeq()[snp->getPos()] << "\t" << snp->getPos() << "\n\n";
                    cerr << snp->getSeq() << "\n\n";
                    if (Global::instance()->isDebug3()) {
                        cerr << "SNP is not in the chromosome position provided" << endl;
                        exit(-1);
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
                            try {
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
                            } catch (exceptions::NotFoundException ex) {
                                cerr << ex.what() << endl;
                            }

                        }
                        mean[3] += snp->getDescriptors()[3];
                        mean[4] += snp->getDescriptors()[4];
                    }
                    snps.push_back(snp);
                }
            }
        }
    } catch (exceptions::FileNotFoundException ex) {
        cerr << ex.what() << endl;
        cerr << "Error parsing file" << endl;
        exit(-1);
    } catch (exceptions::ErrorReadingFromFileException ex) {
        cerr << ex.what() << endl;
        cerr << "Error parsing file" << endl;
        exit(-1);
    } catch (exceptions::NotFoundException ex) {
        cerr << ex.what() << endl;
        cerr << "Error retrieving sequences" << endl;
        exit(-1);
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
            featuresFile << s->getId() << "\t";
            zscoreFile << s->getId() << "\t";
        }
        for (i = 0; i < featNumber; i++) {
            x[i].index = i + 1;
            x[i].value = (s->getDescriptors()[i] - mean[i]) / sd[i];
            if (Global::instance()->isDebug3()) {
                featuresFile << s->getDescriptors()[i] << "\t";
                zscoreFile << x[i].value << "\t";
            }
        }
        if (Global::instance()->isDebug3()) {
            featuresFile << endl;
            zscoreFile << endl;
        }
        svmPredict.svmPredictCalulation(x, target_label);
        s->setProbPos(svmPredict.getProbEstimates()[0]);
    }

    if (Global::instance()->isDebug3()) {
        featuresFile.close();
        zscoreFile.close();
    }
    free(x);
    free(mean);
    free(sd);
    cout.precision(2);

    return count;
}


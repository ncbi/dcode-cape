/* 
 * File:   SNPFactory.cpp
 * Author: veraalva
 * 
 * Created on February 25, 2016, 9:22 AM
 */

#include <math.h>

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <vector>
#include <algorithm>
#include <fstream>

#include "berror.h"
#include "bmemory.h"
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
    this->pos = 0;
    this->chrPos = 0;
    this->ref = 0;
    this->alt = 0;
    this->probPos = 0.0;
}

SNP::~SNP() {
}

void SNP::calculateKmerDescriptors(kmers::KmersFactory& kmersFactory, unsigned long int featNumber) {
    unsigned long int i, j, k, startPos, endPos;

    for (i = 0; i < featNumber; i++) {
        descriptors.push_back(0.0000);
    }

    k = 0;
    for (unsigned int kmerIndex = 0; kmerIndex < kmersFactory.getKmerNumber(); kmerIndex++) {
        for (auto oIt = Global::instance()->getOrders().begin(); oIt != Global::instance()->getOrders().end(); ++oIt) {
            unsigned long int order = *oIt;
            double overlapMutated = 0.0;
            double overlapping = 0.0;
            double neighboring = 0.0;
            if (static_cast<int> (pos - order - 1) <= 0) {
                startPos = 0;
            } else {
                startPos = pos - order + 1;
            }
            if (pos + order >= length) {
                endPos = length - order;
            } else {
                endPos = pos;
            }

            j = order - 1;
            for (i = 0; i <= length - order; i++) {
                string sub(seq.c_str() + i, order);

                if (i >= startPos && i <= endPos) {
                    // Overlapping region
                    overlapping += kmersFactory.getKmerSig(sub, kmerIndex);
                    sub[j] = alt;
                    overlapMutated += kmersFactory.getKmerSig(sub, kmerIndex);
                    j--;
                } else {
                    // Neighboring region
                    neighboring += kmersFactory.getKmerSig(sub, kmerIndex);
                }
            }
            descriptors[k] = std::fabs(overlapping - overlapMutated);
            descriptors[k + 1] = overlapping + neighboring;
            k += 2;
        }
    }
}

SNPFactory::SNPFactory() {
}

SNPFactory::SNPFactory(const SNPFactory& orig) {
}

SNPFactory::~SNPFactory() {
}

void SNPFactory::parseSNPFile(std::string snpFileName, unsigned long int neighbors, FastaFactory &chrFactory) {
    bool process = true;
    FileParserFactory fParser;
    std::shared_ptr<Seq> f;
    shared_ptr<SNP> snp;
    int snpPos, startPos, endPos;

    try {
        f = chrFactory.getFirstSequence();
    } catch (exceptions::NotFoundException) {
        cerr << "Not chromosome sequences loaded" << endl;
        exit(-1);
    }

    try {
        fParser.setFileToParse(snpFileName);
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() != 5) {
                cerr << "Input SNP file with a wrong format " << endl;
                exit(-1);
            }
            if (f->getId().compare(fParser.getWords()[0]) != 0) {
                try {
                    f = chrFactory.getSequenceFromID(fParser.getWords()[0]);
                    process = true;
                } catch (exceptions::NotFoundException) {
                    process = false;
                }
            }
            if (process) {
                snp = make_shared<SNP>();
                snp->setId(fParser.getWords()[2]);
                snp->setChr(fParser.getWords()[0]);

                snp->setRef(fParser.getWords()[3].at(0));
                snp->setAlt(fParser.getWords()[4].at(0));
                snpPos = atoi((fParser.getWords()[1]).c_str()) - 1;
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
                    snp->setSeq(f->getSubStr(startPos, endPos - startPos + 1));
                    snp->setLength(endPos - startPos + 1);

                    if (Global::instance()->isDebug3()) {
                        cout << snp->getRef() << " " << snp->getSeq()[snp->getPos()]
                                << " " << snp->getPos()
                                << " " << snp->getLength()
                                << " " << snp->getSeq()
                                << endl;
                    }

                    if (snp->getSeq().at(snp->getPos()) != snp->getRef()) {
                        cerr << "\nERROR1:\n\n" << snp->getSeq().at(snp->getPos()) << " != " << fParser.getWords()[3].at(0) << "\n\n";
                        cerr << snp->getRef() << "\t" << snp->getSeq().at(snp->getPos()) << "\t" << snp->getPos() << "\n\n";
                        cerr << snp->getSeq() << "\n\n";
                        if (Global::instance()->isDebug3()) {
                            cerr << "SNP is not in the chromosome position provided" << endl;
                            exit(-1);
                        }
                    } else {
                        snps.push_back(snp);
                    }
                } catch (std::out_of_range) {
                    cerr << "Out of range coordinates for sequence." << endl;
                    exit(-1);
                }

            }
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file: " << snpFileName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file: " << snpFileName << endl;
        exit(-1);
    }
}

void SNPFactory::writeEnhansersFastaFile(std::string fastaFile, bool binary) {
    FastaFactory fastaFactory;
    shared_ptr<SNP> snp;
    string id;

    for (auto it = snps.begin(); it != snps.end(); ++it) {
        snp = *it;
        shared_ptr<Seq> f = make_shared<Seq>();
        id = snp->getId();
        f->setId(id);
        f->setSeq(snp->getSeq());
        fastaFactory.getSequenceContainter().insert(pair<string, shared_ptr < Seq >> (id, f));
    }

    fastaFactory.writeSequencesToFile(fastaFile, binary);
}

int SNPFactory::processSNPFromFile(std::string snpFileName, unsigned long int neighbors, FastaFactory &chrFactory, KmersFactory& kmersFactory, SVMPredict& svmPredict, FimoFactory & fimoFactory, TFBSFactory & tFBSFactory, std::string outputFileName) {
    bool process = true;
    FileParserFactory fParser;
    double overlapValue, neighborSum;
    unsigned long int i, count = 0;

    std::shared_ptr<Seq> f;
    shared_ptr<SNP> snp;
    int snpPos, startPos, endPos;

    ofstream featuresFile;
    ofstream zscoreFile;

    std::pair<std::string, double> cPair;

    unsigned long int baseFeatNumber = 2 * Global::instance()->getOrders().size() * kmersFactory.getKmerNumber();
    unsigned long int featNumber = baseFeatNumber;
    double target_label = 0.0;

    if (!fimoFactory.getSnpIDContainer().empty() || tFBSFactory.isReady()) {
        featNumber += 2;
    }

    struct svm_node *x = (struct svm_node *) allocate(sizeof (struct svm_node) * (featNumber + 1), __FILE__, __LINE__);
    x[featNumber].index = -1;

    vector<double> mean(featNumber);
    vector<double> sd(featNumber);
    for (i = 0; i < featNumber; i++) {
        mean[i] = 0.0;
        sd[i] = 0.0;
    }

    try {
        f = chrFactory.getFirstSequence();
    } catch (exceptions::NotFoundException) {
        cerr << "Not chromosome sequences loaded" << endl;
        exit(-1);
    }

    try {
        fParser.setFileToParse(snpFileName);
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() != 5) {
                cerr << "Input SNP file with a wrong format" << endl;
                exit(-1);
            }

            if (f->getId().compare(fParser.getWords()[0]) != 0) {
                try {
                    f = chrFactory.getSequenceFromID(fParser.getWords()[0]);
                    process = true;
                } catch (exceptions::NotFoundException) {
                    process = false;
                }
            }

            if (process) {
                snp = make_shared<SNP>();
                snp->setId(fParser.getWords()[2]);
                snp->setChr(fParser.getWords()[0]);
                snp->setRef(fParser.getWords()[3].at(0));
                snp->setAlt(fParser.getWords()[4].at(0));
                snpPos = atoi((fParser.getWords()[1]).c_str()) - 1;

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
                    snp->setSeq(f->getSubStr(startPos, endPos - startPos + 1));
                    snp->setLength(endPos - startPos + 1);

                    if (snp->getSeq().at(snp->getPos()) != snp->getRef()) {
                        cerr << "\nERROR1:\n\n" << snp->getSeq().at(snp->getPos()) << " != " << fParser.getWords()[3].at(0) << "\n\n";
                        cerr << snp->getRef() << "\t" << snp->getSeq().at(snp->getPos()) << "\t" << snp->getPos() << "\n\n";
                        cerr << snp->getSeq() << "\n\n";
                        if (Global::instance()->isDebug3()) {
                            cerr << "SNP is not in the chromosome position provided" << endl;
                            exit(-1);
                        }
                    } else {
                        count++;
                        snp->calculateKmerDescriptors(kmersFactory, featNumber);

                        /*
                         * Calculating overall sum for mean and sd
                         */
                        for (i = 0; i < baseFeatNumber; i++) {
                            mean[i] += snp->getDescriptors()[i];
                        }

                        if (featNumber == baseFeatNumber + 2) {
                            if (!fimoFactory.getSnpIDContainer().empty()) {
                                auto fimoMapIt = fimoFactory.getSnpIDContainer().find(snp->getId());
                                if (fimoMapIt != fimoFactory.getSnpIDContainer().end()) {
                                    snp->getDescriptors()[baseFeatNumber] = fimoMapIt->second[0];
                                    snp->getDescriptors()[baseFeatNumber + 1] = fimoMapIt->second[1];
                                }
                            } else {
                                overlapValue = neighborSum = 0;
                                try {
                                    tFBSFactory.extractTFBSFromFile(startPos, endPos, f);
                                    for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
                                        shared_ptr<TFBS> t = *it;
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
                                    snp->getDescriptors()[baseFeatNumber] = overlapValue;
                                    snp->getDescriptors()[baseFeatNumber + 1] = neighborSum;
                                } catch (exceptions::NotFoundException) {
                                    cerr << "Not indexes available for " << f->getId() << ". Ignoring" << endl;
                                }

                            }
                            mean[baseFeatNumber] += snp->getDescriptors()[baseFeatNumber];
                            mean[baseFeatNumber + 1] += snp->getDescriptors()[baseFeatNumber + 1];
                        }
                        snps.push_back(snp);
                    }
                } catch (std::out_of_range) {
                    cerr << "Out of range coordinates for sequence." << endl;
                    exit(-1);
                }
            }
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file:" << snpFileName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file:" << snpFileName << endl;
        exit(-1);
    }

    if (!snps.empty()) {

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
            shared_ptr<SNP> s = *it;
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

        if (Global::instance()->isDebug3()) {
            featuresFile.open(outputFileName + "_CAPE_features.txt");
            if (!featuresFile.is_open()) {
                cerr << "Can't open features debug file named: " + outputFileName + "_CAPE_features.txt" << endl;
                exit(-1);
            }
            featuresFile.precision(8);
            zscoreFile.open(outputFileName + "_CAPE_zscores.txt");
            if (!zscoreFile.is_open()) {
                cerr << "Can't open zscores debug file named: " + outputFileName + "_CAPE_zscores.txt" << endl;
                exit(-1);
            }
            zscoreFile.precision(8);
        }

        /*
         * Calculating ZScore terms
         */
        for (auto it = snps.begin(); it != snps.end(); ++it) {
            shared_ptr<SNP>s = *it;

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
    }
    free(x);
    return count;
}

int SNPFactory::createSVMModelFromSNPFile(std::string snpFileName, unsigned long int neighbors, FastaFactory &chrFactory, KmersFactory& kmersFactory, svm::SVMPredict& svmPredict, FimoFactory & fimoFactory, TFBSFactory & tFBSFactory, std::string &outputFileName) {
    bool process = true;
    FileParserFactory fParser;
    double overlapValue, neighborSum;
    unsigned long int i, j, count = 0;

    std::shared_ptr<Seq> f;
    shared_ptr<SNP> snp;
    int snpPos, startPos, endPos;

    ofstream featuresFile;
    ofstream zscoreFile;

    std::pair<std::string, double> cPair;

    unsigned long int baseFeatNumber = 2 * Global::instance()->getOrders().size() * kmersFactory.getKmerNumber();
    unsigned long int featNumber = baseFeatNumber;

    if (!fimoFactory.getSnpIDContainer().empty() || tFBSFactory.isReady()) {
        featNumber += 2;
    }

    vector<double> y;
    struct svm_node **x = NULL;

    vector<double> mean(featNumber);
    vector<double> sd(featNumber);
    for (i = 0; i < featNumber; i++) {
        mean[i] = 0.0;
        sd[i] = 0.0;
    }

    try {
        f = chrFactory.getFirstSequence();
    } catch (exceptions::NotFoundException) {
        cerr << "Not chromosome sequences loaded" << endl;
        exit(-1);
    }

    try {
        fParser.setFileToParse(snpFileName);
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() != 6) {
                cerr << "Input SNP file with a wrong format. It has " << fParser.getWords().size() << " words ans should be 6" << endl;
                cerr << fParser.getLine() << endl;
                exit(-1);
            }

            if (f->getId().compare(fParser.getWords()[0]) != 0) {
                try {
                    f = chrFactory.getSequenceFromID(fParser.getWords()[0]);
                    process = true;
                } catch (exceptions::NotFoundException) {
                    process = false;
                }
            }

            if (process) {
                snp = make_shared<SNP>();
                snp->setId(fParser.getWords()[2]);
                snp->setChr(fParser.getWords()[0]);
                snp->setRef(fParser.getWords()[3].at(0));
                snp->setAlt(fParser.getWords()[4].at(0));
                snpPos = atoi((fParser.getWords()[1]).c_str()) - 1;

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
                    snp->setSeq(f->getSubStr(startPos, endPos - startPos + 1));
                    snp->setLength(endPos - startPos + 1);

                    if (snp->getSeq().at(snp->getPos()) != snp->getRef()) {
                        cerr << "\nERROR1:\n\n" << snp->getSeq().at(snp->getPos()) << " != " << fParser.getWords()[3].at(0) << "\n\n";
                        cerr << snp->getRef() << "\t" << snp->getSeq().at(snp->getPos()) << "\t" << snp->getPos() << "\n\n";
                        cerr << snp->getSeq() << "\n\n";
                        if (Global::instance()->isDebug3()) {
                            cerr << "SNP is not in the chromosome position provided" << endl;
                            exit(-1);
                        }
                    } else {
                        count++;
                        snp->calculateKmerDescriptors(kmersFactory, featNumber);

                        /*
                         * Calculating overall sum for mean and sd
                         */
                        for (i = 0; i < baseFeatNumber; i++) {
                            mean[i] += snp->getDescriptors()[i];
                        }

                        if (featNumber == baseFeatNumber + 2) {
                            if (!fimoFactory.getSnpIDContainer().empty()) {
                                auto fimoMapIt = fimoFactory.getSnpIDContainer().find(snp->getId());
                                if (fimoMapIt != fimoFactory.getSnpIDContainer().end()) {
                                    snp->getDescriptors()[baseFeatNumber] = fimoMapIt->second[0];
                                    snp->getDescriptors()[baseFeatNumber + 1] = fimoMapIt->second[1];
                                }
                            } else {
                                overlapValue = neighborSum = 0;
                                try {
                                    tFBSFactory.extractTFBSFromFile(startPos, endPos, f);
                                    if (Global::instance()->isDebug3()) {
                                        cerr << "startPos: " << startPos << "\tendPos: " << endPos << endl;
                                        cerr << "TBSF: " << tFBSFactory.getTfbs().size() << endl;
                                        cerr << "ID: " << f->getId() << endl;
                                        cerr << "Seq: " << f->getSeq() << endl;
                                    }
                                    for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
                                        shared_ptr<TFBS> t = *it;
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
                                    snp->getDescriptors()[baseFeatNumber] = overlapValue;
                                    snp->getDescriptors()[baseFeatNumber + 1] = neighborSum;
                                } catch (exceptions::NotFoundException) {
                                    cerr << "Not indexes available for " << f->getId() << ". Ignoring" << endl;
                                }

                            }
                            mean[baseFeatNumber] += snp->getDescriptors()[baseFeatNumber];
                            mean[baseFeatNumber + 1] += snp->getDescriptors()[baseFeatNumber + 1];
                        }
                        y.push_back(atof(fParser.getWords()[5].c_str()));
                        snps.push_back(snp);
                    }
                } catch (std::out_of_range) {
                    cerr << "Out of range coordinates for sequence." << endl;
                    exit(-1);
                }
            }
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file:" << snpFileName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file:" << snpFileName << endl;
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
        shared_ptr<SNP> s = *it;
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

    if (Global::instance()->isDebug3()) {
        featuresFile.open(outputFileName + "_CAPE_features.txt");
        if (!featuresFile.is_open()) {
            cerr << "Can't open features debug file named: " + outputFileName + "_CAPE_features.txt" << endl;
            exit(-1);
        }
        featuresFile.precision(8);
        zscoreFile.open(outputFileName + "_CAPE_zscores.txt");
        if (!zscoreFile.is_open()) {
            cerr << "Can't open zscores debug file named: " + outputFileName + "_CAPE_zscores.txt" << endl;
            exit(-1);
        }
        zscoreFile.precision(8);
    }

    /*
     * Calculating ZScore terms and creating the SVM prob data
     */
    j = 0;
    x = (struct svm_node **) allocate(sizeof (struct svm_node *) * snps.size(), __FILE__, __LINE__);
    for (auto it = snps.begin(); it != snps.end(); ++it) {
        shared_ptr<SNP>s = *it;

        x[j] = (struct svm_node *) allocate(sizeof (struct svm_node) * (featNumber + 1), __FILE__, __LINE__);
        x[j][featNumber].index = -1;

        if (Global::instance()->isDebug3()) {
            featuresFile << s->getId() << "\t";
            zscoreFile << s->getId() << "\t";
        }
        for (i = 0; i < featNumber; i++) {
            x[j][i].index = i + 1;
            x[j][i].value = (s->getDescriptors()[i] - mean[i]) / sd[i];
            if (Global::instance()->isDebug3()) {
                featuresFile << s->getDescriptors()[i] << "\t";
                zscoreFile << x[j][i].value << "\t";
            }
        }
        if (Global::instance()->isDebug3()) {
            featuresFile << endl;
            zscoreFile << endl;
        }
        j++;
    }

    if (x != NULL) {
        svmPredict.svmTrainModel(j, &y[0], x, outputFileName.c_str());
        for (i = 0; i < j; i++) {
            free(x[i]);
        }
        free(x);
    }
    if (Global::instance()->isDebug3()) {
        featuresFile.close();
        zscoreFile.close();
    }
    return count;
}
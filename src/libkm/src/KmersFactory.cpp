/* 
 * File:   Kmers.cpp
 * Author: veraalva
 * 
 * Created on February 17, 2016, 10:30 AM
 */

#include <math.h>
#include <inttypes.h>

#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include <set>

#include "bmath.h"
#include "Global.h"
#include "cstring.h"
#include "FileParserFactory.h"
#include "KmersFactory.h"

using namespace std;
using namespace parsers;
using namespace kmers;

Kmer::Kmer() {
    this->controlFreq = 0;
    this->negativePeak = 0;
    this->negativeControl = 0;
    this->peakFreq = 0;
    this->pValue = INFINITY;
    this->sig = 0.0;
    this->pf = 0.0;
}

Kmer::~Kmer() {
}

KmersFactory::KmersFactory() {
    this->kmerNumber = 0;
    this->totalNRnt_control = 0;
    this->totalNRnt_peak = 0;
}

KmersFactory::~KmersFactory() {
}

void KmersFactory::createGenomeWideKmers() {
    unsigned long int i, count = 1;
    set<string> newWords;
    set<string> kmersTemp;

    this->kmersGenome.clear();
    kmersTemp.clear();
    kmersTemp.insert("A");
    kmersTemp.insert("C");
    kmersTemp.insert("G");
    kmersTemp.insert("T");

    for (auto oIt = Global::instance()->getOrders().begin(); oIt != Global::instance()->getOrders().end(); ++oIt) {
        unsigned long int order = *oIt;
        for (i = count; i < order; i++) {
            newWords.clear();
            for (auto it = kmersTemp.begin(); it != kmersTemp.end(); ++it) {
                string key = (*it);
                newWords.insert(key + "A");
                newWords.insert(key + "C");
                newWords.insert(key + "G");
                newWords.insert(key + "T");
            }
            kmersTemp.clear();
            kmersTemp.insert(newWords.begin(), newWords.end());
        }
        for (auto it = newWords.begin(); it != newWords.end(); ++it) {
            string key = (*it);
            string rcKey = cstring::reverseComplement(key);
            if (this->kmersGenome.find(key) == this->kmersGenome.end() &&
                    this->kmersGenome.find(rcKey) == this->kmersGenome.end()) {
                this->kmersGenome.insert(key);
            }
        }
        count = i;
    }
    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> totalNumber of k-mers: " << this->kmersGenome.size() << endl;
    }
}

void Kmer::calculatePValue(double totalNRnt_peak, double totalNRnt_control) {
    /*
     * 2x2 Matrix
     * 
     *      peakFreq                        controlFreq
     *  totalNRnt_peak - peakFreq     totalNRnt_control - controlFreq
     * 
     * m is the sum over column 0
     * n is the sum over column 1
     * k is the sum over row 0
     * 
     * phyper(matrix[0][0], m, n, k, false, false);
     */
    this->pValue = phyper(this->peakFreq - 1, totalNRnt_peak, totalNRnt_control,
            this->peakFreq + this->controlFreq, false, false);
    this->sig = -1.0 * log10(this->pValue);
    this->pf = (static_cast<double> (this->peakFreq) / static_cast<double> (this->controlFreq)) /
            (static_cast<double> (totalNRnt_peak) / static_cast<double> (totalNRnt_control));
}

void KmersFactory::buildKmers() {
    unordered_map<string, std::vector<std::shared_ptr < Kmer>>>::iterator it;

    std::unordered_map<std::string, unsigned long int>::iterator controlFreq_it;
    for (auto peakFreq_it = this->kmer2peakFreq.begin(); peakFreq_it != this->kmer2peakFreq.end(); ++peakFreq_it) {
        string kmer = peakFreq_it->first;
        unsigned long int kmerPeakFreq = peakFreq_it->second;
        unsigned long int kmerControlFreq = 0;
        if ((controlFreq_it = this->kmer2controlFreq.find(kmer)) != this->kmer2controlFreq.end()) {
            kmerControlFreq = controlFreq_it->second;
        }
        shared_ptr<Kmer> k = make_shared<Kmer>();
        k->setPeakFreq(kmerPeakFreq);
        k->setControlFreq(kmerControlFreq);
        k->setNegativePeak(this->totalNRnt_peak - kmerPeakFreq);
        k->setNegativeControl(this->totalNRnt_control - kmerControlFreq);
        k->calculatePValue(this->totalNRnt_peak, this->totalNRnt_control);
        it = this->kmers.find(kmer);
        if (it != this->kmers.end()) {
            it->second.push_back(k);
        } else {
            std::vector<std::shared_ptr < Kmer>> kmerVec;
            kmerVec.push_back(k);
            this->kmers.insert(make_pair(kmer, kmerVec));
        }
    }
}

void KmersFactory::clearKmerControlData() {
    this->kmer2controlFreq.clear();
    this->totalNRnt_control = 0;
}

void KmersFactory::clearKmerPeakData() {
    this->kmer2peakFreq.clear();
    this->totalNRnt_peak = 0;
}

void KmersFactory::scanSequences(string inputSeq, bool control) {
    unsigned long int i;
    string kmer, rc_kmer, dealedKmer;
    for (auto oIt = Global::instance()->getOrders().begin(); oIt != Global::instance()->getOrders().end(); ++oIt) {
        int order = *oIt;
        for (i = 0; i < inputSeq.length() - order + 1; i++) {
            kmer = inputSeq.substr(i, order);
            rc_kmer = cstring::reverseComplement(kmer);
            dealedKmer.clear();
            if (this->kmersGenome.find(kmer) != this->kmersGenome.end()) {
                dealedKmer = kmer;
            } else if (this->kmersGenome.find(rc_kmer) != this->kmersGenome.end()) {
                dealedKmer = rc_kmer;
            }
            if (!dealedKmer.empty()) {
                if (control) {
                    if (this->kmer2controlFreq.find(dealedKmer) != this->kmer2controlFreq.end()) {
                        this->kmer2controlFreq[dealedKmer]++;
                    } else {
                        this->kmer2controlFreq[dealedKmer] = 1;
                    }
                } else {
                    if (this->kmer2peakFreq.find(dealedKmer) != this->kmer2peakFreq.end()) {
                        this->kmer2peakFreq[dealedKmer]++;
                    } else {
                        this->kmer2peakFreq[dealedKmer] = 1;
                    }
                }
            }
        }
    }

    if (control) {
        this->totalNRnt_control += (inputSeq.length() - cstring::countCharacter(inputSeq, "Nn"));
    } else {
        this->totalNRnt_peak += (inputSeq.length() - cstring::countCharacter(inputSeq, "Nn"));
    }
}

void KmersFactory::readKmersFromFile(std::string fileName) {
    shared_ptr<Kmer> k, kr;
    string rc_kmer;
    vector<double> maxSig;
    set<string> infSig;
    unordered_map<string, std::vector<std::shared_ptr < Kmer>>>::iterator it;

    FileParserFactory fParser;

    try {
        fParser.setFileToParse(fileName);
        while (fParser.iterate("#", "\t")) {
            if (this->kmerNumber == 0) {
                this->kmerNumber = (fParser.getWords().size() - 1) / 3;
                printf("kmerNumber %d\n", this->kmerNumber);
            }
            if (maxSig.size() == 0) {
                maxSig.resize(this->kmerNumber, NAN);
            }
            rc_kmer = cstring::reverseComplement(fParser.getWords()[0]);
            for (int i = 0; i < this->kmerNumber; i++) {
                k = make_shared<Kmer>();
                kr = make_shared<Kmer>();

                k->setValue(atof((fParser.getWords()[3 * i + 1]).c_str()));
                kr->setValue(k->getValue());

                k->setSig(atof((fParser.getWords()[3 * i + 2]).c_str()));
                kr->setSig(k->getSig());
                if (std::isinf(k->getSig())) {
                    infSig.insert(fParser.getWords()[0]);
                } else {
                    if (std::isnan(maxSig[i]) || maxSig[i] < std::fabs(k->getSig())) maxSig[i] = std::fabs(k->getSig());
                }

                k->setPf(atof((fParser.getWords()[3 * i + 3]).c_str()));
                kr->setPf(k->getPf());

                it = this->kmers.find(fParser.getWords()[0]);
                if (it != this->kmers.end()) {
                    it->second.push_back(k);
                } else {
                    std::vector<std::shared_ptr < Kmer>> kmerVec;
                    kmerVec.push_back(k);
                    this->kmers.insert(make_pair(fParser.getWords()[0], kmerVec));
                }

                if (rc_kmer.compare(fParser.getWords()[0]) != 0) {
                    it = this->kmers.find(rc_kmer);
                    if (it != this->kmers.end()) {
                        it->second.push_back(k);
                    } else {
                        std::vector<std::shared_ptr < Kmer>> kmerVec;
                        kmerVec.push_back(k);
                        this->kmers.insert(make_pair(rc_kmer, kmerVec));
                    }
                }
            }
            if (fParser.getWords().size() != 4) {
                cerr << "Input kmer weight file with a wrong format " << endl;
                exit(-1);
            }

        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file: " << fileName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file: " << fileName << endl;
        exit(-1);
    }

    if (!infSig.empty()) {
        if (Global::instance()->isInfo()) {
            cout << "\tINFO ==> There are " << infSig.size() << " infinities to fix" << endl;
        }
        for (auto it1 = infSig.begin(); it1 != infSig.end(); ++it1) {
            rc_kmer = cstring::reverseComplement(*it1);
            it = this->kmers.find(*it1);
            int i = 0;
            for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                if (std::isinf(it2->get()->getSig())) {
                    it2->get()->setSig(maxSig[i] + 10);
                }
                i++;
            }
            it = this->kmers.find(rc_kmer);
            i = 0;
            for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                if (std::isinf(it2->get()->getSig())) {
                    it2->get()->setSig(maxSig[i] + 10);
                }
                i++;
            }
        }
    }
}

void KmersFactory::writeKmersToFile(std::string fileName) {
    ofstream outputFile(fileName);
    if (!outputFile.is_open()) {
        cerr << "Can't open output file " << fileName << endl;
        exit(-1);
    }
    outputFile.precision(12);
    for (auto it = this->kmers.begin(); it != this->kmers.end(); ++it) {
        string kmer = it->first;
        outputFile << kmer;
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            shared_ptr<Kmer> k = *it2;
            outputFile << "\t" << k->getValue() << "\t" << k->getSig() << "\t" << k->getPf();
        }
        outputFile << endl;

    }
    outputFile.close();
}

double KmersFactory::getKmerSig(std::string kmer, unsigned int index) {
    std::unordered_map<std::string, std::vector<std::shared_ptr < Kmer>>>::iterator it;
    it = this->kmers.find(kmer);
    if (it == this->kmers.end()) {
        return 0.0;
    }
    if (it->second.size() <= index) {
        cerr << "Index number is bigger that the vector size" << endl;
        exit(-1);
    }
    return it->second[index]->getSig();
}

void KmersFactory::mergeKmers(KmersFactory& kmersFactory) {
    std::unordered_map<std::string, std::vector<std::shared_ptr < Kmer>>>::iterator destIt;
    for (auto srcIt = kmersFactory.getKmers().begin(); srcIt != kmersFactory.getKmers().end(); ++srcIt) {
        destIt = kmers.find(srcIt->first);
        if (destIt != kmers.end()) {
            destIt->second.insert(std::end(destIt->second), std::begin(srcIt->second), std::end(srcIt->second));
        } else {
            std::vector<std::shared_ptr < Kmer>> kmerVec;
            for (int i = 0; i < kmerNumber; i++) {
                shared_ptr<Kmer> k = make_shared<Kmer>();
                kmerVec.push_back(k);
            }
            kmerVec.insert(std::end(kmerVec), std::begin(srcIt->second), std::end(srcIt->second));
            kmers.insert(make_pair(srcIt->first, kmerVec));
        }
    }
    kmerNumber += kmersFactory.getKmerNumber();
}

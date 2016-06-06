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
        this->getKmers().insert(make_pair(kmer, k));
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
    double maxSig = NAN;
    vector<string> infSig;
    unordered_map<string, shared_ptr < Kmer>>::iterator it;

    FileParserFactory fParser;

    try {
        fParser.setFileToParse(fileName);
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() != 4) {
                cerr << "Input kmer weight file with a wrong format " << endl;
                exit(-1);
            }
            k = make_shared<Kmer>();
            kr = make_shared<Kmer>();
            rc_kmer = cstring::reverseComplement(fParser.getWords()[0]);

            k->setValue(atof((fParser.getWords()[1]).c_str()));
            kr->setValue(k->getValue());

            k->setSig(atof((fParser.getWords()[2]).c_str()));
            kr->setSig(k->getSig());
            if (std::isinf(k->getSig())) {
                infSig.push_back(fParser.getWords()[0]);
            } else {
                if (std::isnan(maxSig) || maxSig < std::fabs(k->getSig())) maxSig = std::fabs(k->getSig());
            }

            k->setPf(atof((fParser.getWords()[3]).c_str()));
            kr->setPf(k->getPf());

            this->kmers.insert(make_pair(fParser.getWords()[0], k));

            if (rc_kmer.compare(fParser.getWords()[0]) != 0) {
                this->kmers.insert(make_pair(rc_kmer, kr));
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
            cout << "\tINFO ==> There are " << infSig.size() << " infinities to fix using value: " << maxSig + 10 << endl;
        }
        for (auto it1 = infSig.begin(); it1 != infSig.end(); ++it1) {
            rc_kmer = cstring::reverseComplement(*it1);
            it = this->kmers.find(*it1);
            it->second->setSig(maxSig + 10);
            it = this->kmers.find(rc_kmer);
            it->second->setSig(maxSig + 10);
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
        shared_ptr<Kmer> k = it->second;
        outputFile << kmer << "\t" << k->getValue() << "\t" << k->getSig() << "\t" << k->getPf() << endl;
    }
    outputFile.close();
}

double KmersFactory::getKmerSig(std::string kmer) {
    std::unordered_map<std::string, shared_ptr < Kmer>>::iterator it;
    it = this->kmers.find(kmer);
    if (it == this->kmers.end()) {
        return 0.0;
    }
    return it->second->getSig();
}


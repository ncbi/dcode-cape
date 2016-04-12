/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Kmers.cpp
 * Author: veraalva
 * 
 * Created on February 17, 2016, 10:30 AM
 */

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>
#include <memory>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <vector>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "bmath.h"
#include "Global.h"
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
    for (auto it = this->getKmers().begin(); it != this->getKmers().end(); ++it) {
        delete (it->second);
    }
}

void KmersFactory::createGenomeWideKmers() {
    char *comp;
    string key, key1, rcKey1;
    set<string> newWords;
    set<string> ::iterator it, it1;

    this->kmersGenome.clear();
    this->kmersGenome.insert("A");
    this->kmersGenome.insert("C");
    this->kmersGenome.insert("G");
    this->kmersGenome.insert("T");

    for (unsigned long int i = 1; i < Global::instance()->getOrder(); i += 1) {
        newWords.clear();
        for (it = this->kmersGenome.begin(); it != this->kmersGenome.end(); ++it) {
            key = (*it);
            newWords.insert(key + "A");
            newWords.insert(key + "C");
            newWords.insert(key + "G");
            newWords.insert(key + "T");
        }
        this->kmersGenome.clear();
        this->kmersGenome.insert(newWords.begin(), newWords.end());
    }
    this->kmersGenome.clear();
    this->kmersGenome.insert(newWords.begin(), newWords.end());
    newWords.clear();

    for (it1 = this->kmersGenome.begin(); it1 != this->kmersGenome.end(); ++it1) {//this->kmers now are redundant
        key1 = (*it1);
        comp = complement(key1.c_str());
        rcKey1 = comp;
        free(comp);
        reverse(rcKey1.begin(), rcKey1.end());
        if (newWords.find(key1) != newWords.end()) {
            continue;
        }
        if (newWords.find(rcKey1) != newWords.end()) {
            continue;
        }
        newWords.insert(key1);
    }
    this->kmersGenome.clear();
    this->kmersGenome.insert(newWords.begin(), newWords.end());
    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> totalNumber of " << Global::instance()->getOrder() << "-mers: " << this->kmersGenome.size() << endl;
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

    std::map<std::string, unsigned long int>::iterator controlFreq_it;
    for (auto peakFreq_it = this->kmer2peakFreq.begin(); peakFreq_it != this->kmer2peakFreq.end(); ++peakFreq_it) {
        string kmer = peakFreq_it->first;
        unsigned long int kmerPeakFreq = peakFreq_it->second;
        unsigned long int kmerControlFreq = 0;
        if ((controlFreq_it = this->kmer2controlFreq.find(kmer)) != this->kmer2controlFreq.end()) {
            kmerControlFreq = controlFreq_it->second;
        }
        Kmer *k = new Kmer();
        k->setPeakFreq(kmerPeakFreq);
        k->setControlFreq(kmerControlFreq);
        k->setNegativePeak(this->totalNRnt_peak - kmerPeakFreq);
        k->setNegativeControl(this->totalNRnt_control - kmerControlFreq);
        k->calculatePValue(this->totalNRnt_peak, this->totalNRnt_control);
        this->getKmers().insert(pair<std::string, Kmer *>(kmer, move(k)));
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
    char *comp;
    string kmer, rc_kmer, dealedKmer;
    for (i = 0; i < inputSeq.length() - Global::instance()->getOrder() + 1; i++) {
        kmer = inputSeq.substr(i, Global::instance()->getOrder());
        comp = complement(kmer.c_str());
        rc_kmer = comp;
        free(comp);
        reverse(rc_kmer.begin(), rc_kmer.end());
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

    if (control) {
        this->totalNRnt_control += (inputSeq.length() - countCharacter(inputSeq.c_str(), string("Nn").c_str()));
    } else {
        this->totalNRnt_peak += (inputSeq.length() - countCharacter(inputSeq.c_str(), string("Nn").c_str()));
    }
}

void KmersFactory::readKmersFromFile(std::string fileName, bool binary) {
    unsigned long int i, index;
    double value;
    FILE *poutputFile = NULL;
    Kmer *k, *kr;
    char *seq;
    char *comp;
    string rc_kmer;
    double maxSig = NAN;
    char **infSig = NULL;
    size_t infSigSize = 0;
    map<string, Kmer *>::iterator it;

    if (binary) {
        poutputFile = (FILE *) checkPointerError(fopen(fileName.c_str(), "rb"), "Can't open output file", __FILE__, __LINE__, -1);
        fread(&i, sizeof (unsigned long int), 1, poutputFile);
        Global::instance()->setOrder(i);
        fread(&index, sizeof (unsigned long int), 1, poutputFile);
        while (index) {
            seq = (char *) allocate(sizeof (char) * (Global::instance()->getOrder() + 1), __FILE__, __LINE__);
            k = new Kmer();
            kr = new Kmer();

            fread(seq, sizeof (char), Global::instance()->getOrder() + 1, poutputFile);
            comp = complement(seq);
            rc_kmer = comp;
            free(comp);
            reverse(rc_kmer.begin(), rc_kmer.end());

            fread(&value, sizeof (double), 1, poutputFile);
            k->setValue(value);
            kr->setValue(value);

            fread(&value, sizeof (double), 1, poutputFile);
            k->setSig(value);
            kr->setSig(value);
            if (std::isinf(k->getSig())) {
                infSig = (char **) reallocate(infSig, sizeof (char **) * (infSigSize + 1), __FILE__, __LINE__);
                infSig[infSigSize] = strdup(seq);
                infSigSize++;
            } else {
                if (std::isnan(maxSig) || maxSig < std::fabs(k->getSig())) maxSig = std::fabs(k->getSig());
            }

            fread(&value, sizeof (double), 1, poutputFile);
            k->setPf(value);
            kr->setPf(value);

            fread(&i, sizeof (unsigned long int), 1, poutputFile);
            k->setControlFreq(i);
            kr->setControlFreq(i);

            fread(&i, sizeof (unsigned long int), 1, poutputFile);
            k->setNegativeControl(i);
            kr->setNegativeControl(i);

            fread(&i, sizeof (unsigned long int), 1, poutputFile);
            k->setNegativePeak(i);
            kr->setNegativePeak(i);

            fread(&i, sizeof (unsigned long int), 1, poutputFile);
            k->setPeakFreq(i);
            kr->setPeakFreq(i);

            this->kmers.insert(pair<string, Kmer *>(seq, move(k)));
            if (rc_kmer.compare(seq) != 0) {
                this->kmers.insert(pair<string, Kmer *>(rc_kmer, move(kr)));
            }
            index--;
            free(seq);
        }
        fclose(poutputFile);
    } else {
        FileParserFactory fParser(fileName);

        while (fParser.iterate('#', "\t")) {
            if (fParser.getNWords() != 4) {
                printLog(stderr, "Input kmer weight file with a wrong format ", __FILE__, __LINE__, -1);
            }
            k = new Kmer();
            kr = new Kmer();
            seq = strdup(fParser.getWords()[0]);
            comp = complement(seq);
            rc_kmer = comp;
            free(comp);
            reverse(rc_kmer.begin(), rc_kmer.end());

            k->setValue(strtod(fParser.getWords()[1], NULL));
            kr->setValue(k->getValue());

            k->setSig(strtod(fParser.getWords()[2], NULL));
            kr->setSig(k->getSig());
            if (std::isinf(k->getSig())) {
                infSig = (char **) reallocate(infSig, sizeof (char **) * (infSigSize + 1), __FILE__, __LINE__);
                infSig[infSigSize] = strdup(seq);
                infSigSize++;
            } else {
                if (std::isnan(maxSig) || maxSig < std::fabs(k->getSig())) maxSig = std::fabs(k->getSig());
            }

            k->setPf(strtod(fParser.getWords()[3], NULL));
            kr->setPf(k->getPf());

            this->kmers.insert(pair<string, Kmer *>(seq, move(k)));

            if (rc_kmer.compare(seq) != 0) {
                this->kmers.insert(pair<string, Kmer *>(rc_kmer, move(kr)));
            }
            free(seq);
        }
    }


    if (infSig) {
        if (Global::instance()->isInfo()) {
            streamsize initPrecision = std::cout.precision();
            std::cout.precision(12);
            cout << "\tINFO ==> There are " << infSigSize << " infinities to fix using value: " << maxSig + 10 << endl;
            std::cout.precision(initPrecision);
        }
        for (i = 0; i < infSigSize; i++) {
            comp = complement(infSig[i]);
            rc_kmer = comp;
            free(comp);
            reverse(rc_kmer.begin(), rc_kmer.end());

            if (Global::instance()->isDebug3()) {
                cout << "\tDEBUG3 ==> " << infSig[i] << endl;
                cout << "\tDEBUG3 ==> " << rc_kmer << endl;
            }

            it = this->kmers.find(infSig[i]);
            it->second->setSig(maxSig + 10);
            it = this->kmers.find(rc_kmer);
            it->second->setSig(maxSig + 10);

            free(infSig[i]);
        }
        free(infSig);
    }
}

void KmersFactory::writeKmersToFile(std::string fileName, bool binary) {
    unsigned long int i;
    double value;
    FILE *poutputFile = NULL;

    if (binary) {
        poutputFile = (FILE *) checkPointerError(fopen(fileName.c_str(), "wb"), "Can't open output file", __FILE__, __LINE__, -1);
        i = Global::instance()->getOrder();
        fwrite(&i, sizeof (unsigned long int), 1, poutputFile);
        i = this->kmers.size();
        fwrite(&i, sizeof (unsigned long int), 1, poutputFile);
        for (auto it = this->kmers.begin(); it != this->kmers.end(); ++it) {
            string kmer = it->first;
            Kmer *k = it->second;
            fwrite(kmer.c_str(), sizeof (char), Global::instance()->getOrder() + 1, poutputFile);
            value = k->getValue();
            fwrite(&value, sizeof (double), 1, poutputFile);
            value = k->getSig();
            fwrite(&value, sizeof (double), 1, poutputFile);
            value = k->getPf();
            fwrite(&value, sizeof (double), 1, poutputFile);
            i = k->getControlFreq();
            fwrite(&i, sizeof (unsigned long int), 1, poutputFile);
            i = k->getNegativeControl();
            fwrite(&i, sizeof (unsigned long int), 1, poutputFile);
            i = k->getNegativePeak();
            fwrite(&i, sizeof (unsigned long int), 1, poutputFile);
            i = k->getPeakFreq();
            fwrite(&i, sizeof (unsigned long int), 1, poutputFile);
        }
    } else {
        poutputFile = (FILE *) checkPointerError(fopen(fileName.c_str(), "w"), "Can't open output file", __FILE__, __LINE__, -1);
        streamsize initPrecision = std::cout.precision();
        std::cout.precision(12);
        for (auto it = this->kmers.begin(); it != this->kmers.end(); ++it) {
            string kmer = it->first;
            Kmer *k = it->second;
            fprintf(poutputFile, "%s\t%.14e\t%.14f\t%.12f\n", kmer.c_str(), k->getValue(), k->getSig(), k->getPf());
        }
        std::cout.precision(initPrecision);
    }

    fclose(poutputFile);
}

double KmersFactory::getKmerSig(std::string kmer) {
    std::map<std::string, Kmer *>::iterator it;
    it = this->kmers.find(kmer);
    if (it == this->kmers.end()) {
        return 0.0;
    }
    return it->second->getSig();
}


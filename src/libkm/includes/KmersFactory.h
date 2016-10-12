/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Kmers.h
 * Author: veraalva
 *
 * Created on February 17, 2016, 10:30 AM
 */

#ifndef KMERS_H
#define KMERS_H

namespace kmers {

    class Kmer {
    public:
        Kmer();
        virtual ~Kmer();

        unsigned long int getControlFreq() const {
            return controlFreq;
        }

        void setControlFreq(unsigned long int controlFreq) {
            this->controlFreq = controlFreq;
        }

        unsigned long int getNegativeControl() const {
            return negativeControl;
        }

        void setNegativeControl(unsigned long int negativeControl) {
            this->negativeControl = negativeControl;
        }

        unsigned long int getNegativePeak() const {
            return negativePeak;
        }

        void setNegativePeak(unsigned long int negativePeak) {
            this->negativePeak = negativePeak;
        }

        unsigned long int getPeakFreq() const {
            return peakFreq;
        }

        void setPeakFreq(unsigned long int peakFreq) {
            this->peakFreq = peakFreq;
        }

        double getValue() const {
            return pValue;
        }

        double getPf() const {
            return pf;
        }

        double getSig() const {
            return sig;
        }

        void setValue(double Value) {
            pValue = Value;
        }

        void setPf(double pf) {
            this->pf = pf;
        }

        void setSig(double sig) {
            this->sig = sig;
        }

        void calculatePValue(double totalNRnt_peak, double totalNRnt_control);

    private:
        unsigned long int peakFreq;
        unsigned long int controlFreq;
        unsigned long int negativePeak;
        unsigned long int negativeControl;
        double pValue;
        double sig;
        double pf;
    };

    class KmersFactory {
    public:
        KmersFactory();
        virtual ~KmersFactory();

        std::unordered_set<std::string>& getKmersGenome() {
            return kmersGenome;
        }

        void setKmersGenome(std::unordered_set<std::string> kmersGenome) {
            this->kmersGenome = kmersGenome;
        }

        std::unordered_map<std::string, std::vector<std::shared_ptr<Kmer>>>& getKmers() {
            return kmers;
        }

        std::unordered_map<std::string, unsigned long int>& getKmer2controlFreq() {
            return kmer2controlFreq;
        }

        std::unordered_map<std::string, unsigned long int>& getKmer2peakFreq() {
            return kmer2peakFreq;
        }

        unsigned long int getTotalNRNTControl() const {
            return totalNRnt_control;
        }

        unsigned long int getTotalNRNTPeak() const {
            return totalNRnt_peak;
        }
        
        int getKmerNumber() const {
            return kmerNumber;
        }

        void setKmerNumber(int kmerNumber) {
            this->kmerNumber = kmerNumber;
        }

        void createGenomeWideKmers();

        void clearKmerControlData();

        void clearKmerPeakData();

        void scanSequences(std::string inputSeq, bool control);

        void buildKmers();

        void writeKmersToFile(std::string fileName);

        void readKmersFromFile(std::string fileName);

        double dropSigKmerSigForRefandAlt(char *seq, unsigned long int refPos, char ref, char alt);

        double getKmerSig(std::string kmer, int index);
        
        void mergeKmers(KmersFactory &kmersFactory);

    private:
        int kmerNumber;
        std::unordered_set<std::string> kmersGenome;
        std::unordered_map<std::string, std::vector<std::shared_ptr<Kmer>>> kmers;
        std::unordered_map<std::string, unsigned long int> kmer2controlFreq;
        std::unordered_map<std::string, unsigned long int> kmer2peakFreq;
        unsigned long int totalNRnt_control;
        unsigned long int totalNRnt_peak;

    };
}

#endif /* KMERS_H */


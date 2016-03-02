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

        unsigned long int GetControlFreq() const {
            return controlFreq;
        }

        void SetControlFreq(unsigned long int controlFreq) {
            this->controlFreq = controlFreq;
        }

        unsigned long int GetNegativeControl() const {
            return negativeControl;
        }

        void SetNegativeControl(unsigned long int negativeControl) {
            this->negativeControl = negativeControl;
        }

        unsigned long int GetNegativePeak() const {
            return negativePeak;
        }

        void SetNegativePeak(unsigned long int negativePeak) {
            this->negativePeak = negativePeak;
        }

        unsigned long int GetPeakFreq() const {
            return peakFreq;
        }

        void SetPeakFreq(unsigned long int peakFreq) {
            this->peakFreq = peakFreq;
        }

        double GetValue() const {
            return pValue;
        }

        double GetPf() const {
            return pf;
        }

        double GetSig() const {
            return sig;
        }

        void SetValue(double Value) {
            pValue = Value;
        }

        void SetPf(double pf) {
            this->pf = pf;
        }

        void SetSig(double sig) {
            this->sig = sig;
        }

        void CalculatePValue(double totalNRnt_peak, double totalNRnt_control);

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

        std::set<std::string>& GetKmersGenome() {
            return kmersGenome;
        }

        void SetKmersGenome(std::set<std::string> kmersGenome) {
            this->kmersGenome = kmersGenome;
        }

        std::map<std::string, Kmer *>& GetKmers() {
            return kmers;
        }

        std::map<std::string, unsigned long int>& GetKmer2controlFreq() {
            return kmer2controlFreq;
        }

        std::map<std::string, unsigned long int>& GetKmer2peakFreq() {
            return kmer2peakFreq;
        }

        unsigned long int GetTotalNRnt_control() const {
            return totalNRnt_control;
        }

        unsigned long int GetTotalNRnt_peak() const {
            return totalNRnt_peak;
        }

        void CreateGenomeWideKmers();

        void ClearKmerControlData();

        void ClearKmerPeakData();

        void scanSequences(std::string inputSeq, bool control);

        void BuildKmers();

        void WriteKmersToFile(char *fileName, bool binary);

        void ReadKmersFromFile(char *fileName, bool binary);

        double DropSigKmerSigForRefandAlt(char *seq, unsigned long int refPos, char ref, char alt);

        double GetKmerSig(std::string kmer);

    private:
        std::set<std::string> kmersGenome;
        std::map<std::string, Kmer *> kmers;
        std::map<std::string, unsigned long int> kmer2controlFreq;
        std::map<std::string, unsigned long int> kmer2peakFreq;
        unsigned long int totalNRnt_control;
        unsigned long int totalNRnt_peak;

    };
}

#endif /* KMERS_H */


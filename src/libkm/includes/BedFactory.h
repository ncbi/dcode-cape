/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BedFactory.h
 * Author: veraalva
 *
 * Created on February 11, 2016, 3:52 PM
 */

#ifndef BEDFACTORY_H
#define BEDFACTORY_H

#include "FastaFactory.h"


namespace peak {

    class Peak {
    public:
        Peak();
        virtual ~Peak();

        void CalculateContent();

        std::string GetChr() const {
            return chr;
        }

        void SetChr(std::string chr) {
            this->chr = chr;
        }

        unsigned long int GetGCCount() const {
            return GCCount;
        }

        void SetGCCount(unsigned long int GCCount) {
            this->GCCount = GCCount;
        }

        unsigned long int GetNCount() const {
            return NCount;
        }

        void SetNCount(unsigned long int NCount) {
            this->NCount = NCount;
        }

        unsigned long int GetNRCount() const {
            return NRCount;
        }

        double GetNPercent() const {
            return NPercent;
        }

        unsigned long int GetEnd() const {
            return end;
        }

        void SetEnd(unsigned long int end) {
            this->end = end;
        }

        char *GetSeq() {
            return seq;
        }

        void SetSeq(char **seq) {
            this->seq = *seq;
        }

        unsigned long int GetStart() const {
            return start;
        }

        void SetStart(unsigned long int start) {
            this->start = start;
        }

        unsigned long int GetLength() const {
            return end - start + 1;
        }

        std::pair<int, int> GetGCNcontentBin();

        char *RandomizePeakSeq();

    private:
        std::string chr;
        unsigned long int start;
        unsigned long int end;
        unsigned long int GCCount;
        unsigned long int NCount;
        unsigned long int NRCount;
        double NPercent;
        char *seq;
    };

    class BedFactory {
    public:
        BedFactory();
        virtual ~BedFactory();

        std::string GetGenomeName() const {
            return genomeName;
        }

        std::vector<Peak*>& GetPeaks() {
            return peaks;
        }

        std::map<std::string, std::map<int, std::map<std::pair<int, int>, int>>>& GetGCNcontentBin() {
            return GCNcontentBin;
        }

        void CreatePeaksFromBedFile(fasta::FastaFactory& chrFactory, std::string genomeName, std::string bedFileName, double maxNPercent, kmers::KmersFactory &kmersFactory);

        void GeneratingControlsFromChromosomes(fasta::FastaFactory &chrFactory, unsigned long int hit_Num, kmers::KmersFactory &kmersFactory);

        void GeneratingControlsFromShufflingPeaks(unsigned long int hit_Num, kmers::KmersFactory &kmersFactory);

        void ReadingControlsFromFile(char *controlFileName, fasta::FastaFactory &chrFactory, kmers::KmersFactory &kmersFactory);

    private:
        unsigned long int NCount;
        unsigned long int GCCount;
        std::string genomeName;
        std::vector<Peak *> peaks;
        std::map<std::string, std::map<int, std::map<std::pair<int, int>, int>>> GCNcontentBin;

        void plus_ocur(char nt);
        void minus_ocur(char nt);

    };

}

#endif /* BEDFACTORY_H */


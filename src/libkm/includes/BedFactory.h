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

        void calculateContent();

        std::string getChr() const {
            return chr;
        }

        void setChr(std::string chr) {
            this->chr = chr;
        }

        unsigned long int getGCCount() const {
            return GCCount;
        }

        void setGCCount(unsigned long int GCCount) {
            this->GCCount = GCCount;
        }

        unsigned long int getNCount() const {
            return NCount;
        }

        void setNCount(unsigned long int NCount) {
            this->NCount = NCount;
        }

        unsigned long int getNRCount() const {
            return NRCount;
        }

        double getNPercent() const {
            return NPercent;
        }

        unsigned long int getEnd() const {
            return end;
        }

        void setEnd(unsigned long int end) {
            this->end = end;
        }

        char *getSeq() {
            return seq;
        }

        void setSeq(char **seq) {
            this->seq = *seq;
        }

        unsigned long int getStart() const {
            return start;
        }

        void setStart(unsigned long int start) {
            this->start = start;
        }

        unsigned long int getLength() const {
            return end - start + 1;
        }

        std::pair<int, int> getGCNcontentBin();

        std::string randomizePeakSeq();

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

        std::vector<Peak*>& getPeaks() {
            return peaks;
        }

        std::map<std::string, std::map<int, std::map<std::pair<int, int>, int>>>& getGCNcontentBin() {
            return GCNcontentBin;
        }

        void createPeaksFromBedFile(sequence::FastaFactory& chrFactory, std::string bedFileName, double maxNPercent, kmers::KmersFactory &kmersFactory);

        void generatingControlsFromChromosomes(sequence::FastaFactory &chrFactory, unsigned long int hit_Num, kmers::KmersFactory &kmersFactory);

        void generatingControlsFromShufflingPeaks(unsigned long int hit_Num, kmers::KmersFactory &kmersFactory);

        void readControlsFromFile(std::string controlFileName, sequence::FastaFactory &chrFactory, kmers::KmersFactory &kmersFactory);

    private:
        unsigned long int NCount;
        unsigned long int GCCount;
        std::vector<Peak *> peaks;
        std::map<std::string, std::map<int, std::map<std::pair<int, int>, int>>> GCNcontentBin;

        void plus_ocur(char nt);
        void minus_ocur(char nt);

    };

}

#endif /* BEDFACTORY_H */


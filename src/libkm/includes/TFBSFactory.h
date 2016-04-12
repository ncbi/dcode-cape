/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TFBSFactory.h
 * Author: veraalva
 *
 * Created on March 24, 2016, 12:46 PM
 */

#ifndef TFBSFACTORY_H
#define TFBSFACTORY_H

namespace tfbs {

    class Tib {
    public:
        Tib();
        virtual ~Tib();

        long int getLen() const {
            return len;
        }

        std::string getName() const {
            return name;
        }

        void setLen(long int len) {
            this->len = len;
        }

        void setName(std::string name) {
            this->name = name;
        }

    private:
        std::string name;
        long int len;
    };

    class TFBS {
    public:
        TFBS(unsigned long int d, short int i);
        virtual ~TFBS();

        unsigned long int getDelta() const {
            return delta;
        }

        short int getIndex() const {
            return index;
        }

        char getStrand() const {
            return strand;
        }

        long int getEnd() const {
            return end;
        }

        void setEnd(long int end) {
            this->end = end;
        }

        long int getStart() const {
            return start;
        }

        void setStart(long int start) {
            this->start = start;
        }

    private:
        unsigned long int delta;
        short int index;
        char strand;
        long int start;
        long int end;
    };

    class TFBSFactory {
    public:
        TFBSFactory();
        virtual ~TFBSFactory();

        void createTFBSFileIndexMap(std::string dirName, std::string prefix, std::string idxExtension, std::string tibExtension);
        void createPWMIndexFromTibInfoFile(std::string tibInfoFileName);
        void extractTFBSFromFile(long int from, long int to, sequence::Seq *chr);

        std::vector<Tib *>& getPwmIndex() {
            return pwmIndex;
        }

        std::vector<TFBS *>& getTfbs() {
            return tfbs;
        }

        long int getLongestPWM() const {
            return longestPWM;
        }

        void setLongestPWM(long int longestPWM) {
            this->longestPWM = longestPWM;
        }

        bool isReady() {
            return !this->tfbsFileIndex.empty();
        }

    private:
        std::vector<Tib *> pwmIndex;
        long int longestPWM;
        std::vector<TFBS *> tfbs;
        std::unordered_map<std::string, std::pair<FILE *, FILE *>> tfbsFileIndex;
        std::string currentChr;
        FILE *chrIdxFile;
        FILE *chrTibFile;
    };
}

#endif /* TFBSFACTORY_H */


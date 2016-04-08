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

        long int GetLen() const {
            return len;
        }

        std::string GetName() const {
            return name;
        }

        void SetLen(long int len) {
            this->len = len;
        }

        void SetName(std::string name) {
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

        unsigned long int GetDelta() const {
            return delta;
        }

        short int GetIndex() const {
            return index;
        }

        char GetStrand() const {
            return strand;
        }

        long int GetEnd() const {
            return end;
        }

        void SetEnd(long int end) {
            this->end = end;
        }

        long int GetStart() const {
            return start;
        }

        void SetStart(long int start) {
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

        void CreateTFBSFileIndexMap(char *dirName, const char *prefix, const char *idxExtension, const char *tibExtension);
        void CreatePWMIndexFromTibInfoFile(const char *tibInfoFileName);
        void ExtractTFBSFromFile(long int from, long int to, fasta::Fasta *chr);

        std::vector<Tib *>& GetPwmIndex() {
            return pwmIndex;
        }

        std::vector<TFBS *>& GetTfbs() {
            return tfbs;
        }

        long int GetLongestPWM() const {
            return longestPWM;
        }

        void SetLongestPWM(long int longestPWM) {
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


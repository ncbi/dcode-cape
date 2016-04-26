/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SNPFactory.h
 * Author: veraalva
 *
 * Created on February 25, 2016, 9:22 AM
 */

#ifndef SNPFACTORY_H
#define SNPFACTORY_H

#include "KmersFactory.h"
#include "SVMPredict.h"


namespace snp {

    class SNP {
    public:
        SNP();
        virtual ~SNP();

        char getAlt() const {
            return alt;
        }

        void setAlt(char alt) {
            this->alt = alt;
        }

        std::string getId() const {
            return id;
        }

        void setId(std::string id) {
            this->id = id;
        }

        std::string getChr() const {
            return chr;
        }

        void setChr(std::string chr) {
            this->chr = chr;
        }

        unsigned long int getLength() const {
            return length;
        }

        void setLength(unsigned long int length) {
            this->length = length;
        }

        unsigned long int getPos() const {
            return pos;
        }

        void setPos(unsigned long int pos) {
            this->pos = pos;
        }

        unsigned long int getChrPos() const {
            return chrPos;
        }

        void setChrPos(unsigned long int chrPos) {
            this->chrPos = chrPos;
        }

        char getRef() const {
            return ref;
        }

        void setRef(char ref) {
            this->ref = ref;
        }

        std::string getSeq() {
            return seq;
        }

        void setSeq(std::string seq) {
            this->seq = seq;
        }

        std::vector<double>& getDescriptors() {
            return descriptors;
        }

        double getProbPos() const {
            return probPos;
        }

        void setProbPos(double ProbPos) {
            probPos = ProbPos;
        }

        void calculateKmerDescriptors(kmers::KmersFactory& kmersFactory, unsigned long int featNumber);

    private:
        std::string id;
        std::string chr;
        unsigned long int length;
        std::string seq;
        unsigned long int pos;
        unsigned long int chrPos;
        char ref;
        char alt;
        std::vector<double> descriptors;
        double probPos;


    };

    class SNPFactory {
    public:
        SNPFactory();
        SNPFactory(const SNPFactory& orig);
        virtual ~SNPFactory();

        std::vector<std::shared_ptr<SNP>>& getSnps() {
            return snps;
        }
        
        std::string getExpressionCode() const {
            return expressionCode;
        }

        void setExpressionCode(std::string expressionCode) {
            this->expressionCode = expressionCode;
        }

        void parseSNPFile(std::string snpFileName, unsigned long int neighbors, sequence::FastaFactory &chrFactory);
        void writeEnhansersFastaFile(std::string fastaFile, bool binary);
        int processSNPFromFile(std::string snpFileName, unsigned long int neighbors, sequence::FastaFactory &chrFactory, kmers::KmersFactory& kmersFactory, svm::SVMPredict& svmPredict, fimo::FimoFactory & fimoFactory, tfbs::TFBSFactory & tFBSFactory);
    private:
        std::string expressionCode;
        std::vector<std::shared_ptr<SNP>> snps;

    };
}

#endif /* SNPFACTORY_H */


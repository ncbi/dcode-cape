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

        char GetAlt() const {
            return alt;
        }

        void SetAlt(char alt) {
            this->alt = alt;
        }

        std::string GetId() const {
            return id;
        }

        void SetId(std::string id) {
            this->id = id;
        }

        std::string GetChr() const {
            return chr;
        }

        void SetChr(std::string chr) {
            this->chr = chr;
        }

        unsigned long int GetLength() const {
            return length;
        }

        void SetLength(unsigned long int length) {
            this->length = length;
        }

        unsigned long int GetPos() const {
            return pos;
        }

        void SetPos(unsigned long int pos) {
            this->pos = pos;
        }

        unsigned long int GetChrPos() const {
            return chrPos;
        }

        void SetChrPos(unsigned long int chrPos) {
            this->chrPos = chrPos;
        }

        char GetRef() const {
            return ref;
        }

        void SetRef(char ref) {
            this->ref = ref;
        }

        char* GetSeq() {
            return seq;
        }

        void SetSeq(char** seq) {
            this->seq = *seq;
            this->length = strlen(this->seq);
        }

        std::vector<double>& GetDescriptors() {
            return descriptors;
        }

        double GetProbPos() const {
            return probPos;
        }

        void SetProbPos(double ProbPos) {
            probPos = ProbPos;
        }

        void CalculateKmerDescriptors(kmers::KmersFactory& kmersFactory, unsigned long int featNumber);

    private:
        std::string id;
        std::string chr;
        unsigned long int length;
        char *seq;
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

        std::vector<SNP*>& GetSnps() {
            return snps;
        }

        void ReadSNPFromFile(char* snpFileName, unsigned long int neighbors, fasta::FastaFactory &chrFactory);
        void WriteEnhansersFastaFile(char* fastaFile, bool binary);
        int ProcessSNPFromFile(char* snpFileName, unsigned long int neighbors, fasta::FastaFactory &chrFactory, kmers::KmersFactory& kmersFactory, svm::SVMPredict& svmPredict, fimo::FimoFactory & fimoFactory);
    private:
        std::vector<SNP *> snps;

    };
}

#endif /* SNPFACTORY_H */


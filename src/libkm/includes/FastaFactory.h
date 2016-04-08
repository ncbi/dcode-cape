/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FastaFactory.h
 * Author: veraalva
 *
 * Created on February 10, 2016, 3:41 PM
 */

#ifndef FASTAFACTORY_H
#define FASTAFACTORY_H

namespace fasta {

    class Fasta {
    public:
        Fasta();

        Fasta(const std::string& i, const std::string& d) : id(i), description(d) {
            this->length = 0;
            this->seq = NULL;
        }
        Fasta(const Fasta& orig);
        virtual ~Fasta();

        char *GetSubStr(int pos, int length) {
            return strndup(seq + pos, length);
        }

        std::string GetId() const {
            return id;
        }

        void SetId(std::string id) {
            this->id = id;
        }

        char *GetSeq() {
            return seq;
        }

        void SetSeq(char **seq) {
            this->seq = *seq;
        }

        std::string GetDescription() const {
            return description;
        }

        void SetDescription(std::string desc) {
            description = desc;
        }

        unsigned long int GetLength() {
            return length;
        }

        void SetLength(unsigned long int len) {
            length = len;
        }

    private:
        std::string id;
        std::string description;
        char *seq;
        unsigned long int length;
    };

    class FastaFactory {
    public:
        FastaFactory();
        FastaFactory(const FastaFactory& orig);
        virtual ~FastaFactory();

        void LoadFastaInDirectory(char *dirName, const char *prefix, const char *sufix, bool binary);
        long unsigned int ParseFastaFile(FILE *fName, int numberSeqTotalRead, bool cleanContainers, bool binary);
        Fasta *GetFastaFromID(std::string id);

        void WriteSequencesToFile(char *fileName, bool binary);

        std::unordered_map<std::string, Fasta*>& GetFastaMap() {
            return fastaMap;
        }
    private:
        std::unordered_map<std::string, Fasta *> fastaMap;
    };

}

#endif /* FASTAFACTORY_H */


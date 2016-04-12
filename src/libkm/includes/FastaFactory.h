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

namespace sequence {

    class Seq {
    public:
        Seq();

        Seq(const std::string& i, const std::string& d) : id(i), description(d) {
            this->length = 0;
            this->seq = NULL;
        }
        Seq(const Seq& orig);
        virtual ~Seq();

        char *getSubStr(int pos, int length) {
            return strndup(seq + pos, length);
        }

        std::string getId() const {
            return id;
        }

        void setId(std::string id) {
            this->id = id;
        }

        char *getSeq() {
            return seq;
        }

        void setSeq(char **seq) {
            this->seq = *seq;
        }

        std::string getDescription() const {
            return description;
        }

        void setDescription(std::string desc) {
            description = desc;
        }

        unsigned long int getLength() {
            return length;
        }

        void setLength(unsigned long int len) {
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

        std::unordered_map<std::string, Seq*>& getSequenceContainter() {
            return sequenceContainer;
        }

        Seq *getSequenceFromID(std::string id) {
            std::unordered_map<std::string, Seq *>::iterator it = sequenceContainer.find(id);
            if (it == sequenceContainer.end()) {
                throw exceptions::NotFoundException("Id " + id + " was not found in the sequence container");
            }
            return it->second;
        }

        void parseFastaInDirectory(std::string dirName, std::string prefix, std::string sufix, bool binary);
        long unsigned int parseFastaFile(FILE *fName, int numberSeqTotalRead, bool cleanContainers, bool binary);
        void writeSequencesToFile(std::string fileName, bool binary);
    private:
        std::unordered_map<std::string, Seq *> sequenceContainer;
    };

}

#endif /* FASTAFACTORY_H */


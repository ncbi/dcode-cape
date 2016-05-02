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

#include "Exceptions.h"


namespace sequence {

    class Seq {
    public:
        Seq();

        virtual ~Seq();

        std::string getSubStr(unsigned long int pos, unsigned long int length) {
            return seq.substr(pos, length);
        }

        std::string &getId() {
            return id;
        }

        void setId(std::string id) {
            this->id = id;
        }

        std::string &getSeq() {
            return seq;
        }

        void setSeq(std::string seq) {
            this->seq = seq;
        }

        std::string &getDescription() {
            return description;
        }

        void setDescription(std::string desc) {
            description = desc;
        }

        unsigned long int getLength() {
            return seq.size();
        }

    private:
        std::string id;
        std::string description;
        std::string seq;
    };

    class FastaFactory {
    public:
        FastaFactory();
        FastaFactory(const FastaFactory& orig);
        virtual ~FastaFactory();

        std::unordered_map<std::string, std::shared_ptr<Seq>>&getSequenceContainter() {
            return sequenceContainer;
        }

        std::shared_ptr<Seq> getFirstSequence() {
            std::unordered_map<std::string, std::shared_ptr < Seq>>::iterator it = sequenceContainer.begin();
            if (it == sequenceContainer.end()) {
                throw exceptions::NotFoundException("Not sequences on the container");
            }
            return it->second;
        }

        std::shared_ptr<Seq> getSequenceFromID(std::string id) {
            std::unordered_map<std::string, std::shared_ptr < Seq>>::iterator it = sequenceContainer.find(id);
            if (it == sequenceContainer.end()) {
                throw exceptions::NotFoundException("Id " + id + " was not found in the sequence container");
            }
            return it->second;
        }

        void parseFastaInDirectory(std::string dirName, std::string prefix, std::string sufix, bool binary);
        long unsigned int parseFastaFile(std::string fName, bool binary);
        void writeSequencesToFile(std::string fileName, bool binary);
    private:
        std::unordered_map<std::string, std::shared_ptr<Seq>> sequenceContainer;
    };

}

#endif /* FASTAFACTORY_H */


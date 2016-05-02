/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FileParserFactory.h
 * Author: veraalva
 *
 * Created on April 11, 2016, 12:50 PM
 */

#ifndef FILEPARSERFACTORY_H
#define FILEPARSERFACTORY_H

#include "Exceptions.h"


namespace parsers {

    class FileParserFactory {
    public:
        FileParserFactory();
        virtual ~FileParserFactory();

        void setFileToParse(std::string fileToParseName) {
            clean();
            fileToParse.open(fileToParseName, std::ios::in | std::ios::binary);
            if (!fileToParse) {
                std::cerr << "Error opening file: " << fileToParseName << std::endl;
                exit(-1);
            }
            closeFile = true;
            fileToParse.seekg(0, fileToParse.end);
            if (static_cast<unsigned long int> (fileToParse.tellg()) < bufferSize) {
                bufferSize = static_cast<unsigned long int> (fileToParse.tellg());
            }
            buffer.resize(bufferSize + 1);
            fileToParse.seekg(0, fileToParse.beg);
        }

        std::vector<std::string>& getWords() {
            return words;
        }

        std::string& getLine() {
            return line;
        }

        bool lineStartWith(std::string s) {
            if (line.compare(0, s.size(), s) == 0) return true;
            return false;
        }

        bool iterate(std::string dontStartWith);
        bool iterate(std::string dontStartWith, std::string delimiter);

    private:
        bool closeFile;
        bool backup;
        std::ifstream fileToParse;
        std::string buffer;
        std::string line;
        std::vector<std::string> words;
        unsigned long int bufferSize;
        unsigned long int currPosition;

        void clean();
    };

}

#endif /* FILEPARSERFACTORY_H */


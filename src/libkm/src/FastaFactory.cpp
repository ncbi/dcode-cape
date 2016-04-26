/* 
 * File:   FastaFactory.cpp
 * Author: veraalva
 * 
 * Created on February 10, 2016, 3:41 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <dirent.h>
#include <inttypes.h>

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <utility>
#include <fstream>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"

#include "Global.h"
#include "TimeUtils.h"
#include "FileParserFactory.h"
#include "FastaFactory.h"

using namespace std;
using namespace parsers;
using namespace sequence;

Seq::Seq() {
    seq = NULL;
    length = 0;
}

Seq::~Seq() {
    if (seq) free(seq);
}

FastaFactory::FastaFactory() {
}

FastaFactory::FastaFactory(const FastaFactory& orig) {
}

FastaFactory::~FastaFactory() {
}

long unsigned int FastaFactory::parseFastaFile(std::string fName, bool binary) {

    int numberSeqCurrentRead = 0;
    uint64_t i, len;
    char *buffer, *seq, *seq_end;
    shared_ptr<Seq> fasta = nullptr;

    sequenceContainer.clear();
    if (binary) {
        ifstream  inFile(fName, std::ifstream::binary);
        inFile.read((char *) &i, sizeof (uint64_t));
        for (unsigned long int j = 0; j < i; j++) {
            fasta = make_shared<Seq>();
            inFile.read((char *) &len, sizeof (uint64_t));
            buffer = (char *) allocate(sizeof (char) * len, __FILE__, __LINE__);
            inFile.read(buffer, sizeof (char) * len);
            fasta->setId(buffer);
            inFile.read((char *) &len, sizeof (uint64_t));
            buffer = (char *) reallocate(buffer, sizeof (char) * len, __FILE__, __LINE__);
            inFile.read(buffer, sizeof (char) * len);
            fasta->setSeq(&buffer);
            fasta->setLength(len - 1);
            if (Global::instance()->isDebug3()) cerr << "\tDEBUG3 ==> Adding seq: " << fasta->getId() << " with " << fasta->getLength() << " bp" << endl;
            sequenceContainer.insert(pair<string, shared_ptr < Seq >> (fasta->getId(), fasta));
            numberSeqCurrentRead++;
        }
        inFile.close();
    } else {
        FileParserFactory fParser;
        seq = seq_end = NULL;

        size_t str_length;
        size_t readTotal = 0;
        size_t numLines = 0;
        size_t seq_length = 0;

        try {
            fParser.setFileToParse(fName);
            while (fParser.iterate('#')) {
                char *line = fParser.getLine();
                str_length = fParser.getLineLength();
                if (*line == '>') {
                    if (seq_end != seq) {
                        seq = (char *) reallocate(seq, sizeof (char) * (seq_length + 1), __FILE__, __LINE__);
                        *(seq + seq_length) = 0;
                        fasta->setSeq(&seq);
                        fasta->setLength(strlen(seq));
                        if (Global::instance()->isDebug3()) cout << "\tDEBUG3 ==> Adding seq: " << fasta->getId().c_str() << " with " << fasta->getLength() << " bp" << endl;
                        sequenceContainer.insert(pair<string, shared_ptr < Seq >> (fasta->getId(), fasta));
                        numberSeqCurrentRead++;
                    }
                    seq_length = 0;
                    seq = NULL;
                    readTotal += str_length;
                    fasta = make_shared<Seq>();
                    fasta->setId(string(line + 1, str_length - 1));
                } else {
                    if (fasta == nullptr){
                        cerr << "Fasta file does not start with the header (>)" << endl;
                        exit(-1);
                    }
                    readTotal += str_length;
                    seq_length += str_length;

                    seq = (char *) reallocate(seq, sizeof (char) * seq_length, __FILE__, __LINE__);
                    seq_end = seq + (seq_length - str_length);
                    memcpy(seq_end, line, str_length);
                    seq_end += str_length;
                }
                numLines++;
            }

            if (seq) {
                numberSeqCurrentRead++;
                seq = (char *) reallocate(seq, sizeof (char) * (seq_length + 1), __FILE__, __LINE__);
                *(seq + seq_length) = 0;
                fasta->setSeq(&seq);
                fasta->setLength(seq_length);
                if (Global::instance()->isDebug3()) cout << "\tDEBUG3 ==> Adding seq: " << fasta->getId().c_str() << " with " << fasta->getLength() << " bp" << endl;
                sequenceContainer.insert(pair<string, shared_ptr < Seq >> (fasta->getId(), fasta));
            }
        } catch (exceptions::FileNotFoundException ex) {
            cerr << ex.what() << endl;
            cerr << "Error parsing file" << endl;
            exit(-1);
        } catch (exceptions::ErrorReadingFromFileException ex) {
            cerr << ex.what() << endl;
            cerr << "Error parsing file" << endl;
            exit(-1);
        }
    }

    return numberSeqCurrentRead;
}

void FastaFactory::parseFastaInDirectory(std::string dirName, std::string prefix, std::string sufix, bool binary) {
    struct dirent *dp;
    DIR *dirp = (DIR *) checkPointerError(opendir(dirName.c_str()), "Can't open input directory", __FILE__, __LINE__, -1);

    while ((dp = readdir(dirp)) != NULL) {
        bool read = false;
        string fName(dp->d_name);
        if (Global::instance()->isDebug3()) {
            cout << "\tDEBUG3 ==> Found file: " << fName << endl;
        }
        if (fName[0] != '.') {
            if (prefix.empty() && sufix.empty()) {
                read = true;
            } else {
                if (!prefix.empty() && sufix.empty()) {
                    if (prefix.size() <= fName.size() &&
                            fName.compare(0, prefix.size(), prefix) == 0) read = true;
                } else if (prefix.empty() && !sufix.empty()) {
                    if (sufix.size() <= fName.size() &&
                            fName.compare(fName.size() - sufix.size(), sufix.size(), sufix) == 0) read = true;
                } else if (!prefix.empty() && !sufix.empty()) {
                    if (prefix.size() <= fName.size() &&
                            sufix.size() <= fName.size() &&
                            fName.compare(0, prefix.size(), prefix) == 0 &&
                            fName.compare(fName.size() - sufix.size(), sufix.size(), sufix) == 0) read = true;
                }
            }
        }
        if (read) {
            if (Global::instance()->isInfo()) {
                TimeUtils::instance()->setClock();
                cout << "\tINFO ==> Parsing file: " << dirName + "/" + fName << endl;
            }
            int seqs = parseFastaFile(dirName + "/" + fName, binary);
            if (Global::instance()->isInfo()) cout << "\tINFO ==> " << seqs << " sequences read in " << TimeUtils::instance()->getTimeSec() << " sec" << endl;
        }
    }
    closedir(dirp);
}

void FastaFactory::writeSequencesToFile(std::string fileName, bool binary) {
    uint64_t i, len;
    shared_ptr<Seq> f;

    if (binary) {
        std::ofstream outputFile(fileName, std::ofstream::binary);
        if (!outputFile.is_open()) {
            cerr << "Can't open output file " << fileName << endl;
            exit(-1);
        }
        i = sequenceContainer.size();
        outputFile.write((char *) &i, sizeof (uint64_t));
        for (auto it = sequenceContainer.begin(); it != sequenceContainer.end(); ++it) {
            f = it->second;
            len = f->getId().size() + 1;
            outputFile.write((char *) &len, sizeof (uint64_t));
            outputFile.write(f->getId().c_str(), sizeof(char) * len);
            len = f->getLength() + 1;
            outputFile.write((char *) &len, sizeof (uint64_t));
            outputFile.write(f->getSeq(), sizeof(char) * len);
        }
        outputFile.close();
    } else {
        ofstream outputFile(fileName);
        if (!outputFile.is_open()) {
            cerr << "Can't open output file " << fileName << endl;
            exit(-1);
        }
        for (auto it = sequenceContainer.begin(); it != sequenceContainer.end(); ++it) {
            f = it->second;
            outputFile << ">" << f->getId() << endl;
            for (i = 0; i < f->getLength(); i += 50) {
                char t = 0;
                if (i + 50 < f->getLength()) {
                    t = f->getSeq()[i + 50];
                    f->getSeq()[i + 50] = 0;
                }
                outputFile << (f->getSeq() + i) << endl;
                if (i + 50 < f->getLength()) {
                    f->getSeq()[i + 50] = t;
                }
            }
        }
        outputFile.close();
    }

}

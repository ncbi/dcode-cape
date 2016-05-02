/* 
 * File:   FastaFactory.cpp
 * Author: veraalva
 * 
 * Created on February 10, 2016, 3:41 PM
 */

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

#include "Global.h"
#include "TimeUtils.h"
#include "FileParserFactory.h"
#include "FastaFactory.h"

using namespace std;
using namespace parsers;
using namespace sequence;

Seq::Seq() {
}

Seq::~Seq() {
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
    string id;
    shared_ptr<Seq> fasta;
    pair < unordered_map<string, shared_ptr < Seq>>::iterator, bool> result;

    if (binary) {
        ifstream inFile(fName, std::ifstream::binary);
        inFile.read((char *) &i, sizeof (uint64_t));
        for (unsigned long int j = 0; j < i; j++) {
            inFile.read((char *) &len, sizeof (uint64_t));
            id.resize(len);
            inFile.read(&(id[0]), len);
            result = sequenceContainer.insert(make_pair(id, make_shared<Seq>()));
            if (result.second) {
                fasta = result.first->second;
            } else {
                cerr << "Duplicated sequence ID " << id << endl;
                exit(-1);
            }
            fasta->setId(id);
            inFile.read((char *) &len, sizeof (uint64_t));
            fasta->getSeq().resize(len);
            inFile.read(&(fasta->getSeq()[0]), len);
            numberSeqCurrentRead++;
        }
        inFile.close();
    } else {
        FileParserFactory fParser;
        try {
            fParser.setFileToParse(fName);
            while (fParser.iterate("#")) {
                string line = fParser.getLine();
                if (fParser.lineStartWith(">")) {
                    id = line.substr(1, fParser.getLine().size() - 1);
                    result = sequenceContainer.insert(make_pair(id, make_shared<Seq>()));
                    if (result.second) {
                        fasta = result.first->second;
                    } else {
                        cerr << "Duplicated sequence ID " << id << endl;
                        exit(-1);
                    }
                    fasta->setId(id);
                    numberSeqCurrentRead++;
                } else {
                    if (fasta == nullptr) {
                        cerr << "Fasta file does not start with the header (>)" << endl;
                        exit(-1);
                    }
                    fasta->getSeq().append(line);
                }
            }
        } catch (exceptions::FileNotFoundException) {
            cerr << "Error parsing file: " << fName << endl;
            exit(-1);
        } catch (ios::failure) {
            cerr << "Error parsing file: " << fName << endl;
            exit(-1);
        }
    }

    if (Global::instance()->isDebug3()) {
        cout << "\tDEBUG3 ==> " << sequenceContainer.size() << " sequences in the container."<< endl;
        for (auto it = sequenceContainer.begin(); it != sequenceContainer.end(); ++it) {
            fasta = it->second;
            cout << "\tDEBUG3 ==>\t\t" << fasta->getId() << " with " << fasta->getLength() << " bp" << endl;
        }
    }

    return numberSeqCurrentRead;
}

void FastaFactory::parseFastaInDirectory(std::string dirName, std::string prefix, std::string sufix, bool binary) {
    struct dirent *dp;
    DIR *dirp = (DIR *) opendir(dirName.c_str());
    if (!dirp) {
        cerr << "Can't open directory: " << dirName << endl;
        exit(-1);
    }

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
    
    if (Global::instance()->isDebug3()) {
        cout << "\tDEBUG3 ==> " << sequenceContainer.size() << " sequences in the container to write"<< endl;
    }

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
            if (Global::instance()->isDebug3()) {
                cout << "\tDEBUG3 ==> Writing sequence: " << f->getId() << " with length " << f->getLength() << endl;
            }
            len = f->getId().size();
            outputFile.write((char *) &len, sizeof (uint64_t));
            outputFile.write(f->getId().c_str(), len);
            len = f->getLength();
            outputFile.write((char *) &len, sizeof (uint64_t));
            outputFile.write(f->getSeq().c_str(), len);
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
            if (Global::instance()->isDebug3()) {
                cout << "\tDEBUG3 ==> Writing sequence: " << f->getId() << " with length " << f->getLength() << endl;
            }
            outputFile << ">" << f->getId() << endl;
            for (i = 0; i < f->getLength(); i += 50) {
                char t = 0;
                if (i + 50 < f->getLength()) {
                    t = f->getSeq()[i + 50];
                    f->getSeq()[i + 50] = 0;
                }
                outputFile << (f->getSeq().c_str() + i) << endl;
                if (i + 50 < f->getLength()) {
                    f->getSeq()[i + 50] = t;
                }
            }
        }
        outputFile.close();
    }

}

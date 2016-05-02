/* 
 * File:   TFBSFactory.cpp
 * Author: veraalva
 * 
 * Created on March 24, 2016, 12:46 PM
 */

#include <dirent.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>

#include "Global.h"
#include "FileParserFactory.h"
#include "Exceptions.h"
#include "FastaFactory.h"
#include "TFBSFactory.h"

using namespace std;
using namespace parsers;
using namespace exceptions;
using namespace sequence;
using namespace tfbs;

Tib::Tib() {
    this->len = 0;
}

Tib::~Tib() {

}

TFBS::TFBS(unsigned long int d, short int i) {
    this->start = 0;
    this->end = 0;
    delta = d;
    if (i < 0) {
        index = -1 * i;
        strand = '-';
    } else {
        index = i;
        strand = '+';
    }
}

TFBS::~TFBS() {

}

TFBSFactory::TFBSFactory() {
    longestPWM = 0;
}

TFBSFactory::~TFBSFactory() {
    if (chrIdxFile.is_open()) chrIdxFile.close();
    if (chrTibFile.is_open()) chrTibFile.close();
}

void TFBSFactory::createTFBSFileIndexMap(string dirName, string prefix, string idxExtension, string tibExtension) {
    struct dirent *dp;
    DIR *dirp = (DIR *) opendir(dirName.c_str());
    if (!dirp) {
        cerr << "Can't open directory: " << dirName << endl;
        exit(-1);
    }

    while ((dp = readdir(dirp)) != NULL) {
        bool read = false;
        string fName(dp->d_name);
        if (fName[0] != '.') {
            if (prefix.empty() && idxExtension.empty()) {
                read = true;
            } else {
                if (!prefix.empty() && idxExtension.empty()) {
                    if (prefix.size() <= fName.size() &&
                            fName.compare(0, prefix.size(), prefix) == 0) read = true;
                } else if (prefix.empty() && !idxExtension.empty()) {
                    if (idxExtension.size() <= fName.size() &&
                            fName.compare(fName.size() - idxExtension.size(), idxExtension.size(), idxExtension) == 0) read = true;
                } else if (!prefix.empty() && !idxExtension.empty()) {
                    if (prefix.size() <= fName.size() &&
                            idxExtension.size() <= fName.size() &&
                            fName.compare(0, prefix.size(), prefix) == 0 &&
                            fName.compare(fName.size() - idxExtension.size(), idxExtension.size(), idxExtension) == 0) read = true;
                }
            }
        }
        if (read) {
            fName.replace(fName.size() - idxExtension.size(), idxExtension.size(), "");
            chrIdxFile.open(dirName + "/" + fName + idxExtension, ifstream::binary);
            chrTibFile.open(dirName + "/" + fName + tibExtension, ifstream::binary);
            if (!chrIdxFile.is_open() || !chrTibFile.is_open()) {
                cerr << "ERROR!!" << endl;
                cerr << "Can't open TFBS index files" << endl;
                exit(-1);
            }
            
            chrIdxFile.close();
            chrTibFile.close();
            tfbsFileIndex.insert(make_pair(fName, make_pair(dirName + "/" + fName + idxExtension, dirName + "/" + fName + tibExtension)));
        }
    }
    closedir(dirp);
}

void TFBSFactory::createPWMIndexFromTibInfoFile(string tibInfoFileName) {
    FileParserFactory fParser;

    try {
        fParser.setFileToParse(tibInfoFileName);
        while (fParser.iterate("#", " ")) {
            if (fParser.getWords().size() != 3) {
                cerr << "Tib info file with a wrong format " << endl;
                exit(-1);
            }
            shared_ptr<Tib> tib = make_shared<Tib>();
            tib->setName(fParser.getWords()[2]);
            tib->setLen(atol((fParser.getWords()[1]).c_str()));
            if (longestPWM < tib->getLen()) {
                longestPWM = tib->getLen();
            }
            pwmIndex.push_back(tib);
        }
    } catch (exceptions::FileNotFoundException ex) {
        cerr << "Error parsing file: " << tibInfoFileName << endl;
        exit(-1);
    } catch (ios::failure ex) {
        cerr << "Error parsing file: " << tibInfoFileName << endl;
        exit(-1);
    }
}

class TFBSList {
public:
    TFBSList();
    virtual ~TFBSList();

    vector<shared_ptr<TFBS>>&GetElements() {
        return elements;
    }

private:
    vector<shared_ptr<TFBS>> elements;
};

TFBSList::TFBSList() {

}

TFBSList::~TFBSList() {

}

void TFBSFactory::extractTFBSFromFile(long int from, long int to, shared_ptr<Seq> chr) {
    long int i, j, k;
    uint32_t b;
    uint16_t intFromByte;
    bool seek_performed = false;
    long int start, end, rFrom, rTo;
    shared_ptr<TFBS> tfbsElement;
    shared_ptr<TFBSList> tfbsPtr;
    vector<shared_ptr < TFBSList>> tfbsList;
    unordered_map<string, pair<string, string>>::iterator tfbsFileIndexIt;

    if (chr->getId().compare(currentChr) != 0) {
        tfbsFileIndexIt = tfbsFileIndex.find(chr->getId());
        if (tfbsFileIndexIt == tfbsFileIndex.end()) {
            throw NotFoundException("Indexes files are not available");
        }

        currentChr = chr->getId();
        if (chrIdxFile.is_open()) chrIdxFile.close();
        if (chrTibFile.is_open()) chrTibFile.close();

        chrIdxFile.open(tfbsFileIndexIt->second.first);
        chrTibFile.open(tfbsFileIndexIt->second.second);
    }

    if (!chrIdxFile.is_open() || !chrTibFile.is_open()) {
        throw NotFoundException("Current indexes files are not available");
    }

    rFrom = from - longestPWM - 1;
    rTo = to + longestPWM;
    if (rFrom < 1) {
        rFrom = 1;
    }
    if (rTo > static_cast<long int> (chr->getLength())) {
        rTo = chr->getLength();
    }

    if (rFrom > rTo) {
        throw std::out_of_range("Positions out of range in CreateTFBSFromFile");
    }

    vector<unsigned long int> offset((rTo - rFrom + 1));
    try {
        chrIdxFile.seekg((rFrom - 1)*4);
    } catch (ios::failure ex) {
        cerr << "Error seeking position in index file" << endl;
        exit(-1);
    }

    for (i = rFrom; i <= rTo; i++) {
        chrIdxFile.read((char *) &b, 4);
        offset[i - rFrom] = b;
    }

    for (i = rFrom; i <= rTo; i++) {
        if (offset[ i - rFrom ] > 0) {
            if (!seek_performed) {
                try {
                    chrTibFile.seekg(2 * offset[ i - rFrom ]);
                } catch (ios::failure ex) {
                    cerr << "Error seeking position in index file" << endl;
                    exit(-1);
                }
            }
            for (j = i + 1; j <= rTo; j++) {
                if (offset[ j - rFrom ] > 0) {
                    tfbsPtr = make_shared<TFBSList>();
                    for (k = 0; k < static_cast<long int> (offset[ j - rFrom ] - offset[ i - rFrom ]); k++) {
                        chrTibFile.read((char *) &intFromByte, 2);
                        tfbsElement = make_shared<TFBS>(i - rFrom, intFromByte);
                        tfbsPtr->GetElements().push_back(tfbsElement);
                    }

                    tfbsList.push_back(tfbsPtr);
                    break;
                }
            }
            i = j - 1;
        }
    }


    tfbs.clear();
    for (auto tfbsListIt = tfbsList.begin(); tfbsListIt != tfbsList.end(); ++tfbsListIt) {
        tfbsPtr = *tfbsListIt;

        for (auto tfbsPtrIt = tfbsPtr->GetElements().begin(); tfbsPtrIt != tfbsPtr->GetElements().end(); ++tfbsPtrIt) {
            tfbsElement = *tfbsPtrIt;

            start = rFrom + tfbsElement->getDelta();
            end = start + pwmIndex[tfbsElement->getIndex() - 1]->getLen() - 1;
            if ((start >= from) && (end <= to)) {
                tfbsElement->setStart(start);
                tfbsElement->setEnd(end);
                tfbs.push_back(tfbsElement);
            }
        }
    }
}



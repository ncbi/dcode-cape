/* 
 * File:   TFBSFactory.cpp
 * Author: veraalva
 * 
 * Created on March 24, 2016, 12:46 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <dirent.h>

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"

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
    this->chrIdxFile = NULL;
    this->chrTibFile = NULL;
}

TFBSFactory::~TFBSFactory() {
    for (auto it = pwmIndex.begin(); it != pwmIndex.end(); ++it) {
        delete (*it);
    }
    for (auto it = tfbs.begin(); it != tfbs.end(); ++it) {
        delete (*it);
    }
    for (auto it = tfbsFileIndex.begin(); it != tfbsFileIndex.end(); ++it) {
        pair<FILE *, FILE*> cPair = it->second;
        if (cPair.first) fclose(cPair.first);
        if (cPair.second) fclose(cPair.second);
    }
}

void TFBSFactory::createTFBSFileIndexMap(std::string dirName, std::string prefix, std::string idxExtension, std::string tibExtension) {
    struct dirent *dp;
    DIR *dirp = (DIR *) checkPointerError(opendir(dirName.c_str()), "Can't open input directory", __FILE__, __LINE__, -1);

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
            pair<FILE *, FILE*> cPair;
            fName.replace(fName.size() - idxExtension.size(), idxExtension.size(), "");
            cPair.first = (FILE *) checkPointerError(fopen((dirName + "/" + fName + idxExtension).c_str(), "rb"), "Can't open input file", __FILE__, __LINE__, -1);
            cPair.second = (FILE *) checkPointerError(fopen((dirName + "/" + fName + tibExtension).c_str(), "rb"), "Can't open input file", __FILE__, __LINE__, -1);
            if (!cPair.first || !cPair.second) {
                cerr << "ERROR!!" << endl;
                cerr << "Can't open TFBS index files" << endl;
                exit(-1);
            }
            tfbsFileIndex.insert(pair<string, pair < FILE *, FILE*>>(fName, cPair));
        }
    }
    closedir(dirp);
}

void TFBSFactory::createPWMIndexFromTibInfoFile(std::string tibInfoFileName) {
    FileParserFactory fParser(tibInfoFileName);

    try {
        while (fParser.iterate('#', " ")) {
            if (fParser.getNWords() != 3) {
                cerr << "Tib info file with a wrong format " << endl;
                exit(-1);
            }
            Tib *tib = new Tib();
            tib->setName(fParser.getWords()[2]);
            tib->setLen(atol(fParser.getWords()[1]));
            if (longestPWM < tib->getLen()) {
                longestPWM = tib->getLen();
            }
            pwmIndex.push_back(tib);
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

class TFBSList {
public:
    TFBSList();
    virtual ~TFBSList();

    std::vector<TFBS*>& GetElements() {
        return elements;
    }

private:
    std::vector<TFBS *> elements;
};

TFBSList::TFBSList() {

}

TFBSList::~TFBSList() {

}

void TFBSFactory::extractTFBSFromFile(long int from, long int to, Seq *chr) {
    long int i, j, k;
    unsigned long int *offset = NULL;
    uint32_t b;
    uint16_t intFromByte;
    bool seek_performed = false;
    long int start, end, rFrom, rTo;
    TFBS *tfbsElement;
    TFBSList *tfbsPtr;
    std::vector<TFBSList *> tfbsList;
    unordered_map<string, pair < FILE *, FILE *>>::iterator tfbsFileIndexIt;

    if (chr->getId().compare(currentChr) != 0) {
        tfbsFileIndexIt = tfbsFileIndex.find(chr->getId());
        if (tfbsFileIndexIt == tfbsFileIndex.end()) {
            throw NotFoundException("Indexes files are not available");
        }

        currentChr = chr->getId();
        chrIdxFile = tfbsFileIndexIt->second.first;
        chrTibFile = tfbsFileIndexIt->second.second;
    }

    if (!chrIdxFile || !chrTibFile) {
        throw NotFoundException("Current indexes pointers are NULL");
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
        throw OutOfRangeException("Positions out of range in CreateTFBSFromFile");
    }

    offset = (unsigned long int *) allocate(sizeof (unsigned long int) * (rTo - rFrom + 1), __FILE__, __LINE__);
    if (fseek(chrIdxFile, (rFrom - 1)*4, SEEK_SET) == -1) {
        checkPointerError(NULL, "Error trying to seek the offset position in the index file", __FILE__, __LINE__, -1);
    }

    for (i = rFrom; i <= rTo; i++) {
        fread(&b, 4, 1, chrIdxFile);
        if (ferror(chrIdxFile)) {
            checkPointerError(NULL, "Error while reading the index file", __FILE__, __LINE__, -1);
        }
        offset[i - rFrom] = b;
    }

    for (i = rFrom; i <= rTo; i++) {
        if (offset[ i - rFrom ] > 0) {
            if (!seek_performed) {
                if (fseek(chrTibFile, 2 * offset[ i - rFrom ], SEEK_SET) == -1) {
                    checkPointerError(NULL, "Error trying to seek the offset position in the tib file", __FILE__, __LINE__, -1);
                }
            }
            for (j = i + 1; j <= rTo; j++) {
                if (offset[ j - rFrom ] > 0) {
                    tfbsPtr = new TFBSList();
                    for (k = 0; k < static_cast<long int> (offset[ j - rFrom ] - offset[ i - rFrom ]); k++) {
                        fread(&intFromByte, 2, 1, chrTibFile);
                        if (ferror(chrTibFile)) {
                            checkPointerError(NULL, "Error while reading the tib file", __FILE__, __LINE__, -1);
                        }

                        tfbsElement = new TFBS(i - rFrom, intFromByte);
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
            } else {
                delete (*tfbsPtrIt);
            }
        }
        delete (*tfbsListIt);
    }

    free(offset);
}



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

#include "Exceptions.h"
#include "FastaFactory.h"
#include "TFBSFactory.h"
#include "Global.h"

using namespace std;
using namespace exceptions;
using namespace fasta;
using namespace tfbs;

Tib::Tib(const char* n, long int l) {
    name = n;
    len = l;
}

Tib::~Tib() {

}

TFBS::TFBS(unsigned long int d, short int i) {
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

void TFBSFactory::CreateTFBSFileIndexMap(char* dirName, const char *prefix, const char *idxExtension, const char *tibExtension) {
    struct dirent *dp;
    bool read = false;

    pair<FILE *, FILE*> cPair;
    char *index;
    size_t len = strlen(dirName);
    char *fileName = (char *) allocate(sizeof (char) * (len + 1), __FILE__, __LINE__);

    DIR *dirp = (DIR *) checkPointerError(opendir(dirName), "Can't open input directory", __FILE__, __LINE__, -1);

    while ((dp = readdir(dirp)) != NULL) {
        if (Global::instance()->isDebug3()) {
            cout << "\tDEBUG3 ==> Found file: " << dp->d_name << endl;
        }
        if (dp->d_name[0] != '.') {
            if (!prefix && !idxExtension) {
                read = true;
            } else {
                if (prefix && !idxExtension) {
                    if (strncmp(dp->d_name, prefix, strlen(prefix)) == 0) read = true;
                } else if (!prefix && idxExtension) {
                    if (strbcmp(dp->d_name, idxExtension) == 0) read = true;
                } else if (prefix && idxExtension) {
                    if (strncmp(dp->d_name, prefix, strlen(prefix)) == 0 && strbcmp(dp->d_name, idxExtension) == 0) read = true;
                }
            }
        }
        if (read) {
            if (len < strlen(dirName) + strlen(dp->d_name) + strlen(tibExtension) + 2) {
                len = strlen(dirName) + strlen(dp->d_name) + strlen(tibExtension) + 2;
                fileName = (char *) reallocate(fileName, sizeof (char) * len, __FILE__, __LINE__);
            }
            sprintf(fileName, "%s/%s", dirName, dp->d_name);
            if (Global::instance()->isDebug3()) {
                cout << "\tDEBUG3 ==> Opening file: " << fileName << endl;
            }
            cPair.first = (FILE *) checkPointerError(fopen(fileName, "rb"), "Can't open input file", __FILE__, __LINE__, -1);
            index = strstr(dp->d_name, idxExtension);
            strncpy(index, tibExtension, strlen(tibExtension));
            sprintf(fileName, "%s/%s", dirName, dp->d_name);
            if (Global::instance()->isDebug3()) {
                cout << "\tDEBUG3 ==> Opening file: " << fileName << endl;
            }
            cPair.second = (FILE *) checkPointerError(fopen(fileName, "rb"), "Can't open input file", __FILE__, __LINE__, -1);

            *index = 0;
            if (Global::instance()->isInfo()) {
                cout << "\tDEBUG3 ==> Chr file name: " << dp->d_name << endl;
            }

            if (!cPair.first || !cPair.second) {
                cerr << "ERROR!!" << endl;
                cerr << "Can't open TFBS index files" << endl;
                exit(-1);
            }

            tfbsFileIndex.insert(pair<string, pair < FILE *, FILE*>>(dp->d_name, cPair));
            
            strncpy(index, idxExtension, strlen(idxExtension));
            if (Global::instance()->isDebug3()) {
                cout << "\tDEBUG3 ==> Setting index file name to the initial state: " << dp->d_name << endl;
            }

            read = false;
        }
    }
    closedir(dirp);
    free(fileName);
}

void TFBSFactory::CreatePWMIndexFromTibInfoFile(const char* tibInfoFileName) {
    FILE *tibInfoFile = (FILE *) checkPointerError(fopen(tibInfoFileName, "r"), "Can't open Tib info file", __FILE__, __LINE__, -1);

    size_t bufferSize, read, backupLineSize;
    char *buffer, *newLine, *str, *backupLine, *completeLine;
    char **fields = NULL;
    size_t fieldsSize = 0;
    Tib *tib;

    backupLineSize = bufferSize = 100000000;
    buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
    backupLine = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);

    *backupLine = 0;
    while (!feof(tibInfoFile)) {
        read = fread(buffer, sizeof (char), bufferSize, tibInfoFile);
        buffer[read] = 0;
        if (feof(tibInfoFile)) {
            if (buffer[read - 1] != '\n') {
                buffer[read] = '\n';
                buffer[read + 1] = 0;
            }
        }
        str = buffer;
        while ((newLine = strchr(str, '\n')) != NULL) {
            *newLine = 0;
            if (*backupLine != 0) {
                if (strlen(backupLine) + strlen(str) + 1 > backupLineSize) {
                    backupLineSize += backupLineSize;
                    backupLine = (char *) reallocate(backupLine, sizeof (char) * (backupLineSize + 1), __FILE__, __LINE__);
                }
                strcat(backupLine, str);
                completeLine = backupLine;
            } else {
                completeLine = str;
            }

            if (*completeLine != '#') {
                fieldsSize = splitString(&fields, completeLine, " ");
                if (fieldsSize != 3) {
                    fprintf(stderr, "%s\n", completeLine);
                    printLog(stderr, "Tib info file with a wrong format ", __FILE__, __LINE__, -1);
                }
                tib = new Tib(fields[2], atol(fields[1]));
                if (longestPWM < tib->GetLen()) {
                    longestPWM = tib->GetLen();
                }
                pwmIndex.push_back(move(tib));

                freeArrayofPointers((void **) fields, fieldsSize);
            }
            *backupLine = 0;
            str = newLine + 1;
        }

        if (strlen(str) > 0) {
            if (strlen(backupLine) + strlen(str) + 1 > backupLineSize) {
                backupLineSize += backupLineSize;
                backupLine = (char *) reallocate(backupLine, sizeof (char) * (backupLineSize + 1), __FILE__, __LINE__);
            }
            strcat(backupLine, str);
        }
    }

    free(buffer);
    free(backupLine);
    fclose(tibInfoFile);
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

void TFBSFactory::ExtractTFBSFromFile(long int from, long int to, Fasta *chr) {
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

    if (chr->GetId().compare(currentChr) != 0) {
        tfbsFileIndexIt = tfbsFileIndex.find(chr->GetId());
        if (tfbsFileIndexIt == tfbsFileIndex.end()) {
            throw OutOfRangeException("Indexes files are not available");
        }

        currentChr = chr->GetId();
        chrIdxFile = tfbsFileIndexIt->second.first;
        chrTibFile = tfbsFileIndexIt->second.second;
    }

    if (!chrIdxFile || !chrTibFile) {
        throw OutOfRangeException("Current indexes pointers are NULL");
    }

    rFrom = from - longestPWM - 1;
    rTo = to + longestPWM;
    if (rFrom < 1) {
        rFrom = 1;
    }
    if (rTo > static_cast<long int>(chr->GetLength())) {
        rTo = chr->GetLength();
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
                    for (k = 0; k < static_cast<long int>(offset[ j - rFrom ] - offset[ i - rFrom ]); k++) {
                        fread(&intFromByte, 2, 1, chrTibFile);
                        if (ferror(chrTibFile)) {
                            checkPointerError(NULL, "Error while reading the tib file", __FILE__, __LINE__, -1);
                        }

                        tfbsElement = new TFBS(i - rFrom, intFromByte);
                        tfbsPtr->GetElements().push_back(move(tfbsElement));
                    }

                    tfbsList.push_back(move(tfbsPtr));
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

            start = rFrom + tfbsElement->GetDelta();
            end = start + pwmIndex[tfbsElement->GetIndex() - 1]->GetLen() - 1;
            if ((start >= from) && (end <= to)) {
                tfbsElement->SetStart(start);
                tfbsElement->SetEnd(end);
                tfbs.push_back(move(tfbsElement));
                if (Global::instance()->isDebug3()) {
                    cerr << "\tDEBUG3 ==> \tTFBS\t" << pwmIndex[tfbsElement->GetIndex() - 1]->GetName() << "\t" << start << "\t" << end << "\t" << tfbsElement->GetStrand() << endl;
                }
            } else {
                delete (*tfbsPtrIt);
            }
        }
        delete (*tfbsListIt);
    }

    free(offset);
}



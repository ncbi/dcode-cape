/* 
 * File:   FastaFactory.cpp
 * Author: veraalva
 * 
 * Created on February 10, 2016, 3:41 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <utility>
#include <stdbool.h>
#include <dirent.h>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "FastaFactory.h"
#include "Global.h"
#include "TimeUtils.h"

using namespace std;
using namespace fasta;

Fasta::Fasta() {
    this->seq = NULL;
    this->length = 0;
}

Fasta::Fasta(const Fasta& orig) {
    this->seq = NULL;
    this->length = 0;
}

Fasta::~Fasta() {
    if (this->seq) free(this->seq);
}

FastaFactory::FastaFactory() {
}

FastaFactory::FastaFactory(const FastaFactory& orig) {
}

FastaFactory::~FastaFactory() {
}

Fasta *FastaFactory::GetFastaFromID(string id) {
    unordered_map<string, long unsigned int>::iterator it = this->fastaUnorderedMap.find(id);
    if (it != this->fastaUnorderedMap.end()) {
        if (this->fastaVector.size() > it->second) {
            return this->fastaVector[it->second];
        }
    }
    return nullptr;
}

long unsigned int FastaFactory::ParseFastaFile(FILE *fName, int numberSeqTotalRead, bool cleanContainers, bool binary) {
    int numberSeqCurrentRead;
    unsigned long int i, j, len;
    off_t pos;
    size_t bufferSize, read, readTotal, seq_length, numLines, str_length;
    char *buffer, *line, *s, *str, *seq, *seq_end;
    Fasta *fasta;

    if (cleanContainers) {
        this->fastaVector.clear();
    }

    pos = ftello(fName);
    numberSeqCurrentRead = 0;
    if (binary) {
        fread(&i, sizeof (unsigned long int), 1, fName);
        for (j = 0; j < i; j++) {
            fasta = new Fasta();

            fread(&len, sizeof (unsigned long int), 1, fName);
            buffer = (char *) allocate(sizeof (char) * len, __FILE__, __LINE__);
            fread(buffer, sizeof (char), len, fName);
            fasta->SetId(buffer);
            fread(&len, sizeof (unsigned long int), 1, fName);
            buffer = (char *) reallocate(buffer, sizeof (char) * len, __FILE__, __LINE__);
            fread(buffer, sizeof (char), len, fName);
            fasta->SetSeq(&buffer);
            fasta->SetLength(len - 1);
            if (Global::instance()->isDebug3()) cerr << "\tDEBUG3 ==> Adding seq: " << fasta->GetId() << " with " << fasta->GetLength() << " bp" << endl;
            this->fastaUnorderedMap.insert(pair<string, long unsigned int> (fasta->GetId(), this->fastaVector.size()));
            this->fastaVector.push_back(move(fasta));
            if (++numberSeqCurrentRead == numberSeqTotalRead) {
                break;
            }
        }
    } else {
        bufferSize = 100000000;
        buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
        seq = seq_end = NULL;

        readTotal = numLines = seq_length = 0;
        while (!feof(fName)) {
            read = fread(buffer, sizeof (char), bufferSize, fName);
            buffer[read] = '\0';
            str = buffer;
            while (1) {
                line = strchr(str, '\n');
                if (line) *line = '\0';
                if (*str != '\0' && *str != '\n') {
                    str_length = strlen(str);
                    if (*str == '>') {
                        if (seq_end != seq) {
                            seq = (char *) reallocate(seq, sizeof (char) * (seq_length + 1), __FILE__, __LINE__);
                            *(seq + seq_length) = '\0';
                            fasta->SetSeq(&seq);
                            fasta->SetLength(strlen(seq));
                            if (Global::instance()->isDebug3()) cerr << "\tDEBUG3 ==> Adding seq: " << fasta->GetId().c_str() << " with " << fasta->GetLength() << " bp" << endl;
                            this->fastaUnorderedMap.insert(pair<string, long unsigned int> (fasta->GetId(), this->fastaVector.size()));
                            this->fastaVector.push_back(move(fasta));
                            if (++numberSeqCurrentRead == numberSeqTotalRead) {
                                pos = pos + sizeof (char) * (readTotal + numLines);
                                fseeko(fName, pos, SEEK_SET);
                                break;
                            }
                        }
                        seq_length = 0;
                        seq = NULL;
                        readTotal += str_length;
                        fasta = new Fasta();
                        s = strndup(str + 1, str_length - 1);
                        fasta->SetId(s);
                        free(s);
                    } else {
                        checkPointerError(fasta, "Fasta file does not start with the header (>)", __FILE__, __LINE__, -1);
                        readTotal += str_length;
                        seq_length += str_length;

                        seq = (char *) reallocate(seq, sizeof (char) * seq_length, __FILE__, __LINE__);
                        seq_end = seq + (seq_length - str_length);
                        memcpy(seq_end, str, str_length);
                        seq_end += str_length;
                    }
                }
                if (!line) break;
                str = line + 1;
                numLines++;
            }
            if (*str == '>') break;
        }

        if (numberSeqCurrentRead != numberSeqTotalRead && seq && strlen(seq) > 0) {
            numberSeqCurrentRead++;
            seq = (char *) reallocate(seq, sizeof (char) * (seq_length + 1), __FILE__, __LINE__);
            *(seq + seq_length) = '\0';
            fasta->SetSeq(&seq);
            fasta->SetLength(seq_length);
            if (Global::instance()->isDebug3()) cerr << "\tDEBUG3 ==> Adding seq: " << fasta->GetId().c_str() << " with " << fasta->GetLength() << " bp" << endl;
            this->fastaUnorderedMap.insert(pair<string, long unsigned int> (fasta->GetId(), this->fastaVector.size()));
            this->fastaVector.push_back(move(fasta));
        }
        free(buffer);
    }

    return numberSeqCurrentRead;
}

void FastaFactory::LoadFastaInDirectory(char* dirName, const char *prefix, const char *sufix, bool binary) {
    int seqs = 0;
    FILE *fName;
    struct dirent *dp;
    bool read = false;

    size_t len = strlen(dirName);
    char *fileName = (char *) allocate(sizeof (char) * (len + 1), __FILE__, __LINE__);

    DIR *dirp = (DIR *) checkPointerError(opendir(dirName), "Can't open input directory", __FILE__, __LINE__, -1);

    while ((dp = readdir(dirp)) != NULL) {
        if (Global::instance()->isDebug3()) {
            cout << "\tDEBUG3 ==> Found file: " << dp->d_name << endl;
        }
        if (dp->d_name[0] != '.') {
            if (!prefix && !sufix) {
                read = true;
            } else {
                if (prefix && !sufix) {
                    if (strncmp(dp->d_name, prefix, strlen(prefix)) == 0) read = true;
                } else if (!prefix && sufix) {
                    if (strbcmp(dp->d_name, sufix) == 0) read = true;
                } else if (prefix && sufix) {
                    if (strncmp(dp->d_name, prefix, strlen(prefix)) == 0 && strbcmp(dp->d_name, sufix) == 0) read = true;
                }
            }
        }
        if (read) {
            if (len < strlen(dirName) + strlen(dp->d_name) + 2) {
                len = strlen(dirName) + strlen(dp->d_name) + 2;
                fileName = (char *) reallocate(fileName, sizeof (char) * len, __FILE__, __LINE__);
            }
            sprintf(fileName, "%s/%s", dirName, dp->d_name);
            if (Global::instance()->isInfo()) {
                TimeUtils::instance()->SetClock();
                cout << "\tINFO ==> Parsing file: " << fileName << endl;
            }
            if (binary) {
                fName = (FILE *) checkPointerError(fopen(fileName, "rb"), "Can't open input file", __FILE__, __LINE__, -1);
            } else {
                fName = (FILE *) checkPointerError(fopen(fileName, "r"), "Can't open input file", __FILE__, __LINE__, -1);
            }
            seqs = this->ParseFastaFile(fName, -1, false, binary);
            fclose(fName);
            read = false;
            if (Global::instance()->isInfo()) cout << "\tINFO ==> " << seqs << " sequences read in " << TimeUtils::instance()->GetTimeSec() << " sec" << endl;
        }
    }
    closedir(dirp);
    free(fileName);
}

void FastaFactory::WriteSequencesToFile(char* fileName, bool binary) {
    unsigned long int i, len;
    FILE *outputFile = NULL;
    Fasta *f;
    char t = '\0';

    if (binary) {
        outputFile = (FILE *) checkPointerError(fopen(fileName, "wb"), "Can't open output file", __FILE__, __LINE__, -1);
        i = this->fastaVector.size();
        fwrite(&i, sizeof (unsigned long int), 1, outputFile);
        for (auto it = this->fastaVector.begin(); it != this->fastaVector.end(); ++it) {
            f = *it;
            len = f->GetId().size() + 1;
            fwrite(&len, sizeof (unsigned long int), 1, outputFile);
            fwrite(f->GetId().c_str(), sizeof (char), len, outputFile);
            len = f->GetLength() + 1;
            fwrite(&len, sizeof (unsigned long int), 1, outputFile);
            fwrite(f->GetSeq(), sizeof (char), len, outputFile);
        }
    } else {
        outputFile = (FILE *) checkPointerError(fopen(fileName, "w"), "Can't open output file", __FILE__, __LINE__, -1);
        for (auto it = this->fastaVector.begin(); it != this->fastaVector.end(); ++it) {
            f = *it;
            fprintf(outputFile, ">%s\n", f->GetId().c_str());
            for (i = 0; i < f->GetLength(); i += 50) {
                if (i + 50 < f->GetLength()) {
                    t = f->GetSeq()[i + 50];
                    f->GetSeq()[i + 50] = '\0';
                }
                fprintf(outputFile, "%s\n", f->GetSeq() + i);
                if (i + 50 < f->GetLength()) {
                    f->GetSeq()[i + 50] = t;
                }
            }
        }
    }

    fclose(outputFile);
}

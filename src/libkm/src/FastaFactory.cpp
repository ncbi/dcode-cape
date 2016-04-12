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

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <utility>

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

Seq::Seq(const Seq& orig) {
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
    for (auto it = sequenceContainer.begin(); it != sequenceContainer.end(); ++it) {
        delete it->second;
    }
}

long unsigned int FastaFactory::parseFastaFile(FILE *fName, int numberSeqTotalRead, bool cleanContainers, bool binary) {

    int numberSeqCurrentRead;
    unsigned long int i, len;
    off_t pos;
    char *buffer, *seq, *seq_end;
    Seq *fasta = NULL;

    if (cleanContainers) {
        sequenceContainer.clear();
    }

    pos = ftello(fName);
    numberSeqCurrentRead = 0;
    if (binary) {
        fread(&i, sizeof (unsigned long int), 1, fName);
        for (unsigned long int j = 0; j < i; j++) {
            fasta = new Seq();

            fread(&len, sizeof (unsigned long int), 1, fName);
            buffer = (char *) allocate(sizeof (char) * len, __FILE__, __LINE__);
            fread(buffer, sizeof (char), len, fName);
            fasta->setId(buffer);
            fread(&len, sizeof (unsigned long int), 1, fName);
            buffer = (char *) reallocate(buffer, sizeof (char) * len, __FILE__, __LINE__);
            fread(buffer, sizeof (char), len, fName);
            fasta->setSeq(&buffer);
            fasta->setLength(len - 1);
            if (Global::instance()->isDebug3()) cerr << "\tDEBUG3 ==> Adding seq: " << fasta->getId() << " with " << fasta->getLength() << " bp" << endl;
            sequenceContainer.insert(pair<string, Seq *> (fasta->getId(), fasta));
            if (++numberSeqCurrentRead == numberSeqTotalRead) {
                break;
            }
        }
    } else {
        FileParserFactory fParser(fName);
        seq = seq_end = NULL;

        size_t str_length;
        size_t readTotal = 0;
        size_t numLines = 0;
        size_t seq_length = 0;

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
                    sequenceContainer.insert(pair<string, Seq *> (fasta->getId(),fasta));
                    if (++numberSeqCurrentRead == numberSeqTotalRead) {
                        pos = pos + sizeof (char) * (readTotal + numLines);
                        fseeko(fName, pos, SEEK_SET);
                        break;
                    }
                }
                seq_length = 0;
                seq = NULL;
                readTotal += str_length;
                fasta = new Seq();
                char *s = strndup(line + 1, str_length - 1);
                fasta->setId(s);
                free(s);
            } else {
                checkPointerError(fasta, "Fasta file does not start with the header (>)", __FILE__, __LINE__, -1);
                readTotal += str_length;
                seq_length += str_length;

                seq = (char *) reallocate(seq, sizeof (char) * seq_length, __FILE__, __LINE__);
                seq_end = seq + (seq_length - str_length);
                memcpy(seq_end, line, str_length);
                seq_end += str_length;
            }
            numLines++;
        }

        if (numberSeqCurrentRead != numberSeqTotalRead && seq && strlen(seq) > 0) {
            numberSeqCurrentRead++;
            seq = (char *) reallocate(seq, sizeof (char) * (seq_length + 1), __FILE__, __LINE__);
            *(seq + seq_length) = 0;
            fasta->setSeq(&seq);
            fasta->setLength(seq_length);
            if (Global::instance()->isDebug3()) cout << "\tDEBUG3 ==> Adding seq: " << fasta->getId().c_str() << " with " << fasta->getLength() << " bp" << endl;
            sequenceContainer.insert(pair<string, Seq *> (fasta->getId(), fasta));
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
            FILE *file;
            if (binary) {
                file = (FILE *) checkPointerError(fopen((dirName + "/" + fName).c_str(), "rb"), "Can't open input file", __FILE__, __LINE__, -1);
            } else {
                file = (FILE *) checkPointerError(fopen((dirName + "/" + fName).c_str(), "r"), "Can't open input file", __FILE__, __LINE__, -1);
            }
            int seqs = parseFastaFile(file, -1, false, binary);
            fclose(file);
            if (Global::instance()->isInfo()) cout << "\tINFO ==> " << seqs << " sequences read in " << TimeUtils::instance()->getTimeSec() << " sec" << endl;
        }
    }
    closedir(dirp);
}

void FastaFactory::writeSequencesToFile(std::string fileName, bool binary) {
    unsigned long int i, len;
    FILE *outputFile = NULL;
    Seq *f;

    if (binary) {
        outputFile = (FILE *) checkPointerError(fopen(fileName.c_str(), "wb"), "Can't open output file", __FILE__, __LINE__, -1);
        i = sequenceContainer.size();
        fwrite(&i, sizeof (unsigned long int), 1, outputFile);
        for (auto it = sequenceContainer.begin(); it != sequenceContainer.end(); ++it) {
            f = it->second;
            len = f->getId().size() + 1;
            fwrite(&len, sizeof (unsigned long int), 1, outputFile);
            fwrite(f->getId().c_str(), sizeof (char), len, outputFile);
            len = f->getLength() + 1;
            fwrite(&len, sizeof (unsigned long int), 1, outputFile);
            fwrite(f->getSeq(), sizeof (char), len, outputFile);
        }
    } else {
        outputFile = (FILE *) checkPointerError(fopen(fileName.c_str(), "w"), "Can't open output file", __FILE__, __LINE__, -1);
        for (auto it = sequenceContainer.begin(); it != sequenceContainer.end(); ++it) {
            f = it->second;
            fprintf(outputFile, ">%s\n", f->getId().c_str());
            for (i = 0; i < f->getLength(); i += 50) {
                char t = 0;
                if (i + 50 < f->getLength()) {
                    t = f->getSeq()[i + 50];
                    f->getSeq()[i + 50] = 0;
                }
                fprintf(outputFile, "%s\n", f->getSeq() + i);
                if (i + 50 < f->getLength()) {
                    f->getSeq()[i + 50] = t;
                }
            }
        }
    }

    fclose(outputFile);
}

/* 
 * File:   FimoFactory.cpp
 * Author: veraalva
 * 
 * Created on March 4, 2016, 3:48 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <locale>
#include <unordered_map>
#include <set>
#include <vector>
#include <algorithm>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "Global.h"
#include "TimeUtils.h"
#include "FimoFactory.h"

using namespace std;
using namespace fimo;

FimoFactory::FimoFactory() {
}

FimoFactory::FimoFactory(const FimoFactory& orig) {
}

FimoFactory::~FimoFactory() {

}

void FimoFactory::CreateTissueIndexFromFiles(char *pwm_EnsembleID, char *tissue_file) {
    FILE *pwmEnsembleIDFile = (FILE *) checkPointerError(fopen(pwm_EnsembleID, "r"), "Can't open pwm_tFName file", __FILE__, __LINE__, -1);
    FILE *tissueFile = (FILE *) checkPointerError(fopen(tissue_file, "r"), "Can't open tissue file", __FILE__, __LINE__, -1);

    char *line = NULL;
    size_t len = 0;
    char **fields = NULL;
    size_t fieldsSize = 0;
    char **fields2 = NULL;
    size_t fields2Size = 0;
    unordered_map<string, set < string>> tFNameReverseMap;
    unordered_map<string, set < string>>::iterator tFNameReverseMapIt;

    while (getline(&line, &len, pwmEnsembleIDFile) != -1) {
        if (*line != '#') {
            if (*(line + (strlen(line) - 1)) == '\n') *(line + (strlen(line) - 1)) = '\0';
            fieldsSize = splitString(&fields, line, "\t");
            if (fieldsSize < 2) {
                printLog(stderr, "Input pwm_tFName file with a wrong format ", __FILE__, __LINE__, -1);
            }
            fields2Size = splitString(&fields2, fields[fieldsSize - 1], ";");
            set<string> a;
            for (int i = 0; i < fields2Size; i++) {
                a.insert(fields2[i]);

                tFNameReverseMapIt = tFNameReverseMap.find(fields2[i]);
                if (tFNameReverseMapIt == tFNameReverseMap.end()) {
                    set<string> s;
                    s.insert(fields[0]);
                    tFNameReverseMap.insert(pair<string, set < string >> (fields2[i], s));
                } else {
                    tFNameReverseMapIt->second.insert(fields[0]);
                }
            }
            freeArrayofPointers((void **) fields2, fields2Size);
            freeArrayofPointers((void **) fields, fieldsSize);
        }
    }
    fclose(pwmEnsembleIDFile);

    getline(&line, &len, tissueFile);
    if (*(line + (strlen(line) - 1)) == '\n') *(line + (strlen(line) - 1)) = '\0';
    fieldsSize = splitString(&fields, line, "\t");

    while (getline(&line, &len, tissueFile) != -1) {
        if (*line != '#') {
            if (*(line + (strlen(line) - 1)) == '\n') *(line + (strlen(line) - 1)) = '\0';
            fields2Size = splitString(&fields2, line, "\t");
            if (fields2Size != fieldsSize) {
                printLog(stderr, "Tissue file with a wrong format ", __FILE__, __LINE__, -1);
            }

            tFNameReverseMapIt = tFNameReverseMap.find(fields2[0]);
            if (tFNameReverseMapIt != tFNameReverseMap.end()) {
                for (auto it1 = tFNameReverseMapIt->second.begin(); it1 != tFNameReverseMapIt->second.end(); ++it1) {

                    auto tissueMapIt = tissueIndex.find(*it1);

                    for (int i = 1; i < fields2Size; i++) {
                        if (tissueMapIt == tissueIndex.end()) {
                            pair<string, double> tissueEnsemblePair(fields2[0], strtod(fields2[i], NULL));
                            unordered_map<string, pair<string, double>> tissues;
                            tissues.insert(pair<string, pair<string, double>>(fields[i], tissueEnsemblePair));
                            tissueIndex.insert(pair<string, unordered_map<string, pair<string, double>>>(*it1, tissues));
                            tissueMapIt = tissueIndex.find(*it1);
                        } else {
                            auto tissueIt = tissueMapIt->second.find(fields[i]);
                            if (tissueIt == tissueMapIt->second.end()) {
                                pair<string, double> tissueEnsemblePair(fields2[0], strtod(fields2[i], NULL));
                                tissueMapIt->second.insert(pair<string, pair<string, double>>(fields[i], tissueEnsemblePair));
                            } else {
                                if (tissueIt->second.second < strtod(fields2[i], NULL)) {
                                    pair<string, double> tissueEnsemblePair(fields2[0], strtod(fields2[i], NULL));
                                    tissueIt->second = pair<string, double>(fields2[0], strtod(fields2[i], NULL));
                                }
                            }
                        }
                    }
                }
            }
            freeArrayofPointers((void **) fields2, fields2Size);
        }
    }
    fclose(tissueFile);
    freeArrayofPointers((void **) fields, fieldsSize);

    if (line) free(line);
}

void FimoFactory::CreateCutoffIndexFromFile(char* cutoffFileName, int column) {
    FILE *cutoffFile = (FILE *) checkPointerError(fopen(cutoffFileName, "r"), "Can't open pwm_tFName file", __FILE__, __LINE__, -1);

    char *line = NULL;
    size_t len = 0;
    char **fields = NULL;
    size_t fieldsSize = 0;

    while (getline(&line, &len, cutoffFile) != -1) {
        if (*line != '#') {
            if (*(line + (strlen(line) - 1)) == '\n') *(line + (strlen(line) - 1)) = '\0';
            fieldsSize = splitString(&fields, line, " ");
            if (fieldsSize <= column) {
                fprintf(stderr, "%s\n", line);
                printLog(stderr, "Input cutoff file with less columns of required ", __FILE__, __LINE__, -1);
            }
            cutoffIndex.insert(pair<string, double>(fields[0], strtod(fields[column], NULL)));
            freeArrayofPointers((void **) fields, fieldsSize);
        }
    }

    if (line) free(line);
    fclose(cutoffFile);
}

void FimoFactory::ParseFimoOutput(char* fimoOuputName, char *tissueCode, unsigned long int snpPos) {
    int i;
    FILE *fimoOuputFile = (FILE *) checkPointerError(fopen(fimoOuputName, "r"), "Can't open pwm_tFName file", __FILE__, __LINE__, -1);
    size_t bufferSize, read, backupLineSize;
    char *buffer, *newLine, *str, *backupLine, *completeLine;
    char **fields = NULL;
    size_t fieldsSize = 0;
    char **fields2 = NULL;
    size_t fields2Size = 0;
    Fimo *f;
    unordered_map<string, set < Fimo *, PointerCompare>> fimoMap;
    unordered_map <string, set < Fimo *, PointerCompare>>::iterator fimoIt;
    unordered_map<string, set < string>> ensemblMap;
    unordered_map<string, set < string>>::iterator ensemblMapIt;
    pair < unordered_map<string, set < string>>::iterator, bool> res;

    backupLineSize = bufferSize = 100000000;
    buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
    backupLine = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);

    *backupLine = '\0';
    while (!feof(fimoOuputFile)) {
        read = fread(buffer, sizeof (char), bufferSize, fimoOuputFile);
        buffer[read] = '\0';
        str = buffer;
        while (*str != '\0') {
            newLine = strchr(str, '\n');
            if (newLine) *newLine = '\0';
            if (*backupLine != '\0') {
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
                fieldsSize = splitString(&fields, completeLine, "\t");
                if (fieldsSize != 8) {
                    printLog(stderr, "FIMO output file with a wrong format ", __FILE__, __LINE__, -1);
                }

                res = ensemblMap.insert(pair<string, set < string >> (fields[1], set < string>()));
                ensemblMapIt = res.first;

                fields2Size = splitString(&fields2, fields[0], ";");
                for (i = 0; i < fields2Size; i++) {
                    pair<string, double> motifExpression = GetTissueValue(fields2[i], tissueCode);
                    if (fabs(motifExpression.second) >= 1.0) {
                        double pValue = strtod(fields[6], NULL);
                        if (pValue < GetCutoffValue(fields2[i])) {

                            fimoIt = fimoMap.find(fields[1]);
                            f = new Fimo();
                            f->SetMotif(fields2[i]);
                            f->SetId(fields[1]);
                            f->SetStart(static_cast<unsigned long int> (atoi(fields[2]) - 1));
                            f->SetEnd(static_cast<unsigned long int> (atoi(fields[3]) - 1));
                            f->SetStrand(fields[4][0]);
                            f->SetScore(strtod(fields[5], NULL));
                            f->SetValue(pValue);
                            f->SetSeq(fields[7]);
                            f->SetExpression(motifExpression.second);
                            f->SetExpEnsembl(motifExpression.first);

                            if (fimoIt == fimoMap.end()) {
                                vector<double> v;
                                if (f->GetStart() <= snpPos && f->GetEnd() >= snpPos) {
                                    v.push_back(motifExpression.second);
                                    v.push_back(10E-10);
                                } else {
                                    v.push_back(10E-10);
                                    v.push_back(motifExpression.second);
                                    ensemblMapIt->second.insert(motifExpression.first);
                                    if (f->GetId().compare("rs10000101-A-C") == 0) {
                                        cout << "First " << *f << endl;
                                    }
                                }
                                snpIDMap.insert(pair<string, vector<double>>(fields[1], v));

                                set<Fimo *, PointerCompare> s;
                                s.insert(move(f));
                                fimoMap.insert(pair<string, set < Fimo *, PointerCompare >> (fields[1], s));
                                f = NULL;
                            } else {
                                pair < set < Fimo *, PointerCompare>::iterator, bool> iter;

                                iter = fimoIt->second.insert(move(f));
                                if (iter.second) {
                                    f = *iter.first;
                                    auto snpMapIt = snpIDMap.find(fields[1]);

                                    if (f->GetStart() <= snpPos && f->GetEnd() >= snpPos) {
                                        if (snpMapIt->second[0] < motifExpression.second) {
                                            snpMapIt->second[0] = motifExpression.second;
                                        }
                                    } else {
                                        if (ensemblMapIt->second.find(motifExpression.first) == ensemblMapIt->second.end()) {
                                            ensemblMapIt->second.insert(motifExpression.first);
                                            snpMapIt->second[1] += motifExpression.second;
                                            if (f->GetId().compare("rs10000101-A-C") == 0) {
                                                cout << "Second " << *f << endl;
                                            }
                                        }
                                    }
                                    f = NULL;
                                }
                            }
                            if (f) delete f;
                        }
                    }
                }

                freeArrayofPointers((void **) fields2, fields2Size);
                freeArrayofPointers((void **) fields, fieldsSize);
            }
            *backupLine = '\0';
            if (newLine) str = newLine + 1;
            else str += strlen(str);
        }

        if (strlen(str) > 0) {
            if (strlen(backupLine) + strlen(str) + 1 > backupLineSize) {
                backupLineSize += backupLineSize;
                backupLine = (char *) reallocate(backupLine, sizeof (char) * (backupLineSize + 1), __FILE__, __LINE__);
            }
            strcat(backupLine, str);
        }
    }

    for (auto it = fimoMap.begin(); it != fimoMap.end(); ++it) {
        for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
            delete (*it1);
        }
    }

    fclose(fimoOuputFile);
    free(buffer);
    free(backupLine);
}

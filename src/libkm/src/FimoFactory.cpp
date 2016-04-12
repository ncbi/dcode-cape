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
#include "FileParserFactory.h"
#include "FimoFactory.h"

using namespace std;
using namespace parsers;
using namespace fimo;

FimoFactory::FimoFactory() {
}

FimoFactory::FimoFactory(const FimoFactory& orig) {
}

FimoFactory::~FimoFactory() {

}

void FimoFactory::createTissueIndexFromFiles(std::string pwm_EnsembleID, std::string tissue_file) {
    char **words = NULL;
    size_t nWords = 0;
    size_t wordsSize = 0;
    FileParserFactory fParser(pwm_EnsembleID);
    unordered_map<string, set < string>> tFNameReverseMap;
    unordered_map<string, set < string>>::iterator tFNameReverseMapIt;
    vector<string> header;

    while (fParser.iterate('#', "\t")) {
        if (fParser.getNWords() < 2) {
            printLog(stderr, "Input pwm_tFName file with a wrong format ", __FILE__, __LINE__, -1);
        }
        nWords = strsep_ptr(&words, &wordsSize, fParser.getWords()[fParser.getNWords() - 1], ";");
        set<string> a;
        for (size_t i = 0; i < nWords; i++) {
            a.insert(words[i]);

            tFNameReverseMapIt = tFNameReverseMap.find(words[i]);
            if (tFNameReverseMapIt == tFNameReverseMap.end()) {
                set<string> s;
                s.insert(fParser.getWords()[0]);
                tFNameReverseMap.insert(pair<string, set < string >> (words[i], s));
            } else {
                tFNameReverseMapIt->second.insert(fParser.getWords()[0]);
            }
        }
    }

    fParser.clean();
    fParser.setFileToParse(tissue_file);

    fParser.iterate('#', "\t");
    fParser.wordsToVector(header);

    while (fParser.iterate('#', "\t")) {
        if (fParser.getNWords() != header.size()) {
            printLog(stderr, "Tissue file with a wrong format ", __FILE__, __LINE__, -1);
        }

        tFNameReverseMapIt = tFNameReverseMap.find(fParser.getWords()[0]);
        if (tFNameReverseMapIt != tFNameReverseMap.end()) {
            for (auto it1 = tFNameReverseMapIt->second.begin(); it1 != tFNameReverseMapIt->second.end(); ++it1) {

                auto tissueMapIt = tissueIndex.find(*it1);

                for (size_t i = 1; i < fParser.getNWords(); i++) {
                    if (tissueMapIt == tissueIndex.end()) {
                        pair<string, double> tissueEnsemblePair(fParser.getWords()[0], strtod(fParser.getWords()[i], NULL));
                        unordered_map<string, pair<string, double>> tissues;
                        tissues.insert(pair<string, pair<string, double>>(header[i], tissueEnsemblePair));
                        tissueIndex.insert(pair<string, unordered_map<string, pair<string, double>>>(*it1, tissues));
                        tissueMapIt = tissueIndex.find(*it1);
                    } else {
                        auto tissueIt = tissueMapIt->second.find(header[i]);
                        if (tissueIt == tissueMapIt->second.end()) {
                            pair<string, double> tissueEnsemblePair(fParser.getWords()[0], strtod(fParser.getWords()[i], NULL));
                            tissueMapIt->second.insert(pair<string, pair<string, double>>(header[i], tissueEnsemblePair));
                        } else {
                            if (tissueIt->second.second < strtod(fParser.getWords()[i], NULL)) {
                                pair<string, double> tissueEnsemblePair(fParser.getWords()[0], strtod(fParser.getWords()[i], NULL));
                                tissueIt->second = pair<string, double>(fParser.getWords()[0], strtod(fParser.getWords()[i], NULL));
                            }
                        }
                    }
                }
            }
        }
    }

    if (words) free(words);
}

void FimoFactory::createCutoffIndexFromFile(std::string cutoffFileName, size_t column) {
    FileParserFactory fParser(cutoffFileName);
    while (fParser.iterate('#', " ")) {
        if (fParser.getNWords() <= column) {
            printLog(stderr, "Input cutoff file with less columns of required ", __FILE__, __LINE__, -1);
        }
        cutoffIndex.insert(pair<string, double>(fParser.getWords()[0], strtod(fParser.getWords()[column], NULL)));
    }
}

void FimoFactory::parseFimoOutput(std::string fimoOuputName, std::string tissueCode, unsigned long int snpPos) {
    FileParserFactory fParser(fimoOuputName);
    unordered_map<string, set < Fimo *, PointerCompare>> fimoMap;
    unordered_map <string, set < Fimo *, PointerCompare>>::iterator fimoIt;
    unordered_map<string, set < string>> ensemblMap;
    unordered_map<string, set < string>>::iterator ensemblMapIt;
    pair < unordered_map<string, set < string>>::iterator, bool> res;
    char **words = NULL;
    size_t nWords = 0;
    size_t wordsSize = 0;

    while (fParser.iterate('#', "\t")) {
        if (fParser.getNWords() != 8) {
            printLog(stderr, "FIMO output file with a wrong format ", __FILE__, __LINE__, -1);
        }

        res = ensemblMap.insert(pair<string, set < string >> (fParser.getWords()[1], set < string>()));
        ensemblMapIt = res.first;

        nWords = strsep_ptr(&words, &wordsSize, fParser.getWords()[0], ";");

        for (size_t i = 0; i < nWords; i++) {
            pair<string, double> motifExpression = getTissueValue(words[i], tissueCode);
            if (fabs(motifExpression.second) >= 1.0) {
                double pValue = strtod(fParser.getWords()[6], NULL);
                if (pValue < getCutoffValue(words[i])) {

                    fimoIt = fimoMap.find(fParser.getWords()[1]);
                    Fimo *f = new Fimo();
                    f->setMotif(words[i]);
                    f->setId(fParser.getWords()[1]);
                    f->setStart(static_cast<unsigned long int> (atoi(fParser.getWords()[2]) - 1));
                    f->setEnd(static_cast<unsigned long int> (atoi(fParser.getWords()[3]) - 1));
                    f->setStrand(fParser.getWords()[4][0]);
                    f->setScore(strtod(fParser.getWords()[5], NULL));
                    f->setValue(pValue);
                    f->setSeq(fParser.getWords()[7]);
                    f->setExpression(motifExpression.second);
                    f->setExpEnsembl(motifExpression.first);

                    if (fimoIt == fimoMap.end()) {
                        vector<double> v;
                        if (f->getStart() <= snpPos && f->getEnd() >= snpPos) {
                            v.push_back(motifExpression.second);
                            v.push_back(10E-10);
                        } else {
                            v.push_back(10E-10);
                            v.push_back(motifExpression.second);
                            ensemblMapIt->second.insert(motifExpression.first);
                        }
                        snpIDContainer.insert(pair<string, vector<double>>(fParser.getWords()[1], v));

                        set<Fimo *, PointerCompare> s;
                        s.insert(move(f));
                        fimoMap.insert(pair<string, set < Fimo *, PointerCompare >> (fParser.getWords()[1], s));
                        f = NULL;
                    } else {
                        pair < set < Fimo *, PointerCompare>::iterator, bool> iter;

                        iter = fimoIt->second.insert(move(f));
                        if (iter.second) {
                            f = *iter.first;
                            auto snpMapIt = snpIDContainer.find(fParser.getWords()[1]);

                            if (f->getStart() <= snpPos && f->getEnd() >= snpPos) {
                                if (snpMapIt->second[0] < motifExpression.second) {
                                    snpMapIt->second[0] = motifExpression.second;
                                }
                            } else {
                                if (ensemblMapIt->second.find(motifExpression.first) == ensemblMapIt->second.end()) {
                                    ensemblMapIt->second.insert(motifExpression.first);
                                    snpMapIt->second[1] += motifExpression.second;
                                }
                            }
                            f = NULL;
                        }
                    }
                    if (f) delete f;
                }
            }
        }
    }

    if (words) free(words);

    for (auto it = fimoMap.begin(); it != fimoMap.end(); ++it) {
        for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
            delete (*it1);
        }
    }
}

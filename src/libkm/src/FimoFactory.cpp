/* 
 * File:   FimoFactory.cpp
 * Author: veraalva
 * 
 * Created on March 4, 2016, 3:48 PM
 */

#include <iostream>
#include <memory>
#include <fstream>
#include <cstdlib>
#include <string>
#include <locale>
#include <unordered_map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>

#include "Global.h"
#include "TimeUtils.h"
#include "FileParserFactory.h"
#include "FimoFactory.h"
#include "cstring.h"

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
    vector<string> words;
    FileParserFactory fParser;
    unordered_map<string, set < string>> tFNameReverseMap;
    unordered_map<string, set < string>>::iterator tFNameReverseMapIt;
    vector<string> header;

    try {
        fParser.setFileToParse(pwm_EnsembleID);
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() < 2) {
                cerr << "Input pwm_tFName file with a wrong format " << endl;
                exit(-1);
            }
            cstring::split(fParser.getWords()[fParser.getWords().size() - 1], ";", words);
            for (auto wIt = words.begin(); wIt != words.end(); ++wIt) {
                tFNameReverseMapIt = tFNameReverseMap.find(*wIt);
                if (tFNameReverseMapIt == tFNameReverseMap.end()) {
                    set<string> s;
                    s.insert(fParser.getWords()[0]);
                    tFNameReverseMap.insert(make_pair(*wIt, s));
                } else {
                    tFNameReverseMapIt->second.insert(fParser.getWords()[0]);
                }
            }
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file: " << pwm_EnsembleID << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file: " << pwm_EnsembleID << endl;
        exit(-1);
    }

    try {
        fParser.setFileToParse(tissue_file);
        fParser.iterate("#", "\t");
        header.resize(fParser.getWords().size());
        std::copy(fParser.getWords().begin(), fParser.getWords().end(), header.begin());
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() != header.size()) {
                cerr << "Tissue file with a wrong format " << endl;
                exit(-1);
            }
            tFNameReverseMapIt = tFNameReverseMap.find(fParser.getWords()[0]);
            if (tFNameReverseMapIt != tFNameReverseMap.end()) {
                for (auto it1 = tFNameReverseMapIt->second.begin(); it1 != tFNameReverseMapIt->second.end(); ++it1) {
                    auto tissueMapIt = tissueIndex.find(*it1);

                    for (size_t i = 1; i < fParser.getWords().size(); i++) {
                        double value = std::log2(atof((fParser.getWords()[i]).c_str()) + 1);
                        if (!std::isinf(value) && !std::isnan(value)) {
                            if (tissueMapIt == tissueIndex.end()) {
                                pair<string, double> tissueEnsemblePair(fParser.getWords()[0], value);
                                unordered_map<string, pair<string, double>> tissues;
                                tissues.insert(pair<string, pair<string, double>>(header[i], tissueEnsemblePair));
                                tissueIndex.insert(pair<string, unordered_map<string, pair<string, double>>>(*it1, tissues));
                                tissueMapIt = tissueIndex.find(*it1);
                            } else {
                                auto tissueIt = tissueMapIt->second.find(header[i]);
                                if (tissueIt == tissueMapIt->second.end()) {
                                    pair<string, double> tissueEnsemblePair(fParser.getWords()[0], value);
                                    tissueMapIt->second.insert(make_pair(header[i], tissueEnsemblePair));
                                } else {
                                    if (tissueIt->second.second < value) {
                                        tissueIt->second = make_pair(fParser.getWords()[0], value);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file: " << tissue_file << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file: " << tissue_file << endl;
        exit(-1);
    }
}

void FimoFactory::createCutoffIndexFromFile(std::string cutoffFileName, size_t column) {
    FileParserFactory fParser;
    try {
        fParser.setFileToParse(cutoffFileName);
        while (fParser.iterate("#", " ")) {
            if (fParser.getWords().size() <= column) {
                cerr << "Input cutoff file with less columns of required " << endl;
                exit(-1);
            }
            cutoffIndex.insert(make_pair(fParser.getWords()[0], atof((fParser.getWords()[column]).c_str())));
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file: " << cutoffFileName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file: " << cutoffFileName << endl;
        exit(-1);
    }
}

void FimoFactory::parseFimoOutput(std::string fimoOuputName, std::string tissueCode, unsigned long int snpPos) {
    FileParserFactory fParser;
    unordered_map<string, set < shared_ptr<Fimo>, PointerCompare>> fimoMap;
    unordered_map <string, set < shared_ptr<Fimo>, PointerCompare>>::iterator fimoIt;
    unordered_map<string, set < string>> ensemblMap;
    unordered_map<string, set < string>>::iterator ensemblMapIt;
    pair < unordered_map<string, set < string>>::iterator, bool> res;
    vector<string> words;

    try {
        fParser.setFileToParse(fimoOuputName);
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() != 8) {
                cerr << "FIMO output file with a wrong format " << endl;
                exit(-1);
            }

            res = ensemblMap.insert(make_pair(fParser.getWords()[1], set < string>()));
            ensemblMapIt = res.first;

            cstring::split(fParser.getWords()[0], ";", words);

            for (auto wIt = words.begin(); wIt != words.end(); ++wIt) {
                pair<string, double> motifExpression = getTissueValue(*wIt, tissueCode);
                if (fabs(motifExpression.second) >= 1.0) {
                    double pValue = atof((fParser.getWords()[6]).c_str());
                    if (pValue < getCutoffValue(*wIt)) {

                        fimoIt = fimoMap.find(fParser.getWords()[1]);
                        shared_ptr<Fimo> f = make_shared<Fimo>();
                        f->setMotif(*wIt);
                        f->setId(fParser.getWords()[1]);
                        f->setStart(static_cast<unsigned long int> (atoi((fParser.getWords()[2]).c_str()) - 1));
                        f->setEnd(static_cast<unsigned long int> (atoi((fParser.getWords()[3]).c_str()) - 1));
                        f->setStrand(fParser.getWords()[4].at(0));
                        f->setScore(atof((fParser.getWords()[5]).c_str()));
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
                            snpIDContainer.insert(make_pair(fParser.getWords()[1], v));

                            set<shared_ptr<Fimo>, PointerCompare> s;
                            s.insert(f);
                            fimoMap.insert(make_pair(fParser.getWords()[1], s));
                        } else {
                            pair < set < shared_ptr<Fimo>, PointerCompare>::iterator, bool> iter;

                            iter = fimoIt->second.insert(f);
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
                            }
                        }
                    }
                }
            }
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file:" << fimoOuputName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file:" << fimoOuputName << endl;
        exit(-1);
    }
}

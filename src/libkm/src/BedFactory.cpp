/* 
 * File:   BedFactory.cpp
 * Author: veraalva
 * 
 * Created on February 11, 2016, 3:52 PM
 */

#include <dirent.h>

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <random>
#include <chrono>

#include "TimeUtils.h"
#include "Global.h"
#include "FileParserFactory.h"
#include "cstring.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"

using namespace std;
using namespace parsers;
using namespace sequence;
using namespace peak;
using namespace kmers;

Peak::Peak() {
    this->start = 0;
    this->end = 0;
    this->GCCount = 0;
    this->NCount = 0;
    this->NRCount = 0;
    this->NPercent = 0.0;
}

Peak::~Peak() {
}

void Peak::calculateContent() {
    for (auto it = seq.begin(); it != seq.end(); ++it) {
        if (*it == 'G' || *it == 'C') this->GCCount++;
        else if (*it == 'N') this->NCount++;
    }
    this->NRCount = this->getLength() - this->NCount;
    this->NPercent = double(this->NCount) / double(this->getLength());
}

pair<int, int> Peak::getGCNcontentBin() {
    int n_binID = 0;
    double Ncontent = static_cast<double> (this->NCount) / static_cast<double> (this->getLength());
    double GCcontent = static_cast<double> (this->GCCount) / static_cast<double> (this->getLength() - this->NCount);
    int gc_binID = static_cast<int> (GCcontent / Global::instance()->getBin1());
    if (Ncontent != 0) {
        n_binID = static_cast<int> (Ncontent / Global::instance()->getBin2()) + 1;
    }
    return pair<int, int>(gc_binID, n_binID);
}

std::string Peak::randomizePeakSeq() {
    string::iterator it1;
    for (auto it = seq.begin(); it != seq.end(); ++it) {
        if (*it != 'N' && *it != 'n') {
            for (it1 = it; it1 != seq.end() && *it1 != 'N' && *it1 != 'n'; ++it1);
            if (it1 != it) {
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::shuffle(it, it1, std::default_random_engine(seed));
            }
        }
    }
    return seq;
}

BedFactory::BedFactory() {
    this->NCount = 0;
    this->GCCount = 0;
}

BedFactory::~BedFactory() {
}

void BedFactory::createPeaksFromBedFile(FastaFactory& chrFactory, std::string bedFileName, double maxNPercent, KmersFactory &kmersFactory) {
    bool process = true;
    FileParserFactory fParser;
    std::shared_ptr<Seq> f;

    try {
        f = chrFactory.getFirstSequence();
    } catch (exceptions::NotFoundException) {
        cerr << "Not chromosome sequences loaded" << endl;
        exit(-1);
    }

    try {
        fParser.setFileToParse(bedFileName);
        kmersFactory.clearKmerPeakData();
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() < 3) {
                cerr << "Bed file with a wrong format " << endl;
                exit(-1);
            }
            if (f->getId().compare(fParser.getWords()[0]) != 0) {
                try {
                    f = chrFactory.getSequenceFromID(fParser.getWords()[0]);
                    process = true;
                } catch (exceptions::NotFoundException) {
                    process = false;
                }
            }
            if (process) {
                if (static_cast<unsigned long int> (atoi((fParser.getWords()[2]).c_str())) <= f->getLength()) {
                    shared_ptr<Peak> p = make_shared<Peak>();
                    p->setChr(fParser.getWords()[0]);
                    
                    // For coordinates base see:
                    // http://genome.ucsc.edu/FAQ/FAQformat.html#format12
                    p->setStart(atoi((fParser.getWords()[1]).c_str()));
                    p->setEnd(atoi((fParser.getWords()[2]).c_str()) - 1); 
                    try {
                        p->setSeq(f->getSubStr(p->getStart(), p->getLength()));
                        p->calculateContent();
                        if (p->getNRCount() >= 50 && p->getNPercent() < maxNPercent) {
                            pair<int, int> k = p->getGCNcontentBin();

                            if (this->GCNcontentBin.find(f->getId()) == this->GCNcontentBin.end()) {
                                std::map<int, std::map < std::pair<int, int>, int>> wm;
                                std::map<std::pair<int, int>, int> m;
                                m.insert(make_pair(k, 1));
                                wm.insert(make_pair(p->getLength(), m));
                                this->GCNcontentBin.insert(make_pair(f->getId(), wm));
                            } else {
                                std::map<int, std::map < std::pair<int, int>, int>> *wm = &this->GCNcontentBin.find(f->getId())->second;
                                if (wm->find(p->getLength()) == wm->end()) {
                                    std::map<std::pair<int, int>, int> m;
                                    m.insert(make_pair(k, 1));
                                    wm->insert(make_pair(p->getLength(), m));
                                } else {
                                    std::map<std::pair<int, int>, int> *m = &wm->find(p->getLength())->second;
                                    if (m->find(k) == m->end()) {
                                        m->insert(make_pair(k, 1));
                                    } else {
                                        (*m)[k]++;
                                    }
                                }
                            }
                            kmersFactory.scanSequences(p->getSeq(), false);
                            this->peaks.push_back(p);
                        }
                    } catch (std::out_of_range) {
                        cerr << "Out of range coordinates for sequence. Ignoring peak" << endl;
                    }
                }
            }
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error opening file: " << bedFileName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error retrieving sequences" << endl;
        exit(-1);
    }

    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> " << "Total number of non N for peaks: " << kmersFactory.getTotalNRNTPeak() << endl;
    }
}

void BedFactory::plus_ocur(char nt) {
    switch (nt) {
        case 'N': this->NCount++;
            break;
        case 'G': this->GCCount++;
            break;
        case 'C': this->GCCount++;
            break;
    }
}

void BedFactory::minus_ocur(char nt) {
    switch (nt) {
        case 'N': this->NCount--;
            break;
        case 'G': this->GCCount--;
            break;
        case 'C': this->GCCount--;
            break;
    }
}

void BedFactory::generatingControlsFromChromosomes(FastaFactory &chrFactory, unsigned long int hit_Num, KmersFactory & kmersFactory) {
    string seq;
    clock_t begin = clock();
    unsigned long int i, l, window_len, total_controlNum;
    int gc_binID, n_binID;
    double Ncontent, GCcontent;
    int peakNum;
    std::mt19937 rng;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    rng.seed(seed);

    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> Generating controls from chromosomes" << endl;
    }

    kmersFactory.clearKmerControlData();

    for (auto chr_it = this->GCNcontentBin.begin(); chr_it != this->GCNcontentBin.end(); ++chr_it) {
        try {
            std::shared_ptr<Seq> f = chrFactory.getSequenceFromID(chr_it->first);
            std::map<int, std::map < std::pair<int, int>, int>> wm = chr_it->second;

            if (Global::instance()->isInfo()) {
                cout << "\tINFO ==> \t" << "Chromosome number: " << distance(this->GCNcontentBin.begin(), chr_it) + 1 << ", name: " << f->getId() << " and size: " << f->getLength() << endl;
            }

            for (auto window_it = wm.begin(); window_it != wm.end(); ++window_it) {
                std::map<std::pair<int, int>, std::vector < shared_ptr<Peak> >> GCNcontentControls;
                window_len = window_it->first;
                std::map<std::pair<int, int>, int> contentBin = window_it->second;

                this->GCCount = this->NCount = 0;

                if (Global::instance()->isInfo()) {
                    begin = clock();
                    cout << "\tINFO ==> \t\t" << "Window "
                            << distance(wm.begin(), window_it) + 1 << "/" << wm.size()
                            << ", with length: " << window_len << " and "
                            << contentBin.size() << " bins" << endl;
                }
                for (unsigned long int j = 0; j < window_len; j++) {
                    this->plus_ocur(f->getSeq()[j]);
                }
                if (this->NCount < window_len - 7) {
                    n_binID = 0;
                    Ncontent = static_cast<double> (this->NCount) / static_cast<double> (window_len);
                    GCcontent = static_cast<double> (this->GCCount) / static_cast<double> (window_len - this->NCount);
                    gc_binID = static_cast<int> (GCcontent / Global::instance()->getBin1());
                    if (Ncontent != 0) {
                        n_binID = static_cast<int> (Ncontent / Global::instance()->getBin2()) + 1;
                    }
                    pair<int, int> k(gc_binID, n_binID);

                    if (contentBin.find(k) != contentBin.end()) {
                        shared_ptr<Peak> p = make_shared<Peak>();
                        p->setChr(f->getId());
                        p->setStart(0);
                        p->setEnd(window_len - 1);
                        p->setGCCount(this->GCCount);
                        p->setNCount(this->NCount);
                        if (GCNcontentControls.find(k) == GCNcontentControls.end()) {
                            std::vector <shared_ptr < Peak>> v;
                            v.push_back(p);
                            GCNcontentControls.insert(make_pair(k, v));
                        } else {
                            GCNcontentControls[k].push_back(p);
                        }
                    }
                }

                for (unsigned long int j = 1; j < f->getLength() - window_len + 1; j++) {
                    this->minus_ocur(f->getSeq()[j - 1]);
                    this->plus_ocur(f->getSeq()[j + window_len - 1]);
                    if (this->NCount < window_len - 7) {
                        n_binID = 0;
                        Ncontent = static_cast<double> (this->NCount) / static_cast<double> (window_len);
                        GCcontent = static_cast<double> (this->GCCount) / static_cast<double> (window_len - this->NCount);
                        gc_binID = static_cast<int> (GCcontent / Global::instance()->getBin1());
                        if (Ncontent != 0) {
                            n_binID = static_cast<int> (Ncontent / Global::instance()->getBin2()) + 1;
                        }
                        pair<int, int> k(gc_binID, n_binID);

                        if (contentBin.find(k) != contentBin.end()) {
                            shared_ptr<Peak> p = make_shared<Peak>();
                            p->setChr(f->getId());
                            p->setStart(j);
                            p->setEnd(j + window_len - 1);
                            p->setGCCount(this->GCCount);
                            p->setNCount(this->NCount);
                            if (GCNcontentControls.find(k) == GCNcontentControls.end()) {
                                std::vector <shared_ptr < Peak>> v;
                                v.push_back(p);
                                GCNcontentControls.insert(make_pair(k, v));
                            } else {
                                GCNcontentControls[k].push_back(p);
                            }
                        }
                    }
                }

                if (Global::instance()->isInfo()) {
                    cout << "\tINFO ==> \t\t\t" << GCNcontentControls.size() << " control peaks bin generated in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;
                }

                /*
                 * Generating hit_Num of random sequences for each peak that is not in controls
                 */
                l = 0;
                for (auto bin_it = contentBin.begin(); bin_it != contentBin.end(); ++bin_it) {
                    std::pair<int, int> k = bin_it->first;

                    if (GCNcontentControls.find(k) == GCNcontentControls.end()) {
                        if (Global::instance()->isInfo()) {
                            cout << "\tINFO ==> \t\t\tNot hit for bin: " << k.first << "-" << k.second << endl;
                        }
                        for (auto peak_it = this->peaks.begin(); peak_it != this->peaks.end(); ++peak_it) {
                            shared_ptr<Peak> p = *peak_it;
                            if (p->getGCNcontentBin() == k) {
                                for (i = 0; i < hit_Num; i++) {
                                    kmersFactory.scanSequences(cstring::shuffle(p->getSeq()), true);
                                }
                            }
                        }
                        l++;
                    }
                }

                if (Global::instance()->isInfo()) {
                    if (l != 0) {
                        cout << "\tINFO ==> \t\t\t" << l << " content bin with none hits in controls" << endl;
                    }
                }

                for (auto control_it = GCNcontentControls.begin(); control_it != GCNcontentControls.end(); ++control_it) {
                    pair<int, int> k = control_it->first;
                    std::vector <shared_ptr < Peak>> peaksVector = control_it->second;
                    peakNum = contentBin[k];
                    total_controlNum = hit_Num * peakNum;

                    if (peaksVector.size() < total_controlNum) {
                        if (peaksVector.size() > 0) {
                            if (Global::instance()->isInfo()) {
                                cout << "\tINFO ==> \t\t\t" << total_controlNum << " sequences shuffled for bin " << k.first << "-" << k.second << endl;
                            }
                            try {
                                seq = f->getSubStr(peaksVector[0]->getStart(), window_len);
                                for (i = 0; i < total_controlNum; i++) {
                                    kmersFactory.scanSequences(cstring::shuffle(seq), true);
                                }
                            } catch (std::out_of_range) {
                                cerr << "Out of range sequence coordinates. Ignoring peak" << endl;
                            }
                        } else {
                            /*
                             * This case should not occur, but let's be ready for it printing the flags
                             */
                            cout << "ERROR: No controls. " << k.first << "-" << k.second << " " << peaksVector.size() << " " << total_controlNum << endl;
                        }
                    } else if (peaksVector.size() == total_controlNum) {
                        if (Global::instance()->isInfo()) {
                            cout << "\tINFO ==> \t\t\t" << total_controlNum << " sequences selected for bin " << k.first << "-" << k.second << endl;
                        }
                        for (i = 0; i < total_controlNum; i++) {
                            try {
                                seq = f->getSubStr(peaksVector[i]->getStart(), window_len);
                                if (peaksVector[i]->getStart() + window_len <= f->getLength()) {
                                    kmersFactory.scanSequences(seq, true);
                                }
                            } catch (std::out_of_range) {
                                cerr << "Out of range sequence coordinates. Ignoring peak" << endl;
                            }
                        }
                    } else {
                        if (Global::instance()->isInfo()) {
                            cout << "\tINFO ==> \t\t\t" << total_controlNum << " sequences randomly selected for bin " << k.first << "-" << k.second << endl;
                        }
                        i = 0;
                        vector<uint32_t> index;
                        std::uniform_int_distribution<uint32_t> uint_dist(0,peaksVector.size() - 1);
                        while (1) {
                            uint32_t randomIndex = uint_dist(rng);
                            if (find(index.begin(), index.end(), randomIndex) == index.end()) {
                                index.push_back(randomIndex);
                                if (peaksVector[randomIndex]->getStart() + window_len <= f->getLength()) {
                                    try {
                                        kmersFactory.scanSequences(f->getSubStr(peaksVector[randomIndex]->getStart(), window_len), true);
                                        i++;
                                        if (i == total_controlNum) {
                                            break;
                                        }
                                    } catch (std::out_of_range) {
                                        cerr << "Out of range sequence coordinates. Ignoring peak" << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } catch (exceptions::NotFoundException) {
            cerr << "Not sequence for ID: " << chr_it->first << ". Ignoring entry" << endl;
        }
    }

    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> " << "\tTotal number of non N for controls: " << kmersFactory.getTotalNRNTControl() << endl;
    }
}

void BedFactory::generatingControlsFromShufflingPeaks(unsigned long int hit_Num, kmers::KmersFactory & kmersFactory) {
    unsigned long int i;
    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> Generating controls shuffling input peaks" << endl;
    }

    kmersFactory.clearKmerControlData();

    for (auto peaks_it = this->peaks.begin(); peaks_it != this->peaks.end(); ++peaks_it) {
        shared_ptr<Peak> p = *peaks_it;
        for (i = 0; i < hit_Num; i++) {
            string seq = p->randomizePeakSeq();
            kmersFactory.scanSequences(seq, true);
        }
    }

    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> " << "\tTotal number of non N for controls: " << kmersFactory.getTotalNRNTControl() << endl;
    }
}

void BedFactory::readControlsFromFile(std::string controlFileName, FastaFactory &chrFactory, kmers::KmersFactory & kmersFactory) {
    FileParserFactory fParser;
    bool process = true;
    std::shared_ptr<Seq> f;

    try {
        f = chrFactory.getFirstSequence();
    } catch (exceptions::NotFoundException) {
        cerr << "Not chromosome sequences loaded" << endl;
        exit(-1);
    }

    try {
        fParser.setFileToParse(controlFileName);
        kmersFactory.clearKmerControlData();
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() != 3) {
                cerr << "Bed file with a wrong format " << endl;
                exit(-1);
            }
            if (f->getId().compare(fParser.getWords()[0]) != 0) {
                try {
                    f = chrFactory.getSequenceFromID(fParser.getWords()[0]);
                    process = true;
                } catch (exceptions::NotFoundException) {
                    process = false;
                }
            }
            if (process) {
                try {
                    kmersFactory.scanSequences(f->getSubStr(atoi((fParser.getWords()[1]).c_str()), atoi((fParser.getWords()[2]).c_str())), true);
                } catch (std::out_of_range) {
                    cerr << "Out of range sequence coordinates. Ignoring control sequence";
                }
            }
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file: " << controlFileName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file: " << controlFileName << endl;
        exit(-1);
    }
}


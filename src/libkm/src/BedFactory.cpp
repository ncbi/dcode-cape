/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BedFactory.cpp
 * Author: veraalva
 * 
 * Created on February 11, 2016, 3:52 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <stdbool.h>
#include <dirent.h>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "TimeUtils.h"
#include "Global.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"

using namespace std;
using namespace fasta;
using namespace peak;
using namespace kmers;

Peak::Peak() {
    this->start = 0;
    this->end = 0;
    this->seq = NULL;
    this->GCCount = 0;
    this->NCount = 0;
    this->NRCount = 0;
    this->NPercent = 0.0;
}

Peak::~Peak() {
    if (this->seq) free(this->seq);
}

void Peak::CalculateContent() {
    char *s = this->seq;
    while (*s != '\0') {
        if (*s == 'G' || *s == 'C') this->GCCount++;
        else if (*s == 'N') this->NCount++;
        s++;
    }
    this->NRCount = this->GetLength() - this->NCount;
    this->NPercent = double(this->NCount) / double(this->GetLength());
}

pair<int, int> Peak::GetGCNcontentBin() {
    int n_binID = 0;
    double Ncontent = static_cast<double> (this->NCount) / static_cast<double> (this->GetLength());
    double GCcontent = static_cast<double> (this->GCCount) / static_cast<double> (this->GetLength() - this->NCount);
    int gc_binID = int(GCcontent / Global::instance()->GetBin1());
    if (Ncontent != 0) {
        n_binID = int(Ncontent / Global::instance()->GetBin2()) + 1;
    }
    return pair<int, int>(gc_binID, n_binID);
}

char* Peak::RandomizePeakSeq() {
    char c, *s, *n;
    char *seq = strdup(this->seq);
    s = seq;
    while (*s != '\0') {
        n = strchr(s, 'N');
        if (n == s) s++;
        else {
            if (n != NULL) *n = '\0';
            int i = strlen(s) - 1;
            while (i > 0) {
                int j = rand() % i;
                c = s[i];
                s[i] = s[j];
                s[j] = c;
                i--;
            }
            if (n != NULL) {
                *n = 'N';
                s = n;
            } else {
                break;
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
    for (auto it = this->GetPeaks().begin(); it != this->GetPeaks().end(); ++it) {
        delete (*it);
    }
}

void BedFactory::CreatePeaksFromBedFile(FastaFactory& chrFactory, char *bedFileName, double maxNPercent, KmersFactory &kmersFactory) {
    FILE *bedFile = (FILE *) checkPointerError(fopen(bedFileName, "r"), "Can't open bed file", __FILE__, __LINE__, -1);

    size_t bufferSize, backupLineSize;
    char *buffer, *newLine, *str, *backupLine, *completeLine, *seq;
    char **fields = NULL;
    size_t fieldsSize = 0;
    Fasta *f = NULL;

    backupLineSize = bufferSize = 100000000;
    buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
    backupLine = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);

    kmersFactory.ClearKmerPeakData();

    *backupLine = 0;
    while (!feof(bedFile)) {
        size_t read = fread(buffer, sizeof (char), bufferSize, bedFile);
        buffer[read] = 0;
        if (feof(bedFile)) {
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
                fieldsSize = splitString(&fields, completeLine, "\t");
                if (fieldsSize != 3) {
                    printLog(stderr, "Bed file with a wrong format ", __FILE__, __LINE__, -1);
                }
                if (!f || f->GetId().compare(fields[0]) != 0) {
                    f = chrFactory.GetFastaFromID(fields[0]);
                }
                if (f) {
                    if (static_cast<unsigned long int> (atoi(fields[2])) <= f->GetLength()) {
                        Peak *p = new Peak();
                        p->SetChr(fields[0]);
                        p->SetStart(atoi(fields[1]));
                        p->SetEnd(atoi(fields[2]) - 1);
                        seq = f->GetSubStr(p->GetStart(), p->GetLength());
                        p->SetSeq(&seq);
                        p->CalculateContent();
                        if (p->GetNRCount() >= 50 && p->GetNPercent() < maxNPercent) {
                            pair<int, int> k = p->GetGCNcontentBin();

                            if (this->GCNcontentBin.find(f->GetId()) == this->GCNcontentBin.end()) {
                                std::map<int, std::map < std::pair<int, int>, int>> wm;
                                std::map<std::pair<int, int>, int> m;
                                m.insert(pair<pair<int, int>, int>(k, 1));
                                wm.insert(pair<int, std::map < std::pair<int, int>, int>>(p->GetLength(), m));
                                this->GCNcontentBin.insert(pair<string, std::map<int, std::map < std::pair<int, int>, int>>>(f->GetId(), wm));
                            } else {
                                std::map<int, std::map < std::pair<int, int>, int>> *wm = &this->GCNcontentBin.find(f->GetId())->second;
                                if (wm->find(p->GetLength()) == wm->end()) {
                                    std::map<std::pair<int, int>, int> m;
                                    m.insert(pair<pair<int, int>, int>(k, 1));
                                    wm->insert(pair<int, std::map < std::pair<int, int>, int>>(p->GetLength(), m));
                                } else {
                                    std::map<std::pair<int, int>, int> *m = &wm->find(p->GetLength())->second;
                                    if (m->find(k) == m->end()) {
                                        m->insert(pair<pair<int, int>, int>(k, 1));
                                    } else {
                                        (*m)[k]++;
                                    }
                                }
                            }
                            kmersFactory.scanSequences(p->GetSeq(), false);
                            this->peaks.push_back(std::move(p));
                        } else {
                            delete p;
                        }
                    }
                }
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

    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> " << "Total number of non N for peaks: " << kmersFactory.GetTotalNRnt_peak() << endl;
    }

    free(buffer);
    free(backupLine);
    fclose(bedFile);
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

void BedFactory::GeneratingControlsFromChromosomes(FastaFactory &chrFactory, unsigned long int hit_Num, KmersFactory &kmersFactory) {
    char *seq;
    clock_t begin = clock();
    unsigned long int i, l, window_len, total_controlNum;
    int gc_binID, n_binID;
    double Ncontent, GCcontent;
    int peakNum;

    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> Generating controls from chromosomes" << endl;
    }

    kmersFactory.ClearKmerControlData();

    for (auto chr_it = this->GCNcontentBin.begin(); chr_it != this->GCNcontentBin.end(); ++chr_it) {
        Fasta *f = chrFactory.GetFastaFromID(chr_it->first);
        std::map<int, std::map < std::pair<int, int>, int>> wm = chr_it->second;

        if (Global::instance()->isInfo()) {
            cout << "\tINFO ==> \t" << "Chromosome number: " << distance(this->GCNcontentBin.begin(), chr_it) + 1 << ", name: " << f->GetId() << " and size: " << f->GetLength() << endl;
        }

        for (auto window_it = wm.begin(); window_it != wm.end(); ++window_it) {
            std::map<std::pair<int, int>, std::vector < Peak *>> GCNcontentControls;
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
                this->plus_ocur(f->GetSeq()[j]);
            }
            if (this->NCount < window_len - Global::instance()->GetOrder() - 1) {
                n_binID = 0;
                Ncontent = static_cast<double> (this->NCount) / static_cast<double> (window_len);
                GCcontent = static_cast<double> (this->GCCount) / static_cast<double> (window_len - this->NCount);
                gc_binID = int(GCcontent / Global::instance()->GetBin1());
                if (Ncontent != 0) {
                    n_binID = int(Ncontent / Global::instance()->GetBin2()) + 1;
                }
                pair<int, int> k(gc_binID, n_binID);

                if (contentBin.find(k) != contentBin.end()) {
                    Peak *p = new Peak();
                    p->SetChr(f->GetId());
                    p->SetStart(0);
                    p->SetEnd(window_len - 1);
                    p->SetGCCount(this->GCCount);
                    p->SetNCount(this->NCount);
                    if (GCNcontentControls.find(k) == GCNcontentControls.end()) {
                        std::vector < Peak *> v;
                        v.push_back(std::move(p));
                        GCNcontentControls.insert(pair<std::pair<int, int>, std::vector < Peak *>>(k, v));
                    } else {
                        GCNcontentControls[k].push_back(std::move(p));
                    }
                }
            }

            for (unsigned long int j = 1; j < f->GetLength() - window_len + 1; j++) {
                this->minus_ocur(f->GetSeq()[j - 1]);
                this->plus_ocur(f->GetSeq()[j + window_len - 1]);
                if (this->NCount < window_len - Global::instance()->GetOrder() - 1) {
                    n_binID = 0;
                    Ncontent = static_cast<double> (this->NCount) / static_cast<double> (window_len);
                    GCcontent = static_cast<double> (this->GCCount) / static_cast<double> (window_len - this->NCount);
                    gc_binID = int(GCcontent / Global::instance()->GetBin1());
                    if (Ncontent != 0) {
                        n_binID = int(Ncontent / Global::instance()->GetBin2()) + 1;
                    }
                    pair<int, int> k(gc_binID, n_binID);

                    if (contentBin.find(k) != contentBin.end()) {
                        Peak *p = new Peak();
                        p->SetChr(f->GetId());
                        p->SetStart(j);
                        p->SetEnd(j + window_len - 1);
                        p->SetGCCount(this->GCCount);
                        p->SetNCount(this->NCount);
                        if (GCNcontentControls.find(k) == GCNcontentControls.end()) {
                            std::vector < Peak *> v;
                            v.push_back(std::move(p));
                            GCNcontentControls.insert(pair<std::pair<int, int>, std::vector < Peak *>>(k, v));
                        } else {
                            GCNcontentControls[k].push_back(std::move(p));
                        }
                    }
                }
            }

            if (Global::instance()->isInfo()) {
                cout << "\tINFO ==> \t\t\t" << GCNcontentControls.size() << " control peaks bin generated in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;
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
                        Peak *p = *peak_it;
                        if (p->GetGCNcontentBin() == k) {
                            for (i = 0; i < hit_Num; i++) {
                                char *s = shuffle(p->GetSeq());
                                kmersFactory.scanSequences(s, true);
                                free(s);
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
                std::vector < Peak *> peaksVector = control_it->second;
                peakNum = contentBin[k];
                total_controlNum = hit_Num * peakNum;

                if (peaksVector.size() < total_controlNum) {
                    if (peaksVector.size() > 0) {
                        if (Global::instance()->isInfo()) {
                            cout << "\tINFO ==> \t\t\t" << total_controlNum << " sequences shuffled for bin " << k.first << "-" << k.second << endl;
                        }
                        if (peaksVector[0]->GetStart() + window_len <= f->GetLength()) {
                            seq = f->GetSubStr(peaksVector[0]->GetStart(), window_len);
                            for (i = 0; i < total_controlNum; i++) {
                                char *s = shuffle(seq);
                                kmersFactory.scanSequences(s, true);
                                free(s);
                            }
                            free(seq);
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
                        if (peaksVector[i]->GetStart() + window_len <= f->GetLength()) {
                            seq = f->GetSubStr(peaksVector[i]->GetStart(), window_len);
                            kmersFactory.scanSequences(seq, true);
                            free(seq);
                        }
                    }
                } else {
                    if (Global::instance()->isInfo()) {
                        cout << "\tINFO ==> \t\t\t" << total_controlNum << " sequences randomly selected for bin " << k.first << "-" << k.second << endl;
                    }
                    i = 0;
                    vector<int> index;
                    while (1) {
                        int randomIndex = rand() % peaksVector.size();
                        if (find(index.begin(), index.end(), randomIndex) == index.end()) {
                            index.push_back(randomIndex);
                            if (peaksVector[randomIndex]->GetStart() + window_len <= f->GetLength()) {
                                seq = f->GetSubStr(peaksVector[randomIndex]->GetStart(), window_len);
                                kmersFactory.scanSequences(seq, true);
                                free(seq);
                                i++;
                                if (i == total_controlNum) {
                                    break;
                                }
                            }
                        }
                    }
                }
                for (auto it = peaksVector.begin(); it != peaksVector.end(); ++it) {
                    delete (*it);
                }
            }
        }
    }

    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> " << "\tTotal number of non N for controls: " << kmersFactory.GetTotalNRnt_control() << endl;
    }
}

void BedFactory::GeneratingControlsFromShufflingPeaks(unsigned long int hit_Num, kmers::KmersFactory& kmersFactory) {
    unsigned long int i;
    char *seq;
    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> Generating controls shuffling input peaks" << endl;
    }

    kmersFactory.ClearKmerControlData();

    for (auto peaks_it = this->peaks.begin(); peaks_it != this->peaks.end(); ++peaks_it) {
        Peak *p = *peaks_it;
        for (i = 0; i < hit_Num; i++) {
            seq = p->RandomizePeakSeq();
            kmersFactory.scanSequences(seq, true);
            free(seq);
        }
    }

    if (Global::instance()->isInfo()) {
        cout << "\tINFO ==> " << "\tTotal number of non N for controls: " << kmersFactory.GetTotalNRnt_control() << endl;
    }
}

void BedFactory::ReadControlsFromFile(char* controlFileName, FastaFactory &chrFactory, kmers::KmersFactory& kmersFactory) {
    FILE *bedFile = (FILE *) checkPointerError(fopen(controlFileName, "r"), "Can't open bed file", __FILE__, __LINE__, -1);

    size_t bufferSize, backupLineSize;
    char *buffer, *newLine, *str, *backupLine, *completeLine, *seq;
    char **fields = NULL;
    size_t fieldsSize = 0;
    Fasta *f = NULL;

    backupLineSize = bufferSize = 100000000;
    buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
    backupLine = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);

    kmersFactory.ClearKmerControlData();
    *backupLine = 0;
    while (!feof(bedFile)) {
        size_t read = fread(buffer, sizeof (char), bufferSize, bedFile);
        buffer[read] = 0;
        if (feof(bedFile)) {
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
                fieldsSize = splitString(&fields, completeLine, "\t");
                if (fieldsSize != 3) {
                    printLog(stderr, "Bed file with a wrong format ", __FILE__, __LINE__, -1);
                }
                if (!f || f->GetId().compare(fields[0]) != 0) {
                    f = chrFactory.GetFastaFromID(fields[0]);
                }
                if (f) {
                    if (static_cast<unsigned long int> (atoi(fields[1]) + atoi(fields[2])) < f->GetLength()) {
                        seq = f->GetSubStr(atoi(fields[1]), atoi(fields[2]));
                        if (seq) {
                            kmersFactory.scanSequences(seq, true);
                            free(seq);
                        }
                    }
                }
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
    fclose(bedFile);
}


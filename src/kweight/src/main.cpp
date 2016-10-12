/*  $Id: main.cpp February 10, 2016, 3:18 PM veraalva $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Roberto Vera Alvarez
 *
 * File Description:
 *   
 *
 */

#include <iostream>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>

#include "Global.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "cstring.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"

using namespace std;
using namespace sequence;
using namespace peak;
using namespace kmers;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

void print_usage(char *program_name, int exit_code) {
    cerr << "\n********************************************************************************\n";
    cerr << "This option will read the enhancers coordinates and generate the controls by shuffling the enhancers sequences\n";
    cerr << "kweight -b enhancers_coordinates.bed -m ./directory/hg19/chromosomes/hg19_masked.fa.bin -p weight_file_kmers.txt\n\n";
    cerr << "This option will read the enhancers coordinates and generate the controls from the chromosomes\n";
    cerr << "kweight -g -b enhancers_coordinates.bed -m ./directory/hg19/chromosomes/hg19_masked.fa.bin --p weight_file_kmers.txt\n\n";
    cerr << "This option will read the enhancers coordinates and use the controls provided by the users\n";
    cerr << "kweight -b enhancers_coordinates.bed -m ./directory/hg19/hg19_masked.fa.bin -c users_control.txt -p weight_file_kmers.txt\n";
    cerr << "Control file format:\n";
    cerr << "chr<tab>start_position<tab>length\n\n";
    cerr << "\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr << "-v    Print info\n";
    cerr << "-h    Display this usage information.\n";
    cerr << "-b    Coordination of enhancers or DHS peaks.\n";
    cerr << "-m    Chromosomes masked binary fasta file. Created with formatFasta.\n";
    cerr << "-p    Output file with the p-value for all kmers\n";

    cerr << "-o    Order (default: \"4,6,8,10,12\"). No spaces allowed.\n";
    cerr << "-g    Generate control from chromosomes. Default: NO.\n";
    cerr << "-n    Number of controls per peak (default: 10, if -g set default: 3).\n";
    cerr << "-c    Bed file with the control coordinates. This option will use your own control for the calculations.\n";

    cerr << "********************************************************************************\n";
    cerr << "\n            Shan Li (e-mail: lis11@ncbi.nlm.nih.gov)\n";
    cerr << "            Roberto Vera Alvarez (e-mail: veraalva@ncbi.nlm.nih.gov)\n\n";
    cerr << "********************************************************************************\n";
    exit(exit_code);
}

int main(int argc, char** argv) {
    clock_t begin = clock();
    clock_t start = clock();
    unsigned long int hitNum = 0;
    string chrsBinFileName;
    string bedFileName;
    string controlFileName;
    string poutputFileName;
    string order = "4,6,8,10,12";
    vector<string> orders;
    FastaFactory chrFactory;
    BedFactory bedFactory;
    bool genCtrl = false;

    for (int i = 1; i < argc; i++) {
        string option(argv[i]);
        if (option.compare(0, 1, "-") == 0 && option.compare(1, 1, "-") != 0 && option.size() == 2) {
            if (option.compare(1, 1, "h") == 0) {
                print_usage(argv[0], 0);
            } else if (option.compare(1, 1, "v") == 0) {
                Global::instance()->setVerbose(1);
            } else if (option.compare(1, 1, "d") == 0) {
                Global::instance()->setVerbose(3);
            } else if (option.compare(1, 1, "b") == 0) {
                i++;
                if (i < argc) {
                    bedFileName = argv[i];
                    if (bedFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option b require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option b require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "m") == 0) {
                i++;
                if (i < argc) {
                    chrsBinFileName = argv[i];
                    if (chrsBinFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option m require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option m require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "p") == 0) {
                i++;
                if (i < argc) {
                    poutputFileName = argv[i];
                    if (poutputFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option p require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option p require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "o") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option o require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    order = argv[i];
                } else {
                    cerr << "Option o require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "g") == 0) {
                genCtrl = true;
            } else if (option.compare(1, 1, "n") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option n require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    hitNum = atoi(argv[i]);
                } else {
                    cerr << "Option n require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "c") == 0) {
                i++;
                if (i < argc) {
                    controlFileName = argv[i];
                    if (controlFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option c require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option c require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else {
                cerr << "Unsupported option: " << option << endl;
                print_usage(argv[0], -1);
            }
        } else {
            cerr << "Unsupported option: " << option << endl;
            print_usage(argv[0], -1);
        }
    }

    if (genCtrl) {
        if (hitNum == 0) hitNum = 3;
    } else {
        if (hitNum == 0) hitNum = 10;
    }

    if (poutputFileName.empty()) {
        cerr << "\nKmer contingency table output is required. See -k option" << endl;
        print_usage(argv[0], -1);
    }

    if (chrsBinFileName.empty()) {
        cerr << "\nCan't open chromosomes masked binary fasta file. See -m option" << endl;
        print_usage(argv[0], -1);
    }

    if (bedFileName.empty()) {
        cerr << "\nBed file is required. See -b option" << endl;
        print_usage(argv[0], -1);
    }
    
    if (order.empty()) {
        cerr << "\nOrder is required. See -o option" << endl;
        print_usage(argv[0], -1);
    }
    
    if (order.find(' ', 0) != std::string::npos) {
        cerr << "\nOrder can not have spaces. See -o option" << endl;
        print_usage(argv[0], -1);
    }

    Global::instance()->setBin1(0.005);
    Global::instance()->setBin2(0.01);
    
    cstring::split(order, ",", orders);
    for(auto it = orders.begin(); it != orders.end(); ++it){
        Global::instance()->getOrders().insert(atoi((*it).c_str()));
    }

    begin = clock();
    cout << "Creating kmers genomewide" << endl;
    KmersFactory kmersFactory;
    kmersFactory.createGenomeWideKmers();
    cout << kmersFactory.getKmersGenome().size() << " kmers created in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading chromosome sequences from binary file" << endl;
    chrFactory.parseFastaFile(chrsBinFileName, true);
    cout << chrFactory.getSequenceContainter().size() << " chromosomes loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading peaks from file" << endl;
    bedFactory.createPeaksFromBedFile(chrFactory, bedFileName, 0.7, kmersFactory);
    cout << bedFactory.getPeaks().size() << " peaks loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;
    cout << bedFactory.getGCNcontentBin().size() << " chromosomes to analyze" << endl;

    begin = clock();
    if (controlFileName.empty()) {
        if (genCtrl) {
            cout << "Generating controls from Chromosomes" << endl;
            bedFactory.generatingControlsFromChromosomes(chrFactory, hitNum, kmersFactory);
        } else {
            cout << "Generating controls from shuffling peaks" << endl;
            bedFactory.generatingControlsFromShufflingPeaks(hitNum, kmersFactory);
        }
    } else {
        cout << "Reading controls from file" << endl;
        bedFactory.readControlsFromFile(controlFileName, chrFactory, kmersFactory);
    }
    cout << "Controls generated in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Generating kmers from peaks and controls" << endl;
    kmersFactory.buildKmers();
    cout << kmersFactory.getKmers().size() << " kmers generated in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    kmersFactory.writeKmersToFile(poutputFileName);

    delete Global::instance();
    cout << "Total elapse time: " << TimeUtils::instance()->getTimeMinFrom(start) << " minutes" << endl;
    delete TimeUtils::instance();
    return 0;
}


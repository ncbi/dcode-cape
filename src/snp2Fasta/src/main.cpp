/* 
 * File:   main.cpp
 * Author: veraalva
 *
 * Created on March 10, 2016, 11:46 AM
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>

#include "berror.h"
#include "bmemory.h"
#include "svm.h"
#include "Global.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "FimoFactory.h"
#include "BedFactory.h"
#include "TFBSFactory.h"
#include "SNPFactory.h"
#include "SVMPredict.h"

using namespace std;
using namespace sequence;
using namespace peak;
using namespace kmers;
using namespace snp;
using namespace svm;
using namespace fimo;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

void print_usage(char *program_name, int exit_code) {
    cerr << "\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr << "-v    Print info\n";
    cerr << "-h    Display this usage information.\n";
    cerr << "-i    Input file with SNP coordinates.\n";
    cerr << "-c    Chromosomes binary fasta file. Created with formatFasta.\n";
    cerr << "-o    Fasta file with the sequences\n";
    cerr << "-l    Sequence length to be added before and after the SNP position\n\n";

    cerr << "********************************************************************************\n";
    cerr << "\n            Shan Li (e-mail: lis11@ncbi.nlm.nih.gov)\n";
    cerr << "            Roberto Vera Alvarez (e-mail: veraalva@ncbi.nlm.nih.gov)\n\n";
    cerr << "********************************************************************************\n";
    exit(exit_code);
}

/*
 * 
 */
int main(int argc, char** argv) {
    clock_t begin = clock();
    clock_t start = clock();
    string chrsBinFileName;
    string inName;
    string outName;
    unsigned long int length = 0;
    FastaFactory chrFactory;
    SNPFactory snpFactory;

    for (int i = 1; i < argc; i++) {
        string option(argv[i]);
        if (option.compare(0, 1, "-") == 0 && option.compare(1, 1, "-") != 0 && option.size() == 2) {
            if (option.compare(1, 1, "h") == 0) {
                print_usage(argv[0], 0);
            } else if (option.compare(1, 1, "v") == 0) {
                Global::instance()->setVerbose(1);
            } else if (option.compare(1, 1, "d") == 0) {
                Global::instance()->setVerbose(3);
            } else if (option.compare(1, 1, "i") == 0) {
                i++;
                if (i < argc) {
                    inName = argv[i];
                    if (inName.compare(0, 1, "-") == 0) {
                        cerr << "Option i require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option i require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "c") == 0) {
                i++;
                if (i < argc) {
                    chrsBinFileName = argv[i];
                    if (chrsBinFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option c require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option c require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "o") == 0) {
                i++;
                if (i < argc) {
                    outName = argv[i];
                    if (outName.compare(0, 1, "-") == 0) {
                        cerr << "Option o require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option o require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "l") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option l require a numeric argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    length = static_cast<unsigned long int> (atoi(argv[i]));
                } else {
                    cerr << "Option l require an argument" << endl;
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

    if (inName.empty()) {
        cerr << "\nInput SNP coordinate file is required. See -i option" << endl;
        print_usage(argv[0], -1);
    }

    if (outName.empty()) {
        cerr << "\nOutput file is required. See -i option" << endl;
        print_usage(argv[0], -1);
    }

    if (chrsBinFileName.empty()) {
        cerr << "\nCan't open chromosomes masked binary fasta file. See -c option" << endl;
        print_usage(argv[0], -1);
    }

    if (length <= 0) {
        cerr << "\nSequence length to be added before and after the SNP position is required. See -l option" << endl;
        print_usage(argv[0], -1);
    }

    begin = clock();
    cout << "Reading chromosome sequences from binary file" << endl;
    chrFactory.parseFastaFile(chrsBinFileName, true);
    cout << chrFactory.getSequenceContainter().size() << " chromosomes loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    cout << "Reading input SNP coordinates from files" << endl;
    snpFactory.parseSNPFile(inName, length, chrFactory);

    cout << "Writing fasta file" << endl;
    snpFactory.writeEnhansersFastaFile(outName, false);

    delete Global::instance();
    cout << "Total elapse time: " << TimeUtils::instance()->getTimeMinFrom(start) << " minutes" << endl;
    delete TimeUtils::instance();
    return 0;
}


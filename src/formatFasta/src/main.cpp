/* 
 * File:   main.cpp
 * Author: veraalva
 *
 * Created on February 26, 2016, 2:06 PM
 */

#include <iostream>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>

#include "Global.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "FastaFactory.h"

using namespace std;
using namespace sequence;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

void print_usage(char *program_name, int exit_code) {
    cerr << "\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr << "-v    Print info\n";
    cerr << "-h    Display this usage information.\n";
    cerr << "-r    Read binary and print fasta\n";
    cerr << "-i    Input file.\n";
    cerr << "-o    Output file.\n";
    cerr << "-d    Directory with fasta files. Extension: .fa (all files will be combined in one binary file)\n";
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
    string dirName;
    string inName;
    string outName;
    FastaFactory fastaFactory;
    bool reverse = false;

    for (int i = 1; i < argc; i++) {
        string option(argv[i]);
        if (option.compare(0, 1, "-") == 0 && option.size() == 2) {
            if (option.compare(1, 1, "h") == 0) {
                print_usage(argv[0], 0);
            } else if (option.compare(1, 1, "v") == 0) {
                Global::instance()->setVerbose(3);
            } else if (option.compare(1, 1, "r") == 0) {
                reverse = true;
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
            } else if (option.compare(1, 1, "d") == 0) {
                i++;
                if (i < argc) {
                    dirName = argv[i];
                    if (dirName.compare(0, 1, "-") == 0) {
                        cerr << "Option d require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option d require an argument" << endl;
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

    if (outName.empty()) {
        cerr << "\nOutput file is required. See -o option" << endl;
        print_usage(argv[0], -1);
    }

    if (inName.empty() && dirName.empty()) {
        cerr << "\nInput file or directory name are required. See -i or -d options" << endl;
        print_usage(argv[0], -1);
    }

    if (reverse) {
        if (inName.empty()) {
            cerr << "\nIf reverse the input file is required. See -i option" << endl;
            print_usage(argv[0], -1);
        }
        if (!dirName.empty()) {
            cerr << "\nDirectory is not compatible with reverse option. Use -i for input file" << endl;
            print_usage(argv[0], -1);
        }
    }

    begin = clock();
    if (!reverse) {
        if (!dirName.empty()) {
            cout << "Reading sequences from files" << endl;
            fastaFactory.parseFastaInDirectory(dirName, "", ".fa", false);
        } else if (!inName.empty()) {
            cout << "Reading sequences from file" << endl;
            fastaFactory.parseFastaFile(inName, false);
        }
    } else {
        cout << "Reading sequences from file" << endl;
        fastaFactory.parseFastaFile(inName, true);
    }
    cout << fastaFactory.getSequenceContainter().size() << " sequences loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Writing sequences in binary mode ";
    fastaFactory.writeSequencesToFile(outName, !reverse);
    cout << "in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    cout << "Total elapse time: " << TimeUtils::instance()->getTimeMinFrom(start) << " minutes" << endl;
    delete Global::instance();
    delete TimeUtils::instance();
    return 0;
}


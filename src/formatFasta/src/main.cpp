/* 
 * File:   main.cpp
 * Author: veraalva
 *
 * Created on February 26, 2016, 2:06 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <memory>
#include <getopt.h>
#include <stdbool.h>
#include <time.h>
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

char *program_name;

void print_usage(int exit_code) {
    cerr <<"\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr <<"-v,   --verbose                     Print info\n";
    cerr <<"-h,   --help                        Display this usage information.\n";
    cerr <<"-r,   --reverse                     Read binary and print fasta\n";
    cerr <<"-i,   --in                          Input file.\n";
    cerr <<"-o,   --out                         Output file.\n";
    cerr <<"-d,   --dir                         Directory with fasta files. Extension: .fa (all files will be combined in one binary file)\n";
    cerr <<"********************************************************************************\n";
    cerr <<"\n            Shan Li (e-mail: lis11@ncbi.nlm.nih.gov)\n";
    cerr <<"            Roberto Vera Alvarez (e-mail: veraalva@ncbi.nlm.nih.gov)\n\n";
    cerr <<"********************************************************************************\n";
    exit(exit_code);
}

/*
 * 
 */
int main(int argc, char** argv) {
    clock_t begin = clock();
    clock_t start = clock();
    int next_option;
    const char* const short_options = "vhi:o:d:r";
    string dirName;
    string inName;
    string outName;
    FILE *inFile;
    FastaFactory fastaFactory;
    bool reverse = false;
    
    program_name = argv[0];

    srand(time(NULL));

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "out", 1, NULL, 'o'},
        { "in", 1, NULL, 'i'},
        { "dir", 1, NULL, 'd'},
        { "reverse", 0, NULL, 'r'},
        { NULL, 0, NULL, 0} /* Required at end of array.  */
    };

    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);

        switch (next_option) {
            case 'h':
                print_usage(0);

            case 'v':
                Global::instance()->setVerbose(1);
                break;

            case 'i':
                inName = optarg;
                break;

            case 'o':
                outName = optarg;
                break;

            case 'd':
                dirName = optarg;
                break;

            case 'r':
                reverse = true;
                break;
        }
    } while (next_option != -1);

    if (outName.empty()) {
        cerr << "\nOutput file is required. See -o option" << endl;
        print_usage( -1);
    }

    if (inName.empty() && dirName.empty()) {
        cerr << "\nInput file or directory name are required. See -i or -d options" << endl;
        print_usage(-1);
    }

    if (reverse) {
        if (inName.empty()) {
            cerr << "\nIf reverse the input file is required. See -i option" << endl;
            print_usage(-1);
        }
        if (!dirName.empty()) {
            cerr << "\nDirectory is not compatible with reverse option. Use -i for input file" << endl;
            print_usage(-1);
        }
    }

    begin = clock();
    if (!reverse) {
        if (!dirName.empty()) {
            cout << "Reading sequences from files" << endl;
            fastaFactory.parseFastaInDirectory(dirName, "", ".fa", false);
        } else if (!inName.empty()) {
            cout << "Reading sequences from file" << endl;
            inFile = (FILE *) fopen(inName.c_str(), "r");
            if (!inFile) {
                cerr << "Can't open input file" << endl;
                exit(-1);
            }
            fastaFactory.parseFastaFile(inFile, -1, true, false);
            fclose(inFile);
        }
    } else {
        cout << "Reading sequences from file" << endl;
        inFile = (FILE *) fopen(inName.c_str(), "rb");
        if (!inFile) {
            cerr << "Can't open input file" << endl;
            exit(-1);
        }
        fastaFactory.parseFastaFile(inFile, -1, true, true);
        fclose(inFile);
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


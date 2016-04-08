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
#include "FastaFactory.h"

using namespace std;
using namespace fasta;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s \n", program_name);
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "-v,   --verbose                     Print info\n");
    fprintf(stream, "-h,   --help                        Display this usage information.\n");
    fprintf(stream, "-r,   --reverse                     Read binary and print fasta\n");
    fprintf(stream, "-i,   --in                          Input file.\n");
    fprintf(stream, "-o,   --out                         Output file.\n");
    fprintf(stream, "-d,   --dir                         Directory with fasta files. Extension: .fa (all files will be combined in one binary file)\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Shan Li (e-mail: lis11@ncbi.nlm.nih.gov)\n");
    fprintf(stream, "            Roberto Vera Alvarez (e-mail: veraalva@ncbi.nlm.nih.gov)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(0);
}

/*
 * 
 */
int main(int argc, char** argv) {
    clock_t begin = clock();
    clock_t start = clock();
    int next_option;
    const char* const short_options = "vhi:o:d:r";
    char *dirName = NULL;
    char *inName = NULL;
    char *outName = NULL;
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
                print_usage(stdout, 0);

            case 'v':
                Global::instance()->SetVerbose(1);
                break;

            case 'i':
                inName = strdup(optarg);
                break;

            case 'o':
                outName = strdup(optarg);
                break;

            case 'd':
                dirName = strdup(optarg);
                break;

            case 'r':
                reverse = true;
                break;
        }
    } while (next_option != -1);

    if (!outName) {
        cerr << "\nOutput file is required. See -o option" << endl;
        print_usage(stderr, -1);
    }

    if (!inName && !dirName) {
        cerr << "\nInput file or directory name are required. See -i or -d options" << endl;
        print_usage(stderr, -1);
    }

    if (reverse) {
        if (!inName) {
            cerr << "\nIf reverse the input file is required. See -i option" << endl;
            print_usage(stderr, -1);
        }
        if (dirName) {
            cerr << "\nDirectory is not compatible with reverse option. Use -i for input file" << endl;
            print_usage(stderr, -1);
        }
    }

    begin = clock();
    if (!reverse) {
        if (dirName) {
            cout << "Reading sequences from files" << endl;
            fastaFactory.LoadFastaInDirectory(dirName, NULL, ".fa", false);
        } else if (inName) {
            cout << "Reading sequences from file" << endl;
            inFile = (FILE *) fopen(inName, "r");
            if (!inFile) {
                cerr << "Can't open input file" << endl;
                exit(-1);
            }
            fastaFactory.ParseFastaFile(inFile, -1, true, false);
            fclose(inFile);
        }
    } else {
        cout << "Reading sequences from file" << endl;
        inFile = (FILE *) fopen(inName, "rb");
        if (!inFile) {
            cerr << "Can't open input file" << endl;
            exit(-1);
        }
        fastaFactory.ParseFastaFile(inFile, -1, true, true);
        fclose(inFile);
    }
    cout << fastaFactory.GetFastaMap().size() << " sequences loaded in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Writing sequences in binary mode ";
    fastaFactory.WriteSequencesToFile(outName, !reverse);
    cout << "in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    if (outName) free(outName);
    if (inName) free(inName);
    if (dirName) free(dirName);
    cout << "Total elapse time: " << TimeUtils::instance()->GetTimeMinFrom(start) << " minutes" << endl;
    return 0;
}


/* 
 * File:   main.cpp
 * Author: veraalva
 *
 * Created on March 10, 2016, 11:46 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <time.h>

#include <iostream>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "svm.h"
#include "Global.h"
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "FimoFactory.h"
#include "BedFactory.h"
#include "SNPFactory.h"
#include "SVMPredict.h"

using namespace std;
using namespace fasta;
using namespace peak;
using namespace kmers;
using namespace snp;
using namespace svm;
using namespace fimo;

char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s \n", program_name);
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "-v,   --verbose                     Print info\n");
    fprintf(stream, "-h,   --help                        Display this usage information.\n");
    fprintf(stream, "-i,   --in                          Input file with SNP coordinates.\n");
    fprintf(stream, "-c,   --chrs                        Chromosomes binary fasta file. Created with formatFasta.\n");
    fprintf(stream, "-o,   --output                      Fasta file with the sequences\n");
    fprintf(stream, "-l,   --length                      Sequence length to be added before and after the SNP position\n\n");

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
    const char* const short_options = "vhi:c:o:l:";
    FILE *chrsBinFile = NULL;
    char *inName = NULL;
    char *outName = NULL;
    unsigned long int length = 0;
    FastaFactory chrFactory;
    SNPFactory snpFactory;

    program_name = argv[0];

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "in", 1, NULL, 'i'},
        { "output", 1, NULL, 'o'},
        { "chrs", 1, NULL, 'c'},
        { "length", 1, NULL, 'l'},
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

            case 'c':
                chrsBinFile = fopen(optarg, "r");
                break;

            case 'o':
                outName = strdup(optarg);
                break;

            case 'l':
                length = static_cast<unsigned long int> (atoi(optarg));
                break;
        }
    } while (next_option != -1);

    if (!inName) {
        cerr << "\nInput SNP coordinate file is required. See -i option" << endl;
        print_usage(stderr, -1);
    }

    if (!outName) {
        cerr << "\nOutput file is required. See -i option" << endl;
        print_usage(stderr, -1);
    }

    if (!chrsBinFile) {
        cerr << "\nCan't open chromosomes masked binary fasta file. See -c option" << endl;
        print_usage(stderr, -1);
    }

    if (length == 0) {
        cerr << "\nSequence length to be added before and after the SNP position is required. See -l option" << endl;
        print_usage(stderr, -1);
    }

    TimeUtils::instance()->SetStartTime();
    begin = clock();
    cout << "Reading chromosome sequences from binary file" << endl;
    chrFactory.ParseFastaFile(chrsBinFile, -1, true, true);
    cout << chrFactory.GetFastaMap().size() << " chromosomes loaded in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;
    
    cout << "Reading input SNP coordinates from files" << endl;
    snpFactory.ReadSNPFromFile(inName, length, chrFactory);
    
    cout << "Writing fasta file" << endl;
    snpFactory.WriteEnhansersFastaFile(outName, false);

    if (inName) free(inName);
    if (outName) free(outName);
    if (chrsBinFile) fclose(chrsBinFile);
    delete Global::instance();
    cout << "Total elapse time: " << TimeUtils::instance()->GetTimeMinFrom(start) << " minutes" << endl;
    return 0;
}


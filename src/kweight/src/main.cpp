/* 
 * File:   main.cpp
 * Author: veraalva
 *
 * Created on February 10, 2016, 3:18 PM
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

#include "Global.h"
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"

using namespace std;
using namespace fasta;
using namespace peak;
using namespace kmers;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "This option will read the enhancers coordinates and generate the controls by shuffling the enhancers sequences\n");
    fprintf(stream, "kweight --bed enhancers_coordinates.bed -m ./directory/hg19/chromosomes/hg19_masked.fa.bin --poutput weight_file_kmers.txt --eoutput enhancers.fa\n\n");
    fprintf(stream, "This option will read the enhancers coordinates and generate the controls from the chromosomes\n");
    fprintf(stream, "kweight -g --bed enhancers_coordinates.bed -m ./directory/hg19/chromosomes/hg19_masked.fa.bin --poutput weight_file_kmers.txt --eoutput enhancers.fa\n\n");
    fprintf(stream, "This option will read the enhancers coordinates and use the controls provided by the users\n");
    fprintf(stream, "kweight --bed enhancers_coordinates.bed -m ./directory/hg19/hg19_masked.fa.bin --control users_control.txt --poutput weight_file_kmers.txt --eoutput enhancers.fa\n");
    fprintf(stream, "Control file format:\n");
    fprintf(stream, "chr<tab>start_position<tab>length\n\n");
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s \n", program_name);
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "-v,   --verbose                     Print info\n");
    fprintf(stream, "-h,   --help                        Display this usage information.\n");
    fprintf(stream, "-b,   --bed                         Coordination of enhancers or DHS peaks.\n");
    fprintf(stream, "-m,   --masked                      Chromosomes masked binary fasta file. Created with formatFasta.\n");
    fprintf(stream, "-p,   --poutput                     Output file with the p-value for all kmers\n");
    
    fprintf(stream, "-o,   --order                       Order (default: 10).\n");    
    fprintf(stream, "-g,   --genCtrl                     Generate control from chromosomes. Default: NO.\n");
    fprintf(stream, "-n,   --hitNum                      Number of controls per peak (default: 10, if -g set default: 3).\n");
    fprintf(stream, "-c,   --control                     Bed file with the control coordinates. This option will use your own control for the calculations.\n");
    
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Shan Li (e-mail: lis11@ncbi.nlm.nih.gov)\n");
    fprintf(stream, "            Roberto Vera Alvarez (e-mail: veraalva@ncbi.nlm.nih.gov)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(0);
}

int main(int argc, char** argv) {
    clock_t begin = clock();
    clock_t start = clock();
    int next_option;
    unsigned long int hitNum;
    const char* const short_options = "vhm:o:b:n:gc:p:ir:u:";
    FILE *chrsBinFile = NULL;
    char *bedFileName = NULL;
    char *controlFileName = NULL;
    char *poutputFileName = NULL;
    char *prefix = NULL;
    char *subfix = NULL;
    FastaFactory chrFactory;
    BedFactory bedFactory;
    bool genCtrl = false;
    bool binary = false;

    program_name = argv[0];

    srand(time(NULL));
    Global::instance()->SetOrder(10);

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "masked", 1, NULL, 'm'},
        { "order", 1, NULL, 'o'},
        { "bed", 1, NULL, 'b'},
        { "hitNum", 1, NULL, 'n'},
        { "genCtrl", 0, NULL, 'g'},
        { "control", 1, NULL, 'c'},
        { "poutput", 1, NULL, 'p'},
        { "prefix", 1, NULL, 'r'},
        { "subfix", 1, NULL, 'u'},
        { "bin", 0, NULL, 'i'},
        { NULL, 0, NULL, 0} /* Required at end of array.  */
    };

    hitNum = 0;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);

        switch (next_option) {
            case 'h':
                print_usage(stdout, 0);

            case 'v':
                Global::instance()->SetVerbose(1);
                break;

            case 'g':
                genCtrl = true;
                break;

            case 'i':
                binary = true;
                break;

            case 'm':
                chrsBinFile = fopen(optarg, "r");
                break;

            case 'r':
                prefix = strdup(optarg);
                break;

            case 'u':
                subfix = strdup(optarg);
                break;

            case 'p':
                poutputFileName = strdup(optarg);
                break;
            
            case 'c':
                controlFileName = strdup(optarg);
                break;

            case 'b':
                bedFileName = strdup(optarg);
                break;

            case 'n':
                hitNum = atoi(optarg);
                break;

            case 'o':
                Global::instance()->SetOrder(atoi(optarg));
                break;
        }
    } while (next_option != -1);

    if (genCtrl) {
        if (hitNum == 0) hitNum = 3;
    } else {
        if (hitNum == 0) hitNum = 10;
    }

    if (!poutputFileName) {
        cerr << "\nKmer contingency table output is required. See -k option" << endl;
        print_usage(stderr, -1);
    }

    if (!chrsBinFile) {
        cerr << "\nCan't open chromosomes masked binary fasta file. See -m option" << endl;
        print_usage(stderr, -1);
    }

    if (!bedFileName) {
        cerr << "\nBed file is required. See -b option" << endl;
        print_usage(stderr, -1);
    }

    TimeUtils::instance()->SetStartTime();
    Global::instance()->SetBin1(0.005);
    Global::instance()->SetBin2(0.01);
    cout.precision(2);

    begin = clock();
    cout << "Creating kmers genomewide" << endl;
    KmersFactory kmersFactory;
    kmersFactory.CreateGenomeWideKmers();
    cout << kmersFactory.GetKmersGenome().size() << " kmers created in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading chromosome sequences from binary file" << endl;
    chrFactory.ParseFastaFile(chrsBinFile, -1, true, true);
    cout << chrFactory.GetFastaMap().size() << " chromosomes loaded in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading peaks from file" << endl;
    bedFactory.CreatePeaksFromBedFile(chrFactory, bedFileName, 0.7, kmersFactory);
    cout << bedFactory.GetPeaks().size() << " peaks loaded in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;
    cout << bedFactory.GetGCNcontentBin().size() << " chromosomes to analyze" << endl;

    begin = clock();
    if (!controlFileName) {
        if (genCtrl) {
            cout << "Generating controls from Chromosomes" << endl;
            bedFactory.GeneratingControlsFromChromosomes(chrFactory, hitNum, kmersFactory);
        } else {
            cout << "Generating controls from shuffling peaks" << endl;
            bedFactory.GeneratingControlsFromShufflingPeaks(hitNum, kmersFactory);
        }
    } else {
        cout << "Reading controls from file" << endl;
        bedFactory.ReadControlsFromFile(controlFileName, chrFactory, kmersFactory);
    }
    cout << "Controls generated in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Generating kmers from peaks and controls" << endl;
    kmersFactory.BuildKmers();
    cout << kmersFactory.GetKmers().size() << " kmers generated in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    kmersFactory.WriteKmersToFile(poutputFileName, binary);

    if (prefix) free(prefix);
    if (subfix) free(subfix);
    if (chrsBinFile) fclose(chrsBinFile);
    if (bedFileName) free(bedFileName);
    if (poutputFileName) free(poutputFileName);
    delete Global::instance();
    cout << "Total elapse time: " << TimeUtils::instance()->GetTimeMinFrom(start) << " minutes" << endl;
    return 0;
}


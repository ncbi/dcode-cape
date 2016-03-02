
/* 
 * File:   main.cpp
 * Author: veraalva
 *
 * Created on February 25, 2016, 12:06 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_statistics_double.h>
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
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "svm.h"
#include "Global.h"
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"
#include "SNPFactory.h"
#include "SVMPredict.h"

using namespace std;
using namespace fasta;
using namespace peak;
using namespace kmers;
using namespace snp;
using namespace svm;

char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s \n", program_name);
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "-v,   --verbose                     Print info\n");
    fprintf(stream, "-h,   --help                        Display this usage information.\n");
    fprintf(stream, "-i    --in                          Input config file\n");

    fprintf(stream, "\nInput conf file format (tab delimited), copy the next 9 lines to your config file:\n\n");
    fprintf(stream, "in\tinput_file_name.txt\t\t\t# Input file with SNP coordinates\n");
    fprintf(stream, "out\toutput_file_name.out\t\t\t# Output file with SNP coordinates and probabilities\n");
    fprintf(stream, "order\t10\t\t\t\t\t# Order (default: 10)\n");
    fprintf(stream, "chrs\t/path-to/hg19.fa.bin\t\t\t# Chromosomes files in binary mode. Format: hg19.fa.bin. Binary files created by formatFasta\n");
    fprintf(stream, "weight\t/path-to/10mers_sigValue_sorted\t\t# Kmers weight file. Generated with kweight\n");
    fprintf(stream, "binary\t0\t\t\t\t\t# 1 if weight file is binary\n");
    fprintf(stream, "neighbors\t100\t\t\t\t# Pb to be added before and after the SNP position. Default 100\n");
    fprintf(stream, "model\t/path-to/svm.model\t\t\t# SVM Model\n");
    fprintf(stream, "probability\t1\t\t\t\t# 1 if the model use probability estimates\n\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Roberto Vera Alvarez (e-mail: veraalva@ncbi.nlm.nih.gov)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(0);
}


using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    clock_t begin = clock();
    clock_t start = clock();
    int next_option, count;
    unsigned long int hitNum;
    const char* const short_options = "vhci:";
    FILE *chrsBinName = NULL;
    char *inFileName = NULL;
    FILE *inFile = NULL;
    char *svmModelName = NULL;
    char *weightFileName = NULL;
    FILE *outputFile = NULL;
    unsigned long int neighbors = 100;
    FastaFactory chrFactory;
    SNPFactory snpFactory;
    KmersFactory kmersFactory;
    bool binary = false;
    SVMPredict svmPredict;
    char *line = NULL;
    size_t len = 0;
    char **fields = NULL;
    size_t fieldsSize = 0;

    program_name = argv[0];

    srand(time(NULL));
    Global::instance()->SetOrder(10);

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "in", 1, NULL, 'i'},
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

            case 'i':
                inFileName = strdup(optarg);
                break;
        }
    } while (next_option != -1);

    if (!inFileName) {
        cerr << "\nInput file with SNP coordinates is required. See -i option" << endl;
        print_usage(stderr, -1);
    }

    inFile = (FILE *) checkPointerError(fopen(inFileName, "r"), "\nCan't open input file. See -i option\n", __FILE__, __LINE__, -1);

    free(inFileName);
    inFileName = NULL;

    while (getline(&line, &len, inFile) != -1) {
        if (*line != '#') {
            if (*(line + strlen(line) - 1) == '\n') *(line + strlen(line) - 1) = '\0';
            fieldsSize = splitString(&fields, line, "\t");
            if (fieldsSize < 2) {
                printLog(stderr, "Input config file with a wrong format", __FILE__, __LINE__, 0);
                print_usage(stderr, -1);
                exit(-1);
            }
            if (strcmp(fields[0], "in") == 0) {
                inFileName = strdup(fields[1]);
            }
            if (strcmp(fields[0], "out") == 0) {
                outputFile = (FILE *) checkPointerError(fopen(fields[1], "w"), "\nCan't open output file. See -i option\n", __FILE__, __LINE__, -1);
            }
            if (strcmp(fields[0], "order") == 0) {
                Global::instance()->SetOrder(static_cast<unsigned long int> (atoi(fields[1])));
            }
            if (strcmp(fields[0], "chrs") == 0) {
                chrsBinName = (FILE *) checkPointerError(fopen(fields[1], "rb"), "\nCan't open chromosome binary file. See -i option\n", __FILE__, __LINE__, -1);
            }
            if (strcmp(fields[0], "weight") == 0) {
                weightFileName = strdup(fields[1]);
            }
            if (strcmp(fields[0], "binary") == 0) {
                if (atoi(fields[1]) == 1)
                    binary = true;
            }
            if (strcmp(fields[0], "neighbors") == 0) {
                neighbors = static_cast<unsigned long int> (atoi(fields[1]));
            }
            if (strcmp(fields[0], "model") == 0) {
                svmModelName = strdup(fields[1]);
            }
            if (strcmp(fields[0], "probability") == 0) {
                svmPredict.SetPredict_probability(atoi(fields[1]));
            }
            freeArrayofPointers((void **) fields, fieldsSize);
        }
    }
    if (line) free(line);

    if (!inFileName) {
        cerr << "\nInput file with SNP coordinates is required in config file." << endl;
        print_usage(stderr, -1);
    }

    if (!weightFileName) {
        cerr << "\nKmer weight file is required in config file" << endl;
        print_usage(stderr, -1);
    }

    if (!chrsBinName) {
        cerr << "\nChromosomes file in binary mode is required in config file" << endl;
        print_usage(stderr, -1);
    }

    if (!outputFile) {
        cerr << "\nOutput file is required in config file." << endl;
        print_usage(stderr, -1);
    }

    TimeUtils::instance()->SetStartTime();
    Global::instance()->SetBin1(0.005);
    Global::instance()->SetBin2(0.01);
    cout.precision(2);

    begin = clock();
    cout << "Reading chromosome sequences from binary file" << endl;
    chrFactory.ParseFastaFile(chrsBinName, -1, true, true);
    cout << chrFactory.size() << " chromosomes loaded in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading SVM model" << endl;
    svmPredict.SVMLoadModel(svmModelName);
    cout << " SVM model processed in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading kmers weight" << endl;
    kmersFactory.ReadKmersFromFile(weightFileName, binary);
    cout << kmersFactory.GetKmers().size() << " kmers loaded in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading input SNP coordinates from files" << endl;
    count = snpFactory.ProcessSNPFromFiles(inFileName, neighbors, chrFactory, kmersFactory, svmPredict);
    cout << count << " SNP processed in " << TimeUtils::instance()->GetTimeSecFrom(begin) << " seconds" << endl;

    fprintf(outputFile, "#chrom\tpos\trsID\trefAle\taltAle\tscore\n");
    for (auto it = snpFactory.GetSnps().begin(); it != snpFactory.GetSnps().end(); ++it) {
        SNP *s = *it;

        fprintf(outputFile, "%s\t%lu\t%s\t%c\t%c\t%.6f\n",
                s->GetChr().c_str(), s->GetChrPos() + 1, s->GetId().c_str(),
                s->GetRef(), s->GetAlt(), s->GetProbPos());
    }

    fclose(outputFile);
    fclose(chrsBinName);
    if (inFileName) free(inFileName);
    if (weightFileName) free(weightFileName);
    if (svmModelName) free(svmModelName);
    delete Global::instance();
    cout << "Total elapse time: " << TimeUtils::instance()->GetTimeSecFrom(start) << " seconds" << endl;
    return 0;
}


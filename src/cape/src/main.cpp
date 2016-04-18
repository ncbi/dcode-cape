
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
#include <getopt.h>
#include <stdbool.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <locale>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "svm.h"

#include "Global.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "FileParserFactory.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "FimoFactory.h"
#include "BedFactory.h"
#include "TFBSFactory.h"
#include "SNPFactory.h"
#include "SVMPredict.h"

using namespace std;
using namespace parsers;
using namespace sequence;
using namespace peak;
using namespace kmers;
using namespace snp;
using namespace svm;
using namespace fimo;
using namespace tfbs;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

char *program_name;

void print_usage(int exit_code) {
    cerr << "\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr << "-v,   --verbose                     Print info\n";
    cerr << "-h,   --help                        Display this usage information.\n";
    cerr << "-i    --in                          Input config file\n";

    cerr << "\nInput conf file format (tab delimited), copy the next 9 lines to your config file:\n\n";
    cerr << "in\tinput_file_name.txt\t\t\t# Input file with SNP coordinates\n";
    cerr << "out\toutput_file_name.out\t\t\t# Output file with SNP coordinates and probabilities\n";
    cerr << "order\t10\t\t\t\t\t# Order (default: 10)\n";
    cerr << "chrs\t/path-to/hg19.fa.bin\t\t\t# Chromosomes files in binary mode. Format: hg19.fa.bin. Binary files created by formatFasta\n";
    cerr << "weight\t/path-to/10mers_sigValue_sorted\t\t# Kmers weight file. Generated with kweight\n";
    cerr << "neighbors\t100\t\t\t\t# Pb to be added before and after the SNP position. Default 100\n";
    cerr << "model\t/path-to/svm.model\t\t\t# SVM Model\n";
    cerr << "probability\t1\t\t\t\t# 1 if the model use probability estimates\n";
    cerr << "fimo\tfimo_output_file.txt\t\t\t# Use FIMO output. Set to: 0 for not using FIMO output\n";
    cerr << "pwm_EnsembleID\tpwm_EnsembleID_mapping\t\t# File mapping TF names with Ensembl IDs. Provided in resources folder\n";
    cerr << "expression\t57epigenomes.RPKM.pc\t\t# Expression file\n";
    cerr << "expression_code\tE116\t\t\t\t# Tissue code used to extract expression data from the expression file.\n\t\t\t\t\t\t  Extract this code from the tissue ID mapping file EG.name.txt\n";
    cerr << "abbrev-mtf-mapped\tabbrev-mtf-mapped-to-whole-label.all.info.renamed\t\t\t\t# TF name cutoff P-Value mapped by our group\n\n";
    cerr << "********************************************************************************\n";
    cerr << "For internal use at NCBI:\n\n";
    cerr << "\tThese parameters can be include into the input conf file to use FIMO\n\toutput through NCBI internal index files\n";
    cerr << "\tPlease, note that \"fimo\" option should be set to \"0\" or simple\n\tdelete it from the input file\n\n";
    cerr << "TibInfoFileName\ttib/hg19/tib.info\n";
    cerr << "TFBSIdxDirName\tcommon/tib/hg19/5\t\t\t# The program reads files with name chrN.idx and chrN.tib\n\n";
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
    int next_option, count;
    const char* const short_options = "vhci:d";
    FILE *chrsBinFile = NULL;
    string inFileName;
    string svmModelName;
    string weightFileName;
    string fimoFileName;
    string pwmEnsembleIDFileName;
    string expressionFileName;
    string expressionCode;
    string abbrevmtfmappedFileName;
    string tibInfoFileName;
    string tFBSIdxDirName;
    ofstream outputFile;
    unsigned long int neighbors = 100;
    FastaFactory chrFactory;
    SNPFactory snpFactory;
    KmersFactory kmersFactory;
    FimoFactory fimoFactory;
    SVMPredict svmPredict;
    TFBSFactory tFBSFactory;
    FileParserFactory fParser;

    program_name = argv[0];

    srand(time(NULL));
    Global::instance()->setOrder(10);

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "in", 1, NULL, 'i'},
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
                inFileName = optarg;
                break;

            case 'd':
                Global::instance()->setVerbose(3);
                break;
        }
    } while (next_option != -1);

    if (inFileName.empty()) {
        cerr << "\nInput file with SNP coordinates is required. See -i option" << endl;
        print_usage(-1);
    }

    try {
        fParser.setFileToParse(inFileName);
        inFileName.clear();

        while (fParser.iterate('#', "\t")) {
            if (fParser.getNWords() >= 2) {
                string field(fParser.getWords()[0]);
                if (field.compare("in") == 0) {
                    inFileName = fParser.getWords()[1];
                }
                if (field.compare("out") == 0) {
                    outputFile.open(fParser.getWords()[1]);
                }
                if (field.compare("order") == 0) {
                    Global::instance()->setOrder(static_cast<unsigned long int> (atoi(fParser.getWords()[1])));
                }
                if (field.compare("chrs") == 0) {
                    chrsBinFile = (FILE *) checkPointerError(fopen(fParser.getWords()[1], "rb"), "\nCan't open chromosome binary file. See -i option\n", __FILE__, __LINE__, -1);
                }
                if (field.compare("weight") == 0) {
                    weightFileName = fParser.getWords()[1];
                }
                if (field.compare("neighbors") == 0) {
                    neighbors = static_cast<unsigned long int> (atoi(fParser.getWords()[1]));
                }
                if (field.compare("model") == 0) {
                    svmModelName = fParser.getWords()[1];
                }
                if (field.compare("probability") == 0) {
                    svmPredict.setPredictProbability(atoi(fParser.getWords()[1]));
                }
                if (field.compare("fimo") == 0) {
                    string fimoTmpName(fParser.getWords()[1]);
                    if (fimoTmpName.compare("0") != 0) {
                        fimoFileName = fParser.getWords()[1];
                    }
                }
                if (field.compare("pwm_EnsembleID") == 0) {
                    pwmEnsembleIDFileName = fParser.getWords()[1];
                }
                if (field.compare("expression") == 0) {
                    expressionFileName = fParser.getWords()[1];
                }
                if (field.compare("expression_code") == 0) {
                    expressionCode = fParser.getWords()[1];
                }
                if (field.compare("abbrev-mtf-mapped") == 0) {
                    abbrevmtfmappedFileName = fParser.getWords()[1];
                }
                if (field.compare("TibInfoFileName") == 0) {
                    tibInfoFileName = fParser.getWords()[1];
                }
                if (field.compare("TFBSIdxDirName") == 0) {
                    tFBSIdxDirName = fParser.getWords()[1];
                }
            }
        }
    } catch (exceptions::FileNotFoundException ex) {
        cerr << ex.what() << endl;
        cerr << "Error parsing file" << endl;
        exit(-1);
    } catch (exceptions::ErrorReadingFromFileException ex) {
        cerr << ex.what() << endl;
        cerr << "Error parsing file" << endl;
        exit(-1);
    }

    if (inFileName.empty()) {
        cerr << "\nInput file with SNP coordinates is required in config file." << endl;
        print_usage(-1);
    }

    if (weightFileName.empty()) {
        cerr << "\nKmer weight file is required in config file" << endl;
        print_usage(-1);
    }

    if (!chrsBinFile) {
        cerr << "\nChromosomes file in binary mode is required in config file" << endl;
        print_usage(-1);
    }

    if (!outputFile.is_open()) {
        cerr << "\nOutput file is required in config file." << endl;
        print_usage(-1);
    }

    if (fimoFileName.empty()) {
        if (pwmEnsembleIDFileName.empty()) {
            cerr << "\npwm_EnsembleID option is required in config file if FIMO ouput is provided." << endl;
            print_usage(-1);
        }
        if (expressionFileName.empty()) {
            cerr << "\nexpression option  is required in config file if FIMO ouput is provided." << endl;
            print_usage(-1);
        }
        if (expressionCode.empty()) {
            cerr << "\nexpression_code option  is required in config file if FIMO ouput is provided." << endl;
            print_usage(-1);
        }
        if (abbrevmtfmappedFileName.empty()) {
            cerr << "\nabbrev-mtf-mapped option  is required in config file if FIMO ouput is provided." << endl;
            print_usage(-1);
        }
    } else if (!tFBSIdxDirName.empty() || !tibInfoFileName.empty()) {
        if (expressionCode.empty()) {
            cerr << "\nexpression_code option  is required in config file if FIMO indexes are provided." << endl;
            print_usage(-1);
        }
        if (tFBSIdxDirName.empty()) {
            cerr << "\tTFBSIdxDirName option  is required in config file if FIMO indexes are provided." << endl;
            print_usage(-1);
        }
        if (tibInfoFileName.empty()) {
            cerr << "\tTibInfoFileName option  is required in config file if FIMO indexes are provided." << endl;
            print_usage(-1);
        }

        snpFactory.setExpressionCode(expressionCode);

        tFBSFactory.createTFBSFileIndexMap(tFBSIdxDirName, "chr", ".idx", ".tib");
        tFBSFactory.createPWMIndexFromTibInfoFile(tibInfoFileName);
    }

    TimeUtils::instance()->setStartTime();
    Global::instance()->setBin1(0.005);
    Global::instance()->setBin2(0.01);
    cout.precision(2);

    begin = clock();
    cout << "Reading chromosome sequences from binary file" << endl;
    chrFactory.parseFastaFile(chrsBinFile, -1, true, true);
    cout << chrFactory.getSequenceContainter().size() << " chromosomes loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    if (!pwmEnsembleIDFileName.empty() && !expressionFileName.empty()) {
        cout << "Creating PWM indexes from files" << endl;
        fimoFactory.createTissueIndexFromFiles(pwmEnsembleIDFileName, expressionFileName);
    }

    if (!fimoFileName.empty()) {
        begin = clock();
        cout << "Parsing FIMO output file" << endl;
        fimoFactory.createCutoffIndexFromFile(abbrevmtfmappedFileName, 4);
        fimoFactory.parseFimoOutput(fimoFileName, expressionCode, neighbors);
        cout << fimoFactory.getSnpIDContainer().size() << " SNP with FIMO expression loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;
    }

    begin = clock();
    cout << "Reading SVM model" << endl;
    svmPredict.svmLoadModel(svmModelName);
    cout << " SVM model processed in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading kmers weight" << endl;
    kmersFactory.readKmersFromFile(weightFileName, false);
    cout << kmersFactory.getKmers().size() << " kmers loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Reading input SNP coordinates from files" << endl;
    count = snpFactory.processSNPFromFile(inFileName, neighbors, chrFactory, kmersFactory, svmPredict, fimoFactory, tFBSFactory);
    cout << count << " SNP processed in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    outputFile << "#chrom\tpos\trsID\trefAle\taltAle\tscore\n";
    for (auto it = snpFactory.getSnps().begin(); it != snpFactory.getSnps().end(); ++it) {
        SNP *s = *it;

        outputFile << s->getChr() << "\t"
                << s->getChrPos() + 1 << "\t"
                << s->getId() << "\t" <<
                s->getRef() << "\t"
                << s->getAlt() << "\t"
                << s->getProbPos() << endl;
    }

    outputFile.close();
    fclose(chrsBinFile);
    delete Global::instance();
    cout << "Total elapse time: " << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << endl;
    delete TimeUtils::instance();
    return 0;
}


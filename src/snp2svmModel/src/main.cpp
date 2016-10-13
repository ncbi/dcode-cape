/*  $Id: main.cpp October 13, 2016 10:03 AM veraalva $
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
 *   Alternative splicing detector
 *
 */
#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <locale>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>

#include "berror.h"
#include "bmemory.h"
#include "svm.h"

#include "Global.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "cstring.h"
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

void print_usage(char *program_name, int exit_code) {
    cerr << "\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr << "-v    Print info\n";
    cerr << "-h    Display this usage information.\n";
    cerr << "-i    Input config file\n";

    cerr << "\nInput conf file format (tab delimited), copy the next 9 lines to your config file:\n\n";
    cerr << "in\tinput_file_name.txt\t\t\t# Input file with SNP coordinates and the (1/-1) values for SNPs\n";
    cerr << "svm_type\t0\t\t\tSet type of SVM (default 0)\n"
            << "\t\t\t\t0 -- C-SVC		(multi-class classification)\n"
            << "\t\t\t\t1 -- nu-SVC		(multi-class classification)\n"
            << "\t\t\t\t2 -- one-class SVM\n"
            << "\t\t\t\t3 -- epsilon-SVR	(regression)\n"
            << "\t\t\t\t4 -- nu-SVR		(regression)\n";
    cerr << "kernel_type\t0\t\t\tSet type of kernel function (default 2)\n"
            << "\t\t\t\t0 -- linear: u'*v\n"
            << "\t\t\t\t1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
            << "\t\t\t\t2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
            << "\t\t\t\t3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
            << "\t\t\t\t4 -- precomputed kernel (kernel values in training_set_file)\n";
    cerr << "w1\t4\t\t\tWeight for the positive set.\n";
    cerr << "w-1\t1\t\t\tWeight for the negative set.\n";
    cerr << "probability\t1\t\t\t\t# Whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 1)\n";
    cerr << "model\t/path-to/svm.model\t\t\t# SVM Model obtained from the data\n";
    cerr << "order\t4,6,8,10,12\t\t\t\t\t# Order (default: 4,6,8,10,12)\n";
    cerr << "chrs\t/path-to/hg19.fa.bin\t\t\t# Chromosomes files in binary mode. Format: hg19.fa.bin. Binary files created by formatFasta\n";
    cerr << "weight\t/path-to/kmers_sigValue_sorted\t\t# Kmers weight file. Generated with kweight o kmerge \n";
    cerr << "neighbors\t100\t\t\t\t# Pb to be added before and after the SNP position. Default 100\n";
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

int main(int argc, char** argv) {

    clock_t begin = clock();
    clock_t start = clock();
    int count;
    string chrsBinFileName;
    string inFileName;
    string weightFileName;
    string fimoFileName;
    string pwmEnsembleIDFileName;
    string expressionFileName;
    string expressionCode;
    string abbrevmtfmappedFileName;
    string tibInfoFileName;
    string tFBSIdxDirName;
    string svmModelName;
    vector<string> orders;
    unsigned long int neighbors = 100;
    FastaFactory chrFactory;
    SNPFactory snpFactory;
    KmersFactory kmersFactory;
    FimoFactory fimoFactory;
    SVMPredict svmPredict;
    TFBSFactory tFBSFactory;
    FileParserFactory fParser;

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
                    inFileName = argv[i];
                    if (inFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option i require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option i require an argument" << endl;
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

    if (inFileName.empty()) {
        cerr << "\nInput file with SNP coordinates is required. See -i option" << endl;
        print_usage(argv[0], -1);
    }

    try {
        fParser.setFileToParse(inFileName);
        while (fParser.iterate("#", "\t")) {
            if (fParser.getWords().size() >= 2) {
                if (fParser.getWords()[0].compare("in") == 0) {
                    inFileName = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("order") == 0) {
                    cstring::split(fParser.getWords()[1], ",", orders);
                    for (auto it = orders.begin(); it != orders.end(); ++it) {
                        Global::instance()->getOrders().insert(atoi((*it).c_str()));
                    }
                }
                if (fParser.getWords()[0].compare("chrs") == 0) {
                    chrsBinFileName = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("weight") == 0) {
                    weightFileName = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("neighbors") == 0) {
                    neighbors = static_cast<unsigned long int> (atoi((fParser.getWords()[1]).c_str()));
                }
                if (fParser.getWords()[0].compare("fimo") == 0) {
                    if (fParser.getWords()[1].compare("0") != 0) {
                        fimoFileName = fParser.getWords()[1];
                    }
                }
                if (fParser.getWords()[0].compare("model") == 0) {
                    svmModelName = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("probability") == 0) {
                    svmPredict.setPredictProbability(atoi((fParser.getWords()[1]).c_str()));
                }
                if (fParser.getWords()[0].compare("svm_type") == 0) {
                    svmPredict.setSVMType(atoi((fParser.getWords()[1]).c_str()));
                }
                if (fParser.getWords()[0].compare("kernel_type") == 0) {
                    svmPredict.setSVMKernelType(atoi((fParser.getWords()[1]).c_str()));
                }
                if (fParser.getWords()[0].compare("w1") == 0) {
                    svmPredict.setSVMw1(atof((fParser.getWords()[1]).c_str()));
                }
                if (fParser.getWords()[0].compare("w-1") == 0) {
                    svmPredict.setSVMwMinus1(atof((fParser.getWords()[1]).c_str()));
                }
                if (fParser.getWords()[0].compare("pwm_EnsembleID") == 0) {
                    pwmEnsembleIDFileName = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("expression") == 0) {
                    expressionFileName = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("expression_code") == 0) {
                    expressionCode = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("abbrev-mtf-mapped") == 0) {
                    abbrevmtfmappedFileName = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("TibInfoFileName") == 0) {
                    tibInfoFileName = fParser.getWords()[1];
                }
                if (fParser.getWords()[0].compare("TFBSIdxDirName") == 0) {
                    tFBSIdxDirName = fParser.getWords()[1];
                }
            }
        }
    } catch (std::exception) {
        cerr << "Error parsing file: " << inFileName << endl;
        exit(-1);
    }

    if (inFileName.empty()) {
        cerr << "\nInput file with SNP coordinates is required in config file." << endl;
        print_usage(argv[0], -1);
    }

    if (weightFileName.empty()) {
        cerr << "\nKmer weight file is required in config file" << endl;
        print_usage(argv[0], -1);
    }

    if (chrsBinFileName.empty()) {
        cerr << "\nChromosomes file in binary mode is required in config file" << endl;
        print_usage(argv[0], -1);
    }

    if (svmModelName.empty()) {
        cerr << "\nSVM Model Output file is required in config file." << endl;
        print_usage(argv[0], -1);
    }

    if (Global::instance()->getOrders().empty()) {
        Global::instance()->getOrders().insert(4);
        Global::instance()->getOrders().insert(6);
        Global::instance()->getOrders().insert(8);
        Global::instance()->getOrders().insert(10);
        Global::instance()->getOrders().insert(12);
    }

    if (!fimoFileName.empty()) {
        if (pwmEnsembleIDFileName.empty()) {
            cerr << "\npwm_EnsembleID option is required in config file if FIMO ouput is provided." << endl;
            print_usage(argv[0], -1);
        }
        if (expressionFileName.empty()) {
            cerr << "\nexpression option  is required in config file if FIMO ouput is provided." << endl;
            print_usage(argv[0], -1);
        }
        if (expressionCode.empty()) {
            cerr << "\nexpression_code option  is required in config file if FIMO ouput is provided." << endl;
            print_usage(argv[0], -1);
        }
        if (abbrevmtfmappedFileName.empty()) {
            cerr << "\nabbrev-mtf-mapped option  is required in config file if FIMO ouput is provided." << endl;
            print_usage(argv[0], -1);
        }
    } else if (!tFBSIdxDirName.empty() || !tibInfoFileName.empty()) {
        if (expressionCode.empty()) {
            cerr << "\nexpression_code option  is required in config file if FIMO indexes are provided." << endl;
            print_usage(argv[0], -1);
        }
        if (tFBSIdxDirName.empty()) {
            cerr << "\tTFBSIdxDirName option  is required in config file if FIMO indexes are provided." << endl;
            print_usage(argv[0], -1);
        }
        if (tibInfoFileName.empty()) {
            cerr << "\tTibInfoFileName option  is required in config file if FIMO indexes are provided." << endl;
            print_usage(argv[0], -1);
        }

        begin = clock();
        cout << "Loading FIMO indexes address" << endl;
        snpFactory.setExpressionCode(expressionCode);
        tFBSFactory.createTFBSFileIndexMap(tFBSIdxDirName, "chr", ".idx", ".tib");
        cout << tFBSFactory.getTfbsFileIndexSize() << " indexes loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

        begin = clock();
        cout << "Loading TIB data" << endl;
        tFBSFactory.createPWMIndexFromTibInfoFile(tibInfoFileName);
        cout << tFBSFactory.getPwmIndex().size() << " TIBs loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;
    }

    Global::instance()->setBin1(0.005);
    Global::instance()->setBin2(0.01);

    begin = clock();
    cout << "Reading chromosome sequences from binary file" << endl;
    chrFactory.parseFastaFile(chrsBinFileName, true);
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
    cout << "Reading kmers weight" << endl;
    kmersFactory.readKmersFromFile(weightFileName);
    cout << kmersFactory.getKmers().size() << " kmers loaded in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    begin = clock();
    cout << "Processing input SNP coordinates from files" << endl;
    count = snpFactory.createSVMModelFromSNPFile(inFileName, neighbors, chrFactory, kmersFactory, svmPredict, fimoFactory, tFBSFactory, svmModelName);
    cout << count << " SNP processed in " << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds" << endl;

    delete Global::instance();
    cout << "Total elapse time: " << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << endl;
    delete TimeUtils::instance();
    return 0;
}



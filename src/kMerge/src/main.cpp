/*  $Id: main.cpp October 12, 2016 2:18 PM veraalva $
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
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <fstream>

#include "Global.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "cstring.h"
#include "FileParserFactory.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"

using namespace std;
using namespace parsers;
using namespace sequence;
using namespace peak;
using namespace kmers;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

void print_usage(char *program_name, int exit_code) {
    cerr << "\n********************************************************************************\n";
    cerr << "This program reads from a text file with a list of kweight file names and concatenate them together   \n";
    cerr << "\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr << "-v    Print info\n";
    cerr << "-h    Display this usage information.\n";
    cerr << "-i    Text file with a list of kweight file names.\n";
    cerr << "-o    Output file.\n";

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
    string inputFileName;
    string outputFileName;
    FileParserFactory fParser;
    KmersFactory kmersFactory;

    for (int i = 1; i < argc; i++) {
        string option(argv[i]);
        if (option.compare(0, 1, "-") == 0 && option.compare(1, 1, "-") != 0 && option.size() == 2) {
            if (option.compare(1, 1, "h") == 0) {
                print_usage(argv[0], 0);
            } else if (option.compare(1, 1, "v") == 0) {
                Global::instance()->setVerbose(1);
            } else if (option.compare(1, 1, "i") == 0) {
                i++;
                if (i < argc) {
                    inputFileName = argv[i];
                    if (inputFileName.compare(0, 1, "-") == 0) {
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
                    outputFileName = argv[i];
                    if (outputFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option o require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option o require an argument" << endl;
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

    if (inputFileName.empty()) {
        cerr << "\nText file with a list of kweight file names is required. See -i option" << endl;
        print_usage(argv[0], -1);
    }

    if (outputFileName.empty()) {
        cerr << "\nOutput file name is required. See -o option" << endl;
        print_usage(argv[0], -1);
    }

    try {
        fParser.setFileToParse(inputFileName);
        while (fParser.iterate("#", "\t")) {
            KmersFactory kmersFac;
            cout << "Reading kmers weight from file " << fParser.getWords()[0] << endl;
            kmersFac.readKmersFromFile(fParser.getWords()[0]);
            kmersFactory.mergeKmers(kmersFac);
            cout << kmersFac.getKmers().size() << " kmers loaded"<< endl;
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file: " << inputFileName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file: " << inputFileName << endl;
        exit(-1);
    }

    
    kmersFactory.writeKmersToFile(outputFileName);

    return 0;
}


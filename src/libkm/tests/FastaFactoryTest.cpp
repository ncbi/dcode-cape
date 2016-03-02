/* 
 * File:   FastaFactoryTest.cpp
 * Author: veraalva
 *
 * Created on February 11, 2016, 9:07 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <memory>
#include <getopt.h>
#include <stdbool.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "Global.h"
#include "KmersFactory.h"
#include "BedFactory.h"

using namespace std;
using namespace kmers;
using namespace fasta;
using namespace peak;

void testParseFastaFile(char *fileName) {
    Fasta *fasta;
    FastaFactory fastaFactory;
    char *seg;
    std::string test_seg("ATCCGACATCAAGTGCCCACCTTGGCTCGTGGCTCTCACTGCAACGGGAAAGCCACAGACTGGGGTGAAGAGTTCAGTCACATGCGACCGGTGACTCCCTGTCCCCACCCCCATGACACTCCCCAGCCCTCCAAGGCCACTGTGTTTCCCAGTTAGCTCAGAGCCTCAGTCGATCCCTGACCCAGCACCGGGCACTGATGAGACAGCGGCTGTTTGAGGAGCCACCTCCCAGCCACCTCGGGGCCAGGGCCAGGGTGTGCAGCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAACCATAGTGCCCAGGGCACTGCCGCTGCAGGCGCAGGCATCGCATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACCATTGCCTTGGACCGCTCTTGGCAGTCGAAGAAGATTCTCCTGTCAGTTTGAGCTGGGTGAGCTTAGAGAGGAAAGCTCCACTATGGCTCCCAAACCAGGAAGGAGCCATAGCCCAGGCAGGAGGGCTGAGGACCTCTGGTGGCGGCCCAGGGCTTCCAGCATGTGCCCTAGGGGAAGCAGGGGCCAGCTGGCAAGAGCAGGGGGTGGGCAGAAAGCACCCGGTGGACTCAGGGCTGGAGGGGAGGAGGCGATCTTGCCCAAGGCCCTCCGACTGCAAGCTCCAGGGCCCGCTCACCTTGCTCCTGCTCCTTCTGCTGCTGCTTCTCCAGCTTTCGCTCCTTCATGCTGCGCAGCTTGGCCTTGCCGATGCCCCCAGCTTGGCGGATGGACTCTAGCAGAGTGGCCAGCCACCGGAGGGGTCAACCACTTCCCTGGGAGCTCCCTGGACTGGAGCCGGGAGGTGGGGAACAGGGCAAGGAGGAAAGGCTGCTCAGGCAGGG");

    FILE *fName = (FILE *) checkPointerError(fopen(fileName, "r"), "Can't open input file", __FILE__, __LINE__, -1);

    long unsigned int seqRead = fastaFactory.ParseFastaFile(fName, -1, true, false);
    if (fastaFactory.GetFastaVector().size() != 1 || seqRead != 1) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=It should read 1 sequence" << endl;
    }
    fasta = (*fastaFactory.GetFastaVector().begin()).get();
    if (fasta->GetId().compare("chr1") != 0) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The sequence Id should be chr1" << endl;
    }
    if (fasta->GetLength() != 1249950 || strlen(fasta->GetSeq()) != fasta->GetLength()) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The sequence should have 1249950 bp and it is " << fasta->GetLength() << endl;
    }
    
    seg = fasta->GetSubStr(15000, 1000);
    if (!seg){
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=Can't get the segment of length 1050 starting at 15000 " << endl;
    }
    if (strcmp(seg, test_seg.c_str()) != 0){
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The segment of length 1050 starting at 15000 is not equal to the test" << endl;
    }

    free(seg);
    fclose(fName);
}

void testParseFastaMultipleFile(char *fileName) {
    long unsigned int seqRead, totalRead; 
    FastaFactory fFactory;
    Fasta *f;

    map<string, long unsigned int> seqSize;

    seqSize.insert(pair<string, long unsigned int>("chr1", 950));
    seqSize.insert(pair<string, long unsigned int>("chr2", 700));
    seqSize.insert(pair<string, long unsigned int>("chr3", 1200));
    seqSize.insert(pair<string, long unsigned int>("chr4", 1450));
    seqSize.insert(pair<string, long unsigned int>("chr5", 1700));

    FILE *fName = (FILE *) checkPointerError(fopen(fileName, "r"), "Can't open input file", __FILE__, __LINE__, -1);

    totalRead = 0;
    while ((seqRead = fFactory.ParseFastaFile(fName, 2, false, false)) != 0) {
        totalRead += seqRead;
        if (fFactory.size() != totalRead) {
            cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=It should read 1 sequence" << endl;
        }
        
        for (auto it = fFactory.GetFastaVector().begin(); it != fFactory.GetFastaVector().end(); ++it) {
            f = (*it).get();
            if (seqSize.find(f->GetId())->second != f->GetLength() || strlen(f->GetSeq()) != f->GetLength()) {
                cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Wrong sequence size "
                        << seqSize.find(f->GetId())->second << " != " << f->GetLength() << endl;
            }
        }
        
    }

    f = fFactory.GetFastaFromID("chr4");
    if (!f) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Cant find chr4 in the factory container " << endl;
    } else {
        if (seqSize.find(f->GetId())->second != f->GetLength() || strlen(f->GetSeq()) != f->GetLength()) {
            cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Wrong sequence size "
                    << seqSize.find(f->GetId())->second << " != " << f->GetLength() << endl;
        }
    }

    fclose(fName);
}

void testParseFastaDir(char *fileName) {
    FastaFactory fFactory;
    string prefix("chr");
    string suxif( ".fa.masked");

    fFactory.LoadFastaInDirectory(fileName, prefix.c_str(), suxif.c_str(), false);

}

int main(int argc, char** argv) {
    char *singleFileName;
    char *multipleFileName;
    char *dirName = NULL;
    cout << "%SUITE_STARTING% FastaFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    Global::instance()->SetVerbose(0);

    if (argc == 4) {
        singleFileName = argv[1];
        multipleFileName = argv[2];
        dirName = argv[3];
    } else {
        singleFileName = strdup("resources/chr1.fa.masked");
        multipleFileName = strdup("resources/multiple.fa.masked");
    }

    cout << "%TEST_STARTED% testParseFastaFile (FastaFactoryTest)" << endl;
    testParseFastaFile(singleFileName);
    cout << "%TEST_FINISHED% time=0 testParseFastaFile (FastaFactoryTest)" << endl;

    cout << "%TEST_STARTED% testParseFastaMultipleFile (FastaFactoryTest)" << endl;
    testParseFastaMultipleFile(multipleFileName);
    cout << "%TEST_FINISHED% time=0 testParseFastaMultipleFile (FastaFactoryTest)" << endl;
    
    if (dirName != NULL) {
        Global::instance()->SetVerbose(3);
        cout << "%TEST_STARTED% testParseFastaDir (FastaFactoryTest)" << endl;
        testParseFastaDir(dirName);
        cout << "%TEST_FINISHED% time=0 testParseFastaDir (FastaFactoryTest)" << endl;
    }

    cout << "%SUITE_FINISHED% time=0" << endl;

    if (argc != 4) {
        free(singleFileName);
        free(multipleFileName);
    }
    
    delete Global::instance();
    return (EXIT_SUCCESS);
}


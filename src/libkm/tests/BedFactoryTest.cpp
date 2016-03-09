/* 
 * File:   BedFactoryTest.cpp
 * Author: veraalva
 *
 * Created on February 11, 2016, 4:14 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>

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
#include "Global.h"
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"

using namespace std;
using namespace kmers;
using namespace fasta;
using namespace peak;


/*
 * Simple C++ Test Suite
 */

void testCreatePeaksFromBedFile(const char *bFName, char *dirName, char *testName) {
    int i = 0;
    fasta::FastaFactory chrFactory;
    string genomeName("hg19");
    string bedFileName(bFName);
    peak::BedFactory bedFactory;
    char *line = NULL;
    size_t len = 0;
    string prefix("chr");
    string sufix(".fa.masked");
    map<unsigned long int, map<string, unsigned long int>> testMap;
    char **fields = NULL;
    size_t fieldsSize = 0;
    KmersFactory kmersFactory;
    kmersFactory.CreateGenomeWideKmers();

    string test("TGAGAGGAAAGCTTTCCCACATTATACAGCTTCTGAAAGGGTTGCTTGACCCACAGATGTGAAGCTGAGGCTGAAGGAGACTGATGTGGTTTCTCCTCAGTTTCTCTGTGCAGCACCAGGTGGCAGCAGAGGTCAGCAAGGCAAACCCGAGCCCGGGGATGCGGAGTGGGGGCAGCTACGTCCTCTCTTG");

    FILE *fd = (FILE *) checkPointerError(fopen(testName, "r"), "Can't open test file", __FILE__, __LINE__, -1);
    while (getline(&line, &len, fd) != -1) {
        fieldsSize = splitString(&fields, line, "\t");
        map<string, unsigned long int> inMap;
        char *n = strchr(fields[2], '-');
        *n = '\0';
        i = atoi(fields[2]);
        inMap.insert(pair<string, unsigned long int>("start", i));
        inMap.insert(pair<string, unsigned long int>("end", atoi(fields[2] + strlen(fields[2]) + 1)));
        inMap.insert(pair<string, unsigned long int>("width", atoi(fields[3])));
        inMap.insert(pair<string, unsigned long int>("gccount", atoi(fields[4])));
        inMap.insert(pair<string, unsigned long int>("ncount", atoi(fields[5])));
        inMap.insert(pair<string, unsigned long int>("nrcount", atoi(fields[6])));

        testMap.insert(pair<unsigned long int, map<string, unsigned long int>>(i, inMap));

        freeArrayofPointers((void **) fields, fieldsSize);
    }
    if (line) free(line);
    fclose(fd);

    chrFactory.LoadFastaInDirectory(dirName, prefix.c_str(), sufix.c_str(), false);
    
    bedFactory.CreatePeaksFromBedFile(chrFactory, genomeName, bedFileName, 0.7, kmersFactory);
    
    if (bedFactory.GetPeaks().size() != 88) {
        cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Peaks to print should be equal to 88 and it is " << i << endl;
    }
    
    
    for (auto it = bedFactory.GetPeaks().begin(); it != bedFactory.GetPeaks().end(); ++it) {
        Peak *p = *it;
        map<string, unsigned long int> inMap = testMap.find(p->GetStart())->second;
        if (inMap.find("width")->second != p->GetLength()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong width: "
                    << p->GetLength() << "!=" << inMap.find("width")->second << endl;
        }
        if (inMap.find("gccount")->second != p->GetGCCount()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong GCCount: "
                    << p->GetGCCount() << "!=" << inMap.find("gccount")->second << endl;
        }
        if (inMap.find("ncount")->second != p->GetNCount()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong NCount: "
                    << p->GetNCount() << "!=" << inMap.find("ncount")->second << endl;
        }
        if (inMap.find("nrcount")->second != p->GetNRCount()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong NRCount: "
                    << p->GetNRCount() << "!=" << inMap.find("nrcount")->second << endl;
        }
        if (strlen(p->GetSeq()) != p->GetLength()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Sequence length not equal to width" << endl;
        }
        if (p->GetStart() == 237640 && p->GetEnd() == 237829) {
            if (strcmp(p->GetSeq(), test.c_str()) != 0) {
                cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=First peak is not equal to the test" << endl;
            }
        }
    }
    
    bedFactory.GeneratingControlsFromChromosomes(chrFactory, 3, kmersFactory);
    kmersFactory.BuildKmers();
    if (kmersFactory.GetKmers().size() != 11848 ){
        cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong number of kmers 11848 != " << kmersFactory.GetKmers().size() << endl;
    }
}

void testReadingControlsFromFile(char *dirName, char *controlName){
    fasta::FastaFactory chrFactory;
    string genomeName("hg19");
    peak::BedFactory bedFactory;
    string prefix("chr");
    string sufix(".fa.masked");
    KmersFactory kmersFactory;
    kmersFactory.CreateGenomeWideKmers();
    
    chrFactory.LoadFastaInDirectory(dirName, prefix.c_str(), sufix.c_str(), false);
    
    bedFactory.ReadingControlsFromFile(controlName, chrFactory, kmersFactory);
    if (kmersFactory.GetTotalNRnt_control() != 1065){
        cout << "%TEST_FAILED% time=0 testname=testReadingControlsFromFile (BedFactoryTest) message=Wrong calculation of Total non N 1065 != " << kmersFactory.GetTotalNRnt_control() << endl;
    }
}

int main(int argc, char** argv) {
    char *bFName = NULL;
    char *dirName = NULL;
    char *testName = NULL;
    char *controlName = NULL;
    cout << "%SUITE_STARTING% BedFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    Global::instance()->SetVerbose(0);
    Global::instance()->SetOrder(10);
    Global::instance()->SetBin1(0.005);
    Global::instance()->SetBin2(0.01);

    if (argc == 3) {
        bFName = argv[1];
        dirName = argv[2];
    } else {
        bFName = strdup("resources/test.bed");
        dirName = strdup("resources/");
        testName = strdup("resources/test.train.control.info.chr1");
        controlName = strdup("resources/randControl.txt");
    }

    cout << "%TEST_STARTED% testCreatePeaksFromBedFile (BedFactoryTest)" << endl;
    testCreatePeaksFromBedFile(bFName, dirName, testName);
    cout << "%TEST_FINISHED% time=0 testCreatePeaksFromBedFile (BedFactoryTest)" << endl;
    
    cout << "%TEST_STARTED% testReadingControlsFromFile (BedFactoryTest)" << endl;
    testReadingControlsFromFile(dirName, controlName);
    cout << "%TEST_FINISHED% time=0 testReadingControlsFromFile (BedFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=0" << endl;

    if (bFName) free(bFName);
    if (dirName) free(dirName);
    if (testName) free(testName);
    if (controlName) free(controlName);
    return (EXIT_SUCCESS);
}


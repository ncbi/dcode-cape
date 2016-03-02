/* 
 * File:   SNPFactoryTest.cpp
 * Author: veraalva
 *
 * Created on February 25, 2016, 10:39 AM
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
#include "Global.h"
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"
#include "SNPFactory.h"

using namespace std;
using namespace kmers;
using namespace fasta;
using namespace snp;

/*
 * Simple C++ Test Suite
 */

void testReadSNPFromFiles(char *dirName, char *snpFileName) {
    fasta::FastaFactory chrFactory;
    string prefix("chr");
    string sufix(".fa.masked");
    KmersFactory kmersFactory;
    SNPFactory snpFactory;
    
    chrFactory.LoadFastaInDirectory(dirName, prefix.c_str(), sufix.c_str(), false);
    
    snpFactory.ProcessSNPFromFiles(snpFileName, 100, chrFactory, kmersFactory);   
    
}

int main(int argc, char** argv) {
    char *dirName = strdup("resources/");
    char *snpFileName = strdup("resources/snp.txt");   
    
    std::cout << "%SUITE_STARTING% SNPFactoryTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;
    
    Global::instance()->SetVerbose(0);
    Global::instance()->SetOrder(10);
    Global::instance()->SetBin1(0.005);
    Global::instance()->SetBin2(0.01);

    std::cout << "%TEST_STARTED% testReadSNPFromFiles (SNPFactoryTest)" << std::endl;
    testReadSNPFromFiles(dirName, snpFileName);
    std::cout << "%TEST_FINISHED% time=0 testReadSNPFromFiles (SNPFactoryTest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    free(snpFileName);
    free(dirName);
    return (EXIT_SUCCESS);
}


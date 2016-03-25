/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TFBSFactoryTest.cpp
 * Author: veraalva
 *
 * Created on March 24, 2016, 1:51 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"

#include "Exceptions.h"
#include "FastaFactory.h"
#include "FimoFactory.h"
#include "TFBSFactory.h"

using namespace std;
using namespace exceptions;
using namespace fasta;
using namespace fimo;
using namespace tfbs;

/*
 * Simple C++ Test Suite
 */

void testExtractTFBSFromFile(char *pwm_tFName, char *tissue_file, char *tibInfoFileName, char* chrIdxFileName, char* chrTibFileName, char *fastaFile) {
    FastaFactory fastaFactory;
    TFBSFactory tFBSFactory;
    FimoFactory fimoFactory;
    double c, v;
    std::pair<std::string, double> cPair;

    FILE *fName = (FILE *) checkPointerError(fopen(fastaFile, "r"), "Can't open input file", __FILE__, __LINE__, -1);
    fastaFactory.ParseFastaFile(fName, -1, true, false);
    
    fimoFactory.CreateTissueIndexFromFiles(pwm_tFName, tissue_file);
    
    cPair = fimoFactory.GetTissueValue("NFYA_f1", "E000");
    c = cPair.second;
    v = 17.511;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    tFBSFactory.CreatePWMIndexFromTibInfoFile(tibInfoFileName);
    
    tFBSFactory.ExtractTFBSFromFile(chrIdxFileName, chrTibFileName, 100871, 100881, fastaFactory.GetFastaFromID("chr1"));
    
    for(auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it){
        TFBS *t  = *it;
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        cout << tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName() << " " << t->GetStart() << "-" << t->GetEnd() << " " << cPair.second << endl;
        
    }

    fclose(fName);
}

int main(int argc, char** argv) {
    char *pwm_tFName= strdup("resources/pwm_EnsembleID_mapping");
    char *tissue_file= strdup("resources/57epigenomes.RPKM.pc");
    char *tibInfoFileName = strdup("/Volumes/devdcode/common/tib/hg19/tib.info");
    char* chrIdxFileName = strdup("/Volumes/devdcode/common/tib/hg19/5/chr1.idx");
    char* chrTibFileName = strdup("/Volumes/devdcode/common/tib/hg19/5/chr1.tib");
    char *fastaFile = strdup("/Users/veraalva/Work/resources/genomes/human/hg19/Chromosomes/fasta/chr1.fa");
    cout << "%SUITE_STARTING% TFBSFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    cout << "%TEST_STARTED% testExtractTFBSFromFile (TFBSFactoryTest)" << endl;
    testExtractTFBSFromFile(pwm_tFName, tissue_file, tibInfoFileName, chrIdxFileName, chrTibFileName, fastaFile);
    cout << "%TEST_FINISHED% time=0 testExtractTFBSFromFile (TFBSFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=0" << endl;

    free(tibInfoFileName);
    free(chrIdxFileName);
    free(chrTibFileName);
    free(fastaFile);
    return (EXIT_SUCCESS);
}


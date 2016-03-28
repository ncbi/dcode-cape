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
#include <math.h>

#include <iostream>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"

#include "Exceptions.h"
#include "Global.h"
#include "TimeUtils.h"
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

void testExtractTFBSFromFile(char *pwm_tFName, char *tissue_file, char *tFBSIdxDirName, char *tibInfoFileName, char *fastaFile) {
    FastaFactory fastaFactory;
    TFBSFactory tFBSFactory;
    FimoFactory fimoFactory;
    double c, v;
    long int start, end;
    string name;
    std::pair<std::string, double> cPair;
    TFBS *t;

    FILE *fName = (FILE *) checkPointerError(fopen(fastaFile, "r"), "Can't open input file", __FILE__, __LINE__, -1);
    fastaFactory.ParseFastaFile(fName, -1, true, false);

    fimoFactory.CreateTissueIndexFromFiles(pwm_tFName, tissue_file);

    cPair = fimoFactory.GetTissueValue("NFYA_f1", "E000");
    c = cPair.second;
    v = 17.511;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    tFBSFactory.CreateTFBSFileIndexMap(tFBSIdxDirName, "chr", ".idx", ".tib");
    tFBSFactory.CreatePWMIndexFromTibInfoFile(tibInfoFileName);

    tFBSFactory.ExtractTFBSFromFile(0, 100, fastaFactory.GetFastaFromID("chr1"));

    // ATF5_si 2-13 23.466
    name = "ATF5_si";
    start = 2;
    end = 13;
    v = 23.466;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }

    // CDC5L_si 2-16 16.729
    name = "CDC5L_si";
    start = 2;
    end = 16;
    v = 16.729;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }

    // CEBPG_si 2-14 19.485
    name = "CEBPG_si";
    start = 2;
    end = 14;
    v = 19.485;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // ESR1_do 2-19 1.535
    name = "ESR1_do";
    start = 2;
    end = 19;
    v = 1.535;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // FOXC1_DBD_1 2-14 11.194
    name = "FOXC1_DBD_1";
    start = 2;
    end = 14;
    v = 11.194;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // FOXD3.p2 2-13 1.304
    name = "FOXD3.p2";
    start = 2;
    end = 13;
    v = 1.304;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // FOXD3_f1 2-10 1.304
    name = "FOXD3_f1";
    start = 2;
    end = 10;
    v = 1.304;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // FOXJ3_dimer 2-13 18.426
    name = "FOXJ3_dimer";
    start = 2;
    end = 13;
    v = 18.426;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // FOXJ3_si 2-14 18.426
    name = "FOXJ3_si";
    start = 2;
    end = 14;
    v = 18.426;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // FOXO1 2-15 0
    name = "FOXO1";
    start = 2;
    end = 15;
    v = 0;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // UP00034_1 79-100 0.0
    name = "UP00034_1";
    start = 79;
    end = 100;
    v = 0.0;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // cell08_Six3_1732 32-41 0.699
    name = "cell08_Six3_1732";
    start = 32;
    end = 41;
    v = 0.699;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }
    
    // EPAS1_si 46-58 30.333
    name = "EPAS1_si";
    start = 46;
    end = 58;
    v = 30.333;
    t = NULL;
    for (auto it = tFBSFactory.GetTfbs().begin(); it != tFBSFactory.GetTfbs().end(); ++it) {
        if (tFBSFactory.GetPwmIndex()[(*it)->GetIndex() - 1]->GetName().compare(name) == 0 && (*it)->GetStart() == start && (*it)->GetEnd() == end) {
            t = *it;
            break;
        }
    }
    if (!t) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
    } else {
        cPair = fimoFactory.GetTissueValue(tFBSFactory.GetPwmIndex()[t->GetIndex() - 1]->GetName(), "E000");
        if (cPair.second != v) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
        }
    }

    fclose(fName);
}

int main(int argc, char** argv) {
    char *pwm_tFName = strdup("resources/pwm_EnsembleID_mapping");
    char *tissue_file = strdup("resources/57epigenomes.RPKM.pc");
    char *tibInfoFileName = strdup("resources/tib.info");
    char* tFBSIdxDirName = strdup("resources/");
    char *fastaFile = strdup("resources/chr1.fa.masked");
    cout << "%SUITE_STARTING% TFBSFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    Global::instance()->SetVerbose(0);

    cout << "%TEST_STARTED% testExtractTFBSFromFile (TFBSFactoryTest)" << endl;
    testExtractTFBSFromFile(pwm_tFName, tissue_file, tFBSIdxDirName, tibInfoFileName, fastaFile);
    cout << "%TEST_FINISHED% time=0 testExtractTFBSFromFile (TFBSFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=0" << endl;

    free(pwm_tFName);
    free(tissue_file);
    free(tibInfoFileName);
    free(tFBSIdxDirName);
    free(fastaFile);
    return (EXIT_SUCCESS);
}


/* 
 * File:   FimoFactoryTest.cpp
 * Author: veraalva
 *
 * Created on March 4, 2016, 8:38 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <locale>
#include <unordered_map>
#include <set>
#include <vector>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "Global.h"
#include "TimeUtils.h"
#include "FimoFactory.h"

using namespace std;
using namespace fimo;

/*
 * Simple C++ Test Suite
 */

void testCreateIndexFromFiles(char *pwm_tFName, char *tissue_file, char* cutoffFileName, char *fimoOutput) {
    FimoFactory fimoFactory;
    double c, v;
    std::pair<std::string, double> cPair;

    fimoFactory.CreateTissueIndexFromFiles(pwm_tFName, tissue_file);
    fimoFactory.CreateCutoffIndexFromFile(cutoffFileName, 4);
    
    cPair = fimoFactory.GetTissueValue("NFYA_f1", "E000");
    c = cPair.second;
    v = 17.511;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.GetTissueValue("ATF3_f1", "E000");
    c = cPair.second;
    v = 7.478;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.GetTissueValue("RXRA_f1", "E006");
    c = cPair.second;
    v = 12.171;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.GetTissueValue("M4681_1.02", "E006");
    c = cPair.second;
    v = 0.174;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.GetTissueValue("UP00066_1", "E000");
    c = cPair.second;
    v = 10.457;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.GetTissueValue("UP00055_1", "E000");
    c = cPair.second;
    v = 68.989;
    /*
     * ENSG00000115677 68.989
     * ENSG00000172232 29.033
     * ENSG00000013583 0.646
     * ENSG00000163950 13.156
     * ENSG00000115145 60.357     * 
     * */
    if (cPair.first.compare("ENSG00000115677") != 0){
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result ENSG00000115677 " << " != " << cPair.first << endl;
    }
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    c = fimoFactory.GetCutoffValue("UP00066_1");
    v = 1E-04;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong cutoff value " << v << " != " << c << endl;
    }
    c = fimoFactory.GetCutoffValue("UP00066_1");
    v = 1E-04;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong cutoff value " << v << " != " << c << endl;
    }

    c = fimoFactory.GetCutoffValue("ATF5_si");
    v = 8.15E-05;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong cutoff value " << v << " != " << c << endl;
    }
    c = fimoFactory.GetCutoffValue("ALX1_si");
    v = 8.33E-05;
    if (fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong cutoff value " << v << " != " << c << endl;
    }

    fimoFactory.ParseFimoOutput(fimoOutput, "E116", 100);
}

int main(int argc, char** argv) {
    cout << "%SUITE_STARTING% FimoFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    char *pwm_tFName = strdup("resources/pwm_EnsembleID_mapping");
    char *tissue_file = strdup("resources/57epigenomes.RPKM.pc");
    char *cutoffFileName = strdup("resources/abbrev-mtf-mapped-to-whole-label.all.info.renamed");
    char *fimoOutput = strdup("resources/chr1.fimo.txt");

    cout << "%TEST_STARTED% test1 (FimoFactoryTest)" << endl;
    testCreateIndexFromFiles(pwm_tFName, tissue_file, cutoffFileName, fimoOutput);
    cout << "%TEST_FINISHED% time=0 test1 (FimoFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=0" << endl;

    if (pwm_tFName) free(pwm_tFName);
    if (tissue_file) free(tissue_file);
    if (cutoffFileName) free(cutoffFileName);
    if (fimoOutput) free(fimoOutput);
    return (EXIT_SUCCESS);
}


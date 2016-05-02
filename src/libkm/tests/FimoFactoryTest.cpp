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
#include <memory>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <locale>
#include <unordered_map>
#include <set>
#include <vector>
#include <ctime>

#include "berror.h"
#include "bmemory.h"
#include "Global.h"
#include "TimeUtils.h"
#include "FimoFactory.h"

using namespace std;
using namespace fimo;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

/*
 * Simple C++ Test Suite
 */

void testCreateIndexFromFiles(string pwm_tFName, string tissue_file, string cutoffFileName, string fimoOutput) {
    FimoFactory fimoFactory;
    double c, v;
    std::pair<std::string, double> cPair;

    fimoFactory.createTissueIndexFromFiles(pwm_tFName, tissue_file);
    fimoFactory.createCutoffIndexFromFile(cutoffFileName, 4);

    cPair = fimoFactory.getTissueValue("NFYA_f1", "E000");
    c = cPair.second;
    v = 17.511;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.getTissueValue("ATF3_f1", "E000");
    c = cPair.second;
    v = 7.478;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.getTissueValue("RXRA_f1", "E006");
    c = cPair.second;
    v = 12.171;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.getTissueValue("M4681_1.02", "E006");
    c = cPair.second;
    v = 0.174;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.getTissueValue("UP00066_1", "E000");
    c = cPair.second;
    v = 10.457;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    cPair = fimoFactory.getTissueValue("UP00055_1", "E000");
    c = cPair.second;
    v = 68.989;
    /*
     * ENSG00000115677 68.989
     * ENSG00000172232 29.033
     * ENSG00000013583 0.646
     * ENSG00000163950 13.156
     * ENSG00000115145 60.357     * 
     * */
    if (cPair.first.compare("ENSG00000115677") != 0) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result ENSG00000115677 " << " != " << cPair.first << endl;
    }
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong result " << v << " != " << c << endl;
    }

    c = fimoFactory.getCutoffValue("UP00066_1");
    v = 1E-04;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong cutoff value " << v << " != " << c << endl;
    }
    c = fimoFactory.getCutoffValue("UP00066_1");
    v = 1E-04;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong cutoff value " << v << " != " << c << endl;
    }

    c = fimoFactory.getCutoffValue("ATF5_si");
    v = 8.15E-05;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong cutoff value " << v << " != " << c << endl;
    }
    c = fimoFactory.getCutoffValue("ALX1_si");
    v = 8.33E-05;
    if (std::fabs(c - v) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testCreateIndexFromFiles (FimoFactoryTest) message=Wrong cutoff value " << v << " != " << c << endl;
    }

    fimoFactory.parseFimoOutput(fimoOutput, "E116", 100);
}

int main(int argc, char** argv) {
    clock_t start = clock();
    clock_t begin;
    cout << "%SUITE_STARTING% FimoFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    string pwm_tFName("resources/pwm_EnsembleID_mapping");
    string tissue_file("resources/57epigenomes.RPKM.pc");
    string cutoffFileName("resources/abbrev-mtf-mapped-to-whole-label.all.info.renamed");
    string fimoOutput("resources/chr1.fimo.txt");

    begin = clock();
    cout << "%TEST_STARTED% testCreateIndexFromFiles (FimoFactoryTest)" << endl;
    testCreateIndexFromFiles(pwm_tFName, tissue_file, cutoffFileName, fimoOutput);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testCreateIndexFromFiles (FimoFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << endl;

    delete Global::instance();
    delete TimeUtils::instance();
    return (EXIT_SUCCESS);
}


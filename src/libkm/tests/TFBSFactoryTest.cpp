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
#include <memory>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <ctime>
#include <set>

#include "berror.h"
#include "bmemory.h"

#include "Exceptions.h"
#include "Global.h"
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "FimoFactory.h"
#include "TFBSFactory.h"

using namespace std;
using namespace exceptions;
using namespace sequence;
using namespace fimo;
using namespace tfbs;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

/*
 * Simple C++ Test Suite
 */

void testExtractTFBSFromFile(string pwm_tFName, string tissue_file, string tFBSIdxDirName, string tibInfoFileName, string fastaFile) {
    FastaFactory fastaFactory;
    TFBSFactory tFBSFactory;
    FimoFactory fimoFactory;
    double c, v;
    long int start, end;
    string name;
    std::pair<std::string, double> cPair;
    shared_ptr<TFBS> t;

    fastaFactory.parseFastaFile(fastaFile, false);

    fimoFactory.createTissueIndexFromFiles(pwm_tFName, tissue_file);

    cPair = fimoFactory.getTissueValue("NFYA_f1", "E000");
    c = cPair.second;
    v = 4.210311;
    if (std::fabs(c - v) >= 10E-6) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message=Wrong result " << v << " != " << c << " " << std::fabs(c - v) << endl;
    }

    tFBSFactory.createTFBSFileIndexMap(tFBSIdxDirName, "chr", ".idx", ".tib");
    tFBSFactory.createPWMIndexFromTibInfoFile(tibInfoFileName);

    try {
        tFBSFactory.extractTFBSFromFile(0, 100, fastaFactory.getSequenceFromID("chr1"));

        // ATF5_si 2-13 23.466
        name = "ATF5_si";
        start = 2;
        end = 13;
        v = 4.612706;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // CDC5L_si 2-16 16.729
        name = "CDC5L_si";
        start = 2;
        end = 16;
        v = 4.148039;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // CEBPG_si 2-14 19.485
        name = "CEBPG_si";
        start = 2;
        end = 14;
        v = 4.356496;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // ESR1_do 2-19 1.535
        name = "ESR1_do";
        start = 2;
        end = 19;
        v = 1.341986;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // FOXC1_DBD_1 2-14 11.194
        name = "FOXC1_DBD_1";
        start = 2;
        end = 14;
        v = 3.6081;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // FOXD3.p2 2-13 1.304
        name = "FOXD3.p2";
        start = 2;
        end = 13;
        v = 1.204141;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // FOXD3_f1 2-10 1.304
        name = "FOXD3_f1";
        start = 2;
        end = 10;
        v = 1.204141;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // FOXJ3_dimer 2-13 18.426
        name = "FOXJ3_dimer";
        start = 2;
        end = 13;
        v = 4.279917;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // FOXJ3_si 2-14 18.426
        name = "FOXJ3_si";
        start = 2;
        end = 14;
        v = 4.279917;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // FOXO1 2-15 0
        name = "FOXO1";
        start = 2;
        end = 15;
        v = 0;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // UP00034_1 79-100 0.0
        name = "UP00034_1";
        start = 79;
        end = 100;
        v = 0.0;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // cell08_Six3_1732 32-41 0.699
        name = "cell08_Six3_1732";
        start = 32;
        end = 41;
        v = 0.7646859;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }

        // EPAS1_si 46-58 30.333
        name = "EPAS1_si";
        start = 46;
        end = 58;
        v = 4.969611;
        t = NULL;
        for (auto it = tFBSFactory.getTfbs().begin(); it != tFBSFactory.getTfbs().end(); ++it) {
            if (tFBSFactory.getPwmIndex()[(*it)->getIndex() - 1]->getName().compare(name) == 0 && (*it)->getStart() == start && (*it)->getEnd() == end) {
                t = *it;
                break;
            }
        }
        if (!t) {
            cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " does not exist at position " << start << "-" << end << endl;
        } else {
            cPair = fimoFactory.getTissueValue(tFBSFactory.getPwmIndex()[t->getIndex() - 1]->getName(), "E000");
            if (std::fabs(cPair.second - v) >= 10E-6) {
                cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message= TFBS " << name << " with wrong tissue value " << v << " != " << cPair.second << endl;
            }
        }
    } catch (exceptions::NotFoundException) {
        cout << "%TEST_FAILED% time=0 testname=testExtractTFBSFromFile (TFBSFactoryTest) message=Not found sequence chr1" << endl;
    }
}

int main() {
    clock_t start = clock();
    clock_t begin;
    string pwm_tFName("resources/pwm_EnsembleID_mapping");
    string tissue_file("resources/57epigenomes.RPKM.pc");
    string tibInfoFileName("resources/tib.info");
    string tFBSIdxDirName("resources/");
    string fastaFile("resources/chr1.fa.masked");
    cout << "%SUITE_STARTING% TFBSFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    Global::instance()->setVerbose(0);

    begin = clock();
    cout << "%TEST_STARTED% testExtractTFBSFromFile (TFBSFactoryTest)" << endl;
    testExtractTFBSFromFile(pwm_tFName, tissue_file, tFBSIdxDirName, tibInfoFileName, fastaFile);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testExtractTFBSFromFile (TFBSFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << endl;

    delete Global::instance();
    delete TimeUtils::instance();
    return (EXIT_SUCCESS);
}


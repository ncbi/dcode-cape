/* 
 * File:   PeakTest.cpp
 * Author: veraalva
 *
 * Created on February 23, 2016, 3:16 PM
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
#include <unordered_set>
#include <map>
#include <set>
#include <ctime>

#include "berror.h"
#include "bmemory.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "FastaFactory.h"
#include "Global.h"
#include "KmersFactory.h"
#include "BedFactory.h"
#include "bmath.h"
#include "cstring.h"

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

using namespace std;

/*
 * Simple C++ Test Suite
 */

void testRandomizePeakSeq() {
    peak::Peak peak;
    peak.setChr("chr1");
    
    string seq = "NNNAAGCAAATGGC";
    peak.setSeq("NNNAAGCAAATGGC");
    string result = peak.randomizePeakSeq();
    if (cstring::countCharacter(seq, "A") != cstring::countCharacter(result, "A")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count A 1" << endl;
    }
    if (cstring::countCharacter(seq, "G") != cstring::countCharacter(result, "G")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count G 1" << endl;
    }
    if (cstring::countCharacter(seq, "C") != cstring::countCharacter(result, "C")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count C 1" << endl;
    }
    if (cstring::countCharacter(seq, "T") != cstring::countCharacter(result, "T")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count T 1" << endl;
    }
    if (result.compare(0, 3, seq, 0, 3) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 1" << endl;
    }
    
    seq = "NNNAAGCAAATGGCNNAAGATGAG";
    peak.setSeq("NNNAAGCAAATGGCNNAAGATGAG");
    result = peak.randomizePeakSeq();
    if (cstring::countCharacter(seq, "A") != cstring::countCharacter(result, "A")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count A 2" << endl;
    }
    if (cstring::countCharacter(seq, "G") != cstring::countCharacter(result, "G")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count G 2" << endl;
    }
    if (cstring::countCharacter(seq, "C") != cstring::countCharacter(result, "C")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count C 2" << endl;
    }
    if (cstring::countCharacter(seq, "T") != cstring::countCharacter(result, "T")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count T 2" << endl;
    }
    if (result.compare(0, 3, seq, 0, 3) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 2" << endl;
    }
    if (result.compare(14, 2, seq, 14, 2) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 4" << endl;
    }

    seq = "AAGCAAATGGCNNAAGATGAGNNN";
    peak.setSeq("AAGCAAATGGCNNAAGATGAGNNN");
    result = peak.randomizePeakSeq();
    if (cstring::countCharacter(seq, "A") != cstring::countCharacter(result, "A")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count A 3" << endl;
    }
    if (cstring::countCharacter(seq, "G") != cstring::countCharacter(result, "G")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count G 3" << endl;
    }
    if (cstring::countCharacter(seq, "C") != cstring::countCharacter(result, "C")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count C 3" << endl;
    }
    if (cstring::countCharacter(seq, "T") != cstring::countCharacter(result, "T")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count T 3" << endl;
    }
    if (result.compare(11, 2, seq, 11, 2) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 5" << endl;
    }
    if (result.compare(22, 3, seq, 22, 3) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 6" << endl;
    }

    seq = "AAGCAAATGGCNNAAGATGAGNNNAAGCAAATGGCNNAAGATGAGNNN";
    peak.setSeq("AAGCAAATGGCNNAAGATGAGNNNAAGCAAATGGCNNAAGATGAGNNN");
    result = peak.randomizePeakSeq();
    if (cstring::countCharacter(seq, "A") != cstring::countCharacter(result, "A")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count A 4" << endl;
    }
    if (cstring::countCharacter(seq, "G") != cstring::countCharacter(result, "G")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count G 4" << endl;
    }
    if (cstring::countCharacter(seq, "C") != cstring::countCharacter(result, "C")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count C 4" << endl;
    }
    if (cstring::countCharacter(seq, "T") != cstring::countCharacter(result, "T")) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test count T 4" << endl;
    }
    if (result.compare(11, 2, seq, 11, 2) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 7" << endl;
    }
    if (result.compare(22, 2, seq, 22, 2) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 8" << endl;
    }
    if (result.compare(35, 2, seq, 35, 2) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 9" << endl;
    }
    if (result.compare(45, 3, seq, 45, 3) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testRandomizePeakSeq (PeakTest) message=Wrong test 10" << endl;
    }
}

int main(int argc, char** argv) {
    clock_t start = clock();
    clock_t begin;
    std::cout << "%SUITE_STARTING% PeakTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    begin = clock();
    std::cout << "%TEST_STARTED% testRandomizePeakSeq (PeakTest)" << std::endl;
    testRandomizePeakSeq();
    std::cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testRandomizePeakSeq (PeakTest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << std::endl;

    delete Global::instance();
    delete TimeUtils::instance();
    return (EXIT_SUCCESS);
}


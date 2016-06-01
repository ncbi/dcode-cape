/* 
 * File:   PhyperTest.cpp
 * Author: veraalva
 *
 * Created on April 8, 2016, 11:29 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "bmath.h"

#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <set>

#include "TimeUtils.h"
#include "Global.h"

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

/*
 * Simple C++ Test Suite
 */

void testPhyperTest() {
    double x = 15;
    double m = 16885050;
    double n = 50655226;
    double k = 79;

    double pvalue = phyper(x, m, n, k, false, false);
    if (fabs(pvalue - 0.866574039066821) >= 10e-15) {
        std::cout << "%%TEST_FAILED%% time=0 testname=testPhyperTest (PhyperTest) message=Wrong pValue. 0.866574 != " << pvalue << std::endl;
    }
}

int main(int argc, char** argv) {
    clock_t start = clock();
    clock_t begin;
    std::cout << "%SUITE_STARTING% PhyperTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    begin = clock();
    std::cout << "%TEST_STARTED% testPhyperTest (PhyperTest)" << std::endl;
    testPhyperTest();
    std::cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testPhyperTest (PhyperTest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << std::endl;

    delete Global::instance();
    delete TimeUtils::instance();
    return (EXIT_SUCCESS);
}

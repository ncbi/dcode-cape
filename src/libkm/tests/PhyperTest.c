/*
 * File:   PhyperTest.c
 * Author: veraalva
 *
 * Created on February 18, 2016, 3:00 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "bmath.h"

/*
 * Simple C Test Suite
 */

void testPhyperTest() {
    double x = 15;
    double m = 16885050;
    double n = 50655226;
    double k = 79;

    double pvalue = phyper(x, m, n, k, false, false);
    if (fabs(pvalue - 0.866574039066821) >= 10e-15) {
        printf("%%TEST_FAILED%% time=0 testname=testPhyperTest (PhyperTest) message=Wrong pValue. 0.866574 != %f\n", pvalue);
    }
}

int main(int argc, char** argv) {
    printf("%%SUITE_STARTING%% PhyperTest\n");
    printf("%%SUITE_STARTED%%\n");

    printf("%%TEST_STARTED%% test1 (PhyperTest)\n");
    testPhyperTest();
    printf("%%TEST_FINISHED%% time=0 test1 (PhyperTest) \n");

    printf("%%SUITE_FINISHED%% time=0\n");

    return (EXIT_SUCCESS);
}

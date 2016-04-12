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
#include <map>
#include <set>
#include <ctime>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "FastaFactory.h"
#include "Global.h"
#include "KmersFactory.h"
#include "BedFactory.h"
#include "bmath.h"

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

/*
 * Simple C++ Test Suite
 */

void testRandomizePeakSeq() {
    peak::Peak peak;
    peak.setChr("chr1");
    char *seq = strdup("NNNAAGCAAATGGC");
    peak.setSeq(&seq);
    char* result = peak.randomizePeakSeq();   
    free(result); 
    free(seq);
    
    seq = strdup("NNNAAGCAAATGGCNNAAGATGAG");
    peak.setSeq(&seq);
    result = peak.randomizePeakSeq();
    free(result);
    free(seq);
    
    seq = strdup("AAGCAAATGGCNNAAGATGAGNNN");
    peak.setSeq(&seq);
    result = peak.randomizePeakSeq();
    free(result);
    free(seq);
    
    seq = strdup("AAGCAAATGGCNNAAGATGAGNNNAAGCAAATGGCNNAAGATGAGNNN");
    peak.setSeq(&seq);
    result = peak.randomizePeakSeq();
    free(result);
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
    return (EXIT_SUCCESS);
}


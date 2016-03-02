/* 
 * File:   PeakTest.cpp
 * Author: veraalva
 *
 * Created on February 23, 2016, 3:16 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
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
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "Global.h"
#include "KmersFactory.h"
#include "BedFactory.h"
#include "bmath.h"

/*
 * Simple C++ Test Suite
 */

void testRandomizePeakSeq() {
    peak::Peak peak;
    peak.SetChr("chr1");
    char *seq = strdup("NNNAAGCAAATGGC");
    peak.SetSeq(&seq);
    char* result = peak.RandomizePeakSeq();
    printf("%s\n%s\n",peak.GetSeq(), result);    
    free(result);    
    
    seq = strdup("NNNAAGCAAATGGCNNAAGATGAG");
    peak.SetSeq(&seq);
    result = peak.RandomizePeakSeq();
    printf("%s\n%s\n",peak.GetSeq(), result);    
    free(result);
    
    
    seq = strdup("AAGCAAATGGCNNAAGATGAGNNN");
    peak.SetSeq(&seq);
    result = peak.RandomizePeakSeq();
    printf("%s\n%s\n",peak.GetSeq(), result);    
    free(result);
    
    seq = strdup("AAGCAAATGGCNNAAGATGAGNNNAAGCAAATGGCNNAAGATGAGNNN");
    peak.SetSeq(&seq);
    result = peak.RandomizePeakSeq();
    printf("%s\n%s\n",peak.GetSeq(), result);    
    free(result);
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% PeakTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% testRandomizePeakSeq (PeakTest)" << std::endl;
    testRandomizePeakSeq();
    std::cout << "%TEST_FINISHED% time=0 testRandomizePeakSeq (PeakTest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}


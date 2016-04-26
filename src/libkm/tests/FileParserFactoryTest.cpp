/* 
 * File:   FileParserFactoryTest.cpp
 * Author: veraalva
 *
 * Created on April 11, 2016, 1:55 PM
 */

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <map>
#include <ctime>
#include <vector>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"

#include "Global.h"
#include "TimeUtils.h"
#include "FileParserFactory.h"

using namespace std;
using namespace parsers;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

/*
 * Simple C++ Test Suite
 */

void testIterate() {
    parsers::FileParserFactory fileParserFactory;
    fileParserFactory.setFileToParse("resources/57epigenomes.RPKM.pc");
    int i = 0;
    try {
        while (fileParserFactory.iterate('#', "\t")) {
            if (fileParserFactory.getNWords() != 58) {
                cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                break;
            }
            if (i == 0) {
                if (strncmp(fileParserFactory.getWords()[0], "gene_id", 7) != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
                if (strncmp(fileParserFactory.getWords()[57], "E128", 4) != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
            } else if (i == 1) {
                if (strncmp(fileParserFactory.getWords()[0], "ENSG00000000003", 15) != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
                if (strncmp(fileParserFactory.getWords()[57], "15.818", 6) != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
            } else if (i == 19794) {
                if (strncmp(fileParserFactory.getWords()[0], "ENSG00000259765", 15) != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
                if (strncmp(fileParserFactory.getWords()[16], "0.129", 5) != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
            } else if (i == 19795) {
                if (strncmp(fileParserFactory.getWords()[0], "ENSG00000259766", 15) != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
                if (strncmp(fileParserFactory.getWords()[57], "0.000", 5) != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
            }

            i++;
        }
    } catch (exceptions::FileNotFoundException ex) {
        cerr << ex.what() << endl;
        cerr << "Error parsing file" << endl;
        exit(-1);
    } catch (exceptions::ErrorReadingFromFileException ex) {
        cerr << ex.what() << endl;
        cerr << "Error parsing file" << endl;
        exit(-1);
    }
    if (i != 19796) {
        cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file. Read lines " << i << endl;
    }

}

int main() {
    clock_t start = clock();
    clock_t begin;
    cout << "%SUITE_STARTING% FileParserFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    begin = clock();
    cout << "%TEST_STARTED% testIterate (FileParserFactoryTest)" << endl;
    testIterate();
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testIterate (FileParserFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << endl;

    delete Global::instance();
    delete TimeUtils::instance();
    return (EXIT_SUCCESS);
}


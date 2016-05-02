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
    int i = 0;
    try {
        fileParserFactory.setFileToParse("resources/57epigenomes.RPKM.pc");
        while (fileParserFactory.iterate("#", "\t")) {
            if (fileParserFactory.getWords().size() != 58) {
                cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                break;
            }
            if (i == 0) {
                if (fileParserFactory.getWords()[0].compare("gene_id") != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
                if (fileParserFactory.getWords()[57].compare("E128") != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
            } else if (i == 1) {
                if (fileParserFactory.getWords()[0].compare("ENSG00000000003") != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
                if (fileParserFactory.getWords()[57].compare("15.818") != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
            } else if (i == 19794) {
                if (fileParserFactory.getWords()[0].compare("ENSG00000259765") != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
                if (fileParserFactory.getWords()[16].compare("0.129") != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
            } else if (i == 19795) {
                if (fileParserFactory.getWords()[0].compare("ENSG00000259766") != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
                if (fileParserFactory.getWords()[57].compare("0.000") != 0) {
                    cout << "word: " << fileParserFactory.getWords()[0] << endl;
                    cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file" << endl;
                }
            }

            i++;
        }
    } catch (exceptions::FileNotFoundException ex) {
        cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file: resources/57epigenomes.RPKM.pc " << endl;
        return;
    } catch (ios::failure ex) {
        cout << "%TEST_FAILED% time=0 testname=testIterate (FileParserFactoryTest) message=Error parsing file: resources/57epigenomes.RPKM.pc " << endl;
        return;
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


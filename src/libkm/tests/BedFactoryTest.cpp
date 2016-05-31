/* 
 * File:   BedFactoryTest.cpp
 * Author: veraalva
 *
 * Created on February 11, 2016, 4:14 PM
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <ctime>

#include "Global.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "FileParserFactory.h"
#include "FastaFactory.h"
#include "KmersFactory.h"
#include "BedFactory.h"
#include "cstring.h"

using namespace std;
using namespace parsers;
using namespace kmers;
using namespace sequence;
using namespace peak;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

/*
 * Simple C++ Test Suite
 */

void testCreatePeaksFromBedFile(string bFName, string dirName, string testName) {
    int i = 0;
    sequence::FastaFactory chrFactory;
    peak::BedFactory bedFactory;
    FileParserFactory fParser;
    string prefix("chr");
    string sufix(".fa.masked");
    map<unsigned long int, map<string, unsigned long int>> testMap;
    KmersFactory kmersFactory;
    kmersFactory.createGenomeWideKmers();
    vector<string> w;

    string test("TGAGAGGAAAGCTTTCCCACATTATACAGCTTCTGAAAGGGTTGCTTGACCCACAGATGTGAAGCTGAGGCTGAAGGAGACTGATGTGGTTTCTCCTCAGTTTCTCTGTGCAGCACCAGGTGGCAGCAGAGGTCAGCAAGGCAAACCCGAGCCCGGGGATGCGGAGTGGGGGCAGCTACGTCCTCTCTTG");

    try {
        fParser.setFileToParse(testName);
        while (fParser.iterate("#", "\t")) {
            map<string, unsigned long int> inMap;
            cstring::split(fParser.getWords()[2], "-", w);
            i = atoi(w[0].c_str());
            inMap.insert(make_pair("start", i));
            inMap.insert(make_pair("end", atoi(w[1].c_str())));
            inMap.insert(make_pair("width", atoi((fParser.getWords()[3]).c_str())));
            inMap.insert(make_pair("gccount", atoi((fParser.getWords()[4]).c_str())));
            inMap.insert(make_pair("ncount", atoi((fParser.getWords()[5]).c_str())));
            inMap.insert(make_pair("nrcount", atoi((fParser.getWords()[6]).c_str())));

            testMap.insert(pair<unsigned long int, map<string, unsigned long int>>(i, inMap));
        }
    } catch (exceptions::FileNotFoundException) {
        cerr << "Error parsing file: " << testName << endl;
        exit(-1);
    } catch (ios::failure) {
        cerr << "Error parsing file: " << testName << endl;
        exit(-1);
    }
    
    chrFactory.parseFastaInDirectory(dirName, prefix, sufix, false);
    bedFactory.createPeaksFromBedFile(chrFactory, bFName, 0.7, kmersFactory);
    if (bedFactory.getPeaks().size() != 88) {
        cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Peaks to print should be equal to 88 and it is " << i << endl;
    }

    for (auto it = bedFactory.getPeaks().begin(); it != bedFactory.getPeaks().end(); ++it) {
        shared_ptr<Peak> p = *it;
        
        map<string, unsigned long int> inMap = testMap.find(p->getStart())->second;
        
        if (inMap.find("width")->second != p->getLength()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong width: "
                    << p->getLength() << "!=" << inMap.find("width")->second << endl;
        }
        if (inMap.find("gccount")->second != p->getGCCount()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong GCCount: "
                    << p->getGCCount() << "!=" << inMap.find("gccount")->second << endl;
        }
        if (inMap.find("ncount")->second != p->getNCount()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong NCount: "
                    << p->getNCount() << "!=" << inMap.find("ncount")->second << endl;
        }
        if (inMap.find("nrcount")->second != p->getNRCount()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong NRCount: "
                    << p->getNRCount() << "!=" << inMap.find("nrcount")->second << endl;
        }
        if (p->getSeq().size() != p->getLength()) {
            cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Sequence length not equal to width" << endl;
        }
        if (p->getStart() == 237640 && p->getEnd() == 237829) {
            if (p->getSeq().compare(test) != 0) {
                cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=First peak is not equal to the test" << endl;
            }
        }
    }

    bedFactory.generatingControlsFromChromosomes(chrFactory, 3, kmersFactory);
    kmersFactory.buildKmers();
    if (kmersFactory.getKmers().size() != 12122) {
        cout << "%TEST_FAILED% time=0 testname=testCreatePeaksFromBedFile (BedFactoryTest) message=Wrong number of kmers 12122 != " << kmersFactory.getKmers().size() << endl;
    }
}

void testReadingControlsFromFile(string dirName, string controlName) {
    sequence::FastaFactory chrFactory;
    peak::BedFactory bedFactory;
    string prefix("chr");
    string sufix(".fa.masked");
    KmersFactory kmersFactory;
    kmersFactory.createGenomeWideKmers();

    chrFactory.parseFastaInDirectory(dirName, prefix, sufix, false);

    bedFactory.readControlsFromFile(controlName, chrFactory, kmersFactory);

    if (kmersFactory.getTotalNRNTControl() != 1065) {
        cout << "%TEST_FAILED% time=0 testname=testReadingControlsFromFile (BedFactoryTest) message=Wrong calculation of Total non N 1065 != " << kmersFactory.getTotalNRNTControl() << endl;
    }
}

int main() {
    clock_t start = clock();
    clock_t begin;
    string bFName("resources/test.bed");
    string dirName("resources/");
    string testName("resources/test.train.control.info.chr1");
    string controlName("resources/randControl.txt");
    cout << "%SUITE_STARTING% BedFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    Global::instance()->setVerbose(0);
    Global::instance()->setOrder(10);
    Global::instance()->setBin1(0.005);
    Global::instance()->setBin2(0.01);

    begin = clock();
    cout << "%TEST_STARTED% testCreatePeaksFromBedFile (BedFactoryTest)" << endl;
    testCreatePeaksFromBedFile(bFName, dirName, testName);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " seconds testCreatePeaksFromBedFile (BedFactoryTest)" << endl;

    begin = clock();
    cout << "%TEST_STARTED% testReadingControlsFromFile (BedFactoryTest)" << endl;
    testReadingControlsFromFile(dirName, controlName);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testReadingControlsFromFile (BedFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << endl;

    delete Global::instance();
    delete TimeUtils::instance();
    return (EXIT_SUCCESS);
}


/* 
 * File:   FastaFactoryTest.cpp
 * Author: veraalva
 *
 * Created on February 11, 2016, 9:07 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <time.h>
#include <dirent.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <ctime>
#include <set>

#include "berror.h"
#include "bmemory.h"
#include "Global.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "FastaFactory.h"

using namespace std;
using namespace sequence;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

void testParseFastaFile(string fileName, bool binary) {
    shared_ptr<Seq> fasta;
    FastaFactory fastaFactory;
    std::string seg;
    std::string test_seg("ATCCGACATCAAGTGCCCACCTTGGCTCGTGGCTCTCACTGCAACGGGAAAGCCACAGACTGGGGTGAAGAGTTCAGTCACATGCGACCGGTGACTCCCTGTCCCCACCCCCATGACACTCCCCAGCCCTCCAAGGCCACTGTGTTTCCCAGTTAGCTCAGAGCCTCAGTCGATCCCTGACCCAGCACCGGGCACTGATGAGACAGCGGCTGTTTGAGGAGCCACCTCCCAGCCACCTCGGGGCCAGGGCCAGGGTGTGCAGCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAACCATAGTGCCCAGGGCACTGCCGCTGCAGGCGCAGGCATCGCATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACCATTGCCTTGGACCGCTCTTGGCAGTCGAAGAAGATTCTCCTGTCAGTTTGAGCTGGGTGAGCTTAGAGAGGAAAGCTCCACTATGGCTCCCAAACCAGGAAGGAGCCATAGCCCAGGCAGGAGGGCTGAGGACCTCTGGTGGCGGCCCAGGGCTTCCAGCATGTGCCCTAGGGGAAGCAGGGGCCAGCTGGCAAGAGCAGGGGGTGGGCAGAAAGCACCCGGTGGACTCAGGGCTGGAGGGGAGGAGGCGATCTTGCCCAAGGCCCTCCGACTGCAAGCTCCAGGGCCCGCTCACCTTGCTCCTGCTCCTTCTGCTGCTGCTTCTCCAGCTTTCGCTCCTTCATGCTGCGCAGCTTGGCCTTGCCGATGCCCCCAGCTTGGCGGATGGACTCTAGCAGAGTGGCCAGCCACCGGAGGGGTCAACCACTTCCCTGGGAGCTCCCTGGACTGGAGCCGGGAGGTGGGGAACAGGGCAAGGAGGAAAGGCTGCTCAGGCAGGG");

    long unsigned int seqRead = fastaFactory.parseFastaFile(fileName, binary);

    if (fastaFactory.getSequenceContainter().size() != 1 || seqRead != 1) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=It should read 1 sequence" << endl;
    }
    fasta = fastaFactory.getSequenceContainter().begin()->second;
    if (fasta->getId().compare("chr1") != 0) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The sequence Id should be chr1 != " << fasta->getId() << endl;
    }
    if (fasta->getLength() != 1249950) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The sequence should have 1249950 bp and it is " << fasta->getLength() << endl;
    }

    if (fasta->getSubStr(15000, 1000).compare(test_seg) != 0) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The segment of length 1050 starting at 15000 is not equal to the test" << endl;
    }

    //remove(fileName.c_str());
}

void testParseFastaMultipleFile(string fileName) {
    long unsigned int seqRead;
    FastaFactory fFactory;
    shared_ptr<Seq> fasta;
    shared_ptr<Seq> f;

    map<string, long unsigned int> seqSize;

    seqSize.insert(pair<string, long unsigned int>("chr1", 950));
    seqSize.insert(pair<string, long unsigned int>("chr2", 700));
    seqSize.insert(pair<string, long unsigned int>("chr3", 1200));
    seqSize.insert(pair<string, long unsigned int>("chr4", 1450));
    seqSize.insert(pair<string, long unsigned int>("chr5", 1700));

    seqRead = fFactory.parseFastaFile(fileName, false);
    if (fFactory.getSequenceContainter().size() != 5 && seqRead != 5) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=It should read 5 sequence" << endl;
    }
    for (auto it = fFactory.getSequenceContainter().begin(); it != fFactory.getSequenceContainter().end(); ++it) {
        fasta = it->second;
        if (seqSize.find(fasta->getId())->second != fasta->getLength()) {
            cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Wrong sequence size (test 1) "
                    << seqSize.find(fasta->getId())->second << " != " << fasta->getLength() << endl;
        }
    }

    try {
        f = fFactory.getSequenceFromID("chr4");
        if (!f) {
            cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Cant find chr4 in the factory container " << endl;
        } else {
            if (seqSize.find(f->getId())->second != f->getLength()) {
                cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Wrong sequence size (test 2) "
                        << seqSize.find(f->getId())->second << " != " << f->getLength() << endl;
            }
        }

        f = fFactory.getSequenceFromID("chr5");
        if (!f) {
            cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Cant find chr4 in the factory container " << endl;
        } else {
            if (seqSize.find(f->getId())->second != f->getLength() || f->getLength() != 1700) {
                cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Wrong sequence size (test 3)"
                        << seqSize.find(f->getId())->second << " != " << f->getLength() << endl;
            }
        }
    } catch (exceptions::NotFoundException) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaMultipleFile (FastaFactoryTest) message=Not found sequence chr4" << endl;
    }
}

void testParseFastaDir(string fileName) {
    FastaFactory fFactory;
    string prefix("chr");
    string suxif(".fa.masked");
    std::string test_seg("ATCCGACATCAAGTGCCCACCTTGGCTCGTGGCTCTCACTGCAACGGGAAAGCCACAGACTGGGGTGAAGAGTTCAGTCACATGCGACCGGTGACTCCCTGTCCCCACCCCCATGACACTCCCCAGCCCTCCAAGGCCACTGTGTTTCCCAGTTAGCTCAGAGCCTCAGTCGATCCCTGACCCAGCACCGGGCACTGATGAGACAGCGGCTGTTTGAGGAGCCACCTCCCAGCCACCTCGGGGCCAGGGCCAGGGTGTGCAGCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAACCATAGTGCCCAGGGCACTGCCGCTGCAGGCGCAGGCATCGCATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACCATTGCCTTGGACCGCTCTTGGCAGTCGAAGAAGATTCTCCTGTCAGTTTGAGCTGGGTGAGCTTAGAGAGGAAAGCTCCACTATGGCTCCCAAACCAGGAAGGAGCCATAGCCCAGGCAGGAGGGCTGAGGACCTCTGGTGGCGGCCCAGGGCTTCCAGCATGTGCCCTAGGGGAAGCAGGGGCCAGCTGGCAAGAGCAGGGGGTGGGCAGAAAGCACCCGGTGGACTCAGGGCTGGAGGGGAGGAGGCGATCTTGCCCAAGGCCCTCCGACTGCAAGCTCCAGGGCCCGCTCACCTTGCTCCTGCTCCTTCTGCTGCTGCTTCTCCAGCTTTCGCTCCTTCATGCTGCGCAGCTTGGCCTTGCCGATGCCCCCAGCTTGGCGGATGGACTCTAGCAGAGTGGCCAGCCACCGGAGGGGTCAACCACTTCCCTGGGAGCTCCCTGGACTGGAGCCGGGAGGTGGGGAACAGGGCAAGGAGGAAAGGCTGCTCAGGCAGGG");

    try {
        fFactory.parseFastaInDirectory(fileName, prefix, suxif, false);
        shared_ptr<Seq> seq = fFactory.getSequenceFromID("chr1");
        if (seq->getId().compare("chr1") != 0) {
            cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The sequence Id should be chr1" << endl;
        }
        if (seq->getLength() != 1249950) {
            cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The sequence should have 1249950 bp and it is " << seq->getLength() << endl;
        }

        if (seq->getSubStr(15000, 1000).compare(test_seg) != 0) {
            cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The segment of length 1050 starting at 15000 is not equal to the test" << endl;
        }
    } catch (exceptions::NotFoundException) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=Not found sequence chr1" << endl;
    }


}

void testWriteSequencesToFile(std::string inFileName, std::string outFileName, bool binary) {
    FastaFactory fFactory;
    fFactory.parseFastaFile(inFileName, false);
    if (fFactory.getSequenceContainter().size() != 1) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=It should read 1 sequence" << endl;
    }
    shared_ptr<Seq> fasta = fFactory.getSequenceContainter().begin()->second;
    if (fasta->getId().compare("chr1") != 0) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The sequence Id should be chr1 != " << fasta->getId() << endl;
    }
    if (fasta->getLength() != 1249950) {
        cout << "%TEST_FAILED% time=0 testname=testParseFastaFile (FastaFactoryTest) message=The sequence should have 1249950 bp and it is " << fasta->getLength() << endl;
    }
    fFactory.writeSequencesToFile(outFileName, binary);
}

int main() {
    clock_t start = clock();
    clock_t begin;
    string singleFileName("resources/chr1.fa.masked");
    string outFileName("resources/chr1.fa.masked.w");
    string multipleFileName("resources/multiple.fa.masked");
    string dirName("resources/");
    cout << "%SUITE_STARTING% FastaFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    Global::instance()->setVerbose(0);

    begin = clock();
    cout << "%TEST_STARTED% testWriteSequencesToFile (FastaFactoryTest)" << endl;
    testWriteSequencesToFile(singleFileName, outFileName, true);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testWriteSequencesToFile (FastaFactoryTest)" << endl;

    begin = clock();
    cout << "%TEST_STARTED% testParseFastaFile (FastaFactoryTest)" << endl;
    testParseFastaFile(outFileName, true);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testParseFastaFile (FastaFactoryTest)" << endl;

    begin = clock();
    cout << "%TEST_STARTED% testParseFastaMultipleFile (FastaFactoryTest)" << endl;
    testParseFastaMultipleFile(multipleFileName);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testParseFastaMultipleFile (FastaFactoryTest)" << endl;

    begin = clock();
    cout << "%TEST_STARTED% testParseFastaDir (FastaFactoryTest)" << endl;
    testParseFastaDir(dirName);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testParseFastaDir (FastaFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << endl;

    delete Global::instance();
    delete TimeUtils::instance();
    return (EXIT_SUCCESS);
}


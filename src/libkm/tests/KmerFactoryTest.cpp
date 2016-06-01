/* 
 * File:   KmerFactoryTest.cpp
 * Author: veraalva
 *
 * Created on February 17, 2016, 10:44 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <math.h>

#include <iostream>
#include <cmath>
#include <memory>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <algorithm>
#include <ctime>

#include "berror.h"
#include "bmemory.h"
#include "TimeUtils.h"
#include "Exceptions.h"
#include "FastaFactory.h"
#include "Global.h"
#include "cstring.h"
#include "KmersFactory.h"
#include "BedFactory.h"
#include "bmath.h"

using namespace std;
using namespace kmers;
using namespace sequence;
using namespace peak;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

/*
 * Simple C++ Test Suite
 */

void testKmers() {
    KmersFactory kmersFactory;
    kmersFactory.createGenomeWideKmers();
    if (kmersFactory.getKmersGenome().size() != 524800) {
        cout << "%TEST_FAILED% time=0 testname=testKmers (KmerFactoryTest) message=It should generate 524800 kmers not " << kmersFactory.getKmersGenome().size() << endl;
    }
}

void testKmersPValue() {
    double totalNRnt_peak = 16885050;
    double totalNRnt_control = 50655226;

    Kmer kmer;
    kmer.setPeakFreq(16);
    kmer.setControlFreq(63);
    kmer.calculatePValue(totalNRnt_peak, totalNRnt_control);
    if (std::fabs(kmer.getValue() - 0.866574039066821) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pValue. 0.866574 != " << kmer.getValue() << endl;
    }
    if (std::fabs(kmer.getSig() - 0.0621943257550334) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong sig. 0.0621943257550334 != " << kmer.getSig() << endl;
    }
    if (std::fabs(kmer.getPf() - 0.761905905021738) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pf. 0.761905905021738 != " << kmer.getPf() << endl;
    }

    kmer.setPeakFreq(30);
    kmer.setControlFreq(123);
    kmer.calculatePValue(totalNRnt_peak, totalNRnt_control);
    if (std::fabs(kmer.getValue() - 0.951928135758493) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pValue. 0.951928 != " << kmer.getValue() << endl;
    }
    if (std::fabs(kmer.getSig() - 0.0213958367222335) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong sig. 0.0213958367222335 != " << kmer.getSig() << endl;
    }
    if (std::fabs(kmer.getPf() - 0.731708414883682) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pf. 0.731708414883682 != " << kmer.getPf() << endl;
    }

    kmer.setPeakFreq(3);
    kmer.setControlFreq(13);
    kmer.calculatePValue(totalNRnt_peak, totalNRnt_control);
    if (std::fabs(kmer.getValue() - 0.80288828067583) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pValue. 0.802888 != " << kmer.getValue() << endl;
    }
    if (std::fabs(kmer.getSig() - 0.0953448811988908) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong sig. 0.0953448811988908 != " << kmer.getSig() << endl;
    }
    if (std::fabs(kmer.getPf() - 0.69230873100533) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pf. 0.69230873100533 != " << kmer.getPf() << endl;
    }

    kmer.setPeakFreq(40);
    kmer.setControlFreq(152);
    kmer.calculatePValue(totalNRnt_peak, totalNRnt_control);
    if (std::fabs(kmer.getValue() - 0.924001786259078) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pValue. 0.924001 != " << kmer.getValue() << endl;
    }
    if (std::fabs(kmer.getSig() - 0.0343271892109419) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong sig. 0.0343271892109419 != " << kmer.getSig() << endl;
    }
    if (std::fabs(kmer.getPf() - 0.789474868690288) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pf. 0.789474868690288 != " << kmer.getPf() << endl;
    }
}

void testWriteKmersToFile() {
    KmersFactory kmersFactory;
    shared_ptr<Kmer> k;

    //            AAAAAAAAAA 3.25573332278254e-21 20.48735117540051 1.49479032886669
    k = make_shared<Kmer>();
    k->setControlFreq(11);
    k->setNegativeControl(12);
    k->setNegativePeak(13);
    k->setPeakFreq(14);
    k->setValue(3.25573332278254e-21);
    k->setSig(20.48735117540051);
    k->setPf(1.49479032886669);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAAAA", k));

    //            AAAAAAAAAC 3.66936367681041e-05 4.43540924245372 1.62318695308535
    k = make_shared<Kmer>();
    k->setControlFreq(21);
    k->setNegativeControl(22);
    k->setNegativePeak(23);
    k->setPeakFreq(24);
    k->setValue(3.66936367681041e-05);
    k->setSig(4.43540924245372);
    k->setPf(1.62318695308535);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAAAC", k));

    //            AAAAAAAAAG 1.56360980495561e-03 2.80587161489318 1.33064397039723
    k = make_shared<Kmer>();
    k->setControlFreq(31);
    k->setNegativeControl(32);
    k->setNegativePeak(33);
    k->setPeakFreq(34);
    k->setValue(1.56360980495561e-03);
    k->setSig(2.80587161489318);
    k->setPf(1.33064397039723);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAAAG", k));

    //            AAAAAAAAAT 2.19631497914574e-01 0.65830537647355 1.08933620086384
    k = make_shared<Kmer>();
    k->setControlFreq(41);
    k->setNegativeControl(42);
    k->setNegativePeak(43);
    k->setPeakFreq(44);
    k->setValue(2.19631497914574e-01);
    k->setSig(0.65830537647355);
    k->setPf(1.08933620086384);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAAAT", k));

    //            AAAAAAAACA 1.07873928249526e-02 1.96708350606356 1.41176344238936
    k = make_shared<Kmer>();
    k->setControlFreq(51);
    k->setNegativeControl(52);
    k->setNegativePeak(53);
    k->setPeakFreq(54);
    k->setValue(1.07873928249526e-02);
    k->setSig(1.96708350606356);
    k->setPf(1.41176344238936);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAACA", k));

    //            AAAAAAAACC 1.05511751710985e-03 2.97669916672645 1.74193392488365
    k = make_shared<Kmer>();
    k->setControlFreq(61);
    k->setNegativeControl(62);
    k->setNegativePeak(63);
    k->setPeakFreq(64);
    k->setValue(1.05511751710985e-03);
    k->setSig(2.97669916672645);
    k->setPf(1.74193392488365);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAACC", k));

    //            AAAAAAAACG 3.26401549435001e-01 0.48624778830827 1.28571313503317
    k = make_shared<Kmer>();
    k->setControlFreq(71);
    k->setNegativeControl(72);
    k->setNegativePeak(73);
    k->setPeakFreq(74);
    k->setValue(3.26401549435001e-01);
    k->setSig(0.48624778830827);
    k->setPf(1.28571313503317);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAACG", k));

    //            AAAAAAAACT 1.48116148621570e-02 1.82939758918367 1.64515981794567
    k = make_shared<Kmer>();
    k->setControlFreq(81);
    k->setNegativeControl(82);
    k->setNegativePeak(83);
    k->setPeakFreq(84);
    k->setValue(1.48116148621570e-02);
    k->setSig(1.82939758918367);
    k->setPf(1.64515981794567);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAACT", k));

    //            AAAAAAAAGA 2.70816098140979e-01 0.56732552342117 1.09793716170358
    k = make_shared<Kmer>();
    k->setControlFreq(91);
    k->setNegativeControl(92);
    k->setNegativePeak(93);
    k->setPeakFreq(94);
    k->setValue(2.70816098140979e-01);
    k->setSig(0.56732552342117);
    k->setPf(1.09793716170358);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAAGA", k));
    //            AAAAAAAAGC 1.88687646220172e-03 2.72425653308383 1.63888742212561
    k = make_shared<Kmer>();
    k->setControlFreq(101);
    k->setNegativeControl(102);
    k->setNegativePeak(103);
    k->setPeakFreq(104);
    k->setValue(1.88687646220172e-03);
    //k->SetSig(2.72425653308383);
    k->setSig(INFINITY);
    k->setPf(1.63888742212561);

    kmersFactory.getKmers().insert(make_pair("AAAAAAAAGC", k));

    kmersFactory.writeKmersToFile("/tmp/kmers.bin", true);
}

void testReadKmersFromFile(std::string fileName, bool binary) {
    KmersFactory kmersFactory;
    string rc_kmer;

    kmersFactory.readKmersFromFile(fileName, binary);

    if (std::fabs(kmersFactory.getKmerSig("AAAAAAAAAA") - 20.48735117540051) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig. for kmer  != AAAAAAAAAA " << kmersFactory.getKmerSig("AAAAAAAAAA") << endl;
    }

    if (std::fabs(kmersFactory.getKmerSig("AAAAAAAAGC") - (20.48735117540051 + 10)) >= 10E-15) {
        cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig. for kmer  != AAAAAAAAGC " << kmersFactory.getKmerSig("AAAAAAAAGC") << endl;
    }

    for (auto it = kmersFactory.getKmers().begin(); it != kmersFactory.getKmers().end(); ++it) {
        bool right = false;
        string kmer = it->first;
        shared_ptr<Kmer> k = it->second;

        rc_kmer = cstring::reverseComplement(kmer);

        //                AAAAAAAAAA 3.25573332278254e-21 20.48735117540051 1.49479032886669
        if (kmer.compare("AAAAAAAAAA") == 0 || rc_kmer.compare("AAAAAAAAAA") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 3.25573332278254e-21) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 20.48735117540051) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.49479032886669) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 11) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 12) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 13) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 14) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }

        //                AAAAAAAAAC 3.66936367681041e-05 4.43540924245372 1.62318695308535
        if (kmer.compare("AAAAAAAAAC") == 0 || rc_kmer.compare("AAAAAAAAAC") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 3.66936367681041e-05) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 4.43540924245372) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.62318695308535) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 21) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 22) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 23) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 24) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }
        //                AAAAAAAAAG 1.56360980495561e-03 2.80587161489318 1.33064397039723
        if (kmer.compare("AAAAAAAAAG") == 0 || rc_kmer.compare("AAAAAAAAAG") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 1.56360980495561e-03) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 2.80587161489318) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.33064397039723) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 31) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 32) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 33) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 34) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }
        //                AAAAAAAAAT 2.19631497914574e-01 0.65830537647355 1.08933620086384
        if (kmer.compare("AAAAAAAAAT") == 0 || rc_kmer.compare("AAAAAAAAAT") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 2.19631497914574e-01) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 0.65830537647355) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.08933620086384) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 41) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 42) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 43) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 44) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }
        //                AAAAAAAACA 1.07873928249526e-02 1.96708350606356 1.41176344238936
        if (kmer.compare("AAAAAAAACA") == 0 || rc_kmer.compare("AAAAAAAACA") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 1.07873928249526e-02) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 1.96708350606356) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.41176344238936) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 51) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 52) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 53) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 54) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }
        //                AAAAAAAACC 1.05511751710985e-03 2.97669916672645 1.74193392488365
        if (kmer.compare("AAAAAAAACC") == 0 || rc_kmer.compare("AAAAAAAACC") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 1.05511751710985e-03) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 2.97669916672645) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.74193392488365) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 61) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 62) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 63) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 64) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }
        //                AAAAAAAACG 3.26401549435001e-01 0.48624778830827 1.28571313503317
        if (kmer.compare("AAAAAAAACG") == 0 || rc_kmer.compare("AAAAAAAACG") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 3.26401549435001e-01) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 0.48624778830827) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.28571313503317) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 71) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 72) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 73) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 74) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }
        //                AAAAAAAACT 1.48116148621570e-02 1.82939758918367 1.64515981794567
        if (kmer.compare("AAAAAAAACT") == 0 || rc_kmer.compare("AAAAAAAACT") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 1.48116148621570e-02) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 1.82939758918367) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.64515981794567) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 81) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 82) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 83) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 84) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }
        //                AAAAAAAAGA 2.70816098140979e-01 0.56732552342117 1.09793716170358
        if (kmer.compare("AAAAAAAAGA") == 0 || rc_kmer.compare("AAAAAAAAGA") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 2.70816098140979e-01) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - 0.56732552342117) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.09793716170358) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 91) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 92) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 93) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 94) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }
        //                AAAAAAAAGC 1.88687646220172e-03 2.72425653308383 1.63888742212561
        if (kmer.compare("AAAAAAAAGC") == 0 || rc_kmer.compare("AAAAAAAAGC") == 0) {
            right = true;
            if (std::fabs(k->getValue() - 1.88687646220172e-03) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->getValue() << endl;
            }
            if (std::fabs(k->getSig() - (20.48735117540051 + 10)) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->getSig() << endl;
            }
            if (std::fabs(k->getPf() - 1.63888742212561) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->getPf() << endl;
            }
            if (k->getControlFreq() != 101) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->getControlFreq() << endl;
            }
            if (k->getNegativeControl() != 102) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->getNegativeControl() << endl;
            }
            if (k->getNegativePeak() != 103) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->getNegativePeak() << endl;
            }
            if (k->getPeakFreq() != 104) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->getPeakFreq() << endl;
            }
        }

        if (!right) {
            cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong Peak.  != " << kmer << endl;
        }
    }

    remove(fileName.c_str());
}

int main(int argc, char** argv) {
    clock_t start = clock();
    clock_t begin;
    cout << "%SUITE_STARTING% KmerFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    set<int> orders;

    orders.insert(4);
    orders.insert(6);
    orders.insert(8);
//    orders.insert(10);

    Global::instance()->setVerbose(0);
    Global::instance()->setOrder(10);
    Global::instance()->setOrders(orders);
    Global::instance()->setBin1(0.005);
    Global::instance()->setBin2(0.01);

    begin = clock();
    cout << "%TEST_STARTED% testKmers (KmerFactoryTest)" << endl;
    testKmers();
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testKmers (KmerFactoryTest)" << endl;

    begin = clock();
    cout << "%TEST_STARTED% testKmersPValue (KmerFactoryTest)" << endl;
    testKmersPValue();
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testKmersPValue (KmerFactoryTest)" << endl;

    begin = clock();
    cout << "%TEST_STARTED% testWriteKmersToFile (KmerFactoryTest)" << endl;
    testWriteKmersToFile();
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testWriteKmersToFile (KmerFactoryTest)" << endl;

    begin = clock();
    cout << "%TEST_STARTED% testReadKmersFromFile (KmerFactoryTest)" << endl;
    testReadKmersFromFile("/tmp/kmers.bin", true);
    cout << "%TEST_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(begin) << " second testReadKmersFromFile (KmerFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=" << TimeUtils::instance()->getTimeSecFrom(start) << " seconds" << endl;

    delete Global::instance();
    delete TimeUtils::instance();
    return (EXIT_SUCCESS);
}


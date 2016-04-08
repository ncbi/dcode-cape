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
#include <map>
#include <set>
#include <algorithm>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"
#include "TimeUtils.h"
#include "FastaFactory.h"
#include "Global.h"
#include "KmersFactory.h"
#include "BedFactory.h"
#include "bmath.h"

using namespace std;
using namespace kmers;
using namespace fasta;
using namespace peak;

Global *Global::s_instance = 0;
TimeUtils *TimeUtils::s_instance = 0;

/*
 * Simple C++ Test Suite
 */

void testKmers() {
    KmersFactory kmersFactory;
    kmersFactory.CreateGenomeWideKmers();
    if (kmersFactory.GetKmersGenome().size() != 524800) {
        cout << "%TEST_FAILED% time=0 testname=testKmers (KmerFactoryTest) message=It should generate 524800 kmers not " << kmersFactory.GetKmersGenome().size() << endl;
    }
}

void testKmersPValue() {
    double totalNRnt_peak = 16885050;
    double totalNRnt_control = 50655226;

    Kmer kmer;
    kmer.SetPeakFreq(16);
    kmer.SetControlFreq(63);
    kmer.CalculatePValue(totalNRnt_peak, totalNRnt_control);
    if (std::fabs(kmer.GetValue() - 0.866574039066821) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pValue. 0.866574 != " << kmer.GetValue() << endl;
    }
    if (std::fabs(kmer.GetSig() - 0.0621943257550334) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong sig. 0.0621943257550334 != " << kmer.GetSig() << endl;
    }
    if (std::fabs(kmer.GetPf() - 0.761905905021738) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pf. 0.761905905021738 != " << kmer.GetPf() << endl;
    }

    kmer.SetPeakFreq(30);
    kmer.SetControlFreq(123);
    kmer.CalculatePValue(totalNRnt_peak, totalNRnt_control);
    if (std::fabs(kmer.GetValue() - 0.951928135758493) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pValue. 0.951928 != " << kmer.GetValue() << endl;
    }
    if (std::fabs(kmer.GetSig() - 0.0213958367222335) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong sig. 0.0213958367222335 != " << kmer.GetSig() << endl;
    }
    if (std::fabs(kmer.GetPf() - 0.731708414883682) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pf. 0.731708414883682 != " << kmer.GetPf() << endl;
    }

    kmer.SetPeakFreq(3);
    kmer.SetControlFreq(13);
    kmer.CalculatePValue(totalNRnt_peak, totalNRnt_control);
    if (std::fabs(kmer.GetValue() - 0.80288828067583) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pValue. 0.802888 != " << kmer.GetValue() << endl;
    }
    if (std::fabs(kmer.GetSig() - 0.0953448811988908) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong sig. 0.0953448811988908 != " << kmer.GetSig() << endl;
    }
    if (std::fabs(kmer.GetPf() - 0.69230873100533) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pf. 0.69230873100533 != " << kmer.GetPf() << endl;
    }

    kmer.SetPeakFreq(40);
    kmer.SetControlFreq(152);
    kmer.CalculatePValue(totalNRnt_peak, totalNRnt_control);
    if (std::fabs(kmer.GetValue() - 0.924001786259078) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pValue. 0.924001 != " << kmer.GetValue() << endl;
    }
    if (std::fabs(kmer.GetSig() - 0.0343271892109419) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong sig. 0.0343271892109419 != " << kmer.GetSig() << endl;
    }
    if (std::fabs(kmer.GetPf() - 0.789474868690288) >= 10e-15) {
        cout << "%TEST_FAILED% time=0 testname=testKmersPValue (KmerFactoryTest) message=Wrong pf. 0.789474868690288 != " << kmer.GetPf() << endl;
    }
}

void testWriteKmersToFile() {
    KmersFactory kmersFactory;
    Kmer *k;

    //            AAAAAAAAAA 3.25573332278254e-21 20.48735117540051 1.49479032886669
    k = new Kmer();
    k->SetControlFreq(11);
    k->SetNegativeControl(12);
    k->SetNegativePeak(13);
    k->SetPeakFreq(14);
    k->SetValue(3.25573332278254e-21);
    k->SetSig(20.48735117540051);
    k->SetPf(1.49479032886669);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAAAA", move(k)));

    //            AAAAAAAAAC 3.66936367681041e-05 4.43540924245372 1.62318695308535
    k = new Kmer();
    k->SetControlFreq(21);
    k->SetNegativeControl(22);
    k->SetNegativePeak(23);
    k->SetPeakFreq(24);
    k->SetValue(3.66936367681041e-05);
    k->SetSig(4.43540924245372);
    k->SetPf(1.62318695308535);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAAAC", move(k)));

    //            AAAAAAAAAG 1.56360980495561e-03 2.80587161489318 1.33064397039723
    k = new Kmer();
    k->SetControlFreq(31);
    k->SetNegativeControl(32);
    k->SetNegativePeak(33);
    k->SetPeakFreq(34);
    k->SetValue(1.56360980495561e-03);
    k->SetSig(2.80587161489318);
    k->SetPf(1.33064397039723);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAAAG", move(k)));

    //            AAAAAAAAAT 2.19631497914574e-01 0.65830537647355 1.08933620086384
    k = new Kmer();
    k->SetControlFreq(41);
    k->SetNegativeControl(42);
    k->SetNegativePeak(43);
    k->SetPeakFreq(44);
    k->SetValue(2.19631497914574e-01);
    k->SetSig(0.65830537647355);
    k->SetPf(1.08933620086384);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAAAT", move(k)));

    //            AAAAAAAACA 1.07873928249526e-02 1.96708350606356 1.41176344238936
    k = new Kmer();
    k->SetControlFreq(51);
    k->SetNegativeControl(52);
    k->SetNegativePeak(53);
    k->SetPeakFreq(54);
    k->SetValue(1.07873928249526e-02);
    k->SetSig(1.96708350606356);
    k->SetPf(1.41176344238936);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAACA", move(k)));

    //            AAAAAAAACC 1.05511751710985e-03 2.97669916672645 1.74193392488365
    k = new Kmer();
    k->SetControlFreq(61);
    k->SetNegativeControl(62);
    k->SetNegativePeak(63);
    k->SetPeakFreq(64);
    k->SetValue(1.05511751710985e-03);
    k->SetSig(2.97669916672645);
    k->SetPf(1.74193392488365);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAACC", move(k)));

    //            AAAAAAAACG 3.26401549435001e-01 0.48624778830827 1.28571313503317
    k = new Kmer();
    k->SetControlFreq(71);
    k->SetNegativeControl(72);
    k->SetNegativePeak(73);
    k->SetPeakFreq(74);
    k->SetValue(3.26401549435001e-01);
    k->SetSig(0.48624778830827);
    k->SetPf(1.28571313503317);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAACG", move(k)));

    //            AAAAAAAACT 1.48116148621570e-02 1.82939758918367 1.64515981794567
    k = new Kmer();
    k->SetControlFreq(81);
    k->SetNegativeControl(82);
    k->SetNegativePeak(83);
    k->SetPeakFreq(84);
    k->SetValue(1.48116148621570e-02);
    k->SetSig(1.82939758918367);
    k->SetPf(1.64515981794567);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAACT", move(k)));

    //            AAAAAAAAGA 2.70816098140979e-01 0.56732552342117 1.09793716170358
    k = new Kmer();
    k->SetControlFreq(91);
    k->SetNegativeControl(92);
    k->SetNegativePeak(93);
    k->SetPeakFreq(94);
    k->SetValue(2.70816098140979e-01);
    k->SetSig(0.56732552342117);
    k->SetPf(1.09793716170358);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAAGA", move(k)));
    //            AAAAAAAAGC 1.88687646220172e-03 2.72425653308383 1.63888742212561
    k = new Kmer();
    k->SetControlFreq(101);
    k->SetNegativeControl(102);
    k->SetNegativePeak(103);
    k->SetPeakFreq(104);
    k->SetValue(1.88687646220172e-03);
    //k->SetSig(2.72425653308383);
    k->SetSig(INFINITY);
    k->SetPf(1.63888742212561);

    kmersFactory.GetKmers().insert(pair<string, Kmer*>("AAAAAAAAGC", move(k)));

    kmersFactory.WriteKmersToFile("/tmp/kmers.bin", true);
}

void testReadKmersFromFile() {
    KmersFactory kmersFactory;
    char *comp;
    string rc_kmer;

    kmersFactory.ReadKmersFromFile("/tmp/kmers.bin", true);
    
    if (std::fabs(kmersFactory.GetKmerSig("AAAAAAAAAA") - 20.48735117540051) >= 10E-15){
        cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig. for kmer  != AAAAAAAAAA" << endl;
    }
    
    if (std::fabs(kmersFactory.GetKmerSig("AAAAAAAAGC") - (20.48735117540051 + 10)) >= 10E-15){
        cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig. for kmer  != AAAAAAAAAA" << endl;
    }

    for (auto it = kmersFactory.GetKmers().begin(); it != kmersFactory.GetKmers().end(); ++it) {
        bool right = false;
        string kmer = it->first;
        Kmer *k = it->second;

        comp = complement(kmer.c_str());
        rc_kmer = comp;
        free(comp);
        reverse(rc_kmer.begin(), rc_kmer.end());

        //                AAAAAAAAAA 3.25573332278254e-21 20.48735117540051 1.49479032886669
        if (kmer.compare("AAAAAAAAAA") == 0|| rc_kmer.compare("AAAAAAAAAA") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 3.25573332278254e-21) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 20.48735117540051) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.49479032886669) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 11) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 12) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 13) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 14) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }

        //                AAAAAAAAAC 3.66936367681041e-05 4.43540924245372 1.62318695308535
        if (kmer.compare("AAAAAAAAAC") == 0 || rc_kmer.compare("AAAAAAAAAC") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 3.66936367681041e-05) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 4.43540924245372) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.62318695308535) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 21) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 22) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 23) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 24) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }
        //                AAAAAAAAAG 1.56360980495561e-03 2.80587161489318 1.33064397039723
        if (kmer.compare("AAAAAAAAAG") == 0 || rc_kmer.compare("AAAAAAAAAG") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 1.56360980495561e-03) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 2.80587161489318) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.33064397039723) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 31) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 32) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 33) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 34) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }
        //                AAAAAAAAAT 2.19631497914574e-01 0.65830537647355 1.08933620086384
        if (kmer.compare("AAAAAAAAAT") == 0 || rc_kmer.compare("AAAAAAAAAT") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 2.19631497914574e-01) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 0.65830537647355) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.08933620086384) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 41) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 42) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 43) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 44) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }
        //                AAAAAAAACA 1.07873928249526e-02 1.96708350606356 1.41176344238936
        if (kmer.compare("AAAAAAAACA") == 0 || rc_kmer.compare("AAAAAAAACA") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 1.07873928249526e-02) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 1.96708350606356) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.41176344238936) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 51) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 52) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 53) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 54) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }
        //                AAAAAAAACC 1.05511751710985e-03 2.97669916672645 1.74193392488365
        if (kmer.compare("AAAAAAAACC") == 0 || rc_kmer.compare("AAAAAAAACC") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 1.05511751710985e-03) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 2.97669916672645) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.74193392488365) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 61) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 62) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 63) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 64) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }
        //                AAAAAAAACG 3.26401549435001e-01 0.48624778830827 1.28571313503317
        if (kmer.compare("AAAAAAAACG") == 0 || rc_kmer.compare("AAAAAAAACG") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 3.26401549435001e-01) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 0.48624778830827) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.28571313503317) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 71) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 72) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 73) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 74) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }
        //                AAAAAAAACT 1.48116148621570e-02 1.82939758918367 1.64515981794567
        if (kmer.compare("AAAAAAAACT") == 0 || rc_kmer.compare("AAAAAAAACT") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 1.48116148621570e-02) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 1.82939758918367) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.64515981794567) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 81) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 82) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 83) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 84) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }
        //                AAAAAAAAGA 2.70816098140979e-01 0.56732552342117 1.09793716170358
        if (kmer.compare("AAAAAAAAGA") == 0 || rc_kmer.compare("AAAAAAAAGA") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 2.70816098140979e-01) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - 0.56732552342117) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.09793716170358) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 91) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 92) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 93) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 94) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }
        //                AAAAAAAAGC 1.88687646220172e-03 2.72425653308383 1.63888742212561
        if (kmer.compare("AAAAAAAAGC") == 0 || rc_kmer.compare("AAAAAAAAGC") == 0) {
            right = true;
            if (std::fabs(k->GetValue() - 1.88687646220172e-03) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong p-value.  != " << k->GetValue() << endl;
            }
            if (std::fabs(k->GetSig() - (20.48735117540051 + 10)) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong sig.  != " << k->GetSig() << endl;
            }
            if (std::fabs(k->GetPf() - 1.63888742212561) >= 10E-15) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong pf.  != " << k->GetPf() << endl;
            }
            if (k->GetControlFreq() != 101) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong ControlFreq.  != " << k->GetControlFreq() << endl;
            }
            if (k->GetNegativeControl() != 102) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativeControl.  != " << k->GetNegativeControl() << endl;
            }
            if (k->GetNegativePeak() != 103) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong NegativePeak.  != " << k->GetNegativePeak() << endl;
            }
            if (k->GetPeakFreq() != 104) {
                cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong PeakFreq.  != " << k->GetPeakFreq() << endl;
            }
        }

        if (!right) {
            cout << "%TEST_FAILED% time=0 testname=testReadKmersFromFile (KmerFactoryTest) message=Wrong Peak.  != " << kmer << endl;
        }
    }

    remove("/tmp/kmers.bin");
}

int main(int argc, char** argv) {
    cout << "%SUITE_STARTING% KmerFactoryTest" << endl;
    cout << "%SUITE_STARTED%" << endl;

    Global::instance()->SetVerbose(3);
    Global::instance()->SetOrder(10);
    Global::instance()->SetBin1(0.005);
    Global::instance()->SetBin2(0.01);

    cout << "%TEST_STARTED% testKmers (KmerFactoryTest)" << endl;
    testKmers();
    cout << "%TEST_FINISHED% time=0 testKmers (KmerFactoryTest)" << endl;

    cout << "%TEST_STARTED% testKmersPValue (KmerFactoryTest)" << endl;
    testKmersPValue();
    cout << "%TEST_FINISHED% time=0 testKmersPValue (KmerFactoryTest)" << endl;

    cout << "%TEST_STARTED% testWriteKmersToFile (KmerFactoryTest)" << endl;
    testWriteKmersToFile();
    cout << "%TEST_FINISHED% time=0 testWriteKmersToFile (KmerFactoryTest)" << endl;

    cout << "%TEST_STARTED% testReadKmersFromFile (KmerFactoryTest)" << endl;
    testReadKmersFromFile();
    cout << "%TEST_FINISHED% time=0 testReadKmersFromFile (KmerFactoryTest)" << endl;

    cout << "%SUITE_FINISHED% time=0" << endl;

    delete Global::instance();
    return (EXIT_SUCCESS);
}


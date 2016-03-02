/* 
 * File:   Global.h
 * Author: veraalva
 *
 * Created on February 11, 2016, 1:11 PM
 */

#ifndef GLOBAL_H
#define GLOBAL_H

class Global {
    int verbose;
    unsigned long int order;
    double bin1;
    double bin2;
    static Global *s_instance;

    Global(int v = false) {
        verbose = v;
    }
public:

    bool GetVerbose() {
        return verbose;
    }

    bool isInfo() {
        if (verbose >= 1) return true;
        return false;
    }

    bool isDebug2() {
        if (verbose >= 2) return true;
        return false;
    }

    bool isDebug3() {
        if (verbose >= 3) return true;
        return false;
    }

    void SetVerbose(int v) {
        verbose = v;
    }

    unsigned long int GetOrder() const {
        return order;
    }

    void SetOrder(unsigned long int order) {
        this->order = order;
    }
    
    double GetBin1() const {
        return bin1;
    }

    void SetBin1(double bin1) {
        this->bin1 = bin1;
    }

    double GetBin2() const {
        return bin2;
    }

    void SetBin2(double bin2) {
        this->bin2 = bin2;
    }

    static Global *instance() {
        if (!s_instance)
            s_instance = new Global;
        return s_instance;
    }
};
#endif /* GLOBAL_H */


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
    std::set<int> orders;
    double bin1;
    double bin2;
    static Global *s_instance;

    Global() {
        verbose = 0;
    }
public:

    bool getVerbose() {
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

    void setVerbose(int v) {
        verbose = v;
    }

    unsigned long int getOrder() const {
        return order;
    }

    void setOrder(unsigned long int order) {
        this->order = order;
    }
    
    std::set<int>& getOrders() {
        return orders;
    }

    void setOrders(std::set<int> orders) {
        this->orders = orders;
    }

    double getBin1() const {
        return bin1;
    }

    void setBin1(double bin1) {
        this->bin1 = bin1;
    }

    double getBin2() const {
        return bin2;
    }

    void setBin2(double bin2) {
        this->bin2 = bin2;
    }

    static Global *instance() {
        if (!s_instance)
            s_instance = new Global;
        return s_instance;
    }
};
#endif /* GLOBAL_H */


/* 
 * File:   TimeUtils.h
 * Author: veraalva
 *
 * Created on February 16, 2016, 8:49 AM
 */

#include <ctime>

#ifndef TIMEUTILS_H
#define TIMEUTILS_H

class TimeUtils {
    clock_t begin;
    clock_t startTime;
    static TimeUtils *s_instance;
    
    TimeUtils(){}
public:

    static TimeUtils *instance() {
        if (!s_instance)
            s_instance = new TimeUtils;
        return s_instance;
    }

    void setStartTime() {
        startTime = clock();
    }

    void setClock() {
        begin = clock();
    }

    double getTimeSec() {
        return double(clock() - begin) / CLOCKS_PER_SEC;
    }

    double getTimeMin() {
        return getTimeSec() / 60;
    }

    double getTimeHour() {
        return getTimeSec() / 3600;
    }

    double getElapseTimeSec() {
        return double(clock() - startTime) / CLOCKS_PER_SEC;
    }

    double getElapseTimeMin() {
        return getElapseTimeSec() / 60;
    }

    double getElapseTimeHour() {
        return getElapseTimeSec() / 3600;
    }

    double getTimeSecFrom(clock_t b) {
        return double(clock() - b) / CLOCKS_PER_SEC;
    }

    double getTimeMinFrom(clock_t b) {
        return getTimeSecFrom(b) / 60;
    }

    double getTimeHourFrom(clock_t b) {
        return getTimeSecFrom(b) / 3600;
    }
};

#endif /* TIMEUTILS_H */


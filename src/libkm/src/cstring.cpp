/* 
 * File:   cString.cpp
 * Author: veraalva
 * 
 * Created on April 13, 2016, 3:47 PM
 */
#include <stdlib.h>

#include <string>
#include <algorithm>
#include <random>
#include <chrono>

#include "cstring.h"

using namespace std;

cstring::cstring() {
}

cstring::cstring(const cstring& orig) {
}

cstring::~cstring() {
}

std::string cstring::shuffle(std::string str) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(str.begin(), str.end(), std::default_random_engine(seed));
    return str;
}

std::string cstring::complement(std::string str) {
    for (auto it = str.begin(); it < str.end(); ++it) {
        switch (*it) {
            case 'A': *it = 'T';
                break;
            case 'a': *it = 'T';
                break;
            case 'T': *it = 'A';
                break;
            case 't': *it = 'A';
                break;
            case 'C': *it = 'G';
                break;
            case 'c': *it = 'G';
                break;
            case 'G': *it = 'C';
                break;
            case 'g': *it = 'C';
                break;
            case 'N': *it = 'N';
                break;
            case 'n': *it = 'N';
                break;
        }
    }
    return str;
}

std::string cstring::reverseComplement(std::string str) {
    str = cstring::complement(str);
    std::reverse(str.begin(), str.end());
    return str;
}

/**
 * Count the number of occurrences of characters in c in the string str
 * 
 * @param str the string to count on
 * @param c the characters to be counted
 * @return the number of occurrences 
 */
int cstring::countCharacter(std::string str, std::string characters) {
    int count = 0;
    for (auto it = str.begin(); it < str.end(); ++it) {
        for (auto it1 = characters.begin(); it1 < characters.end(); ++it1) {
            if (*it == *it1) count++;
        }
    }
    return count;
}
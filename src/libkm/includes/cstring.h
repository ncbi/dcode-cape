/* 
 * File:   cString.h
 * Author: veraalva
 *
 * Created on April 13, 2016, 3:47 PM
 */

#ifndef CSTRING_H
#define CSTRING_H

class cstring {
public:
    cstring();
    cstring(const cstring& orig);
    virtual ~cstring();

    static std::string shuffle(std::string str);
    static std::string complement(std::string str);
    static std::string reverseComplement(std::string str);

    /**
     * Count the number of occurrences of characters in c in the string str
     * 
     * @param str the string to count on
     * @param c the characters to be counted
     * @return the number of occurrences 
     */
    static int countCharacter(std::string str, std::string characters);
private:

};

#endif /* CSTRING_H */


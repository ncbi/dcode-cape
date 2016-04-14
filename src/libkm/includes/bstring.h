/* 
 * File:   bstring.h
 * Author: roberto
 *
 * Created on April 15, 2014, 3:05 PM
 */

#ifndef BSTRING_H
#define BSTRING_H

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * Split a string into an array of string using the delimiter
     * 
     * @param dest the destination pointer
     * @param src the source pointer
     * @param delimiter the delimiter
     * @return a array with the splited strings
     */
    extern size_t strsep_ptr(char ***tokens, size_t *len, char *src, const char *delimiter);

    /**
     * Change the string to upper characters
     * 
     * @param src the string to be changed
     */
    void toUpper(char **src);


#ifdef __cplusplus
}
#endif

#endif /* BSTRING_H */


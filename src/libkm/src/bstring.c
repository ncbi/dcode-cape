/* 
 * File:   bstring.c
 * Author: roberto
 *
 * Created on April 15, 2014, 3:05 PM
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include "berror.h"
#include "bmemory.h"

/**
 * Split a string into a array of pointers string using the delimiter
 * 
 * @param dest the destination pointer
 * @param len number of elements of the destination pointer
 * @param src the source pointer
 * @param delimiter the delimiter
 * @return a array with the split strings
 */
size_t strsep_ptr(char ***tokens, size_t *len, char *src, const char *delimiter) {
    size_t nWords = 0;
    char *token;

    while ((token = strsep(&src, delimiter)) != NULL) {
        if (*token != 0) {
            if (*len <= nWords) {
                *len = nWords + 1;
                *tokens = (char **) reallocate(*tokens, sizeof (char **) * (*len), __FILE__, __LINE__);
            }
            *(*tokens + nWords++) = token;
        }
    }
    return nWords;
}

/**
 * Change the string to upper characters
 * 
 * @param src the string to be changed
 */
void toUpper(char **src) {
    char *s = *src;
    while (*s != 0) {
        *s = toupper(*s);
        s++;
    }
}
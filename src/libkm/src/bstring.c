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
                *tokens = (char **) realloc(*tokens, sizeof (char **) * (*len));
            }
            *(*tokens + nWords++) = token;
        }
    }
    return nWords;
}

/**
 * Count the number of occurrences of characters in c in the string str
 * 
 * @param str the string to count on
 * @param c the character to be counted ended by \0
 * @return the number of occurrences 
 */
int countCharacter(const char *str, const char *chracters) {
    int count = 0;
    const char *s = str;
    while (*s != 0) {
        const char *c = chracters;
        while (*c != 0) {
            if (*s == *c) count++;
            c++;
        }
        s++;
    }
    return count;
}

/**
 * Create a new string with the complement sequence of input
 * 
 * @param input the sequence 
 * @return a new char * with the complement sequence
 */
char *complement(const char *input) {
    char *output = allocate(sizeof (char) * (strlen(input) + 1), __FILE__, __LINE__);
    char temp, *o;
    const char *i;
    i = input;
    o = output;
    while (*i != 0) {
        switch (*i) {
            case 'A': temp = 'T';
                break;
            case 'a': temp = 'T';
                break;
            case 'T': temp = 'A';
                break;
            case 't': temp = 'A';
                break;
            case 'C': temp = 'G';
                break;
            case 'c': temp = 'G';
                break;
            case 'G': temp = 'C';
                break;
            case 'g': temp = 'C';
                break;
            case 'N': temp = 'N';
                break;
            case 'n': temp = 'N';
                break;
            default:
                temp = *i;
        }
        *o++ = temp;
        i++;
    }
    *o = 0;
    return output;
}

/**
 * Shuffle the string
 * 
 * @param input the input string to be shuffle 
 * @return a new char * with the shuffle string
 */
char *shuffle(const char *input) {
    char *output = strdup(input);
    int i = strlen(output) - 1;
    while (i > 0) {
        int j = rand() % i;
        char c = output[i];
        output[i] = output[j];
        output[j] = c;
        i--;
    }
    return output;
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
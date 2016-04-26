/* 
 * File:   FileParserFactory.cpp
 * Author: veraalva
 * 
 * Created on April 11, 2016, 12:50 PM
 */
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <map>
#include <vector>
#include <fstream>

#include "berror.h"
#include "bmemory.h"
#include "bstring.h"

#include "FileParserFactory.h"

using namespace std;
using namespace parsers;

FileParserFactory::FileParserFactory() {
    this->closeFile = false;
    this->words = NULL;
    this->wordsSize = 0;
    this->nWords = 0;
    this->buffer = NULL;
    this->bufferEndPtr = NULL;
    this->str = NULL;
    this->line = NULL;
    this->backup = NULL;
    this->bufferSize = 10000000;
    this->backupTotalSize = 0;
    this->backupSize = 0;
    this->read = 0;
    this->lineLength = 0;
}

FileParserFactory::~FileParserFactory() {
    if (this->closeFile && this->fileToParse.is_open()) this->fileToParse.close();
    if (this->words) free(this->words);
    if (this->buffer) free(this->buffer);
    if (this->backup) free(this->backup);
}

void FileParserFactory::clean() {
    if (this->closeFile && this->fileToParse.is_open()) this->fileToParse.close();
    if (this->words) free(this->words);
    if (this->buffer) free(this->buffer);
    if (this->backup) free(this->backup);
    this->closeFile = false;
    this->words = NULL;
    this->wordsSize = 0;
    this->nWords = 0;
    this->buffer = NULL;
    this->bufferEndPtr = NULL;
    this->str = NULL;
    this->line = NULL;
    this->backup = NULL;
    this->bufferSize = 10000000;
    this->backupTotalSize = 0;
    this->backupSize = 0;
    this->read = 0;
    this->lineLength = 0;
}

bool FileParserFactory::iterate(const char dontStartWith, const char *delimiter) {
    size_t i;
    char *newLinePtr;

    if (!fileToParse.is_open()) {
        throw exceptions::FileNotFoundException("Please, open the file correctly");
    }

    while (1) {
        if (!str || *str == 0) {
            if (fileToParse.eof()) return false;

            if (!buffer) {
                buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
            }

            fileToParse.read(buffer, sizeof (char) * bufferSize);
            read = fileToParse.gcount();
            buffer[read] = 0;
            bufferEndPtr = buffer + read;

            if (fileToParse.eof()) {
                if (buffer[read - 1] != '\n') {
                    buffer[read] = '\n';
                    buffer[read + 1] = 0;
                    bufferEndPtr = buffer + read + 1;
                }
            }
            str = buffer;
        }

        if ((newLinePtr = strchr(str, '\n')) != NULL) {
            *newLinePtr = 0;
            if (backup && *backup != 0) {
                i = backupSize;
                backupSize += (bufferEndPtr - str);
                if (backupSize > backupTotalSize) {
                    backupTotalSize = backupSize;
                    backup = (char *) reallocate(backup, sizeof (char) * (backupTotalSize + 1), __FILE__, __LINE__);
                }
                for (size_t j = 0; i < backupTotalSize; i++) {
                    backup[i] = str[j];
                    if (str[j++] == 0) break;
                }
                line = backup;
                lineLength = i;
            } else {
                line = str;
            }
            backupSize = 0;
            str = newLinePtr + 1;
            if (line != backup) {
                lineLength = str - line - 1;
            }
            if (*line != dontStartWith) {
                nWords = strsep_ptr(&words, &wordsSize, line, delimiter);
                return true;
            }
        } else {
            if (str && *str != 0) {
                i = backupSize;
                backupSize += (bufferEndPtr - str);
                if (backupSize > backupTotalSize) {
                    backup = (char *) reallocate(backup, sizeof (char) * (backupSize + 1), __FILE__, __LINE__);
                    *(backup + backupTotalSize) = 0;
                    backupTotalSize = backupSize;
                }
                for (size_t j = 0; i < backupTotalSize; i++) {
                    backup[i] = str[j];
                    if (str[j++] == 0) break;
                }
                *str = 0;
            }
        }
    }
    return false;
}

bool FileParserFactory::iterate(const char dontStartWith) {
    size_t i;
    char *newLinePtr;

    if (!fileToParse.is_open()) {
        throw exceptions::FileNotFoundException("Can't do an iteration in a NULL file");
    }

    while (1) {
        if (!str || *str == 0) {
            if (fileToParse.eof()) return false;

            if (!buffer) {
                buffer = (char *) allocate(sizeof (char) * (bufferSize + 1), __FILE__, __LINE__);
            }

            fileToParse.read(buffer, sizeof (char) * bufferSize);
            read = fileToParse.gcount();
            buffer[read] = 0;
            bufferEndPtr = buffer + read;

            if (fileToParse.eof()) {
                if (buffer[read - 1] != '\n') {
                    buffer[read] = '\n';
                    buffer[read + 1] = 0;
                    bufferEndPtr = buffer + read + 1;
                }
            }
            str = buffer;
        }

        if ((newLinePtr = strchr(str, '\n')) != NULL) {
            *newLinePtr = 0;
            lineLength = 0;
            if (backup && *backup != 0) {
                i = backupSize;
                backupSize += (bufferEndPtr - str);
                if (backupSize > backupTotalSize) {
                    backupTotalSize = backupSize;
                    backup = (char *) reallocate(backup, sizeof (char) * (backupTotalSize + 1), __FILE__, __LINE__);
                }
                for (size_t j = 0; i < backupTotalSize; i++) {
                    backup[i] = str[j];
                    if (str[j++] == 0) break;
                }
                lineLength = i;
                line = backup;
            } else {
                line = str;
            }
            backupSize = 0;
            str = newLinePtr + 1;
            if (line != backup) {
                lineLength = str - line - 1;
            }
            if (*line != dontStartWith) return true;
        } else {
            if (str && *str != 0) {
                i = backupSize;
                backupSize += (bufferEndPtr - str);
                if (backupSize > backupTotalSize) {
                    backup = (char *) reallocate(backup, sizeof (char) * (backupSize + 1), __FILE__, __LINE__);
                    *(backup + backupTotalSize) = 0;
                    backupTotalSize = backupSize;
                }
                for (size_t j = 0; i < backupTotalSize; i++) {
                    backup[i] = str[j];
                    if (str[j++] == 0) break;
                }
                *str = 0;
            }
        }
    }
    return false;
}

void FileParserFactory::wordsToVector(std::vector<std::string>& v) {
    for (size_t i = 0; i < this->nWords; i++) {
        v.push_back(words[i]);
    }
}

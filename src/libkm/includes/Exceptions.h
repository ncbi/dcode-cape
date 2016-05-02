/* 
 * File:   Exceptions.h
 * Author: veraalva
 *
 * Created on March 24, 2016, 1:20 PM
 */

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

namespace exceptions {

    class FileNotFoundException : public std::exception {
    public:

        explicit FileNotFoundException(const char* message) : msg(message) {
        }

        explicit FileNotFoundException(const std::string& message) : msg(message) {
        }

        virtual ~FileNotFoundException() {
        }

        virtual const char* what() const throw () {
            return msg.c_str();
        }

    private:
        std::string msg;
    };

    class NotFoundException : public std::exception {
    public:

        explicit NotFoundException(const char* message) : msg(message) {
        }

        explicit NotFoundException(const std::string& message) : msg(message) {
        }

        virtual ~NotFoundException() {
        }

        virtual const char* what() const throw () {
            return msg.c_str();
        }

    private:
        std::string msg;
    };

}

#endif /* EXCEPTIONS_H */


/* 
 * File:   Exceptions.h
 * Author: veraalva
 *
 * Created on March 24, 2016, 1:20 PM
 */

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

namespace exceptions {

    class OutOfRangeException : public std::exception {
    public:

        explicit OutOfRangeException(const char* message) : msg(message) {
        }

        explicit OutOfRangeException(const std::string& message) : msg(message) {
        }

        virtual ~OutOfRangeException() {
        }

        virtual const char* what() const throw () {
            return msg.c_str();
        }

    private:
        std::string msg;
    };

}

#endif /* EXCEPTIONS_H */


/*
 * DNASequenceException.h
 * Created by Tanner Quigley on 1/29/2026
 * Summary:
 *  - Custom exception used to flag composition problems of input sequences.
 */

#ifndef M10EP_TEQUIGLE_DNASEQUENCEERROR_H
#define M10EP_TEQUIGLE_DNASEQUENCEERROR_H

#include <exception>
#include <string>

class DNASequenceException : public std::exception {
private:
    std::string message;
public:
    // Constructor
    DNASequenceException(const std::string& msg) : message(msg) {}

    // Overwriting of what() base method
    const char* what() const noexcept override {
        return message.c_str();
    }
};

#endif //M10EP_TEQUIGLE_DNASEQUENCEERROR_H
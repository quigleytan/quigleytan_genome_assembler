/*
 * DNASequenceException.h
 * Created by Tanner Quigley on 2/15/2026
 * Summary:
 *  - Custom exception used to flag composition problems of input sequences.
 */
#ifndef M10EP_TEQUIGLE_NODENOTFOUNDEXCEPTION_H
#define M10EP_TEQUIGLE_NODENOTFOUNDEXCEPTION_H

#include <exception>
#include <string>
#include <cstdint>

class NodeNotFoundException : public std::exception {
private:
    std::string message;
public:
    explicit NodeNotFoundException(uint64_t node)
        : message("Node not found: " + std::to_string(node)) {}

    // Overwriting of what() base method
    const char* what() const noexcept override {
        return message.c_str();
    }
};
#endif //M10EP_TEQUIGLE_NODENOTFOUNDEXCEPTION_H
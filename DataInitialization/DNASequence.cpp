#include "DNASequence.h"
#include "../CustomExceptions/DNASequenceException.h"

// Struct to allow for easier return of data
struct AnalysisResult {
    std::string reverseSequence;  // Reverse complement sequence
    int gcTotal;                  // Number of G/C bases
    double gcPercent;             // Percentage of G/C
};

// USED ONLY IN CONSTRUCTOR
DNASequence::AnalysisResult DNASequence::analyzeSequence(const std::string& inputSequence) {
    // Ensures no empty sequences are added
    if (inputSequence.empty()) {
        throw DNASequenceException("Sequence is empty");
    }

    std::string reverse;
    int gcCount = 0;

    for (char base : inputSequence) {
        switch(base) {
            case 'A': reverse += 'T'; break;
            case 'T': reverse += 'A'; break;
            case 'C': reverse += 'G'; gcCount++; break;
            case 'G': reverse += 'C'; gcCount++; break;
            // Throws an exception if inputSequence contains an invalid base
            default: throw DNASequenceException("Invalid base: " + std::string(1, base));
        }
    }
    // G and C base percentage calculation
    double gcPercent = double(gcCount) / double(inputSequence.length());

    return {reverse, gcCount, gcPercent};
}

// Constructor
DNASequence::DNASequence(const std::string& name, const std::string& inputSequence)
    : sequence(inputSequence),
    length(inputSequence.length()) {
    this->name = name;
    AnalysisResult analysis = analyzeSequence(inputSequence);
    complementSequence = analysis.complementSequence;
    gcTotal = analysis.gcTotal;
    gcPercent = analysis.gcPercent;
}

DNASequence::DNASequence() {
    sequence = "";
    complementSequence = "";
    gcTotal = 0;
    gcPercent = 0.0;
    length = 0;
}

// Getters
const std::string& DNASequence::getSequence() const {
    return sequence;
}

const std::string& DNASequence::getComplement() const {
    return complementSequence;
}

const std::string& DNASequence::getName() const {
    return name;
}

size_t DNASequence::getLength() const {
    return length;
}

int DNASequence::getGCCount() const {
    return gcTotal;
}

double DNASequence::getGCPercent() const {
    return gcPercent;
}

// Comparison overloaded operators used to compare the composition of two sequences
bool operator == (const DNASequence &sequenceOne, const DNASequence &sequenceTwo) {
    return sequenceOne.getSequence() == sequenceTwo.getSequence();
}

bool operator != (const DNASequence &sequenceOne, const DNASequence &sequenceTwo) {
    return sequenceOne.getSequence() != sequenceTwo.getSequence();
}

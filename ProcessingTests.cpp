#include <iostream>
#include <ostream>

#include "KmerEncoder.h"
#include "KmerTable.h"
#include "DataInitialization/DNASequence.h"

bool kmerEncoderTests();
bool kmerTableTests();

int main() {

    if (kmerEncoderTests()) {
        std::cout << "Kmer Encoder Tests Passed" << std::endl;
    }

    if (kmerTableTests()) {
        std::cout << "Kmer Table Tests Passed" << std::endl;
    }

    return 0;
}

bool kmerEncoderTests() {
    bool passed = true;

    // Should encode as A = 0b00, C = 0b01, G = 0b10, T = 0b11
    if (KmerEncoder::encodeKmer("A") != 0b00) {
        passed = false;
        std::cout << "Improper encoding of A" << std::endl;
    }

    if (KmerEncoder::decodeKmer(0b00, 1) != "A") {
        passed = false;
        std::cout << "Improper decoding of A" << std::endl;
    }

    if (KmerEncoder::encodeKmer("C") != 0b01) {
        passed = false;
        std::cout << "Improper encoding of A" << std::endl;
    }

    if (KmerEncoder::encodeKmer("G") != 0b10) {
        passed = false;
        std::cout << "Improper encoding of A" << std::endl;
    }

    if (KmerEncoder::encodeKmer("T") != 0b11) {
        passed = false;
        std::cout << "Improper encoding of A" << std::endl;
    }

    if (KmerEncoder::encodeKmer("ACGGTGT") != 1723) {
        passed = false;
        std::cout << "Improper encoding of ACGGTGT" << std::endl;
    }

    return passed;
}

bool kmerTableTests() {
    bool passed = true;

    DNASequence genome("Sequence 1", "ACGTACGT");

    KmerTable kTable(genome.getLength(), 3);
    KmerEncoder::encodeSequence(genome.getSequence(), 3, kTable);

    if (kTable.getK() != 3) {
        passed = false;
        std::cout << "Kmer table initialized with incorrect k" << std::endl;
    }

    const uint64_t* p = kTable.find(KmerEncoder::encodeKmer("ACG"));
    uint64_t foundValidValue = p ? *p : 0;

    if (foundValidValue != 2) {
        passed = false;
        std::cout << "Kmer not found in table" << std::endl;
    }

    const uint64_t* p2 = kTable.find(KmerEncoder::encodeKmer("AAA"));
    uint64_t foundInvalidValue = p2 ? *p2 : 0;

    if (foundInvalidValue != 0) {
        passed = false;
        std::cout << "Invalid kmer found in table" << std::endl;
    }

    return passed;
}

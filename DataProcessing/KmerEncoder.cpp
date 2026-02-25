#include "KmerEncoder.h"
#include "KmerTable.h"
#include "CustomExceptions/DNASequenceException.h"

uint64_t KmerEncoder::bitmask(size_t k) {
    // Bitmask to keep last k bases (2 bits per base)
    return (1ULL << (2 * k)) - 1;
}

uint64_t KmerEncoder::encodeBase(char base) {
    // Assigning 2-bit values to each base
    switch (base) {
        case 'A': return 0b00;
        case 'C': return 0b01;
        case 'G': return 0b10;
        case 'T': return 0b11;
        default: throw DNASequenceException("Invalid base");
    }
}

// Rolls the previous k-mer to get the next one
uint64_t KmerEncoder::roll(uint64_t prev, char next, size_t k) {
    // Shift left to make space for the new base
    prev <<= 2;
    prev |= encodeBase(next);
    prev &= bitmask(k); // Trims the sequence to size k-kmers
    return prev;
}

// START OF ENCODING SEQUENCE
uint64_t KmerEncoder::encodeKmer(const std::string& kmer) {
    uint64_t value = 0;
    for (char base : kmer) {
        value <<= 2; // Shift left to include the next base
        value |= encodeBase(base); // Adds new base
    }
    return value;
}

std::string KmerEncoder::decodeKmer(uint64_t encoded, size_t k) {
    std::string result(k, 'A');

    for (size_t i = 0; i < k; ++i) {
        uint64_t bits = encoded & 0b11;

        switch (bits) {
            case 0b00: result[k - 1 - i] = 'A'; break;
            case 0b01: result[k - 1 - i] = 'C'; break;
            case 0b10: result[k - 1 - i] = 'G'; break;
            case 0b11: result[k - 1 - i] = 'T'; break;
        }
        encoded >>= 2;
    }
    return result;
}

void KmerEncoder::encodeSequence(const std::string &dna, size_t k, KmerTable &table){
    if (dna.length() < k) return;

    // Encode first k-mer
    uint64_t kmer = encodeKmer(dna.substr(0, k));
    table.insert(kmer);

    // Roll through the rest of the sequence
    for (size_t i = k; i < dna.length(); ++i) {
        kmer = roll(kmer, dna[i], k);
        table.insert(kmer);
    }
}

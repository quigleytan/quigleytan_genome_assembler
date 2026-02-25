/*
 * KmerEncoder.h
 * Created by Tanner Quigley on 1/30/2026.
 * Summary:
 * - Static class for encoding and rolling k-mers from DNA sequences.
 * - Rolling saves time by building on previous k-mer encoding.
 * - Primary function is to provide utility for encoding a kmer, and also
 *   to encode an entire sequence into a KmerTable.
 * Important notes:
 * - Uses the key: A = 00, C = 01, G = 10, T = 11 for 2-bit encoding.
 * - All encoded kmers are stored as uint64_t.
 * - Used for quick storage and retrieval in hash tables.
 * TODO items for future work;
 * - Eventually will add support for a decoding function to help with reporting
 *   results back to users.
 */

#ifndef M20EP_TEQUIGLE_KMERENCODER_H
#define M20EP_TEQUIGLE_KMERENCODER_H

#include <cstdint>
#include <string>
#include "KmerTable.h"

class KmerEncoder {
private:
    /**
     * @brief Bitmask function to keep the last k bases.
     *
     * Used to ensure that only the last k bases are kept in roll(), which allows
     * for continuous building of k-mers without overflow.
     *
     * @param k Set k-mer length.
     * @return Bitmask for k-mer of length k.
     */
    static uint64_t bitmask(size_t k);

    /**
     * @brief Encodes a single DNA base into its corresponding 2-bit representation.
     * @param base Char input value of a single base.
     * @return 2-bit value as uint64_t.
     */
    static uint64_t encodeBase(char base);

    /**
     * @brief Rolls the previous k-mer to create the next k-mer by adding a new base.
     *
     * Shifts the previous k-mer left by 2 bits, adds the new base and applies a bitmask to ensure
     * only the last k bases are kept.
     *
     * @param prev 2-bit encoded previous kmer as uint64_t.
     * @param next Char input value of the next base to add.
     * @param k Set k-mer length to use when rolling.
     *
     * @return 2-bit encoded next kmer as uint64_t.
     */
    static uint64_t roll(uint64_t prev, char next, size_t k);

public:
    /**
     * @brief Encodes a k-mer string into its 2-bit representation.
     *
     * Encodes the entire k-mer by iterating through each base, shifting left and adding
     * the corresponding 2-bit value. Used for the initial encoding before rolling and
     * for quick encoding of a single kmer.
     *
     * @param kmer String of initial kmer to encode.
     *
     * @return 2-bit encoded kmer as uint64_t.
     */
    static uint64_t encodeKmer(const std::string& kmer);

    /**
     * @brief Decodes a 2-bit encoded k-mer into its string representation.
     *
     * Decodes a k-mer key for results reporting and easier interpretation of results.
     * Only for use with user interaction, keep everything encoded in logical implementation.
     *
     * @param kmer 2-bit encoded kmer as uint64_t.
     * @param kmer The length of the k-mer.
     * @return Equivalent string representation.
     */
    static std::string decodeKmer(uint64_t kmer, size_t k);

    /**
     * @brief Performs k-mer analysis on a DNA sequence and populates a KmerTable.
     *
     * Encodes and rolls through the DNA sequence to extract all k-mers, inserting them into the provided KmerTable.
     * This generates counts and allows for analysis of k-mer frequency in the sequence.
     *
     * @param dna DNA sequence string to extract k-mers from.
     * @param k K-mer size to use for encoding.
     * @param table Table to populate with encoded k-mers.
     */
    static void encodeSequence(const std::string& dna, size_t k, KmerTable& table);
};

#endif //M20EP_TEQUIGLE_KMERENCODER_H
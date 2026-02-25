/*
 * SequenceReader.h
 * Created by Tanner Quigley on 2/15/2026
 * Summary:
 * - Pulls information from FASTA and FASTQ files to create DNASequence objects
 * - Uses std::optional to handle empty files gracefully
 * - Throws std::runtime_error for various format errors in the input files
 * Important notes:
 * - Returns optionals, need to extract the object from the optional
 * TODO items for future work;
 * - Potentially add more support for other files
 */

#ifndef M20EP_TEQUIGLE_SEQUENCE_READER_H
#define M20EP_TEQUIGLE_SEQUENCE_READER_H

#include "DNASequence.h"
#include <optional>

class SequenceReader {

public:

    /**
     * @brief Reads a FASTA format sequence into optional DNASequence.
     * @param in Data file containing sequence information
     * @return Optional DNASequence object if successful, null optional if the file is empty
     */
    static std::optional<DNASequence> readFasta(std::istream& in);

    /**
     * @brief Reads a FASTQ format sequence into optional DNASequence.
     * @param in Data file containing sequence information
     * @return Optional DNASequence object if successful, null optional if the file is empty
     */
    static std::optional<DNASequence> readFastq(std::istream& in);
};

#endif //M20EP_TEQUIGLE_SEQUENCE_READER_H
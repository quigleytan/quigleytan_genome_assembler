/*
 * DNASequence.h
 * Summary:
 * - Represents a DNA sequence and provides analysis of its composition.
 * Features:
 * - Stores a DNA sequence (immutable once constructed).
 * - Computes reverse complement.
 * - Computes GC count and GC percent.
 * - Throws DNASequenceException for an invalid sequence input.
 * Additional:
 * - Internally uses a helper struct AnalysisResult for single-pass computation.
 * - All internal methods and structs are private.
 * TODO items for future work;
 * - Include to_string or summary function.
 */

#ifndef M1OEP_TEQUIGLE_DNA_SEQUENCE_H
#define M1OEP_TEQUIGLE_DNA_SEQUENCE_H

#include <string>
#include <iostream>

class DNASequence {
private:

    // Variables
    std::string name;
    std::string sequence;
    std::string complementSequence;
    int gcTotal;
    double gcPercent;
    size_t length; // # of bases

    // Struct to allow for a single loop data analysis without modifying
    // members outside of the constructor.
    struct AnalysisResult {
        std::string complementSequence; // complement sequence
        int gcTotal;                    // number of G/C bases
        double gcPercent;               // percentage of G/C
    };

    /**
     * @brief Computes basic sequence data
     *
     * Creates the reverse complement, G/C counts and percentage in a single pass through the input sequence.
     * Used only in the constructor.
     *
     * @param inputSequence Sequence passed in through the constructor.
     * @return AnalysisResult Struct containing the reverseSequence, gcTotal, and gcPercent.
     */
    static AnalysisResult analyzeSequence(const std::string& inputSequence);

public:

    /**
     * @brief Constructor for DNASequence.
     *
     * Takes the name and string of bases to create a DNASequence object.
     * The sequence is analyzed for its reverse complement, GC count percentage.
     *
     * @param name Name of the sequence derived from a data file.
     * @param inputSequence String of char's to be analyzed and stored.
     */
    DNASequence(const std::string& name, const std::string& inputSequence);

    /**
     * @brief Basic default constructor, sets all values to "" or 0.
     */
    DNASequence();

    // Getters

    /**
     * @brief Returns the DNA sequence string.
     * @return sequence The DNA sequence string (immutable after construction).
     */
    const std::string& getSequence() const;

    /**
     * @brief Returns the reverse complement of the sequence.
     * @return reverseSequence The reverse complement of the sequence.
     */
    const std::string& getComplement() const;

    /**
     * @brief Returns the associated name with the sequence.
     * @return Name of the sequence derived from a data file.
     */
    const std::string& getName() const;

    /**
     * @brief Returns the length of the sequence.
     * @return size_t length of the sequence (number of bases).
     */
    size_t getLength() const;

    /**
     * @brief Returns the number of instances of G's and C's in the string.
     * @return gcTotal The total count of guanine and cytosine in the sequence.
     */
    int getGCCount() const;

    /**
     * @brief Returns the GC content of the string.
     * @return gcPercent The percentage of guanine and cytosine in the sequence.
     */
    double getGCPercent() const;

    //Overloaded comparisons
    friend bool operator == (const DNASequence &sequenceOne, const DNASequence &sequenceTwo);

    friend bool operator != (const DNASequence &sequenceOne, const DNASequence &sequenceTwo);

};
#endif //M1OEP_TEQUIGLE_DNA_SEQUENCE_H

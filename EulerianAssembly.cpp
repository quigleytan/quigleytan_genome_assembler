//
// Created by quigl on 1/27/2026.
//
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

#include "CustomExceptions/DNASequenceException.h"

#include "DataInitialization/DNASequence.h"
#include "DataInitialization/SequenceReader.h"

#include "DataProcessing/KmerEncoding.h"
#include "DataProcessing/KmerTable.h"

#include "GenomeAssembly/DeBruijnGraph.h"
#include "GenomeAssembly/EulerianTraversal.h"

// -----------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------

int getIntFromUser() {
    while (true) {
        std::string input;
        std::cout << "Enter an integer for a desired k value: ";
        getline(std::cin, input);
        try {
            size_t size = 0;
            int value = std::stoi(input, &size);
            if (size == input.size() && value > 0)
                return value;
            std::cout << "Invalid input. ";
        }
        catch (...) {
            std::cout << "Invalid input. ";
        }
    }
}

// Builds the KMP failure function (LPS array) for a given pattern.
static void constructLps(const std::string& pat, std::vector<int>& lps) {
    int len = 0;
    lps[0] = 0;
    int i = 1;

    while (i < static_cast<int>(pat.length())) {
        if (pat[i] == pat[len]) {
            lps[i++] = ++len;
        } else {
            if (len != 0)
                len = lps[len - 1];
            else
                lps[i++] = 0;
        }
    }
}

// Returns all start indices where pat appears in txt.
static std::vector<int> kmpSearch(const std::string& pat, const std::string& txt) {
    int n = static_cast<int>(txt.length());
    int m = static_cast<int>(pat.length());

    std::vector<int> lps(m);
    std::vector<int> res;
    constructLps(pat, lps);

    int i = 0, j = 0;
    while (i < n) {
        if (txt[i] == pat[j]) {
            i++; j++;
            if (j == m) {
                res.push_back(i - j);
                j = lps[j - 1];
            }
        } else {
            if (j != 0) j = lps[j - 1];
            else        i++;
        }
    }
    return res;
}

static size_t findRotationOffset(const std::string& original, const std::string& assembled) {
    if (original.length() != assembled.length()) return std::string::npos;
    if (original == assembled) return 0;

    size_t n = original.length();
    std::string doubled = assembled + assembled;

    // KMP on 'original' as pattern, 'doubled' as text
    std::vector<int> lps(n, 0);
    int len = 0, i = 1;
    while (i < (int)n) {
        if (original[i] == original[len]) { lps[i++] = ++len; }
        else if (len != 0)                { len = lps[len-1]; }
        else                              { lps[i++] = 0; }
    }

    size_t j = 0;
    for (size_t k = 0; k < doubled.size(); ++k) {
        while (j > 0 && doubled[k] != original[j]) j = lps[j-1];
        if (doubled[k] == original[j]) ++j;
        if (j == n) return k - n + 1; // offset into assembled
    }
    return std::string::npos;
}

// -----------------------------------------------------------------------
// Pipeline stages
// -----------------------------------------------------------------------

/**
 * @brief Stage 1 — Load a FASTA file into a DNASequence.
 *
 * Opens the file at the given path, reads the first FASTA record,
 * and returns the resulting DNASequence. Throws std::runtime_error
 * if the file cannot be opened or the record is missing.
 *
 * @param path Path to the FASTA file.
 * @return DNASequence parsed from the file.
 */
static DNASequence loadGenome(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + path);

    auto genomeOpt = SequenceReader::readFasta(file);
    file.close();

    if (!genomeOpt)
        throw std::runtime_error("File was empty: " + path);

    return *genomeOpt;
}

/**
 * @brief Stage 2 — Build a De Bruijn graph from a DNA sequence.
 *
 * Circularizes the sequence by appending the first k-1 bases, encodes
 * all k-mers into a KmerTable, then inserts them into a DeBruijnGraph
 * pre-sized to avoid rehashing.
 *
 * @param sequence  The raw (non-circularized) DNA sequence string.
 * @param k         K-mer size to use.
 * @return          Populated DeBruijnGraph ready for Eulerian traversal.
 */
static DeBruijnGraph buildGraph(const std::string& sequence, int k) {
    std::string circularSequence = sequence + sequence.substr(0, k - 1);

    KmerTable kTable(circularSequence.length(), k);
    KmerEncoding::encodeSequence(circularSequence, k, kTable);
    std::cout << "Unique k-mers:   " << kTable.getNumItems() << "\n";

    DeBruijnGraph graph(k, kTable.getNumItems() * 2);
    for (const auto& entry : kTable) {
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);
    }
    std::cout << "Graph built:     " << graph.getNodeCount() << " nodes, "
              << graph.getEdgeCount() << " edges\n";

    return graph;
}

/**
 * @brief Stage 3 — Traverse the De Bruijn graph and reconstruct the genome.
 *
 * Runs Hierholzer's algorithm on the graph, then reconstructs and trims
 * the assembled string back to the original sequence length.
 *
 * Since the graph is built from a circularized sequence, the Eulerian path
 * is a circuit and the start node is non-deterministic. The assembled string
 * is normalized back to the original start position by searching for the
 * assembled string as a substring of original+original using KMP. The match
 * position gives the rotation offset needed to align the two sequences.
 *
 * @param graph            The populated De Bruijn graph.
 * @param originalSequence The original (non-circularized) sequence string.
 * @param k                K-mer size used to build the graph (unused directly,
 *                         retained for interface consistency).
 * @return                 Assembled genome string normalized to original start
 *                         position, or the unrotated string if no rotation match
 *                         is found (caller can detect via RECONSTRUCTION: FAILED).
 */
static std::string assembleGenome(DeBruijnGraph& graph,
                                   const std::string& originalSequence,
                                   int k) {
    EulerianTraversal eulerianPath(graph);
    eulerianPath.computePath();

    std::cout << "Path length:     " << eulerianPath.getPath().size() << " nodes\n";

    std::string assembled = eulerianPath.reconstructGenome(true);
    if (assembled.length() > originalSequence.length())
        assembled = assembled.substr(0, originalSequence.length());

    // Find the rotation offset by searching assembled inside original+original.
    // This correctly handles the non-deterministic circuit start without relying
    // on a short anchor that may appear at multiple offsets.

    size_t offset = findRotationOffset(originalSequence, assembled);
    if (offset != std::string::npos && offset != 0) {
        assembled = assembled.substr(offset) + assembled.substr(0, offset);
    }

    std::cout << "Rotation offset: ";
    std::cout << std::scientific << static_cast<double>(offset);
    std::cout << " bases\n";

    // assembled starts offset bases into the original. Rotating assembled
    // left by offset positions aligns it with the original start.
    return assembled;
}

/**
 * @brief Stage 4 — Report assembly results to stdout.
 *
 * Compares the assembled string against the original sequence with a
 * simple exact match. Rotation normalization is handled in assembleGenome(),
 * so no rotation check is needed here.
 *
 * @param original  The original DNASequence loaded from file.
 * @param assembled The assembled string produced by assembleGenome().
 */
static void reportResults(const DNASequence& original, const std::string& assembled) {
    bool match = assembled == original.getSequence();

    std::cout << "ASSEMBLED LENGTH: " << assembled.length() << " bases\n";
    std::cout << "RECONSTRUCTION:   " << (match ? "SUCCESSFUL" : "FAILED") << "\n";
    std::cout << "--------------------------------------\n";
}

// -----------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------

int main() {

    try {
        //const std::string path = "../Data/Escherichia coli str. K-12 substr..fna";
        //const std::string path = "../Data/small_test.fna";
        //const std::string path = "../Data/Escherichia_phage_phiX174.fna"; // Size 5386
        //const std::string path = "../Data/Saccharomyces cerevisiae S288C chromosome I.fna"; // Size 230218
        //const std::string path = "../Data/Mycoplasma_genitalium_G37.fna";
        const std::string path = "../Data/Escherichia_phage_Lambda.fna";

        // Stage 1: Load
        DNASequence genome = loadGenome(path);
        std::cout << genome.getName() << "\n";
        std::cout << "Sequence length: " << genome.getLength() << " bases\n";

        // Stage 2: Build
        //int k = getIntFromUser();
        for (int k = 20; k <= 63; ++k) {
            std::cout << "Kmer size: " << k << "\n";

            DeBruijnGraph graph = buildGraph(genome.getSequence(), k);

            // Stage 3: Assemble
            std::string assembled = assembleGenome(graph, genome.getSequence(), k);

            // Stage 4: Report
            reportResults(genome, assembled);
        }

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
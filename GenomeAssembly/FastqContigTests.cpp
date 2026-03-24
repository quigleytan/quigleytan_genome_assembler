#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

#include "DeBruijnGraph.h"
#include "ContigTraversal.h"
#include "DataInitialization/SequenceReader.h"
#include "DataProcessing/KmerEncoding.h"
#include "DataProcessing/KmerTable.h"

// ─────────────────────────────────────────────
// Flip to false once you've observed expected values
// and filled in the TODOs below.
// ─────────────────────────────────────────────
static constexpr bool DISCOVER_MODE = true;

// ─────────────────────────────────────────────
// EXPECTED VALUES — fill these in after first run
// ─────────────────────────────────────────────
struct Expected {
    size_t minContigs;  // lower bound — at least this many
    size_t maxContigs;  // upper bound — no more than this many
};

// TODO: fill in after first run with DISCOVER_MODE = true
static constexpr Expected EXPECTED_K3 = { 0, 999 };  // placeholder
static constexpr Expected EXPECTED_K4 = { 0, 999 };  // placeholder
static constexpr Expected EXPECTED_K5 = { 0, 999 };  // placeholder

// ─────────────────────────────────────────────
// HELPERS
// ─────────────────────────────────────────────

struct RunResult {
    size_t k;
    size_t contigCount;
    size_t circularCount;
    size_t shortContigCount;  // contigs with length < k
    size_t totalBases;
    size_t shortestContig;
    size_t longestContig;
};

/**
 * @brief Builds a full pipeline from a stream: reads → KmerTable → graph → contigs.
 * Re-opens the stream each call since encodeAllReads() consumes it.
 */
static RunResult runPipeline(const std::string& fastqPath, size_t k, size_t totalBases) {
    std::ifstream file(fastqPath);
    if (!file.is_open())
        throw std::runtime_error("Could not open: " + fastqPath);

    KmerTable kTable(totalBases, k);
    SequenceReader::encodeAllReads(file, k, kTable);
    file.close();

    DeBruijnGraph graph(k, kTable.getNumItems() * 2);
    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    ContigTraversal ct(graph);
    ct.computeContigs();
    const auto& contigs = ct.getContigs();

    RunResult result;
    result.k            = k;
    result.contigCount  = contigs.size();
    result.circularCount  = 0;
    result.shortContigCount = 0;
    result.totalBases   = 0;
    result.shortestContig = SIZE_MAX;
    result.longestContig  = 0;

    for (const auto& contig : contigs) {
        size_t len = contig.sequence.length();
        result.totalBases += len;
        if (contig.isCircular)    ++result.circularCount;
        if (len < k)              ++result.shortContigCount;
        if (len < result.shortestContig) result.shortestContig = len;
        if (len > result.longestContig)  result.longestContig  = len;
    }

    if (contigs.empty()) {
        result.shortestContig = 0;
        result.longestContig  = 0;
    }

    return result;
}

static void printResult(const RunResult& r) {
    std::cout << "  k=" << r.k << "\n";
    std::cout << "    Contigs:         " << r.contigCount     << "\n";
    std::cout << "    Circular:        " << r.circularCount   << "\n";
    std::cout << "    Shorter than k:  " << r.shortContigCount << "\n";
    std::cout << "    Total bases:     " << r.totalBases      << "\n";
    std::cout << "    Shortest contig: " << r.shortestContig  << "\n";
    std::cout << "    Longest contig:  " << r.longestContig   << "\n";
}

// ─────────────────────────────────────────────
// TEST FUNCTIONS
// ─────────────────────────────────────────────

/**
 * @brief Test 1: No contigs shorter than k bases.
 *
 * Provable from first principles — a contig must contain at least one
 * full k-mer, which requires at least k bases. Asserted at all k values.
 */
bool testNoShortContigs(const RunResult& result) {
    bool passed = true;

    if (result.shortContigCount > 0) {
        passed = false;
        std::cout << "  [FAIL] k=" << result.k << ": "
                  << result.shortContigCount
                  << " contig(s) shorter than k="
                  << result.k << "\n";
    }

    return passed;
}

/**
 * @brief Test 2: No circular contigs from linear reads.
 *
 * test_reads.fastq contains only linear reads — no circularized sequences.
 * A circular contig requires a node to appear as both the start and end
 * of a walk, which should not occur from these inputs.
 */
bool testNoCircularContigs(const RunResult& result) {
    bool passed = true;

    if (result.circularCount > 0) {
        passed = false;
        std::cout << "  [FAIL] k=" << result.k << ": "
                  << result.circularCount
                  << " circular contig(s) found from linear reads\n";
    }

    return passed;
}

/**
 * @brief Test 3: Contig count is within expected range.
 *
 * Bounds filled in after first DISCOVER_MODE run.
 * Skipped entirely in DISCOVER_MODE.
 */
bool testContigCountInRange(const RunResult& result, const Expected& expected) {
    bool passed = true;

    if (result.contigCount < expected.minContigs ||
        result.contigCount > expected.maxContigs) {
        passed = false;
        std::cout << "  [FAIL] k=" << result.k
                  << ": contig count " << result.contigCount
                  << " outside expected range ["
                  << expected.minContigs << ", "
                  << expected.maxContigs << "]\n";
    }

    return passed;
}

// ─────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────

int main() {
    const std::string path = "../Data/[7]short reads test.fastq";

    // Total bases across 3 reads of length 8 = 24
    const size_t totalBases = 24;

    const std::vector<size_t> kValues = { 3, 4, 5 };

    // Collect results for all k values up front
    std::vector<RunResult> results;
    for (size_t k : kValues)
        results.push_back(runPipeline(path, k, totalBases));

    // ── DISCOVER MODE ─────────────────────────
    if (DISCOVER_MODE) {
        std::cout << "=== DISCOVER MODE — observed values ===\n";
        std::cout << "Fill these into the Expected structs,\n";
        std::cout << "then set DISCOVER_MODE = false.\n";
        std::cout << "---------------------------------------\n";
        for (const auto& r : results)
            printResult(r);
        std::cout << "---------------------------------------\n";
    }

    // ── PROVABLE ASSERTIONS (always run) ──────
    std::cout << "\n=== Provable assertions ===\n";
    bool allPassed = true;

    for (const auto& r : results) {
        if (!testNoShortContigs(r))    allPassed = false;
        if (!testNoCircularContigs(r)) allPassed = false;
    }

    // ── RANGE ASSERTIONS (only when locked in) ─
    if (!DISCOVER_MODE) {
        std::cout << "\n=== Range assertions ===\n";
        const std::vector<Expected> expectedValues = {
            EXPECTED_K3, EXPECTED_K4, EXPECTED_K5
        };
        for (size_t i = 0; i < results.size(); ++i) {
            if (!testContigCountInRange(results[i], expectedValues[i]))
                allPassed = false;
        }
    }

    // ── SUMMARY ───────────────────────────────
    std::cout << "\n";
    if (allPassed)
        std::cout << "All FastqContig Tests Passed\n";
    else
        std::cout << "Some FastqContig Tests FAILED\n";

    return allPassed ? 0 : 1;
}
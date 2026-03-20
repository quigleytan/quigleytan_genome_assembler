#include <iostream>
#include <ostream>
#include <vector>
#include <stdexcept>

#include "DeBruijnGraph.h"
#include "DataInitialization/DNASequence.h"
#include "DataProcessing/KmerEncoding.h"
#include "DataProcessing/KmerTable.h"
#include "EulerianPath.h"

bool DebruijnGraphTests();
bool EulerianPathTests();

// Copied from CleanSequenceTraversal — checks if assembled is a rotation of original.
// Used for circuit tests where start node is non-deterministic.
static bool isRotation(const std::string& original, const std::string& assembled) {
    if (original.length() != assembled.length()) return false;
    if (original == assembled) return true;

    size_t n = original.length();
    std::string doubled = original + original;

    std::vector<size_t> fail(n, 0);
    size_t len = 0, i = 1;
    while (i < n) {
        if (assembled[i] == assembled[len])      { fail[i++] = ++len; }
        else if (len != 0)                        { len = fail[len - 1]; }
        else                                      { fail[i++] = 0; }
    }

    size_t j = 0;
    for (size_t k = 0; k < doubled.length(); ++k) {
        while (j > 0 && doubled[k] != assembled[j]) j = fail[j - 1];
        if (doubled[k] == assembled[j]) ++j;
        if (j == n) return true;
    }
    return false;
}

// Shared pipeline: encodes a sequence into a KmerTable, builds a DeBruijnGraph,
// and returns it ready for EulerianPath.
static DeBruijnGraph buildGraph(const std::string& sequence, int k) {
    KmerTable kTable(sequence.length(), k);
    KmerEncoding::encodeSequence(sequence, k, kTable);

    DeBruijnGraph graph(k);
    for (const auto& entry : kTable) {
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);
    }
    return graph;
}

int main() {

    if (DebruijnGraphTests()) {
        std::cout << "De Bruijn Graph Tests Passed" << std::endl;
    }

    if (EulerianPathTests()) {
        std::cout << "Eulerian Path Tests Passed" << std::endl;
    }
}

bool DebruijnGraphTests() {
    bool passed = true;

    // Starting the initial sequence
    DNASequence genome("Test Sequence", "AGTGCGTCAGT");
    int k = 3;

    // Entering information into the KmerTable
    KmerTable kTable(genome.getLength(), k);
    KmerEncoding::encodeSequence(genome.getSequence(), k, kTable);

    // Initializing the De Bruijn graph
    DeBruijnGraph graph(k);

    // Inserting kmer information into the table
    for (const auto& entry : kTable) {
        uint64_t graphKmer = entry.key;
        size_t count = entry.value;

        for (size_t i = 0; i < count; ++i)
            graph.addKmer(graphKmer);
    }

    if (graph.getEdgeCount() != 9) {
        passed = false;
        std::cout << "Edge count incorrect" << std::endl;
    }

    if (graph.getNodeCount() != 7) {
        passed = false;
        std::cout << "Node count incorrect" << std::endl;
    }

    // Checking AG node logic
    if (graph.findNode(KmerEncoding::encode("AG"))->getInDegree() != 1) {
        passed = false;
        std::cout << "AG In degree incorrect" << std::endl;
    }

    if (graph.findNode(KmerEncoding::encode("AG"))->getOutDegree() != 2) {
        passed = false;
        std::cout << "AG Out degree incorrect" << std::endl;
    }

    // Checking GT node logic
    if (graph.findNode(KmerEncoding::encode("GT"))->getInDegree() != 3) {
        passed = false;
        std::cout << "GT In degree incorrect" << std::endl;
    }

    if (graph.findNode(KmerEncoding::encode("GT"))->getOutDegree() != 2) {
        passed = false;
        std::cout << "GT Out degree incorrect" << std::endl;
    }

    // The remaining nodes should have in/out degrees of 1
    std::vector<NodeId> remainingNodes = {
        KmerEncoding::encode("TG"), KmerEncoding::encode("GC"),
        KmerEncoding::encode("CG"), KmerEncoding::encode("TC"),
        KmerEncoding::encode("CA")
    };

    for (uint64_t node : remainingNodes) {
        if (graph.findNode(node)->getInDegree() != 1 ||
            graph.findNode(node)->getOutDegree() != 1) {
            passed = false;
            std::cout << "Node " << node << " has incorrect in/out degree" << std::endl;
            }
    }

    return passed;
}

bool EulerianPathTests() {
    bool passed = true;

    // -----------------------------------------------------------------------
    // Test 1: Simple linear path — "ACGTAC", k=3
    //
    // k-mers: ACG, CGT, GTA, TAC
    // Nodes (k-1 mers): AC, CG, GT, TA
    //
    // Degree breakdown:
    //   AC: out=1 (ACG), in=1 (TAC)  -- but ACG gives prefix AC (no in-edge at start)
    //   Start node: AC (out - in == 1), End node: AC after TAC edge
    //
    // Expected: path visits 5 nodes (4 edges + 1), assembled == "ACGTAC"
    // -----------------------------------------------------------------------
    {
        std::string sequence = "ACGTAC";
        int k = 3;

        DeBruijnGraph graph = buildGraph(sequence, k);
        EulerianPath ep(graph);
        ep.computePath();

        // 4 edges → 5 nodes in path
        if (ep.getPath().size() != 5) {
            passed = false;
            std::cout << "[Test 1] Path length incorrect: expected 5, got "
                      << ep.getPath().size() << std::endl;
        }

        std::string assembled = ep.reconstructGenome(false);
        if (assembled != sequence) {
            passed = false;
            std::cout << "[Test 1] Reconstruction incorrect\n"
                      << "  Expected: " << sequence  << "\n"
                      << "  Got:      " << assembled << std::endl;
        }
    }

    // -----------------------------------------------------------------------
    // Test 2: Eulerian circuit via circularization — "ACGT", k=3
    //
    // Circularized: "ACGT" + "AC" = "ACGTAC"
    // All nodes balanced (in == out), so findStartNode() falls back to
    // nodes.front() — start is non-deterministic due to hash table ordering.
    // We therefore accept exact match OR a valid rotation of the original.
    //
    // Expected: path.front() == path.back() (circuit confirmed),
    //           assembled matches original or is a rotation of it
    // -----------------------------------------------------------------------
    {
        std::string original     = "ACGT";
        int k = 3;
        std::string circularized = original + original.substr(0, k - 1); // "ACGTAC"

        DeBruijnGraph graph = buildGraph(circularized, k);
        EulerianPath ep(graph);
        ep.computePath();

        // Circuit: first and last node in path must be identical
        if (ep.getPath().front() != ep.getPath().back()) {
            passed = false;
            std::cout << "[Test 2] Expected a circuit (path.front == path.back)" << std::endl;
        }

        // Trim reconstructed string back to original length
        std::string assembled = ep.reconstructGenome(true);
        if (assembled.length() > original.length())
            assembled = assembled.substr(0, original.length());

        bool exactMatch = (assembled == original);
        bool rotation   = !exactMatch && isRotation(original, assembled);

        if (!exactMatch && !rotation) {
            passed = false;
            std::cout << "[Test 2] Reconstruction incorrect\n"
                      << "  Expected: " << original  << " (or rotation)\n"
                      << "  Got:      " << assembled << std::endl;
        }
    }

    // -----------------------------------------------------------------------
    // Test 3: Repeated k-mers — "AGTGCGTCAGT", k=3
    //
    // Same sequence as DeBruijnGraphTests (graph structure already verified).
    // GT appears as a prefix/suffix multiple times, stressing edge-consumption
    // in Hierholzer's traversal.
    //
    // Graph has 9 edges → path must visit 10 nodes.
    // AG (out=2, in=1) is the unique start node, so reconstruction is deterministic.
    //
    // Expected: path length == 10, assembled == "AGTGCGTCAGT"
    // -----------------------------------------------------------------------
    {
        std::string sequence = "AGTGCGTCAGT";
        int k = 3;

        DeBruijnGraph graph = buildGraph(sequence, k);
        EulerianPath ep(graph);
        ep.computePath();

        // 9 edges → 10 nodes in path
        if (ep.getPath().size() != 10) {
            passed = false;
            std::cout << "[Test 3] Path length incorrect: expected 10, got "
                      << ep.getPath().size() << std::endl;
        }

        std::string assembled = ep.reconstructGenome(false);
        if (assembled != sequence) {
            passed = false;
            std::cout << "[Test 3] Reconstruction incorrect\n"
                      << "  Expected: " << sequence  << "\n"
                      << "  Got:      " << assembled << std::endl;
        }
    }

    return passed;
}
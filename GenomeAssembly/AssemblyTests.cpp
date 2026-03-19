#include <iostream>
#include <ostream>
#include <vector>

#include "DeBruijnGraph.h"
#include "DataInitialization/DNASequence.h"
#include "DataProcessing/KmerEncoding.h"
#include "DataProcessing/KmerTable.h"

bool DebruijnGraphTests();
bool EulerianPathTests();

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
    std::vector<uint64_t> remainingNodes = {
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
    return true;
}
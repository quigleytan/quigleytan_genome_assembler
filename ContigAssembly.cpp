#include <fstream>
#include <iostream>

#include "Data/DataList.h"
#include "DataInitialization/DNASequence.h"
#include "DataInitialization/SequenceReader.h"
#include "DataProcessing/KmerEncoding.h"
#include "DataProcessing/KmerTable.h"
#include "GenomeAssembly/DeBruijnGraph.h"
#include "GenomeAssembly/ContigTraversal.h"

// -----------------------------------------------------------------------
// Pipeline stages
// -----------------------------------------------------------------------

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

static void assembleContigs(DeBruijnGraph& graph, bool overlap) {
    ContigTraversal ct(graph);
    if (overlap)
        ct.setOverlap(true);
    ct.computeContigs();
    ct.printStats();
}

// -----------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------

int main() {
    try {
        const std::string path = "../Data/" + sequence[6];


        // Stage 1: Load
        DNASequence genome = loadGenome(path);
        std::cout << genome.getName() << "\n";
        std::cout << "Sequence length: " << genome.getLength() << " bases\n";

        // Stage 2: Build
        int k = 63;
        std::cout << "Kmer size: " << k << "\n";
        DeBruijnGraph graph = buildGraph(genome.getSequence(), k);

        // Stage 3 & 4: Assemble and report
        assembleContigs(graph, false);

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
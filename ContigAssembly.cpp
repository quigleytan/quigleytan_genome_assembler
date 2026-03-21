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

static DeBruijnGraph buildGraph(const std::string& sequence, int k, bool circular = false) {
    std::cout << "DEBUG: building graph with k=" << k << "\n";


    std::string inputSequence = circular
        ? sequence + sequence.substr(0, k - 1)
        : sequence;
    std::cout << "Input sequence: " << inputSequence << "\n";
    std::cout << "Input length:   " << inputSequence.length() << "\n";

    KmerTable kTable(inputSequence.length(), k);
    KmerEncoding::encodeSequence(inputSequence, k, kTable);
    std::cout << "Unique k-mers:   " << kTable.getNumItems() << "\n";

    DeBruijnGraph graph(k, kTable.getNumItems() * 2);
    for (const auto& entry : kTable) {
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);
    }
    std::cout << "Graph built:     " << graph.getNodeCount() << " nodes, "
              << graph.getEdgeCount() << " edges\n";

    size_t repeated = 0;
    size_t maxCount = 0;
    for (const auto& entry : kTable) {
        if (entry.value > 1) ++repeated;
        if (entry.value > maxCount) maxCount = entry.value;
    }
    std::cout << "Repeated k-mers:  " << repeated << "\n";
    std::cout << "Max k-mer count:  " << maxCount << "\n";


    return graph;
}

static void assembleContigs(DeBruijnGraph& graph, DNASequence sequence, int k) {
    ContigTraversal ct(graph);
    ct.computeContigs();
    // After computeContigs(), concatenate all contig sequences
    // and verify they cover the original
    std::string reconstruction = "";
    for (const auto& contig : ct.getContigs()) {
        reconstruction += contig.sequence;
    }
    std::cout << "Original:       " << sequence.getSequence() << "\n";
    std::cout << "Reconstructed:  " << reconstruction << "\n";
    std::cout << "Length match:   " << (reconstruction == sequence.getSequence() ? "YES" : "NO") << "\n";    ct.printStats();

    for (const auto& contig : ct.getContigs()) {
        std::cout << "  [" << contig.sequence << "] "
                  << "start=" << KmerEncoding::decode(contig.startNode, k-1)
                  << " end=" << KmerEncoding::decode(contig.endNode, k-1)
                  << " circular=" << contig.isCircular << "\n";
    }

    std::cout << "--------------------------------------\n";
    std::cout << "\n";
}

// -----------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------

int main() {
    try {
        const std::string path = "../Data/" + sequence[7];

        DNASequence genome = loadGenome(path);
        std::cout << genome.getName() << "\n";
        std::cout << "Sequence length: " << genome.getLength() << " bases\n";

        for (int k : {4, 8, 16}) {
            std::cout << "Kmer size: " << k << "\n";
            DeBruijnGraph graph = buildGraph(genome.getSequence(), k, false);

            // Debug: print all node degrees
            for (NodeId node : graph.getAllNodes()) {
                const auto* data = graph.findNode(node);
                std::cout << KmerEncoding::decode(node, graph.getK() - 1)
                          << " in=" << data->getInDegree()
                          << " out=" << data->getOutDegree() << "\n";
            }

            assembleContigs(graph, genome, k);
        }

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
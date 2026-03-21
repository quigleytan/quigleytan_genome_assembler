#include "ContigTraversal.h"
#include <algorithm>
#include <iostream>

// PRIVATE

void ContigTraversal::initializeAdjacency() {
    adjCopy = OpenAddressingTable<NodeId, std::vector<NodeId>>(graph.getNodeCount() * 2);

    for (NodeId node : graph.getAllNodes()) {
        const auto* data = graph.findNode(node);
        auto [neighborRef, isNew] = adjCopy.insert(node);
        neighborRef = data->getNeighbors();
    }
}

bool ContigTraversal::isAmbiguous(NodeId node) const {
    const auto* data = graph.findNode(node);
    return data->getInDegree() != 1 || data->getOutDegree() != 1;
}

std::string ContigTraversal::walkContig(NodeId startNode) {
    std::string result = KmerEncoding::decode(startNode, graph.getK() - 1);
    NodeId currentNode = startNode;

    while (true) {
        auto* neighbors = adjCopy.find(currentNode);
        if (!neighbors || neighbors->empty()) break;

        NodeId next = neighbors->back();
        neighbors->pop_back();

        result += KmerEncoding::decode(next, graph.getK() - 1).back();

        if (isAmbiguous(next)) break;
        currentNode = next;
    }

    return result;
}

void ContigTraversal::handleIsolatedCycles() {
    for (NodeId node : graph.getAllNodes()) {
        auto* neighbors = adjCopy.find(node);
        while (neighbors && !neighbors->empty()) {
            contigs.push_back(walkContig(node));
        }
    }
}

// PUBLIC

// Constructor
ContigTraversal::ContigTraversal(DeBruijnGraph& g) : graph(g), adjCopy(g.getNodeCount() * 2) {}

void ContigTraversal::computeContigs() {
    initializeAdjacency();
    for (NodeId node : graph.getAllNodes()) {
        if (isAmbiguous(node)) {
            auto* neighbors = adjCopy.find(node);
            while (neighbors && !neighbors->empty()) {
                contigs.push_back(walkContig(node));
            }
        }
    }
    handleIsolatedCycles();
}

const std::vector<std::string>& ContigTraversal::getContigs() const {
    return contigs;
}

void ContigTraversal::printStats() const {
    if (contigs.empty()) {
        std::cout << "No contigs found\n";
        return;
    }

    size_t totalLength = 0;
    for (const auto& contig : contigs) totalLength += contig.length();

    auto shortest = std::min_element(contigs.begin(), contigs.end(),
        [](const std::string& a, const std::string& b) { return a.length() < b.length(); });
    auto longest = std::max_element(contigs.begin(), contigs.end(),
        [](const std::string& a, const std::string& b) { return a.length() < b.length(); });

    // N50 calculation
    std::vector<size_t> lengths;
    lengths.reserve(contigs.size());
    for (const auto& contig : contigs) lengths.push_back(contig.length());
    std::sort(lengths.rbegin(), lengths.rend());

    size_t half = totalLength / 2;
    size_t accumulated = 0;
    size_t n50 = 0;
    for (size_t len : lengths) {
        accumulated += len;
        if (accumulated >= half) { n50 = len; break; }
    }

    std::cout << "--------------------------------------\n";
    std::cout << "Total contigs:    " << contigs.size() << "\n";
    std::cout << "Total bases:      " << totalLength << "\n";
    std::cout << "Shortest contig:  " << shortest->length() << " bases\n";
    std::cout << "Longest contig:   " << longest->length() << " bases\n";
    std::cout << "N50:              " << n50 << " bases\n";
    std::cout << "--------------------------------------\n";
}
#include "DeBruijnGraph.h"

// Private
std::pair<uint64_t, uint64_t> DeBruijnGraph::chop(uint64_t kmer) const {
    uint64_t prefix = kmer >> 2; // Naturally discards
    uint64_t suffix = kmer & kMask;

    return {prefix, suffix};
}

// Constructor
DeBruijnGraph::DeBruijnGraph(size_t k) {
    if (k < 2)
        throw std::invalid_argument("k must be > 1");
    this->k = k;
    this->kMask = (1ULL << (2 * (k - 1))) - 1;
}

// Getters

size_t DeBruijnGraph::getNodeCount() const {
    return nodeCount;
}

size_t DeBruijnGraph::getEdgeCount() const {
    return edgeCount;
}

OpenAddressingTable<uint64_t, DeBruijnGraph::NodeData>& DeBruijnGraph::getGraph() {
    return this->table;
}

bool DeBruijnGraph::contains(uint64_t node) const {
    return table.find(node) != nullptr;
}

const std::vector<uint64_t>& DeBruijnGraph::getNeighbors(uint64_t node) const {
    auto* nodeData = table.find(node);
    if (!nodeData) {
        throw NodeNotFoundException(node);
    }
    return nodeData->neighbors;
}

void DeBruijnGraph::addKmer(uint64_t inputKmer) {
    // Derives k-1 mers from the input
    auto [prefix, suffix] = chop(inputKmer);

    auto [from, prefixNew] = table.insert(prefix);
    auto [to,   suffixNew] = table.insert(suffix);

    // Updates the number of nodes in the graph if it is a new instance of a k-1 mer.
    if (prefixNew) nodeCount++;
    if (suffixNew) nodeCount++;

    from.neighbors.push_back(suffix);
    to.inDegree++;
    edgeCount++;
}

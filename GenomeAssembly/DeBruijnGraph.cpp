#include "DeBruijnGraph.h"

// Private
std::pair<uint64_t, uint64_t> DeBruijnGraph::chop(uint64_t kmer) const {
    uint64_t prefix = kmer >> 2; // Naturally discards
    uint64_t suffix = kmer & kMask_;

    return {prefix, suffix};
}

// Constructor
DeBruijnGraph::DeBruijnGraph(size_t k)
    : k_(EncodingUtils::validateK(k)),
      kMask_(EncodingUtils::makeMask(k_ - 1)) {}

// Getters

size_t DeBruijnGraph::getNodeCount() const {
    return nodeCount_;
}

size_t DeBruijnGraph::getEdgeCount() const {
    return edgeCount_;
}

const DeBruijnGraph::NodeData *DeBruijnGraph::findNode(NodeId node) const{
    return table.find(node);
}

void DeBruijnGraph::addKmer(uint64_t inputKmer) {
    // Derives k-1 mers from the input
    auto [prefix, suffix] = chop(inputKmer);

    auto [from, prefixNew] = table.insert(prefix);
    auto [to,   suffixNew] = table.insert(suffix);

    // Updates the number of nodes in the graph if it is a new instance of a k-1 mer.
    if (prefixNew) nodeCount_++;
    if (suffixNew) nodeCount_++;

    from.addNeighbor(suffix);
    to.incrementInDegree();
    edgeCount_++;
}

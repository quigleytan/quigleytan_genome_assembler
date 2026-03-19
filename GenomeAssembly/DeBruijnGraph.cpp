#include "DeBruijnGraph.h"
#include "DataProcessing/KmerEncoding.h"

// Private
std::pair<uint64_t, uint64_t> DeBruijnGraph::chop(uint64_t kmer) const {
    uint64_t prefix = kmer >> 2; // Naturally discards
    uint64_t suffix = kmer & kMask_;

    return {prefix, suffix};
}

// Constructor
DeBruijnGraph::DeBruijnGraph(size_t k)
    : k(KmerEncoding::validateK(k)),
      kMask_(KmerEncoding::bitmask(k - 1)) {} // Bitmasks with k-1 for proper k-1 mer sizing.

// Getters

size_t DeBruijnGraph::getNodeCount() const {
    return nodeCount_;
}

size_t DeBruijnGraph::getEdgeCount() const {
    return edgeCount_;
}

const DeBruijnGraph::NodeData *DeBruijnGraph::findNode(NodeId node) const{
    return table_.find(node);
}

void DeBruijnGraph::addKmer(uint64_t inputKmer) {
    // Split k-mer into prefix and suffix (k-1)-mers
    auto [prefix, suffix] = chop(inputKmer);

    auto [from, prefixNew] = table_.insert(prefix);
    auto [to,   suffixNew] = table_.insert(suffix);

    // Updates the number of nodes in the graph if it is a new instance of a k-1 mer.
    if (prefixNew) nodeCount_++;
    if (suffixNew) nodeCount_++;

    from.addNeighbor(suffix);
    to.incrementInDegree();
    edgeCount_++;
}

std::vector<DeBruijnGraph::NodeId> DeBruijnGraph::getAllNodes() const {
    // Initializes the return vector.
    std::vector<NodeId> nodes;
    nodes.reserve(nodeCount_);
    // Populates the vector with all the item keys.
    for (auto it = table_.begin(); it != table_.end(); ++it) {
        nodes.push_back((*it).key);
    }
    return nodes;
}

void DeBruijnGraph::printGraph() const {
    std::cout << "De Bruijn Graph\n";
    std::cout << "--------------------------------------\n";
    auto nodes = getAllNodes();

    for (NodeId node : nodes) {

        std::cout << KmerEncoding::decode(node, k - 1)
                  << " | in: " << findNode(node)->getInDegree()
                  << " | out: " << findNode(node)->getOutDegree()
                  << " | -> ";

        auto neighbors = findNode(node)->getNeighbors();

        for (NodeId neighbor : neighbors) {
            std::cout << KmerEncoding::decode(neighbor, k - 1) << " ";
        }

        std::cout << "\n";
    }

    std::cout << "--------------------------------------\n";
}

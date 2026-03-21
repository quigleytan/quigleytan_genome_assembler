#include "ContigTraversal.h"
#include <algorithm>
#include <iostream>

// PRIVATE

void ContigTraversal::initializeAdjacency() {
    adjCopy_ = OpenAddressingTable<NodeId, std::vector<NodeId>>(graph_.getNodeCount() * 2);

    for (NodeId node : graph_.getAllNodes()) {
        const auto* data = graph_.findNode(node);
        auto [neighborRef, isNew] = adjCopy_.insert(node);
        neighborRef = data->getNeighbors();
    }
}

bool ContigTraversal::isAmbiguous(NodeId node) const {
    const auto* data = graph_.findNode(node);
    return data->getInDegree() != 1 || data->getOutDegree() != 1;
}

std::string ContigTraversal::walkContig(NodeId startNode) {
    auto* neighbors = adjCopy_.find(startNode);
    if (!neighbors || neighbors->empty()) return "";

    NodeId firstStep = neighbors->back();
    neighbors->pop_back();

    // each node contributes exactly 1 new base in non-overlap mode
    std::string result;
    result += KmerEncoding::decode(firstStep, graph_.getK() - 1).back();
    NodeId currentNode = firstStep;

    if (isAmbiguous(firstStep)) return result;

    while (true) {
        auto* nextNeighbors = adjCopy_.find(currentNode);
        if (!nextNeighbors || nextNeighbors->empty()) break;

        NodeId next = nextNeighbors->back();
        nextNeighbors->pop_back();

        if (isAmbiguous(next)) {
            if (overlap_)
                result += KmerEncoding::decode(next, graph_.getK() - 1).back();
            break;
        }

        result += KmerEncoding::decode(next, graph_.getK() - 1).back();
        currentNode = next;
    }

    return result;
}

void ContigTraversal::handleIsolatedCycles() {
    for (NodeId node : graph_.getAllNodes()) {
        auto* neighbors = adjCopy_.find(node);
        while (neighbors && !neighbors->empty()) {
            std::string contig = walkContig(node);
            contigs_.push_back(contig);
        }
    }
}

// PUBLIC

// Constructor
ContigTraversal::ContigTraversal(DeBruijnGraph& g) : graph_(g),
    adjCopy_(g.getNodeCount() * 2), overlap_(false) {}

void ContigTraversal::computeContigs() {
    contigs_.clear();
    initializeAdjacency();

    size_t branchPoints = 0;
    size_t branchContigs = 0;

    for (NodeId node : graph_.getAllNodes()) {
        if (isAmbiguous(node)) {
            branchPoints++;
            auto* neighbors = adjCopy_.find(node);
            while (neighbors && !neighbors->empty()) {
                std::string contig = walkContig(node);
                contigs_.push_back(contig);
                branchContigs++;

            }
        }
    }

    std::cout << "Branch points found: " << branchPoints << "\n";
    std::cout << "Branch contigs:      " << branchContigs << "\n";

    handleIsolatedCycles();

    std::cout << "Total after cycles:  " << contigs_.size() << "\n";
}

void ContigTraversal::setOverlap(bool logical) {
    overlap_ = logical;
}

bool ContigTraversal::getOverlap() const {
    return overlap_;
}

const std::vector<std::string>& ContigTraversal::getContigs() const {
    return contigs_;
}

void ContigTraversal::printStats() const {
    if (contigs_.empty()) {
        std::cout << "No contigs found\n";
        return;
    }

    size_t totalLength = 0;
    std::vector<size_t> lengths;
    lengths.reserve(contigs_.size());

    for (const auto& contig : contigs_) {
        size_t len = contig.length();
        lengths.push_back(len);
        totalLength += len;
    }

    auto [minIt, maxIt] = std::minmax_element(lengths.begin(), lengths.end());
    size_t minLen = *minIt;
    size_t maxLen = *maxIt;
    std::sort(lengths.rbegin(), lengths.rend());

    size_t half = (totalLength + 1) / 2;
    size_t accumulated = 0;
    size_t n50 = 0;

    for (size_t len : lengths) {
        accumulated += len;
        if (accumulated >= half) {
            n50 = len;
            break;
        }
    }

    std::cout << "--------------------------------------\n";
    std::cout << "Overlap:          " << (overlap_ ? "TRUE" : "FALSE") << "\n";
    std::cout << "Total contigs:    " << contigs_.size() << "\n";
    std::cout << "Total bases:      " << totalLength << "\n";
    std::cout << "Shortest contig:  " << minLen << " bases\n";
    std::cout << "Longest contig:   " << maxLen << " bases\n";
    std::cout << "N50:              " << n50 << " bases\n";
    std::cout << "--------------------------------------\n";

    size_t len1 = 0, len2 = 0, lenk = 0, longer = 0;
    for (const auto& c : contigs_) {
        if (c.length() == 1) len1++;
        else if (c.length() == 2) len2++;
        else if (c.length() < graph_.getK()) lenk++;
        else longer++;
    }
    std::cout << "Length 1:    " << len1 << "\n";
    std::cout << "Length 2:    " << len2 << "\n";
    std::cout << "Length 2-63: " << lenk << "\n";
    std::cout << "Length 63+:  " << longer << "\n";

    std::cout << "Genome length:    " << graph_.getEdgeCount() << "\n";
    std::cout << "Total bases:      " << totalLength << "\n";
    std::cout << "Difference:       " << graph_.getEdgeCount() - totalLength << "\n";
}
#include "ContigTraversal.h"

void ContigTraversal::initializeAdjacency() {

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

}

void ContigTraversal::handleIsolatedCycles() {}

explicit ContigTraversal::ContigTraversal(DeBruijnGraph& g) {}

void ContigTraversal::computeContigs() {}

const std::vector<std::string>& ContigTraversal::getContigs() const {}

void ContigTraversal::printStats() const {}
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
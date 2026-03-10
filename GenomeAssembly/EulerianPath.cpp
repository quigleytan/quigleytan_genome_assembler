#include "EulerianPath.h"

#include <algorithm>

void EulerianPath::initializeAdjacency() {

}

[[nodiscard]] uint64_t EulerianPath::findStartNode() const {

    auto nodes = graph.getAllNodes();

    uint64_t startNode = nodes.front(); // Protects in case of cycle.

    int startCount = 0;
    int endCount = 0;

    for (uint64_t node : nodes) {

        const auto* data = graph.findNode(node);

        int diff = static_cast<int>(data->getOutDegree()) -
                   static_cast<int>(data->getInDegree());

        if (diff == 1) {
            startNode = node;
            startCount++;
        }
        else if (diff == -1) {
            endCount++;
        }
        else if (diff != 0) {
            throw std::runtime_error("Graph is not Eulerian");
        }
    }

    // Checking Eulerian conditions
    if (!((startCount == 1 && endCount == 1) || (startCount == 0 && endCount == 0)))
        throw std::runtime_error("Graph does not contain an Eulerian path");

    return startNode; // Returns a valid node in both cycle and path.
}

// Public

EulerianPath::EulerianPath(DeBruijnGraph& g) : graph(g) {

}

void EulerianPath::computePath() {

}

[[nodiscard]] const std::vector<uint64_t>& EulerianPath::getPath() const {

}

[[nodiscard]] std::string EulerianPath::reconstructGenome() const {
}
#include "EulerianPath.h"
#include <stack>
#include <algorithm>

void EulerianPath::initializeAdjacency() {
    for (NodeId node : graph.getAllNodes()) {
        const auto* data = graph.findNode(node);

        auto [neighborRef, isNew] = adjCopy.insert(node);
        neighborRef = data->getNeighbors();

        // Sort neighbors for deterministic traversal order
        std::sort(neighborRef.begin(), neighborRef.end());
    }
}

NodeId EulerianPath::findStartNode() const {

    auto nodes = graph.getAllNodes();

    NodeId startNode = nodes.front(); // Protects in case of cycle.

    int startCount = 0;
    int endCount = 0;

    for (NodeId node : nodes) {

        const auto* data = graph.findNode(node);

        int diff = static_cast<int>(data->getOutDegree()) - static_cast<int>(data->getInDegree());

        if (diff == 1) {
            startNode = node;
            startCount++;
        } else if (diff == -1) {
            endCount++;
        } else if (diff != 0) {
            throw std::runtime_error("Graph is not Eulerian");
        }
    }

    // Checking Eulerian conditions
    if (!((startCount == 1 && endCount == 1) || (startCount == 0 && endCount == 0)))
        throw std::runtime_error("Graph does not contain an Eulerian path");

    return startNode; // Returns a valid node in both cycle and path.
}

// Public methods

EulerianPath::EulerianPath(DeBruijnGraph& g) : graph(g), adjCopy(g.getNodeCount() * 2) {}

void EulerianPath::computePath() {

    path.clear();
    initializeAdjacency();

    std::stack<NodeId> stack;

    NodeId start = findStartNode();
    stack.push(start);

    // Main traversal loop.
    while (!stack.empty()) {

        // Accesses the current node and retrieve its remaining unused edges.
        NodeId currentID = stack.top();
        auto* neighbors = adjCopy.find(currentID);

        // Checks if the current node has neighbors, then takes an edge and moves to the next node.
        if (!neighbors->empty()) {
            NodeId next = neighbors->back();
            neighbors->pop_back();   // Deletes the edge immediately.
            stack.push(next);
        } else { // If no edges remain, adds the node to the path and backtrack.
            path.push_back(currentID);
            stack.pop();
        }
    }

    std::reverse(path.begin(), path.end());
}

const std::vector<NodeId>& EulerianPath::getPath() const {
    return path;
}

std::string EulerianPath::reconstructGenome(bool isCircuit) const {
    if (path.empty())
        throw std::runtime_error("Path is empty — call computePath() first");

    const size_t nodeLen = graph.getK() - 1;
    std::string genome = KmerEncoding::decode(path.front(), nodeLen);

    size_t end = isCircuit ? path.size() - 1 : path.size();

    for (size_t i = 1; i < end; ++i) {
        genome += KmerEncoding::decode(path[i], nodeLen).back();
    }
    return genome;
}
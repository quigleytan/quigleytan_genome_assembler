#include "EulerianPath.h"
#include <utility>

// private
void EulerianPath::initializeAdjacency() {
    auto& tbl = graph.getGraph();

    for (auto& entry : tbl) {
        uint64_t node = entry.key;

        auto* nodeData = tbl.find(node);
        if (!nodeData) continue;

        auto [ref, inserted] = adjCopy.insert(node);
        ref = nodeData->neighbors;  // copy vector
    }
}

bool EulerianPath::isEulerian() const {
    int startNodes = 0;
    int endNodes = 0;

    auto& tbl = graph.getGraph();

    for (auto& entry : tbl) {
        uint64_t node = entry.key;

        auto* nodeData = tbl.find(node);
        if (!nodeData) continue;

        size_t in  = nodeData->inDegree;
        size_t out = nodeData->neighbors.size();

        if (in == out)
            continue;

        if (out == in + 1)
            startNodes++;
        else if (in == out + 1)
            endNodes++;
        else
            return false;
    }

    return (startNodes == 1 && endNodes == 1) ||
           (startNodes == 0 && endNodes == 0);
}

uint64_t EulerianPath::findStartNode() const {
    auto& tbl = graph.getGraph();

    uint64_t anyNodeWithOut = 0;
    bool haveAnyOut = false;

    for (auto& entry : tbl) {
        uint64_t node = entry.key;

        auto* nodeData = tbl.find(node);
        if (!nodeData) continue;

        size_t in  = nodeData->inDegree;
        size_t out = nodeData->neighbors.size();

        if (!haveAnyOut && out > 0) {
            anyNodeWithOut = node;
            haveAnyOut = true;
        }

        if (out == in + 1)
            return node;
    }

    if (haveAnyOut)
        return anyNodeWithOut;

    throw std::runtime_error("No valid start node found.");
}

// public
EulerianPath::EulerianPath(DeBruijnGraph graph)
    : graph(std::move(graph)) {
}

std::vector<uint64_t> EulerianPath::compute() const {


    return path;
}
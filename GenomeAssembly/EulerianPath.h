/*
 * EulerianPath.h
 * Summary:
 * - Implementation of Hierholzer’s Algorithm
 * Features:
 * -
 * Additional:
 * -
 */

#ifndef EULERIAN_PATH_ALGORITHM_H
#define EULERIAN_PATH_ALGORITHM_H

#include "GenomeAssembly/DeBruijnGraph.h"

class EulerianPath {

private:

    // Graph reference
    DeBruijnGraph& graph;

    // Resulting Eulerian path
    std::vector<uint64_t> path;

    // Mutable adjacency list for traversal
    OpenAddressingTable<uint64_t, std::vector<uint64_t>> adjCopy;

    // Number of edges
    size_t edgeCount = 0;

    // Methods
    void initializeAdjacency();

    bool isEulerian() const;

    uint64_t findStartNode() const;

    std::vector<uint64_t> runHierholzer(uint64_t startNode);

public:

    EulerianPath(DeBruijnGraph& graph);

    std::vector<uint64_t> compute();
};

#endif //EULERIAN_PATH_ALGORITHM_H
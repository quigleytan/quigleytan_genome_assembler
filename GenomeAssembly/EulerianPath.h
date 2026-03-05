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

    /**
     *
     */
    void initializeAdjacency();

    /**
     *
     * @return
     */
    bool isEulerian() const;

    /**
     *
     * @return
     */
    uint64_t findStartNode() const;

    /**
     *
     * @param startNode
     * @return
     */
    std::vector<uint64_t> runHierholzer(uint64_t startNode);

public:
    /**
     *
     * @param graph
     */
    EulerianPath(DeBruijnGraph& graph);

    /**
     *
     * @return
     */
    std::vector<uint64_t> compute();
};

#endif //EULERIAN_PATH_ALGORITHM_H
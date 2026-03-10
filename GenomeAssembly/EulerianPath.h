/*
 * EulerianPath.h
 * Summary:
 * - Implementation of modified Hierholzer’s Algorithm
 * Features:
 * -
 * Additional:
 * -
 */

#ifndef EULERIAN_PATH_H
#define EULERIAN_PATH_H

#include <vector>
#include "GenomeAssembly/DeBruijnGraph.h"

using NodeId = uint64_t;

class EulerianPath {

private:

    DeBruijnGraph& graph;

    // copy of adjacency list so edges can be consumed
    OpenAddressingTable<NodeId, std::vector<NodeId>> adjCopy;

    // final path
    std::vector<uint64_t> path;

    // helper functions
    void initializeAdjacency();

    [[nodiscard]] uint64_t findStartNode() const;

public:

    EulerianPath(DeBruijnGraph& g);

    void computePath();

    [[nodiscard]] const std::vector<uint64_t>& getPath() const;

    [[nodiscard]] std::string reconstructGenome() const;
};

#endif
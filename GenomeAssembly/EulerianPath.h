/*
 * EulerianPath.h
 * Summary:
 * - Implementation of modified Hierholzer's Algorithm
 * Features:
 * -
 * Additional:
 * -
 */

#ifndef EULERIAN_PATH_H
#define EULERIAN_PATH_H

#include <vector>
#include "GenomeAssembly/DeBruijnGraph.h"
#include "DataProcessing/OpenAddressingTable.h"
#include "DataProcessing/KmerEncoding.h"

using NodeId = uint64_t;

class EulerianPath {

private:

    DeBruijnGraph& graph;

    // copy of the adjacency list so edges can be consumed
    OpenAddressingTable<NodeId, std::vector<NodeId>> adjCopy;

    // final path
    std::vector<uint64_t> path;

    // helper functions
    void initializeAdjacency();

    [[nodiscard]] uint64_t findStartNode() const;

public:

    EulerianPath(DeBruijnGraph& g);

    /**
     * @breif Computes an Eulerian Path using Hierholzer's algorithm
     *
     * Sets the result of the Eulerian path search to the class member path. Conducts a walk through the
     * de Bruijn graph and uses each edge exactly once. Uses preloaded table of neighbor lists when accessing
     * each node.
     *
     */
    void computePath();

    /**
     *
     * @return
     */
    [[nodiscard]] const std::vector<uint64_t>& getPath() const;

    /**
     *
     * @param isCircular
     * @return
     */
    [[nodiscard]] std::string reconstructGenome(bool isCircular) const;
};

#endif
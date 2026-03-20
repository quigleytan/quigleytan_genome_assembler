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
#include "KmerTypes.h"
#include "GenomeAssembly/DeBruijnGraph.h"
#include "DataProcessing/OpenAddressingTable.h"
#include "DataProcessing/KmerEncoding.h"


class EulerianTraversal {

private:

    DeBruijnGraph& graph;

    OpenAddressingTable<NodeId, std::vector<NodeId>> adjCopy;

    std::vector<NodeId> path;

    void initializeAdjacency();

    [[nodiscard]] NodeId findStartNode() const;

public:

    EulerianTraversal(DeBruijnGraph& g);

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
    [[nodiscard]] const std::vector<NodeId>& getPath() const;

    /**
     *
     * @param isCircular
     * @return
     */
    [[nodiscard]] std::string reconstructGenome(bool isCircular) const;
};

#endif
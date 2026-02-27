/*
 * EulerianPath.h
 * Summary:
 * - Implementation of Hierholzer’s Algorithm
 * Features:
 * -
 * Additional:
 * -
 * TODO items for future work;
 */
#ifndef M30EP_TEQUIGLE_EULERIANPATH_H
#define M30EP_TEQUIGLE_EULERIANPATH_H

#include "GenomeAssembly/DeBruijnGraph.h"

class EulerianPath {

private:

    // Variables
    DeBruijnGraph& graph;
    std::vector<uint64_t> path;
    OpenAddressingTable<uint64_t, std::vector<uint64_t>> adjCopy;

    // Methods

    /**
     *
     * @return Returns a boolean
     */
    void initializeAdjacency();

    /**
     *
     * @return Returns a boolean
     */
    bool isEulerian() const;

    /**
     *
     * @return The start node.
     */
    uint64_t findStartNode() const;

    /**
     * @brief
     * @return
     */
    void runHierholzer();

public:

    /**
     * @brief
     *
     *
     *
     * @param graph The directed edge graph to be traversed through.
     */
    EulerianPath(DeBruijnGraph graph);

    /**
     * @brief
     *
     *
     *
     * @return Returns the path
     */
    std::vector<uint64_t> compute() const;

};

#endif //M30EP_TEQUIGLE_EULERIANPATH_H
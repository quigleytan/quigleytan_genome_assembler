/*
 * DeBruijnGraph.h
 * Created by Tanner Quigley on 2/15/2026
 * Summary:
 * - Representation of a De Bruijn graph to be used for Eularian walks for genome assembly.
 * - Graph nodes are k-1 mers, transitions are represented as full kmers.
 * Important notes:
 * - IMPORTANT: All kmer and node lookups should already be encoded as an u_int64.
 */

#ifndef M20EP_TEQUIGLE_DEBRUIJNGRAPH_H
#define M20EP_TEQUIGLE_DEBRUIJNGRAPH_H
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "CustomExceptions/NodeNotFoundException.h"
#include "DataProcessing/OpenAddressingTable.h"


class DeBruijnGraph {
private:

    // Variables
    size_t k;
    uint64_t kMask;
    size_t nodeCount = 0;
    size_t edgeCount = 0;

    // Note: outdegree is tracked by neighbors.size
    struct NodeData {
        std::vector<uint64_t> neighbors; // Outgoing edges only
        size_t inDegree = 0; // Number of edges entering the node
    };

    OpenAddressingTable<uint64_t, NodeData> table;

    /**
     * @brief Extracts k-1 prefix and suffix to create nodes.
     *
     * "Chops" the kmer into its prefix and suffix, allowing for storage in the table.
     * Applies a bitmask to keep the size of nodes at k-1.
     *
     * @param kmer The kmer in which the prefix and suffix will be pulled from.
     * @return Returns a tuple containing k-1 prefix and suffix.
     */
    std::pair<uint64_t, uint64_t> chop(uint64_t kmer) const;

public:

    /**
     * @brief Constructor for De Bruijn graph class.
     *
     * Initializes a De Bruijn graph for k-1 mers using size k given upon construction.
     * Contains a hash table of nodes containing neighbor and out/indegree information.
     * All kmers to be added to the graph MUST BE ENCODED prior to insertion.
     *
     * @param k The k size for the graph: k must be > 1.
     */
    DeBruijnGraph(size_t k);

    /**
     * @brief Returns the total node counts in the graph.
     * @return Number of nodes (k-1 mers).
     */
    size_t getNodeCount() const;

    /**
     * @brief Returns the total edge counts in the graph.
     * @return Number of edges (kmer transitons).
     */
    size_t getEdgeCount() const;

    /**
     * @brief Checks whether the node is in the table.
     * @param node The node being inquired about.
     * @return Boolean value, true if the node is found in the graph.
     */
    bool contains(uint64_t node) const;

    /**
     * @brief Returns list of neighbors.
     * @param node The k-1 mer to be checked for neighbors.
     * @return The vector list of neighbors to the k-1 mer node specified.
     */
    const std::vector<uint64_t>& getNeighbors(uint64_t node) const;

    /**
     * @brief Getter method for a node's inDegree value.
     * @param node The k-1 mer to get inDegree from.
     * @return The number of edges entering the node.
     */
    size_t getInDegree(uint64_t node) const;

    /**
     * @brief Getter method for a node's outDegree value.
     * @param node The k-1 mer to get outDegree from.
     * @return The number of edges leaving the node.
     */
    size_t getOutDegree(uint64_t node) const;

    /**
     * @brief
     *
     * Adds both the prefix and suffix of inputKmer to the table and adds the suffix to
     * the prefix's list of neighbors (individual adjacency list).
     *
     * @param inputKmer The kmer whose derived k-1mers will be added to the table.
     */
    void addKmer(uint64_t inputKmer);

};

#endif //M20EP_TEQUIGLE_DEBRUIJNGRAPH_H
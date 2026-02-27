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
#include <optional>
#include <vector>
#include <stdexcept>
#include "DataProcessing/KmerEncoding.h"
#include "CustomExceptions/NodeNotFoundException.h"
#include "DataProcessing/OpenAddressingTable.h"


class DeBruijnGraph {
private:

    // Variables
    const size_t k_;
    const uint64_t kMask_;

    size_t nodeCount_ = 0;
    size_t edgeCount_ = 0;

    using NodeId = uint64_t;

    // Note: outdegree is tracked by neighbors.size
    struct NodeData {
        std::vector<uint64_t> neighbors_; // Outgoing edges only
        size_t inDegree_ = 0;             // Number of edges entering the node

        const std::vector<uint64_t>& getNeighbors() const { return neighbors_; }
        std::vector<uint64_t>& getNeighbors() { return neighbors_; }

        size_t getInDegree() const { return inDegree_; }
        size_t getOutDegree() const { return neighbors_.size(); }

        void addNeighbor(uint64_t neighbor) { neighbors_.push_back(neighbor); }
        void incrementInDegree() { ++inDegree_; }
    };

    OpenAddressingTable<NodeId, NodeData> table;

    /**
     * @brief Extracts k-1 prefix and suffix to create nodes.
     *
     * "Chops" the kmer into its prefix and suffix, allowing for storage in the table.
     * Applies a bitmask to keep the size of nodes at k-1.
     *
     * @param kmer The kmer in which the prefix and suffix will be pulled from.
     * @return Returns a tuple containing k-1 prefix and suffix.
     */
    std::pair<NodeId, NodeId> chop(uint64_t kmer) const;

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
     * @return Number of edges (kmer transitions).
     */
    size_t getEdgeCount() const;

    /**
     * @brief Checks table to the desired node.
     * @param queriedNode Encoded node to be found in the table.
     * @return Returns a pointer to the node requested, nullptr if not found.
     */
    const NodeData *findNode(NodeId queriedNode) const;

    /**
     * @brief
     *
     * Adds both the prefix and suffix of inputKmer to the table and adds the suffix to
     * the prefix's list of neighbors (individual adjacency list).
     *
     * @param inputKmer The kmer whose derived k-1mers will be added to the table.
     */
    void addKmer(NodeId inputKmer);

};

#endif //M20EP_TEQUIGLE_DEBRUIJNGRAPH_H
/*
 * ContigScaffolder.h
 * Created by Tanner Quigley on [date]
 * Summary:
 * - Orders contigs into scaffolds using shared boundary node relationships.
 * - Supports multiple branch resolution strategies for comparison.
 * Important notes:
 * - Operates on contigs produced by ContigTraversal.
 * - Scaffolds are ordered lists of contig indices with gap estimates.
 * - Ambiguous branch points are handled according to BranchResolution strategy.
 */

#ifndef CONTIG_SCAFFOLDER_H
#define CONTIG_SCAFFOLDER_H

#include <vector>
#include <string>
#include "KmerTypes.h"
#include "GenomeAssembly/DeBruijnGraph.h"
#include "GenomeAssembly/ContigTraversal.h"
#include "DataProcessing/OpenAddressingTable.h"

// -----------------------------------------------------------------------
// Resolution strategy
// -----------------------------------------------------------------------

enum class BranchResolution {
    Skip,    // Stop at ambiguous branch points, report gap as unknown
    Greedy   // Pick the longest outgoing contig at branch points
};

// -----------------------------------------------------------------------
// Scaffold entry — one contig's position within a scaffold
// -----------------------------------------------------------------------

struct ScaffoldEntry {
    size_t contigIndex; // Index into the contigs vector
    int gapAfter;       // Estimated gap in bases between this contig
                        // and the next. -1 = unknown, 0 = direct overlap
};

// -----------------------------------------------------------------------
// Scaffold — ordered sequence of contigs
// -----------------------------------------------------------------------

struct Scaffold {
    std::vector<ScaffoldEntry> entries;
    bool isCircular = false;

    [[nodiscard]] size_t contigCount() const { return entries.size(); }
};

// -----------------------------------------------------------------------
// Scaffolder class
// -----------------------------------------------------------------------

class ContigScaffolder {

private:

    const std::vector<ContigTraversal::Contig>& contigs_;
    const DeBruijnGraph& graph_;
    BranchResolution strategy_;
    std::vector<Scaffold> scaffolds_;

    // Connection maps built from contig boundary nodes.
    // Both maps use vectors of indices to handle multiple contigs
    // sharing the same boundary node at branch points.
    OpenAddressingTable<NodeId, std::vector<size_t>> endNodeMap_;
    OpenAddressingTable<NodeId, std::vector<size_t>> startNodeMap_;

    /**
     * @brief Populates endNodeMap_ and startNodeMap_ from the contig list.
     *
     * Iterates through all contigs and inserts each contig's index into
     * endNodeMap_ keyed by its endNode, and into startNodeMap_ keyed by
     * its startNode.
     */
    void buildConnectionMap();

    /**
     * @brief Resolves the next contig index from a given boundary node.
     *
     * Looks up which contigs start at boundaryNode and applies the
     * resolution strategy to pick one:
     * - Skip:   returns npos if more than one contig starts at the node
     * - Greedy: returns the index of the longest unvisited outgoing contig
     *
     * @param boundaryNode The endNode of the current contig.
     * @param visited      Tracks which contigs are already assigned.
     * @return Index into contigs_ of the next contig, or npos if unresolvable.
     */
    [[nodiscard]] size_t resolveNext(NodeId boundaryNode,
                                     const std::vector<bool>& visited) const;

    /**
     * @brief Walks the connection map from startIndex to build one scaffold.
     *
     * Follows endNode → startNode links greedily until a dead end,
     * ambiguous branch point, or already-visited contig is reached.
     * Each step appends a ScaffoldEntry with gapAfter = 0 for direct
     * overlaps and gapAfter = -1 for unresolved terminations.
     *
     * @param startIndex Index of the contig to begin the scaffold from.
     * @param visited    Tracks which contigs are already assigned.
     * @return Populated Scaffold.
     */
    [[nodiscard]] Scaffold walkScaffold(size_t startIndex,
                                        std::vector<bool>& visited) const;

    /**
     * @brief Returns true if a contig is a scaffold entry point.
     *
     * A contig is a valid entry point if its startNode has no incoming
     * contigs in endNodeMap_ — i.e. nothing ends where it begins.
     * These are the natural starting points for linear scaffolds.
     *
     * @param contigIndex Index of the contig to check.
     * @return True if no other contig ends at this contig's startNode.
     */
    [[nodiscard]] bool isScaffoldStart(size_t contigIndex) const;

public:

    /**
     * @brief Constructor for ContigScaffolder.
     * @param contigs  Contigs produced by ContigTraversal.
     * @param graph    DeBruijnGraph used to build the contigs.
     * @param strategy Branch resolution strategy to use.
     */
    ContigScaffolder(const std::vector<ContigTraversal::Contig>& contigs,
                     const DeBruijnGraph& graph,
                     BranchResolution strategy = BranchResolution::Skip);

    /**
     * @brief Builds all scaffolds from the contig connection map.
     *
     * Phase 1: Builds a connection map from contig boundary nodes.
     * Phase 2: Walks from all scaffold entry points producing linear scaffolds.
     * Phase 3: Handles any remaining unvisited contigs as isolated scaffolds
     *          or circular scaffolds with no external entry point.
     */
    void buildScaffolds();

    /**
     * @brief Returns the list of scaffolds.
     * @return Const reference to the scaffold vector.
     */
    [[nodiscard]] const std::vector<Scaffold>& getScaffolds() const;

    /**
     * @brief Reports scaffolding statistics to stdout.
     *
     * Strategy:           Resolution strategy used.
     * Total scaffolds:    Number of scaffolds produced.
     * Circular scaffolds: Number of scaffolds forming a closed loop.
     * Longest scaffold:   Length in bases of the longest scaffold.
     * Unresolved gaps:    Number of gaps with unknown size (-1).
     * N50:                Scaffold-level N50.
     */
    void printStats() const;
};

#endif //CONTIG_SCAFFOLDER_H
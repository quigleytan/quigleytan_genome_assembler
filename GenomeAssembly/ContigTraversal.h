#ifndef CONTIG_TRAVERSAL_H
#define CONTIG_TRAVERSAL_H

#include <vector>
#include <string>
#include "KmerTypes.h"
#include "GenomeAssembly/DeBruijnGraph.h"
#include "DataProcessing/OpenAddressingTable.h"
#include "DataProcessing/KmerEncoding.h"

class ContigTraversal {

private:

    DeBruijnGraph& graph;
    OpenAddressingTable<NodeId, std::vector<NodeId>> adjCopy;
    std::vector<std::string> contigs;

    void initializeAdjacency();

    [[nodiscard]] bool isAmbiguous(NodeId node) const;

    std::string walkContig(NodeId startNode);

    void handleIsolatedCycles();

public:

    explicit ContigTraversal(DeBruijnGraph& g);

    void computeContigs();

    [[nodiscard]] const std::vector<std::string>& getContigs() const;

    void printStats() const;

};

#endif //CONTIG_TRAVERSAL_H
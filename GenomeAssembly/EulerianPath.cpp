#include "EulerianPath.h"
#include "GenomeAssembly/DeBruijnGraph.h"

void EulerianPath::initializeAdjacency();

bool EulerianPath::isEulerian() const;

uint64_t EulerianPath::findStartNode() const;


std::vector<uint64_t> EulerianPath::runHierholzer(uint64_t startNode);

//public

EulerianPath::EulerianPath(DeBruijnGraph& graph);

std::vector<uint64_t> compute();

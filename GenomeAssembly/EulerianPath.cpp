#include "EulerianPath.h"
#include <utility>

// private
void EulerianPath::initializeAdjacency() {

}

bool EulerianPath::isEulerian() const {
}

uint64_t EulerianPath::findStartNode() const {

}

// public
EulerianPath::EulerianPath(DeBruijnGraph graph)
    : graph(std::move(graph)) {
}

std::vector<uint64_t> EulerianPath::compute() const {


    return path;
}


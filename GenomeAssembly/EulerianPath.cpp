#include "EulerianPath.h"

void EulerianPath::initializeAdjacency() {

}
[[nodiscard]] bool EulerianPath::isEulerian() const {

}
[[nodiscard]] uint64_t EulerianPath::findStartNode() const {

}

// Public

EulerianPath::EulerianPath(DeBruijnGraph& g) : graph(g) {

}

void EulerianPath::computePath() {

}

[[nodiscard]] const std::vector<uint64_t>& EulerianPath::getPath() const {

}

[[nodiscard]] std::string EulerianPath::reconstructGenome() const {

}
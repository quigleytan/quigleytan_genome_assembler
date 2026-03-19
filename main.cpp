//
// Created by quigl on 1/27/2026.
//
#include <fstream>
#include <iostream>
#include <vector>

#include "CustomExceptions/DNASequenceException.h"

#include "DataInitialization/DNASequence.h"
#include "DataInitialization/SequenceReader.h"

#include "DataProcessing/KmerEncoding.h"
#include "DataProcessing/KmerTable.h"

#include "GenomeAssembly/DeBruijnGraph.h"
#include "GenomeAssembly/EulerianPath.h"

// Function declarations
int getIntFromUser();
bool isRotation(const std::string& original, const std::string& assembled);

int main() {

    std::ifstream file("../Data/genome_sample_ecoli.fna");
    //std::ifstream file("../Data/small_test.fna");
    //std::ifstream file("../Data/genome_small_test.fna");

    auto genomeOpt = SequenceReader::readFasta(file);
    DNASequence genome = *genomeOpt;
    file.close();

    std::cout << genome.getName() << "\n";
    std::cout << "Sequence length: " << genome.getLength() << " bases\n";

    int k = getIntFromUser();

    // Circularize the sequence by appending the first k-1 bases to close the genome loop
    std::string circularSequence = genome.getSequence() + genome.getSequence().substr(0, k - 1);

    // Encode circular sequence into k-mer table
    KmerTable kTable(circularSequence.length(), k);
    KmerEncoding::encodeSequence(circularSequence, k, kTable);
    std::cout << "Unique k-mers:   " << kTable.getNumItems() << "\n";

    // Build De Bruijn graph, pre-sized to avoid rehashing
    DeBruijnGraph dbGraph(k, kTable.getNumItems() * 2);
    for (const auto& entry : kTable) {
        for (size_t i = 0; i < entry.value; ++i)
            dbGraph.addKmer(entry.key);
    }
    std::cout << "Graph built:     " << dbGraph.getNodeCount() << " nodes, "
              << dbGraph.getEdgeCount() << " edges\n";

    // Eulerian path and genome reconstruction
    EulerianPath eulerianPath(dbGraph);
    try {
        eulerianPath.computePath();

        // Trim k-1 bases added by circularization
        std::string assembled = eulerianPath.reconstructGenome();
        if (assembled.length() > genome.getLength())
            assembled = assembled.substr(0, genome.getLength());

        bool match    = assembled == genome.getSequence();
        bool rotation = !match && isRotation(genome.getSequence(), assembled);

        std::cout << "Path length:      " << eulerianPath.getPath().size() << " nodes\n";
        std::cout << "Assembled length: " << assembled.length() << " bases\n";
        std::cout << "Reconstruction:   " << (match ? "YES" : rotation ? "YES (rotation)" : "NO") << "\n";

        std::cout << "--------------------------------------\n";
        /*std::cout << "Original: " << genome.getSequence() << "\n";
        std::cout << "Assembled: " << assembled << "\n";*/
    } catch (const std::runtime_error& e) {
        std::cout << "Eulerian path error: " << e.what() << "\n";
    }
    return 0;
}

int getIntFromUser() {
    while (true) {
        std::string input;
        std::cout << "Enter an integer for a desired k value: ";
        getline(std::cin, input);
        try {
            size_t size = 0;
            int value = std::stoi(input, &size);
            if (size == input.size() && value > 0) {
                return value;
            }
            std::cout << "Invalid input. ";
        }
        catch (...) {
            std::cout << "Invalid input. ";
        }
    }
}

bool isRotation(const std::string& original, const std::string& assembled) {
    if (original.length() != assembled.length()) return false;
    if (original == assembled) return true;

    size_t n = original.length();
    std::string doubled = original + original;

    // Build failure function over assembled (pattern)
    std::vector<size_t> fail(n, 0);
    size_t len = 0;
    size_t i = 1;
    while (i < n) {
        if (assembled[i] == assembled[len]) {
            fail[i++] = ++len;
        } else if (len != 0) {
            len = fail[len - 1];
        } else {
            fail[i++] = 0;
        }
    }

    // KMP search for assembled inside doubled
    size_t j = 0;
    for (size_t k = 0; k < doubled.length(); ++k) {
        while (j > 0 && doubled[k] != assembled[j])
            j = fail[j - 1];
        if (doubled[k] == assembled[j])
            ++j;
        if (j == n)
            return true;
    }
    return false;
}
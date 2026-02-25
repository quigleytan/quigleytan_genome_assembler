//
// Created by quigl on 1/27/2026.
//
#include <fstream>
#include "DataInitialization/DNASequence.h"
#include "DataProcessing/KmerEncoder.h"
#include "CustomExceptions/DNASequenceException.h"
#include "CustomExceptions/NodeNotFoundException.h"
#include "DataProcessing//KmerTable.h"
#include <ostream>

#include "GenomeAssembly/DeBruijnGraph.h"
#include "DataInitialization/SequenceReader.h"

// Functions
int getIntFromUser();
DNASequence getSequenceFromUser(std::string message, int sizeLimit);

int main() {

    std::ifstream file("../Data/genome_small_test.fna");
    //std::ifstream file("../Data/small_test.fna");

    auto genomeOpt = SequenceReader::readFasta(file);

    DNASequence genome = *genomeOpt;

    std::cout << genome.getName() << std::endl;
    file.close();

    int k = getIntFromUser();
    // Encoding sequence into k-mer table
    KmerTable kTable(genome.getLength(), k);
    KmerEncoder::encodeSequence(genome.getSequence(), k, kTable);

    // Getting k-mer from user to check count
    std::string kmer = getSequenceFromUser("Enter a k-mer sequence to check its count in the provided DNA sequence:", k).getSequence();

    const size_t* p = kTable.find(KmerEncoder::encodeKmer(kmer));
    size_t count = p ? *p : 0;
    std::cout << "Found: " << count << " instances of " << kmer << std::endl;

    DeBruijnGraph dbGraph(k);

    for (const auto& entry : kTable) {
        uint64_t graphKmer = entry.key;
        size_t count = entry.value;

        for (size_t i = 0; i < count; ++i)
            dbGraph.addKmer(graphKmer);
    }

    // Asks user for a k-1 (node) to check its carried information and neighbors
    try {
        std::string nodeStr =
            getSequenceFromUser("Enter a (k-1)-mer node to inspect in the De Bruijn graph:", k - 1).getSequence();

        uint64_t node = KmerEncoder::encodeKmer(nodeStr);

        std::cout << "Node: " << nodeStr << "\n";
        std::cout << "Encoded: " << node << "\n";

        if (!dbGraph.contains(node)) {
            std::cout << "Node not found in graph.\n";
        } else {
            size_t inDeg  = dbGraph.getInDegree(node);
            size_t outDeg = dbGraph.getOutDegree(node);
            const auto& nbrs = dbGraph.getNeighbors(node);

            std::cout << "In-degree:  " << inDeg << "\n";
            std::cout << "Out-degree: " << outDeg << "\n";
            std::cout << "Neighbors (encoded): ";
            for (uint64_t n : nbrs) std::cout << KmerEncoder::decodeKmer(n, k - 1) << " ";
            std::cout << "\n";
        }
    } catch (const NodeNotFoundException& e) {
        std::cout << "Graph lookup error: " << e.what() << "\n";
    } catch (const DNASequenceException& e) {
        std::cout << "Input error: " << e.what() << "\n";
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

DNASequence getSequenceFromUser(std::string message, int sizeLimit) {
    std::cout << message << std::endl;
    while (true) {
        std::string input;
        std::getline(std::cin, input);
        if (input.length() != sizeLimit && sizeLimit)  {
            std::cout << "Kmer size is not equal to k, please try again." << std::endl;
            continue;
        }
        try {
            DNASequence sequence("User DNA sequence",input);
            std::cout << "You entered: " << input << std::endl;
            return sequence;   // exit function once valid
        } catch (const DNASequenceException& exception) {
            std::cout << "DNA Sequence Error: " << exception.what() << std::endl;
            std::cout << message << std::endl;
        }
    }
}

#include "DNASequence.h"
#include "SequenceReader.h"
#include <fstream>

bool testDnaSequence();
bool testSequenceReader();

int main() {

    if (testDnaSequence()) {
        std::cout << "DNA Sequence Tests Passed" << std::endl;
    }

    if (testSequenceReader()) {
        std::cout << "Sequence Reader Tests Passed" << std::endl;
    }

    return 0;
}

bool testDnaSequence() {

    bool passed = true;

    DNASequence sequenceOne("Sequence 1", "ACGTACGT");

    if (sequenceOne.getSequence() != "ACGTACGT") {
        passed = false;
        std::cout << "Sequence initialization failed" << std::endl;
    }

    if (sequenceOne.getComplement() != "TGCATGCA") {
        passed = false;
        std::cout << "Reverse complement creation failed" << std::endl;
    }

    if (sequenceOne.getName() != "Sequence 1") {
        passed = false;
        std::cout << "Sequence name initialization failed" << std::endl;
    }

    if (sequenceOne.getLength() != 8) {
        passed = false;
        std::cout << "Sequence length calculation failed" << std::endl;
    }

    if (sequenceOne.getGCCount() != 4) {
        passed = false;
        std::cout << "GC count calculation failed" << std::endl;
    }

    if (sequenceOne.getGCPercent() != 0.5) {
        passed = false;
        std::cout << "GC percentage calculation failed" << std::endl;
    }

    return passed;
};

bool testSequenceReader() {
    bool passed = true;

    std::ifstream testFile("../Data/small_test.fna");

    auto genomeOpt = SequenceReader::readFasta(testFile);
    DNASequence genome = *genomeOpt;

    testFile.close();

    if (genome.getSequence() != "ACGTATTGAC") {
        passed = false;
    }

    if (genome.getName() != "TestSequence") {}

    return passed;
};

cmake_minimum_required(VERSION 3.8)
project(quigleytan_genome_assembler)

set(CMAKE_CXX_STANDARD 17)

add_executable(MainGenomeAssembly
        main.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/EulerianPath.h
        GenomeAssembly/EulerianPath.cpp
        KmerTypes.h
)

add_executable(EulerianTraversal
        EulerianTraversal.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/EulerianPath.h
        GenomeAssembly/EulerianPath.cpp
        KmerTypes.h
)

add_executable(ContigTraversal
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/EulerianPath.h
        GenomeAssembly/EulerianPath.cpp
        KmerTypes.h
)

add_executable(InitializationTests
        DataInitialization/InitializationTests.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        KmerTypes.h
)

add_executable(ProcessingTests
        DataProcessing/KmerEncoding.cpp
        DataProcessing/ProcessingTests.cpp
        DataInitialization/DNASequence.cpp
        KmerTypes.h
)

add_executable(AssemblyTests
        GenomeAssembly/AssemblyTests.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        DataProcessing/KmerEncoding.cpp
        DataInitialization/DNASequence.cpp
        GenomeAssembly/EulerianPath.h
        GenomeAssembly/EulerianPath.cpp
        KmerTypes.h
)

set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR})

target_include_directories(MainGenomeAssembly PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(EulerianTraversal PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(InitializationTests PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(ProcessingTests PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(AssemblyTests PRIVATE ${PROJECT_INCLUDE_DIR})
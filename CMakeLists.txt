cmake_minimum_required(VERSION 3.8)
project(M30EP-tequigle)

set(CMAKE_CXX_STANDARD 17)

add_executable(M3OEP-tequigle
        main.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/EulerianPath.h
        GenomeAssembly/EulerianPath.cpp
)

add_executable(EulerianTraversal
        CleanSequenceTraversal.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/EulerianPath.h
        GenomeAssembly/EulerianPath.cpp
)

add_executable(InitializationTests
        DataInitialization/InitializationTests.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
)

add_executable(ProcessingTests
        DataProcessing/KmerEncoding.cpp
        DataProcessing/ProcessingTests.cpp
        DataInitialization/DNASequence.cpp
)

add_executable(AssemblyTests
        GenomeAssembly/AssemblyTests.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        DataProcessing/KmerEncoding.cpp
        DataInitialization/DNASequence.cpp
        GenomeAssembly/EulerianPath.h
        GenomeAssembly/EulerianPath.cpp
)

set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR})

target_include_directories(M3OEP-tequigle PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(EulerianTraversal PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(InitializationTests PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(ProcessingTests PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(AssemblyTests PRIVATE ${PROJECT_INCLUDE_DIR})
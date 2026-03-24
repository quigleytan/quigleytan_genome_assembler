cmake_minimum_required(VERSION 3.8)
project(quigleytan_genome_assembler)

set(CMAKE_CXX_STANDARD 17)

add_executable(ScaffoldAssembly
        ConstructionPipelines/ScaffoldAssembly.cpp
        DataList.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/ContigTraversal.cpp
        GenomeAssembly/ContigScaffolder.cpp
        ConstructionPipelines/ScaffoldAssembly.cpp

)

add_executable(EulerianAssembly
        ConstructionPipelines/EulerianAssembly.cpp
        DataList.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/EulerianTraversal.cpp
)

add_executable(ContigAssembly
        ConstructionPipelines/ContigAssembly.cpp
        DataList.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/ContigTraversal.cpp
        GenomeAssembly/ContigScaffolder.cpp
        ConstructionPipelines/ScaffoldAssembly.cpp
)

add_executable(ContigTests
        GenomeAssembly/FastqContigTests.cpp
        DataInitialization/DNASequence.cpp
        DataInitialization/SequenceReader.cpp
        DataProcessing/KmerEncoding.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        GenomeAssembly/ContigTraversal.cpp
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
        DataList.cpp
        GenomeAssembly/ConstructionTests.cpp
        GenomeAssembly/DeBruijnGraph.cpp
        DataProcessing/KmerEncoding.cpp
        DataInitialization/DNASequence.cpp
        GenomeAssembly/EulerianTraversal.cpp
        GenomeAssembly/ContigTraversal.cpp
        GenomeAssembly/ContigScaffolder.cpp
)

set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR})

target_include_directories(ScaffoldAssembly PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(EulerianAssembly PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(ContigAssembly PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(InitializationTests PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(ProcessingTests PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(AssemblyTests PRIVATE ${PROJECT_INCLUDE_DIR})
target_include_directories(ContigTests PRIVATE ${PROJECT_INCLUDE_DIR})
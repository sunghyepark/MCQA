cmake_minimum_required(VERSION 3.20)

#User define
set(CMAKE_C_COMPILER "/tool/gcc/install/8.1.0/bin/gcc")
set(CMAKE_CXX_COMPILER "/tool/gcc/install/8.1.0/bin/g++")
set(SPARSEHASH_INCLUDE_DIR "/tool/cpplib/sparsehash/install/include")
set(Boost_INCLUDE_DIR "/tool/cpplib/boost/install/1.79.0/include")

project(MCQA)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-std=c++11 -static-libstdc++ -O3 -fopenmp -fpermissive -DNDEBUG -DIL_STD -pg -Wno-ctor-dtor-privacy")

set(TOOL
        src/qasm-tools/QASMtoken.cpp
        src/qasm-tools/QASMparser.cpp
        src/qasm-tools/QASMscanner.cpp)
set(SRC
        src/parser.cpp
        src/couplinggraph.cpp
        src/graph.cpp
        src/graphFunction.cpp
        src/latency.cpp
        src/layout.cpp
        src/circuit_optimization.cpp
        src/mapping.cpp
        src/mappingPreparation.cpp
        src/mappingFunction.cpp
        src/initial_mapping.cpp
        src/mymeasure.cpp
        src/outputwriter.cpp)

add_executable(MCQA src/main.cpp ${SRC} ${TOOL})

target_include_directories(MCQA PUBLIC ${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIR} ${SPARSEHASH_INCLUDE_DIR})

set_target_properties(MCQA PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/tool")



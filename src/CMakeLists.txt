cmake_minimum_required(VERSION 3.5)
project(FeedbackVertexSet)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_library(BOOST_PROGRAM_OPTIONS boost_program_options)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -pedantic -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_EXE_LINKER_FLAGS -static)
set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)
set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
set(CMAKE_SHARED_LIBRARY_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
set(BUILD_SHARED_LIBS OFF)

set(CMAKE_BUILD_TYPE RELEASE)

# set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(SOURCE_FILES
        Main.cpp
        solutions/HalfIntBranching.h
        solutions/HalfIntBranching.cpp
        solutions/SolutionVariants.h
        solutions/Solution.h
        lib/Graph.h
        lib/Graph.cpp
        lib/Translator.h
        lib/Translator.cpp
        lib/CubicSolver.cpp
        lib/CubicSolver.h
        lib/HalfIntRelax.h
        lib/HalfIntRelax.cpp
     )

set(SOURCE_VERIFIER
	OutputVerifier.cpp)



add_executable(FeedbackVertexSet ${SOURCE_FILES})
add_executable(Verifier ${SOURCE_VERIFIER})


SET_TARGET_PROPERTIES (FeedbackVertexSet PROPERTIES LINK_SEARCH_START_STATIC 1)
SET_TARGET_PROPERTIES (Verifier PROPERTIES LINK_SEARCH_START_STATIC 1)

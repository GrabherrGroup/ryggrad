cmake_minimum_required(VERSION 3.3)
project(ryggrad-base)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14 -O3 -w")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

# ryggrad binaries
set(SOURCE_FILES_TESTCOMPILE src/base/ErrorHandling.cc src/base/FileParser.cc src/base/StringUtil.cc test/TestCompile.cc src/util/mutil.cc src/visual/Eps.cc src/ml/NNet.cc src/ml/NNIO.cc) 

add_executable(TestCompile ${SOURCE_FILES_TESTCOMPILE})





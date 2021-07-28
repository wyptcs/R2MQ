cmake_minimum_required(VERSION 3.17)
project(NEWR2MQ)

set(CMAKE_CXX_STANDARD 11)

add_executable(NEWR2MQ main.cpp eliastest.cpp encodingtest.cpp RunEncoding.cpp RunEncoding.hpp RunEncodingTest.cpp Pawel.cpp huffman.cpp dag.cpp spacetest.cpp auxiliary.cpp lz.cpp newlz.cpp)
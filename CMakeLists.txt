cmake_minimum_required(VERSION 3.9)
project(c_long_read_aligner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O2")
#set(CMAKE_CXX_FLAGS "-mcmodel=medium")

add_executable(c_long_read_aligner main.cpp configs.h hash_utils.cpp hash_utils.h basket_min_hash.cpp basket_min_hash.h bio_utils/sequence.cpp bio_utils/sequence.h bio_utils/bio_reader.cpp bio_utils/bio_reader.h utils/time_profile.cpp utils/time_profile.h utils/logger.cpp utils/logger.h utils/vector_writer.cpp utils/vector_writer.h utils/heap.cpp utils/heap.h)
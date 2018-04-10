cmake_minimum_required(VERSION 3.9)
project(c_long_read_aligner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}  -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fpermissive")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcmodel=medium")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

add_executable(c_long_read_aligner
        main.cpp
        aryana_helper.h aryana_helper.cpp
        configs.h hash_utils.cpp hash_utils.h basket_min_hash.cpp basket_min_hash.h
        bio_utils/sequence.cpp bio_utils/sequence.h bio_utils/bio_reader.cpp bio_utils/bio_reader.h
        utils/time_profile.cpp utils/time_profile.h utils/logger.cpp utils/logger.h utils/vector_writer.cpp
        utils/vector_writer.h utils/heap.cpp utils/heap.h utils/multiproc.cpp utils/multiproc.h
#        aryana/aryana_main.c
        aryana/utils.h aryana/utils.c aryana/bwt.h aryana/bwt.c aryana/bntseq.h aryana/bntseq.c aryana/bwt_gen.c
        aryana/QSufSort.h aryana/QSufSort.c aryana/is.c aryana/bwtindex.c aryana/fa2bin.h aryana/fa2bin.c
        aryana/bwa2.h aryana/bwa2.c aryana/aligner.h aryana/aligner.c aryana/sam.h aryana/sam.c aryana/aryana_args.h
        aryana/hash.h aryana/hash.c aryana/bwtaln.h aryana/bwtaln.c aryana/bwaseqio.c aryana/bwt2.h aryana/bwt2.c
        aryana/bwtgap.h aryana/bwtgap.c aryana/bamlite.h aryana/bamlite.c aryana/bwa.h aryana/bwa.c
        aryana/ksw.h aryana/ksw.c aryana/smith.h aryana/smith.c)

target_link_libraries(c_long_read_aligner -lz -lm -lpthread)
#set(CMAKE_BUILD_TYPE Debug)

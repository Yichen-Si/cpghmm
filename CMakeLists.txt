cmake_minimum_required (VERSION 3.0)
project(cpghmm)


find_package (Eigen3 REQUIRED NO_MODULE HINTS ${CMAKE_SOURCE_DIR}/../eigen/build/)

include_directories(${CMAKE_SOURCE_DIR}/../htslib/include ${CMAKE_SOURCE_DIR}/../optim/header_only_version)

execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)
set(VERSION "1.0.0")

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -no-pie -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2 -DVERSION='\"${VERSION}\"' -DUSER='\"${USER}\"' -DDATE='\"${DATE}\"' -D__ZLIB_AVAILABLE__")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(cpghmm_SOURCES
    src/methy_infer.cpp
    src/cpghmm.cpp
    src/cthmm.cpp
    src/Error.cpp
    src/commands.cpp
    src/params.cpp
    src/tsv_reader.cpp
    src/utils.cpp)

set(cpghmm_HEADERS
    src/cramore.h
    src/brent_obj.h
    src/cthmm.h
    src/Error.h
    src/Constant.h
    src/pException.h
    src/commands.h
    src/params.h
    src/tsv_reader.h
    src/utils.h)

add_executable(cpghmm ${cpghmm_SOURCES} ${cpghmm_HEADERS})
target_link_libraries(cpghmm ${CMAKE_SOURCE_DIR}/../optim/build/lib/liboptim.so ${CMAKE_SOURCE_DIR}/../htslib/libhts.a Eigen3::Eigen dl pthread lzma z bz2 curl crypto)

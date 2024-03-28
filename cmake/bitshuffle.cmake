cmake_minimum_required(VERSION 3.10)
project(bitshuffle C)

set(CMAKE_C_STANDARD 11)

include_directories(lz4)
include_directories(src)

set(
        bitshuffle_headers
        lz4/lz4.h
        src/bitshuffle.h
        src/bitshuffle_core.h
        src/bitshuffle_internals.h
        src/bshuf_h5filter.h
        src/iochain.h
)

set(
        bitshuffle_sources
        lz4/lz4.c
        src/bitshuffle.c
        src/bitshuffle_core.c
        src/bshuf_h5filter.c
        src/bshuf_h5plugin.c
        src/iochain.c
)

add_library(
        bitshuffle
        ${bitshuffle_headers}
        ${bitshuffle_sources}
)

install(TARGETS bitshuffle DESTINATION lib)
install(FILES ${bitshuffle_headers} DESTINATION include/bitshuffle)

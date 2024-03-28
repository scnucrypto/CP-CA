CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(CPLEA)

include_directories(include)

aux_source_directory(src/WBMatrix CPLEA_SRC)
aux_source_directory(src CPLEA_SRC)

add_library(CPLEA_LIB ${CPLEA_SRC})

add_executable(CPLEA test/attack.c)

target_link_libraries(CPLEA CPLEA_LIB m)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")
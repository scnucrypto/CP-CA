CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(CPLEA)

include_directories(include)

aux_source_directory(src/m4ri CPLEA_SRC)

add_library(CPLEA_LIB ${CPLEA_SRC})

add_executable(CPLEA test/white_box_arx.c)

target_link_libraries(CPLEA CPLEA_LIB m)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")
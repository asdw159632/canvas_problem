# sub CMakeLists for tools used in jupyter notebook
project(Canvas_proble)
cmake_minimum_required(VERSION 3.10)

message(STATUS "ipytools tool will also be compiled")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/)


# 找到 chackjkdata 目录下所有 .C 文件
file(GLOB CHACKJKDATA_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/chackjkdata/*.C"
)

# 定义库
add_library(jkdata SHARED ${CHACKJKDATA_SOURCES})
#add_library(jkdataold ${libtype} "${CMAKE_CURRENT_SOURCE_DIR}/src/chackjkdata.C")
set_target_properties(jkdata PROPERTIES
LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
)
#set_target_properties(jkdataold PROPERTIES
#LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
#)

# 链接内部库

# 链接外部库
target_link_libraries(jkdata ${ROOT_LIBS})
#target_link_libraries(jkdataold ${ROOT_LIBS})

cmake_minimum_required(VERSION 3.26)
project(iditm3d
    VERSION 0.1.0 
    DESCRIPTION "An Inductive Dynamic Ionosphere-Thermosphere Model"
    LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

#specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#include source files from a directory
aux_source_directory(${CMAKE_SOURCE_DIR}/src SOURCE_FILES_DIR)
add_executable(iditm3d ${SOURCE_FILES_DIR})

#set directory where executable is produced
#set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY /home/jtu/mycode/iditm3d) doesn't work
set_target_properties(iditm3d PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

#header file directories
target_include_directories(iditm3d PUBLIC
    /usr/local/include
    PRIVATE ${CMAKE_SOURCE_DIR}/src
)

#specify additional libraries. full path required
find_library(convectLib
    NAMES libconvect.so
    HINTS ${CMAKE_SOURCE_DIR}/lib)

target_link_libraries(iditm3d PUBLIC
    /usr/local/lib/libpetsc.so
    /usr/local/lib/libf2cblas.a
    /usr/local/lib/libf2clapack.a
    /usr/local/lib/libmpi.so
    /usr/local/lib/libhdf5.so
    /usr/local/lib/libhdf5_cpp.so
    ${convectLib}
)

#set compile options
target_compile_options(iditm3d PRIVATE
    -Wall
    -Wshadow
    -Wextra
    -Wconversion
    -Wsign-conversion
    -Wuninitialized
    -Wunused
    #$<$<CONFIG:DEBUG>:-g3>
    #$<$<CONFIG:DEBUG>:-Og>
    #$<$<CONFIG:RELEASE>:-O3>  
)

# Modify compile flags to change optimization level from O3 to O2
string(REGEX REPLACE "([\\/\\-]O)3" "\\12"
  CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "([\\/\\-]O)3" "\\12"
  CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
string(REGEX REPLACE "([\\/\\-]O)3" "\\12"
  CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL}")
string(REGEX REPLACE "([\\/\\-]O)3" "\\12"
  CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

  #if want to install executable to a different directory
install(TARGETS iditm3d DESTINATION ${CMAKE_BINARY_DIR})
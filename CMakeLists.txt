cmake_minimum_required(VERSION 3.15)

# set the project name
project(cbns3d_cart_mpi)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find the MPI package
find_package(MPI REQUIRED)

if (MPI_FOUND)
    message(STATUS "MPI found: ${MPI_CXX_COMPILER}")
else()
    message(FATAL_ERROR "MPI not found.")
endif()

message(STATUS "MPI Compile flags: ${MPI_COMPILE_FLAGS}")
message(STATUS "MPI Include path: ${MPI_INCLUDE_PATH}")
message(STATUS "MPI Libraries: ${MPI_LIBRARIES}")

# Bring the headers into the project
include_directories(include)
 
file(GLOB SOURCES "src/*.cpp")

# add the executable
add_executable(cbns3d_cart_mpi ${SOURCES})

target_link_libraries(cbns3d_cart_mpi PRIVATE MPI::MPI_CXX)

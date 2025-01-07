/**
 * @brief High-Order 3D Compressible Navier-Stokes Solver
 *
 * This program implements a high-order finite difference solver for the 
 * 3D compressible Navier-Stokes equations on a Cartesian grid. 
 * It supports parallel execution using MPI, enabling efficient utilization 
 * of multiple cores or processors for large-scale computations.
 *
 * This file contains the main function, which serves as the application
 * entry point. 
 *
 * @author Bei Chen
 */

#include <mpi.h>
#include <iostream>
#include <filesystem>

#include "read_input.h"
#include "process_info.h"


int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // ensure the required command-line arguments is provided
  if (4 > argc) {
    if (0 == world_rank) {
      std::cout << "Usage: program <dimX> <dimY> <dimZ>\n";
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int dims[3] = {std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3])};

  // validate the number of processes
  if (dims[0] * dims[1] * dims[2] != world_size) {
    if (0 == world_rank) {
      std::cerr << "Error: Number of processes must equal dimX * dimY * dimZ.\n";
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  SimulationParams simParam;

  if (0 == world_rank) {
    // load simulation parameters from the input file
    if (0 != read_input(simParam)) MPI_Abort(MPI_COMM_WORLD, 1);

    auto create_dir = [](const char* path) {
      if (!std::filesystem::exists(path)) {
        if (!std::filesystem::create_directory(path)) {
	  std::cerr << "Failed to create directory: " << path << std::endl;
	  MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
    };

    char dirName[20];

    // create directories for storing the output files
    for (size_type i = 0; i <= simParam.num_output; i++) {
      std::snprintf(dirName, 20, "output%.3i", i);
      create_dir(dirName);
    }
  }

  // broadcast simulation parameters to all processes
  MPI_Bcast(&simParam, sizeof(SimulationParams), MPI_BYTE, 0, MPI_COMM_WORLD);

  // create a Cartesian topology
  int periods[3] = {1, 1, 1};
  MPI_Comm cart_comm;
  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cart_comm);

  ProcessInfo proc_info;
  
  MPI_Comm_rank(cart_comm, &proc_info.rank);
  MPI_Cart_coords(cart_comm, proc_info.rank, 3, proc_info.coords);

  proc_info.comm = cart_comm;
  
  // get neighbors along each dimension
  MPI_Cart_shift(cart_comm, 0, 1, &proc_info.left, &proc_info.right); // shift along X-dimension
  MPI_Cart_shift(cart_comm, 1, 1, &proc_info.down, &proc_info.up);    // shift along Y-dimension
  MPI_Cart_shift(cart_comm, 2, 1, &proc_info.back, &proc_info.front); // shift along Z-dimension
  
  MPI_Comm_free(&cart_comm);

  MPI_Finalize();
  
  return 0;

}

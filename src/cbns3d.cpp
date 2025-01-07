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

#include "read_input.h"


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

  // load simulation parameters from the input file
  if (0 == world_rank) {
    if (0 != read_input(simParam)) MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // broadcast simulation parameters to all processes
  MPI_Bcast(&simParam, sizeof(SimulationParams), MPI_BYTE, 0, MPI_COMM_WORLD);

  MPI_Finalize();
  
  return 0;

}

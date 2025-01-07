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
#include <sstream>
#include <chrono>

#include "read_input.h"
#include "process_info.h"
#include "Block3d.h"
#include "output_vtk_xml.h"


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
  
  // --------------------------------------------------------------------------
  // Partition the grid into multiple blocks based on the MPI topology

  size_type IM, JM, KM;

  if (dims[0] - 1 != proc_info.coords[0]) {
    IM = simParam.IM / dims[0];
  } else {
    IM = simParam.IM - (dims[0] - 1) * (simParam.IM / dims[0]);
  }

  if (dims[1] - 1 != proc_info.coords[1]) {
    JM = simParam.JM / dims[1];
  } else {
    JM = simParam.JM - (dims[1] - 1) * (simParam.JM / dims[1]);
  }

  if (dims[2] - 1 != proc_info.coords[2]) {
    KM = simParam.KM / dims[2];
  } else {
    KM = simParam.KM - (dims[2] - 1) * (simParam.KM / dims[2]);
  }

  Block3d block(&proc_info, &simParam, IM, JM, KM);

  block.gen_mesh(proc_info.coords[0] * (simParam.IM / dims[0]),
		 proc_info.coords[1] * (simParam.JM / dims[1]),
		 proc_info.coords[2] * (simParam.KM / dims[2]));

  // --------------------------------------------------------------------------
  // initialize the simulation.
  // It can either start with an initial flow field or resume from a checkpoint.
 
  std::ostringstream oss;
  oss << std::setw(3) << std::setfill('0') << proc_info.rank;
  simParam.checkpoint_fname += oss.str() + ".dat";
  
  if (simParam.is_restart) {
    block.initial_condition();
    block.calc_conservative(block.get_Q());
  } else {
    if (0 != block.read_bin(simParam.checkpoint_fname)) {
      std::cerr << "Unable to open the checkpoint file by process " << world_rank << "\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  // --------------------------------------------------------------------------
  // apply boundary conditions

  boundary_condition(&block, block.get_Q());

  block.calc_primitive(block.get_Q());

  block.init_Q_p();
  
  // --------------------------------------------------------------------------
  // output the initial flow field

  size_type cnt_out = 0;
  oss.str("");
  oss << std::setw(3) << std::setfill('0') << cnt_out;

  if (0 != output_vtk_xml(&block, proc_info.rank, world_size,
			  "output" + oss.str(), "output")) {
      MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  value_type dt_out = simParam.t_end / simParam.num_output;
  value_type t_out = ((int)(simParam.t_cur / dt_out) + 1) * dt_out;
  
  // --------------------------------------------------------------------------
  // perform simulation iterations and track execution time for
  // performance analysis
  
  auto start = std::chrono::high_resolution_clock::now();

  for (size_type n_t = 0; n_t < simParam.nstep_max; n_t++) {

    if (block.is_finished()) {
      if (0 == proc_info.rank) std::cout << "Finished!\n";
      break;
    }

    // ------------------------------------------------------------------------
    // compute the time step

    value_type dt_local = block.calc_dt();
    value_type dt = 0.0;
    
    MPI_Allreduce(&dt_local, &dt, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
    
    if (simParam.t_end - simParam.t_cur < dt) {
      dt = simParam.t_end - simParam.t_cur;
    }
    simParam.t_cur += dt;

    if (0 == proc_info.rank) {
      std::cout << n_t << ", " << dt << ", " << simParam.t_cur << std::endl;
    }

    // ------------------------------------------------------------------------
    
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds =
    std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

  double wall_time_local = elapsed_seconds.count();
  double wall_time = 0.0;
  double cpu_time = 0.0;
  
  MPI_Reduce(&wall_time_local, &wall_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
  MPI_Reduce(&wall_time_local, &cpu_time, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);

  if (0 == proc_info.rank) {
    std::cout << "  Wall time: " << wall_time << " seconds\n"
              << "  CPU time: " << cpu_time << " seconds\n";
  }

  // --------------------------------------------------------------------------

  MPI_Comm_free(&cart_comm);

  MPI_Finalize();
  
  return 0;

}

#include <iostream>
#include <sstream>
#include <iomanip>

#include "output_vtk_xml.h"


/**
 * @brief Outputs the flow field data in VTK XML format.
 *
 * This function writes the flow field data for a 3D block to a file in 
 * VTK XML format, which can be used for visualization and analysis in tools 
 * such as ParaView.
 *
 * @param blk: Pointer to the 3D block containing the flow field data.
 * @param rank: The rank of the current process.
 * @param num_blocks: Total number of blocks in the simulation.
 * @param fname: Name of the output file (excluding path).
 * @param dir_path: Path to the directory where the output file will be saved.
 * @return An integer indicating success (e.g., 0 for success, -1 for failure).
 */
int output_vtk_xml(Block3d* blk, const int rank,
		   const size_type num_blocks,
		   const std::string dir_path,
		   const std::string fname) {

  std::ostringstream oss;

  oss << std::setw(3) << std::setfill('0') << rank;

  std::string output_fname = dir_path + "/" + fname;

  // Each process writes its result to a .vtr file
  if (0 != blk->output_vtr(output_fname + oss.str() + ".vtr")) {
    std::cerr << "Unable to open the .vtr file by process " << rank << "\n";
    return -1;
  }
  
  if (0 == rank) {
    if (0 != blk->output_vtm(output_fname + ".vtm", fname, num_blocks)) {
      std::cerr << "Unable to open the .vtm file\n";
      return -1;
    }
  }

  return 0;
  
}

#ifndef _OUTPUT_VTK_XML_H_
#define _OUTPUT_VTK_XML_H_

#include <string>

#include "Block3d.h"


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
		   const std::string fname);

#endif /* _OUTPUT_VTK_XML_H_ */

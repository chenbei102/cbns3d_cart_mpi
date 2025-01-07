#ifndef _PROCESS_INFO_H_
#define _PROCESS_INFO_H_

#include <mpi.h>

/**
 * @brief Stores process information for use with MPI.
 *
 * This struct encapsulates relevant data for each process 
 * within an MPI communication environment.
 */
struct ProcessInfo {

  MPI_Comm comm;

  int rank;
  int coords[3];
    
  int left;
  int right;
  int down;
  int up;
  int back;
  int front;
  
};

#endif /* _PROCESS_INFO_H_ */

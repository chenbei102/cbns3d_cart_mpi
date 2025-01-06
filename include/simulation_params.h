#ifndef _SIMULATION_PARAMS_H_
#define _SIMULATION_PARAMS_H_

#include <string>

#include "constants.h"


/**
 * @brief Struct to hold simulation parameters for the 3D Navier-Stokes solver.
 *
 * This struct defines the configurable, geometric, and physical parameters 
 * required by the 3D compressible Navier-Stokes solver. 
 */
struct SimulationParams {
  
  std::string checkpoint_fname;

  size_type IM; 
  size_type JM; 
  size_type KM; 

  value_type CFL;

  value_type t_cur;
  value_type t_end;

  size_type nstep_max;
  size_type checkpoint_freq;
  size_type num_output;
  
  bool is_restart;

  value_type L;
  
  value_type Mach;
  value_type Re;
  value_type Mach2;
  value_type gM2;
  value_type Re_inv;

  value_type p_inf;

  value_type gamma;
  value_type Pr;
  value_type gam1;
  value_type gam1_inv;
  value_type Pr_inv;

  value_type C_T_inf;
  value_type C_dt_v;

};

#endif /* _SIMULATION_PARAMS_H_ */

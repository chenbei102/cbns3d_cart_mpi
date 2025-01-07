#include <iostream>
#include <fstream>
#include <cmath>

#include "read_input.h"


/**
 * Reads the input file and initializes the simulation parameters.
 *
 * @param params A reference to a SimulationParams instance where the
 *               input parameters will be stored.
 * @return An integer indicating the success or error code of the operation.
 */
int read_input(SimulationParams& params) {

  std::ifstream fh; 

  try {
    fh.open ("input.txt");
    
    value_type T_inf {0.0};

    fh >> params.checkpoint_fname;
    std::cout << "checkpoint_fname : " << params.checkpoint_fname << std::endl;

    fh >> params.CFL;
    std::cout << "CFL : " << params.CFL << std::endl;

    fh >> params.t_cur;
    std::cout << "t_cur : " << params.t_cur << std::endl;

    fh >> params.t_end;
    std::cout << "t_end : " << params.t_end << std::endl;

    fh >> params.Pr;
    std::cout << "Pr : " << params.Pr << std::endl;

    fh >> params.gamma;
    std::cout << "gamma : " << params.gamma << std::endl;

    fh >> params.Mach;
    std::cout << "Mach : " << params.Mach << std::endl;

    fh >> params.Re;
    std::cout << "Re : " << params.Re << std::endl;

    fh >> T_inf;
    std::cout << "T_inf : " << T_inf << std::endl;

    fh >> params.L;
    std::cout << "L : " << params.L << std::endl;

    fh >> params.IM;
    std::cout << "IM : " << params.IM << std::endl;

    fh >> params.JM;
    std::cout << "JM : " << params.JM << std::endl;

    fh >> params.KM;
    std::cout << "KM : " << params.KM << std::endl;

    fh >> params.nstep_max;
    std::cout << "nstep_max : " << params.nstep_max << std::endl;

    fh >> params.checkpoint_freq;
    std::cout << "checkpoint_freq : " << params.checkpoint_freq << std::endl;

    fh >> params.num_output;
    std::cout << "num_output : " << params.num_output << std::endl;
    
    fh >> params.is_restart;
    std::cout << "is_restart : " << params.is_restart << std::endl;

    // Precompute and store frequently used physical quantities to improve
    // computational efficiency.

    params.Mach2 = params.Mach * params.Mach;

    params.Re_inv = 1.0 / params.Re;
    params.Pr_inv = 1.0 / params.Pr;

    params.gam1 = params.gamma - 1.0;
    params.gM2 = params.gamma * params.Mach2;
    params.gam1_inv = 1.0 / params.gam1;

    params.p_inf = 1.0 / params.gM2;

    params.C_T_inf = 110.5 / T_inf;
    params.C_dt_v = std::max(4.0 / 3.0, params.gamma * params.Pr_inv) * params.Re_inv;

  }
  catch (const std::ifstream::failure& e) {
    std::cout << "Exception opening/reading file\n";
    std::cout << e.what() << std::endl;
    return -1;
  }

  fh.close();

  return 0;

}

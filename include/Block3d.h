#ifndef _BLOCK3D_H_
#define _BLOCK3D_H_

#include "process_info.h"
#include "simulation_params.h"


/**
 * @brief Represents a 3D computational block used for solving the Navier-Stokes equations.
 *
 * This class encapsulates the data structures and algorithms required to perform numerical
 * simulations of the Navier-Stokes equations in three-dimensional space. 
 */
class Block3d {

  ProcessInfo *proc_info;
  
  SimulationParams *sim_pars;
  
  // Constant to define the number of equations to solve
  static const size_type NEQ {constant::NEQ};

  // Number of ghost points (used in boundary conditions)
  static const size_type NG {constant::NG};
  static const size_type NGV {constant::NGV};

  // Number of points in the computational mesh 
  size_type IM; // Number of points in the i-direction
  size_type JM; // Number of points in the j-direction
  size_type KM; // Number of points in the k-direction

  size_type IM_G; 
  size_type JM_G; 
  size_type KM_G; 

  size_type IM_GV; 
  size_type JM_GV; 
  size_type KM_GV; 

  value_type *x;
  value_type *y;
  value_type *z;

  value_type dx;
  value_type dy;
  value_type dz;

  value_type *u;
  value_type *v;
  value_type *w;
  value_type *rho;
  value_type *T;
  value_type *p;

  value_type *mu;

  value_type *Q;
  value_type *Q_p;

  value_type *Ep;
  value_type *Fp;
  value_type *Gp;

  value_type *u_x;
  value_type *v_x;
  value_type *w_x;
  value_type *u_y;
  value_type *v_y;
  value_type *w_y;
  value_type *u_z;
  value_type *v_z;
  value_type *w_z;

  value_type *T_x;
  value_type *T_y;
  value_type *T_z;

  value_type *tau_xx;
  value_type *tau_yy;
  value_type *tau_zz;
  value_type *tau_xy;
  value_type *tau_xz;
  value_type *tau_yz;

  value_type *q_x;
  value_type *q_y;
  value_type *q_z;

  value_type *Ev;
  value_type *Fv;
  value_type *Gv;

  value_type* diff_flux_v;

  value_type* sbuf_left;
  value_type* sbuf_right;
  value_type* sbuf_down;
  value_type* sbuf_up;
  value_type* sbuf_back;
  value_type* sbuf_front;

  value_type* rbuf_left;
  value_type* rbuf_right;
  value_type* rbuf_down;
  value_type* rbuf_up;
  value_type* rbuf_back;
  value_type* rbuf_front;

public:
  Block3d() = delete;
  Block3d(ProcessInfo* proc_info,
	  SimulationParams* sim_pars,
	  const size_type num_x,
	  const size_type num_y,
	  const size_type num_z);

  ~Block3d();
  
  void allocate_mem();
  void free_mem();

  void gen_mesh(const size_type i_begin,
		const size_type j_begin,
		const size_type k_begin);

  void initial_condition();

  void calc_conservative(value_type* Q);

  value_type* get_Q();
  value_type* get_Q_p();

  int read_bin(const std::string fname);

  // --------------------------------------------------------------------------
  // Utility functions for converting multi-dimensional indices to
  // one-dimensional indices.

  inline size_type get_idx(index_type i, index_type j, index_type k) {
    return (NG+i) + IM_G * ((NG+j) + JM_G * (NG+k));
  }

  inline size_type get_idx_ux(index_type i, index_type j, index_type k) {
    return (NGV+i) + IM_GV * ((NGV+j) + JM_GV * (NGV+k));
  }

  inline size_type get_idx_Q(size_type n_eq, index_type i, index_type j, index_type k) {
    return n_eq + NEQ * ((NG+i) + IM_G * ((NG+j) + JM_G * (NG+k)));
  }

  inline size_type get_idx_Ep(size_type n_eq, size_type i, size_type j, size_type k) {
    return n_eq + NEQ * (i + (IM + 1) * (j + JM * k));
  }

  inline size_type get_idx_Fp(size_type n_eq, size_type i, size_type j, size_type k) {
    return n_eq + NEQ * (i + IM * (j + (JM + 1) * k));
  }

  inline size_type get_idx_Gp(size_type n_eq, size_type i, size_type j, size_type k) {
    return n_eq + NEQ * (i + IM * (j + JM * k));
  }

  inline size_type get_idx_Ev(size_type n_eq, index_type i, index_type j, index_type k) {
    return n_eq + (NEQ - 1) * ((NGV+i) + IM_GV * ((NGV+j) + JM_GV * (NGV+k)));
  }

  inline size_type get_idx_dfv(size_type n_eq, size_type i, size_type j, size_type k) {
    return n_eq + (NEQ - 1) * (i + IM * (j + JM * k));
  }
  
};

#endif /* _BLOCK3D_H_ */

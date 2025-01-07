#include <iostream>
#include <fstream>
#include <cmath>

#include "Block3d.h"


Block3d::Block3d(ProcessInfo* proc_info,
		 SimulationParams* sim_pars,
		 const size_type num_x,
		 const size_type num_y,
		 const size_type num_z)
  : proc_info(proc_info)
  , sim_pars(sim_pars)
  , IM(num_x)
  , JM(num_y)
  , KM(num_z)
{

  IM_G = IM + 2*NG;
  JM_G = JM + 2*NG;
  KM_G = KM + 2*NG;

  IM_GV = IM + 2*NGV;
  JM_GV = JM + 2*NGV;
  KM_GV = KM + 2*NGV;

  allocate_mem();

}

Block3d::~Block3d() {

  free_mem();
  
}

void Block3d::allocate_mem() {

  // allocate memory for arrays used in the simulation
  
  x = new value_type[IM + 1];
  y = new value_type[JM + 1];
  z = new value_type[KM + 1];

  size_type array_size = IM_G * JM_G * KM_G;

  u = new value_type[array_size];
  v = new value_type[array_size];
  w = new value_type[array_size];
  rho = new value_type[array_size];
  T = new value_type[array_size];
  p = new value_type[array_size];

  mu = new value_type[array_size];

  array_size = IM_GV * JM_GV * KM_GV;
  
  u_x = new value_type[array_size];
  v_x = new value_type[array_size];
  w_x = new value_type[array_size];
  u_y = new value_type[array_size];
  v_y = new value_type[array_size];
  w_y = new value_type[array_size];
  u_z = new value_type[array_size];
  v_z = new value_type[array_size];
  w_z = new value_type[array_size];

  T_x = new value_type[array_size];
  T_y = new value_type[array_size];
  T_z = new value_type[array_size];

  tau_xx = new value_type[array_size];
  tau_yy = new value_type[array_size];
  tau_zz = new value_type[array_size];
  tau_xy = new value_type[array_size];
  tau_xz = new value_type[array_size];
  tau_yz = new value_type[array_size];

  q_x = new value_type[array_size];
  q_y = new value_type[array_size];
  q_z = new value_type[array_size];
  
  array_size = NEQ * IM_G * JM_G * KM_G;

  Q = new value_type[array_size];
  Q_p = new value_type[array_size];

  array_size = (NEQ - 1) * IM_GV * JM_GV * KM_GV;

  Ev = new value_type[array_size];
  Fv = new value_type[array_size];
  Gv = new value_type[array_size];

  array_size = (NEQ - 1) * IM * JM * KM;
  
  diff_flux_v = new value_type[array_size];

  array_size = NEQ * (IM + 1) * JM * KM;

  Ep = new value_type[array_size];

  array_size = NEQ * IM * (JM + 1) * KM;

  Fp = new value_type[array_size];

  array_size = NEQ * IM * JM * (KM + 1);

  Gp = new value_type[array_size];
  
  array_size = NEQ * NG * JM * KM;

  sbuf_left = new value_type[array_size];
  sbuf_right = new value_type[array_size];
  rbuf_left = new value_type[array_size];
  rbuf_right = new value_type[array_size];

  array_size = NEQ * IM_G * NG * KM;

  sbuf_down = new value_type[array_size];
  sbuf_up = new value_type[array_size];
  rbuf_down = new value_type[array_size];
  rbuf_up = new value_type[array_size];

  array_size = NEQ * IM_G * JM_G * NG;

  sbuf_back = new value_type[array_size];
  sbuf_front = new value_type[array_size];
  rbuf_back = new value_type[array_size];
  rbuf_front = new value_type[array_size];

}

void Block3d::free_mem() {

  // Release memory allocated for dynamic arrays used in the simulation.

  delete[] x;
  delete[] y;
  delete[] z;

  delete[] u;
  delete[] v;
  delete[] w;
  delete[] rho;
  delete[] T;
  delete[] p;

  delete[] mu;

  delete[] Q;
  delete[] Q_p;

  delete[] Ep;
  delete[] Fp;
  delete[] Gp;

  delete[] u_x;
  delete[] v_x;
  delete[] w_x;
  delete[] u_y;
  delete[] v_y;
  delete[] w_y;
  delete[] u_z;
  delete[] v_z;
  delete[] w_z;

  delete[] T_x;
  delete[] T_y;
  delete[] T_z;

  delete[] tau_xx;
  delete[] tau_yy;
  delete[] tau_zz;
  delete[] tau_xy;
  delete[] tau_xz;
  delete[] tau_yz;

  delete[] q_x;
  delete[] q_y;
  delete[] q_z;

  delete[] Ev;
  delete[] Fv;
  delete[] Gv;

  delete[] diff_flux_v;

  delete[] sbuf_left;
  delete[] sbuf_right;
  delete[] sbuf_down;
  delete[] sbuf_up;
  delete[] sbuf_back;
  delete[] sbuf_front;

  delete[] rbuf_left;
  delete[] rbuf_right;
  delete[] rbuf_down;
  delete[] rbuf_up;
  delete[] rbuf_back;
  delete[] rbuf_front;
  
}

void Block3d::gen_mesh(const size_type i_begin,
		       const size_type j_begin,
		       const size_type k_begin) {

  // Generates the mesh data for the 3D block, starting at the given
  // indices (i_begin, j_begin, k_begin).

  value_type tmp = 2.0 * constant::PI;
  
  dx = tmp / sim_pars->IM;
  dy = tmp / sim_pars->JM;
  dz = tmp / sim_pars->KM;

  tmp /= 2.0;
  
  for (size_type i = 0; i < IM+1; i++) {
    x[i] = (i_begin + i) * dx - tmp;
  }

  for (size_type j = 0; j < JM+1; j++) {
    y[j] = (j_begin + j) * dy - tmp;
  }

  for (size_type k = 0; k < KM+1; k++) {
    z[k] = (k_begin + k) * dz - tmp;
  }

}

void Block3d::initial_condition() {

  // initialize the flow field 

  for (size_type k = 0; k < KM; k++) {

    value_type z_l = z[k] + 0.5 * dz;

    for (size_type j = 0; j < JM; j++) {

      value_type y_l = y[j] + 0.5 * dy;

      for (size_type i = 0; i < IM; i++) {

	value_type x_l = x[i] + 0.5 * dx;

	size_type idx1 = get_idx(i, j, k);

	u[idx1] = std::sin(x_l) * std::cos(y_l) * std::cos(z_l); 
	v[idx1] = -std::cos(x_l) * std::sin(y_l) * std::cos(z_l); 
	w[idx1] = 0.0;

	p[idx1] = sim_pars->p_inf + (std::cos(2.0*x_l) + std::cos(2.0*y_l)) *
	  (std::cos(2.0*z_l) + 2.0) / 16.0;

	rho[idx1] = sim_pars->gM2 * p[idx1];

	T[idx1] = 1.0;
	mu[idx1] = 1.0;

      }
    }
  }
  
}

void Block3d::calc_conservative(value_type *Q) {

  // Compute the conservative variables 

  for (size_type k = 0; k < KM; k++) {
    for (size_type j = 0; j < JM; j++) {
      for (size_type i = 0; i < IM; i++) {

	size_type idx1 = get_idx(i, j, k);
	size_type idx2 = get_idx_Q(0, i, j, k);

	value_type rr = rho[idx1];
	value_type uu = u[idx1];
	value_type vv = v[idx1];
	value_type ww = w[idx1];
	value_type pp = p[idx1];

	Q[idx2  ] = rr;
	Q[idx2+1] = rr * uu;
	Q[idx2+2] = rr * vv;
	Q[idx2+3] = rr * ww;
	Q[idx2+4] = pp / sim_pars->gam1 + 0.5 * rr * (uu * uu + vv * vv + ww * ww);

      }
    }
  }
    
}

value_type* Block3d::get_Q() {

  // Returns the private member 'Q'

  return Q;
  
}

value_type* Block3d::get_Q_p() {

  // Returns the private member 'Q_p'

  return Q_p;
  
}

int Block3d::read_bin(const std::string fname) {

  // reads flow field data from a binary file 

  std::ifstream fh;

  try {
    fh.open(fname, std::ios::in | std::ios::binary);

    value_type d_tmp;

    fh.read((char *)&d_tmp, sizeof(value_type));
    sim_pars->t_cur = d_tmp;

    size_type i_tmp;
    fh.read((char *)&i_tmp, sizeof(size_type));
    if (i_tmp != IM_G) throw std::runtime_error("The dimensions do not match.\n");
    fh.read((char *)&i_tmp, sizeof(size_type));
    if (i_tmp != JM_G) throw std::runtime_error("The dimensions do not match.\n");
    fh.read((char *)&i_tmp, sizeof(size_type));
    if (i_tmp != KM_G) throw std::runtime_error("The dimensions do not match.\n");

    size_type array_size = NEQ * IM_G * JM_G * KM_G;

    fh.read((char *)Q, array_size * sizeof(value_type));

  } catch (const std::ifstream::failure& e) {
    std::cout << e.what() << std::endl;
    return -1;
  } catch (const std::runtime_error& e) {
    std::cout << e.what() << std::endl;
    fh.close();
    return -1;
  }
  
  fh.close();

  return 0;

}

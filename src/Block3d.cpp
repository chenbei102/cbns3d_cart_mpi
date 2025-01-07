#include <iostream>
#include <fstream>
#include <cmath>

#include "Block3d.h"
#include "eigenvectors_roe.h"
#include "weno.h"
#include "lf_flux.h"


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

  // release memory allocated for dynamic arrays used in the simulation.

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

  // generate the mesh data for the 3D block, starting at the given
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

  // compute the conservative variables 

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

  // returns the private member 'Q'

  return Q;
  
}

value_type* Block3d::get_Q_p() {

  // returns the private member 'Q_p'

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

void Block3d::calc_primitive(const value_type* Q) {

  // compute the primitive variables 

  static const size_type array_size = IM_G * JM_G * KM_G;
  
  for(size_type idx1 = 0; idx1 < array_size; idx1++) {

    size_type idx2 = NEQ * idx1;

    value_type rr = Q[idx2  ];
    value_type uu = Q[idx2+1] / rr;
    value_type vv = Q[idx2+2] / rr;
    value_type ww = Q[idx2+3] / rr;
    value_type pp = sim_pars->gam1 * (Q[idx2+4] - 0.5 * rr * (uu * uu + vv * vv + ww * ww));

    rho[idx1] = rr;	
    u[idx1] = uu;
    v[idx1] = vv;
    w[idx1] = ww;
    p[idx1] = pp;

    T[idx1] = sim_pars->gM2 * pp / rr;

  }

}

void Block3d::init_Q_p() {
  
  // initialize the Q_p array by copying data from the Q array

  size_type array_size = NEQ * IM_G * JM_G * KM_G;
  for (size_type i = 0; i < array_size; i++) Q_p[i] = Q[i];
  
}

int Block3d::output_vtr(const std::string fname) {

  // output the flow field data of each block to a .vtr file.

  FILE *fh = std::fopen(fname.c_str(), "w");

  if (fh) {
    std::fprintf(fh, "<?xml version=\"1.0\"?>\n");
    std::fprintf(fh, "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
    std::fprintf(fh, "  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",
		 IM, JM, KM);
    std::fprintf(fh, "    <Piece Extent=\"0 %d 0 %d 0 %d\">\n",
		 IM, JM, KM);

    std::fprintf(fh, "    <CellData>\n");

    auto write_cell_data = [this, &fh](const char* name,
				       const value_type* arr) {
      std::fprintf(fh, "      <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", name);
      std::fprintf(fh, "        ");
      size_type cnt = 0;
      for(size_type k = 0; k < KM; k++) {
	for(size_type j = 0; j < JM; j++) {
	  for(size_type i = 0; i < IM; i++) {
	    std::fprintf(fh, "%e ", arr[get_idx(i, j, k)]);
	    if (4 == cnt++ % 5) std::fprintf(fh, "\n        ");
	  }
	}
      }
      if (0 != cnt % 5) std::fprintf(fh, "\n");
      std::fprintf(fh, "      </DataArray>\n");
    };

    write_cell_data("rho", rho);
    write_cell_data("u", u);
    write_cell_data("v", v);
    write_cell_data("w", w);
    write_cell_data("p", p);
    
    std::fprintf(fh, "    </CellData>\n");

    std::fprintf(fh, "    <Coordinates>\n");

    auto write_coord = [&fh](const char* name, const size_type array_size,
			     const value_type* arr) {
      std::fprintf(fh, "      <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", name);
      std::fprintf(fh, "        ");
      for (size_type i = 0; i < array_size; i++) {
	std::fprintf(fh, "%e ", arr[i]);
	if (4 == i % 5) std::fprintf(fh, "\n        ");
      }
      if (0 != array_size % 5) std::fprintf(fh, "\n");
      std::fprintf(fh, "      </DataArray>\n");
    };

    write_coord("x", IM+1, x);
    write_coord("y", JM+1, y);
    write_coord("z", KM+1, z);

    std::fprintf(fh, "    </Coordinates>\n");

    std::fprintf(fh, "    </Piece>\n");

    std::fprintf(fh, "  </RectilinearGrid>\n");
    std::fprintf(fh, "</VTKFile>\n");

    std::fclose(fh);
  } else {
    return -1;
  }

  return 0;
  
}

int Block3d::output_vtm(const std::string fname,
			const std::string vtr_fname,
			const size_type num_blocks) {

  // generate a multiblock (.vtm) file that links to the data files (.vtr)
  // of all blocks.

  FILE *fh = std::fopen(fname.c_str(), "w");

  if (fh) {
    std::fprintf(fh, "<?xml version=\"1.0\"?>\n");
    std::fprintf(fh, "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
    std::fprintf(fh, "  <vtkMultiBlockDataSet>\n");

    for (size_type i = 0; i < num_blocks; i++) {
      std::fprintf(fh, "    <DataSet index=\"%d\" file=\"%s%.3i.vtr\"/>\n",
		   i, vtr_fname.c_str(), i);
    }

    std::fprintf(fh, "  </vtkMultiBlockDataSet>\n");
    std::fprintf(fh, "</VTKFile>\n");

    std::fclose(fh);
  } else {
    return -1;
  }

  return 0;
  
}

bool Block3d::is_finished() {

  // check if the simulation has reached the predefined final time.

  bool flag = false;
  if (sim_pars->t_cur + 1.0e-5 > sim_pars->t_end) {
    flag = true;
  }

  return flag;

}

value_type Block3d::calc_dt() {

  // calculate the time step based on the Courant-Friedrichs-Lewy (CFL) condition

  value_type s_max = 0.0;

  for(size_type k = 0; k < KM; k++) {
    for(size_type j = 0; j < JM; j++) {
      for(size_type i = 0; i < IM; i++) {

	size_type idx1 = get_idx(i, j, k);

	value_type rr = rho[idx1];
	value_type uu = u[idx1];
	value_type vv = v[idx1];
	value_type ww = w[idx1];
	value_type pp = p[idx1];

	value_type c = std::sqrt(sim_pars->gamma * pp / rr);

	value_type factor_v = 2.0 * mu[idx1] * sim_pars->C_dt_v * std::sqrt(uu*uu + vv*vv + ww*ww) / c / rr;

	value_type t1 = std::abs(uu) + c / dx;
	t1 += factor_v / dx / dx;

	value_type t2 = std::abs(vv) + c / dy;
	t2 += factor_v / dy / dy;

	value_type t3 = std::abs(ww) + c / dz;
	t3 += factor_v / dz / dz;

	s_max = ((s_max < t1 + t2 + t3) ? (t1 + t2 + t3) : s_max); 
      }
    }
  }

  return sim_pars->CFL / s_max;

}

void Block3d::rec_riemann_x(const value_type* Q) {

  // compute the numerical flux at all interfaces normal to the x-axis using
  // an approximate Riemann solver.

  value_type R_l[NEQ][NEQ];
  value_type L_l[NEQ][NEQ];
  
  static const value_type gamma = sim_pars->gamma;
  static const value_type gam1 = sim_pars->gam1;

  for(size_type k = 0; k < KM; k++) {
    for(size_type j = 0; j < JM; j++) {
      for(index_type i = 0; i < IM+1; i++) {

	size_type idx1 = get_idx_Ep(0, i, j, k);
	size_type idx2 = get_idx(i-1, j, k);
	size_type idx3 = get_idx(i  , j, k);

	value_type rho_L = rho[idx2];
	value_type u_L = u[idx2];
	value_type v_L = v[idx2];
	value_type w_L = w[idx2];
	value_type p_L = p[idx2];

	value_type rho_R = rho[idx3];
	value_type u_R = u[idx3];
	value_type v_R = v[idx3];
	value_type w_R = w[idx3];
	value_type p_R = p[idx3];

	calc_eigenvectors(gamma,
			  rho_L, u_L, v_L, w_L, p_L,
			  rho_R, u_R, v_R, w_R, p_R,
			  1.0, 0.0, 0.0, R_l, L_l);

	rec_weno5(Q + get_idx_Q(0, i-3, j, k),
		  Q + get_idx_Q(0, i-2, j, k),
		  Q + get_idx_Q(0, i-1, j, k),
		  Q + get_idx_Q(0, i  , j, k),
		  Q + get_idx_Q(0, i+1, j, k),
		  R_l, L_l,
		  gam1,
		  rho_L, u_L, v_L, w_L, p_L);

	rec_weno5(Q + get_idx_Q(0, i+2, j, k),
		  Q + get_idx_Q(0, i+1, j, k),
		  Q + get_idx_Q(0, i  , j, k),
		  Q + get_idx_Q(0, i-1, j, k),
		  Q + get_idx_Q(0, i-2, j, k),
		  R_l, L_l,
		  gam1,
		  rho_R, u_R, v_R, w_R, p_R);

	lf_flux(gamma,
		rho_R, u_R, v_R, w_R, p_R,
		rho_L, u_L, v_L, w_L, p_L,
		1.0, 0.0, 0.0, &Ep[idx1]);

      }
    }
  }

}

void Block3d::rec_riemann_y(const value_type* Q) {

  // compute the numerical flux at all interfaces normal to the y-axis using
  // an approximate Riemann solver.

  value_type R_l[NEQ][NEQ];
  value_type L_l[NEQ][NEQ];
  
  static const value_type gamma = sim_pars->gamma;
  static const value_type gam1 = sim_pars->gam1;

  for(size_type k = 0; k < KM; k++) {
    for(index_type j = 0; j < JM+1; j++) {
      for(size_type i = 0; i < IM; i++) {

	size_type idx1 = get_idx_Fp(0, i, j, k);
	size_type idx2 = get_idx(i, j-1, k);
	size_type idx3 = get_idx(i, j  , k);

	value_type rho_L = rho[idx2];
	value_type u_L = u[idx2];
	value_type v_L = v[idx2];
	value_type w_L = w[idx2];
	value_type p_L = p[idx2];

	value_type rho_R = rho[idx3];
	value_type u_R = u[idx3];
	value_type v_R = v[idx3];
	value_type w_R = w[idx3];
	value_type p_R = p[idx3];

	calc_eigenvectors(gamma,
			  rho_L, u_L, v_L, w_L, p_L,
			  rho_R, u_R, v_R, w_R, p_R,
			  0.0, 1.0, 0.0, R_l, L_l);

	rec_weno5(Q + get_idx_Q(0, i, j-3, k),
		  Q + get_idx_Q(0, i, j-2, k),
		  Q + get_idx_Q(0, i, j-1, k),
		  Q + get_idx_Q(0, i, j  , k),
		  Q + get_idx_Q(0, i, j+1, k),
		  R_l, L_l,
		  gam1,
		  rho_L, u_L, v_L, w_L, p_L);

	rec_weno5(Q + get_idx_Q(0, i, j+2, k),
		  Q + get_idx_Q(0, i, j+1, k),
		  Q + get_idx_Q(0, i, j  , k),
		  Q + get_idx_Q(0, i, j-1, k),
		  Q + get_idx_Q(0, i, j-2, k),
		  R_l, L_l,
		  gam1,
		  rho_R, u_R, v_R, w_R, p_R);

	lf_flux(gamma,
		rho_R, u_R, v_R, w_R, p_R,
		rho_L, u_L, v_L, w_L, p_L,
		0.0, 1.0, 0.0, &Fp[idx1]);

      }
    }
  }

}

void Block3d::rec_riemann_z(const value_type* Q) {

  // compute the numerical flux at all interfaces normal to the z-axis using
  // an approximate Riemann solver.

  value_type R_l[NEQ][NEQ];
  value_type L_l[NEQ][NEQ];
  
  static const value_type gamma = sim_pars->gamma;
  static const value_type gam1 = sim_pars->gam1;

  for(index_type k = 0; k < KM+1; k++) {
    for(size_type j = 0; j < JM; j++) {
      for(size_type i = 0; i < IM; i++) {

	size_type idx1 = get_idx_Gp(0, i, j, k);
	size_type idx2 = get_idx(i, j, k-1);
	size_type idx3 = get_idx(i, j, k  );

	value_type rho_L = rho[idx2];
	value_type u_L = u[idx2];
	value_type v_L = v[idx2];
	value_type w_L = w[idx2];
	value_type p_L = p[idx2];

	value_type rho_R = rho[idx3];
	value_type u_R = u[idx3];
	value_type v_R = v[idx3];
	value_type w_R = w[idx3];
	value_type p_R = p[idx3];

	calc_eigenvectors(gamma,
			  rho_L, u_L, v_L, w_L, p_L,
			  rho_R, u_R, v_R, w_R, p_R,
			  0.0, 0.0, 1.0, R_l, L_l);

	rec_weno5(Q + get_idx_Q(0, i, j, k-3),
		  Q + get_idx_Q(0, i, j, k-2),
		  Q + get_idx_Q(0, i, j, k-1),
		  Q + get_idx_Q(0, i, j, k  ),
		  Q + get_idx_Q(0, i, j, k+1),
		  R_l, L_l,
		  gam1,
		  rho_L, u_L, v_L, w_L, p_L);

	rec_weno5(Q + get_idx_Q(0, i, j, k+2),
		  Q + get_idx_Q(0, i, j, k+1),
		  Q + get_idx_Q(0, i, j, k  ),
		  Q + get_idx_Q(0, i, j, k-1),
		  Q + get_idx_Q(0, i, j, k-2),
		  R_l, L_l,
		  gam1,
		  rho_R, u_R, v_R, w_R, p_R);

	lf_flux(gamma,
		rho_R, u_R, v_R, w_R, p_R,
		rho_L, u_L, v_L, w_L, p_L,
		0.0, 0.0, 1.0, &Gp[idx1]);

      }
    }
  }

}

void Block3d::calc_viscous_terms() {

  // calculate viscosity-related terms, including viscous stress and heat fluxes.

  calc_gradient(this, u, u_x, u_y, u_z);
  calc_gradient(this, v, v_x, v_y, v_z);
  calc_gradient(this, w, w_x, w_y, w_z);
  calc_gradient(this, T, T_x, T_y, T_z);

  static const index_type N_GV = NGV;
  static const index_type IMAX = IM + NGV;
  static const index_type JMAX = JM + NGV;
  static const index_type KMAX = KM + NGV;

  static const value_type C_T_inf = sim_pars->C_T_inf;
  static const value_type Re_inv = sim_pars->Re_inv;
  static const value_type kappa =
    -sim_pars->Pr_inv * sim_pars->gam1_inv / sim_pars->Mach2;

  for(index_type k = -N_GV; k < KMAX; k++) {
    for(index_type j = -N_GV; j < JMAX; j++) {
      for(index_type i = -N_GV; i < IMAX; i++) {

	size_type idx1 = get_idx(i, j, k);
	size_type idx2 = get_idx_ux(i, j, k);

	value_type T_l = T[idx1];
	value_type mu_l = (1.0 + C_T_inf) / (T_l + C_T_inf) * std::sqrt(T_l * T_l * T_l);

	mu[idx1] = mu_l;

	value_type u_x_l = u_x[idx2];
	value_type u_y_l = u_y[idx2];
	value_type u_z_l = u_z[idx2];
	value_type v_x_l = v_x[idx2];
	value_type v_y_l = v_y[idx2];
	value_type v_z_l = v_z[idx2];
	value_type w_x_l = w_x[idx2];
	value_type w_y_l = w_y[idx2];
	value_type w_z_l = w_z[idx2];

	mu_l *= Re_inv;

	tau_xx[idx2] = (2.0/3.0)*mu_l*(2.0*u_x_l - v_y_l - w_z_l);
	tau_yy[idx2] = (2.0/3.0)*mu_l*(2.0*v_y_l - u_x_l - w_z_l);
	tau_zz[idx2] = (2.0/3.0)*mu_l*(2.0*w_z_l - u_x_l - v_y_l);
	tau_xy[idx2] = mu_l*(u_y_l + v_x_l);
	tau_xz[idx2] = mu_l*(u_z_l + w_x_l);
	tau_yz[idx2] = mu_l*(v_z_l + w_y_l);

	mu_l *= kappa;

	q_x[idx2] = mu_l * T_x[idx2];
	q_y[idx2] = mu_l * T_y[idx2];
	q_z[idx2] = mu_l * T_z[idx2];

      }
    }
  }

}

void Block3d::calc_viscous_flux() {
  
  // compute viscous fluxes

  static const index_type N_GV = NGV;
  static const index_type IMAX = IM + NGV;
  static const index_type JMAX = JM + NGV;
  static const index_type KMAX = KM + NGV;

  for(index_type k = -N_GV; k < KMAX; k++) {
    for(index_type j = -N_GV; j < JMAX; j++) {
      for(index_type i = -N_GV; i < IMAX; i++) {

	size_type idx1 = get_idx(i, j, k);
	size_type idx2 = get_idx_ux(i, j, k);
	size_type idx3 = get_idx_Ev(0, i, j, k);

	value_type uu = u[idx1];
	value_type vv = v[idx1];
	value_type ww = w[idx1];

	value_type tau_xx_l = tau_xx[idx2];
	value_type tau_yy_l = tau_yy[idx2];
	value_type tau_zz_l = tau_zz[idx2];
	value_type tau_xy_l = tau_xy[idx2];
	value_type tau_xz_l = tau_xz[idx2];
	value_type tau_yz_l = tau_yz[idx2];

	// ---------------------------------------------------------------------
	// Viscous flux in x-direction 

	Ev[idx3  ] = tau_xx_l;
	Ev[idx3+1] = tau_xy_l;
	Ev[idx3+2] = tau_xz_l;
	Ev[idx3+3] = uu * tau_xx_l + vv * tau_xy_l + ww * tau_xz_l - q_x[idx2]; 

	// ---------------------------------------------------------------------
	// Viscous flux in y-direction

	Fv[idx3  ] = tau_xy_l;
	Fv[idx3+1] = tau_yy_l;
	Fv[idx3+2] = tau_yz_l;
	Fv[idx3+3] = uu * tau_xy_l + vv * tau_yy_l + ww * tau_yz_l - q_y[idx2];
	
	// ---------------------------------------------------------------------
	// Viscous flux in z-direction

	Gv[idx3  ] = tau_xz_l;
	Gv[idx3+1] = tau_yz_l;
	Gv[idx3+2] = tau_zz_l;
	Gv[idx3+3] = uu * tau_xz_l + vv * tau_yz_l + ww * tau_zz_l - q_z[idx2];

      }
    }
  }

}

void Block3d::calc_viscous_flux_contribution() {

  // calculate the total contribution of viscous flux derivatives in the x-, y-,
  // and z-directions.

  calc_viscous_terms();

  calc_viscous_flux();

  static const value_type coeff_c[2] = {1.0/12.0, -2.0/3.0};

  for (index_type k = 0; k < KM; k++) {
    for (index_type j = 0; j < JM; j++) {
      for (index_type i = 0; i < IM; i++) {

	value_type df = 0.0;    

	for(size_type i_eq = 0; i_eq < NEQ-1; i_eq++) {

	  size_type idx1 = get_idx_Ev(i_eq, i-2, j, k);
	  size_type idx2 = get_idx_Ev(i_eq, i-1, j, k);
	  size_type idx3 = get_idx_Ev(i_eq, i+1, j, k);
	  size_type idx4 = get_idx_Ev(i_eq, i+2, j, k);

	  df = (coeff_c[0]*Ev[idx1] + coeff_c[1]*Ev[idx2] - coeff_c[1]*Ev[idx3] - coeff_c[0]*Ev[idx4]) / dx;

	  idx1 = get_idx_Ev(i_eq, i, j-2, k);
	  idx2 = get_idx_Ev(i_eq, i, j-1, k);
	  idx3 = get_idx_Ev(i_eq, i, j+1, k);
	  idx4 = get_idx_Ev(i_eq, i, j+2, k);

	  df += (coeff_c[0]*Fv[idx1] + coeff_c[1]*Fv[idx2] - coeff_c[1]*Fv[idx3] - coeff_c[0]*Fv[idx4]) / dy;

	  idx1 = get_idx_Ev(i_eq, i, j, k-2);
	  idx2 = get_idx_Ev(i_eq, i, j, k-1);
	  idx3 = get_idx_Ev(i_eq, i, j, k+1);
	  idx4 = get_idx_Ev(i_eq, i, j, k+2);

	  df += (coeff_c[0]*Gv[idx1] + coeff_c[1]*Gv[idx2] - coeff_c[1]*Gv[idx3] - coeff_c[0]*Gv[idx4]) / dz;

	  diff_flux_v[get_idx_dfv(i_eq, i, j, k)] = df;

	}

      }
    }
  }

}

void Block3d::update_rk3(const value_type dt, const size_type stage) {

  // 3rd-order Total Variation Diminishing (TVD) Runge-Kutta method for time
  // integration

  for(size_type k = 0; k < KM; k++) {
    for(size_type j = 0; j < JM; j++) {
      for(size_type i = 0; i < IM; i++) {

	for(size_type i_eq = 0; i_eq < NEQ; i_eq++) {
	  size_type idx1 = get_idx_Ep(i_eq, i+1, j, k);
	  size_type idx2 = get_idx_Ep(i_eq, i  , j, k);

	  value_type df = (Ep[idx1] - Ep[idx2]) / dx; 

	  idx1 = get_idx_Fp(i_eq, i, j+1, k);
	  idx2 = get_idx_Fp(i_eq, i, j  , k);

	  df += (Fp[idx1] - Fp[idx2]) / dy;

	  idx1 = get_idx_Gp(i_eq, i, j, k+1);
	  idx2 = get_idx_Gp(i_eq, i, j, k  );

	  df += (Gp[idx1] - Gp[idx2]) / dz; 

	  if (i_eq > 0) df -= diff_flux_v[get_idx_dfv(i_eq-1, i, j, k)];

	  idx1 = get_idx_Q(i_eq, i, j, k);
	    
	  if (1 == stage) {
	    Q_p[idx1] = Q[idx1] - dt * df;
	  } else if (2 == stage) {
	    Q_p[idx1] = 0.75 * Q[idx1] + 0.25 * (Q_p[idx1] - dt * df);
	  } else if (3 == stage) {
	    Q[idx1] = (Q[idx1] + 2.0 * (Q_p[idx1] - dt * df)) / 3.0;
	  }

	}
      }
    }
  }

}

#include "weno.h"

/*
 * @brief 1D WENO Reconstruction (5th-Order Scheme)
 * 
 * This function implements the 5th-order Weighted Essentially Non-Oscillatory
 * (WENO) scheme for reconstructing the numerical flux at a cell interface. The
 * WENO scheme is designed to capture sharp gradients and discontinuities while
 * maintaining high-order accuracy in smooth regions.
 * 
 * The implementation is based on the paper:
 * "Efficient implementation of weighted ENO schemes." Journal of Computational
 * Physics, 126(1), 202-228 (1996).
 *
 * Parameters:
 *   Um2 (value_type): Value of the function at the third cell to the left of the interface.
 *   Um1 (value_type): Value of the function at the second cell to the left of the interface.
 *   Um0 (value_type): Value of the function at the first cell to the left of the interface.
 *   Up1 (value_type): Value of the function at the first cell to the right of the interface.
 *   Up2 (value_type): Value of the function at the second cell to the right of the interface.
 * 
 * Returns:
 *   value_type: Reconstructed value at the cell interface using the 5th-order WENO scheme.
 */
value_type weno5(const value_type Um2,
		 const value_type Um1,
		 const value_type U0,
		 const value_type Up1,
		 const value_type Up2) {

  static const value_type eps = 1.0e-6;
  
  value_type diff1 = Um2 - 2.0 * Um1 + U0;
  value_type diff2 = Um2 - 4.0 * Um1 + 3.0 * U0;
  value_type diff3 = Um1 - 2.0 * U0 + Up1;
  value_type diff4 = Um1 - Up1;
  value_type diff5 = U0 - 2.0 * Up1 + Up2;
  value_type diff6 = 3.0 * U0 - 4.0 * Up1 + Up2;

  value_type IS1 = (13.0/12.0) * diff1 * diff1 + 0.25 * diff2 * diff2 + eps;
  value_type IS2 = (13.0/12.0) * diff3 * diff3 + 0.25 * diff4 * diff4 + eps;
  value_type IS3 = (13.0/12.0) * diff5 * diff5 + 0.25 * diff6 * diff6 + eps;

  value_type w1 = 0.1 / (IS1 * IS1);   
  value_type w2 = 0.6 / (IS2 * IS2);   
  value_type w3 = 0.3 / (IS3 * IS3);   

  value_type sum = w1 + w2 + w3;
	  
  value_type f1 = (2.0 * Um2 - 7.0 * Um1 + 11.0 *U0) / 6.0;
  value_type f2 = (-Um1 + 5.0 * U0 + 2.0 *Up1) / 6.0;
  value_type f3 = (2.0 * U0 + 5.0 * Up1 - Up2) / 6.0;

  return (w1 * f1 + w2 * f2 + w3 * f3) / sum;
  
}

/*
 * @brief 5th-Order WENO Reconstruction for Vector Systems
 *
 * This function performs a characteristic-wise Weighted Essentially Non-Oscillatory (WENO) 
 * reconstruction of a vector system using a 5th-order scheme.
 * 
 * Steps:
 * 1. Transform the conservative variables in the WENO stencil into local characteristic fields.
 * 2. Apply 1D 5th-order WENO reconstruction to each component of the characteristic variables.
 * 3. Transform the reconstructed characteristic variables back into physical space.
 * 
 * @param f1, f2, f3, f4, f5: Pointers to the input conservative variable data in the stencil.
 * @param R_l: Right eigenvector matrix for local characteristic transformation.
 * @param L_l: Left eigenvector matrix for local characteristic transformation.
 * @param gam1: Specific heat ratio minus one (gamma - 1).
 * @param rho_L: Density of the reconstructed state (output).
 * @param u_L, v_L, w_L: Velocity components of the reconstructed state (output).
 * @param p_L: Pressure of the reconstructed state (output).
 */
void rec_weno5(const value_type* f1, const value_type* f2,
	       const value_type* f3, const value_type* f4,
	       const value_type* f5,
	       const value_type R_l[constant::NEQ][constant::NEQ],
	       const value_type L_l[constant::NEQ][constant::NEQ],
	       const value_type gam1,
	       value_type& rho_L,
	       value_type& u_L, value_type& v_L, value_type& w_L,
	       value_type& p_L) {
  
  value_type f_sten[constant::NEQ][5];
  value_type fc_sten[constant::NEQ][5];

  f_sten[0][0] = f1[0];
  f_sten[1][0] = f1[1];
  f_sten[2][0] = f1[2];
  f_sten[3][0] = f1[3];
  f_sten[4][0] = f1[4];

  f_sten[0][1] = f2[0];
  f_sten[1][1] = f2[1];
  f_sten[2][1] = f2[2];
  f_sten[3][1] = f2[3];
  f_sten[4][1] = f2[4];

  f_sten[0][2] = f3[0];
  f_sten[1][2] = f3[1];
  f_sten[2][2] = f3[2];
  f_sten[3][2] = f3[3];
  f_sten[4][2] = f3[4];

  f_sten[0][3] = f4[0];
  f_sten[1][3] = f4[1];
  f_sten[2][3] = f4[2];
  f_sten[3][3] = f4[3];
  f_sten[4][3] = f4[4];

  f_sten[0][4] = f5[0];
  f_sten[1][4] = f5[1];
  f_sten[2][4] = f5[2];
  f_sten[3][4] = f5[3];
  f_sten[4][4] = f5[4];

  for (size_type li = 0; li < constant::NEQ; li++) {
    for (size_type lj = 0; lj < 5; lj++) {

      value_type ss = 0.0;
      for (size_type lk = 0; lk < constant::NEQ; lk++) {
	ss += L_l[li][lk] * f_sten[lk][lj];
      }
      fc_sten[li][lj] = ss;

    }
  }
 
  f_sten[0][0] = weno5(fc_sten[0][0], fc_sten[0][1], fc_sten[0][2], fc_sten[0][3], fc_sten[0][4]);
  f_sten[1][0] = weno5(fc_sten[1][0], fc_sten[1][1], fc_sten[1][2], fc_sten[1][3], fc_sten[1][4]);
  f_sten[2][0] = weno5(fc_sten[2][0], fc_sten[2][1], fc_sten[2][2], fc_sten[2][3], fc_sten[2][4]);
  f_sten[3][0] = weno5(fc_sten[3][0], fc_sten[3][1], fc_sten[3][2], fc_sten[3][3], fc_sten[3][4]);
  f_sten[4][0] = weno5(fc_sten[4][0], fc_sten[4][1], fc_sten[4][2], fc_sten[4][3], fc_sten[4][4]);

  for (size_type li = 0; li < constant::NEQ; li++) {
    value_type ss = 0.0;
    for (size_type lj = 0; lj < constant::NEQ; lj++) {
      ss += R_l[li][lj] * f_sten[lj][0];
    }
    f_sten[li][1] = ss;
  }

  rho_L = f_sten[0][1];
  u_L = f_sten[1][1] / rho_L;
  v_L = f_sten[2][1] / rho_L;
  w_L = f_sten[3][1] / rho_L;
  p_L = gam1 * (f_sten[4][1] - 0.5 * rho_L * (u_L * u_L + v_L * v_L + w_L * w_L));
  
}

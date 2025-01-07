#ifndef _WENO_H_
#define _WENO_H_

#include "constants.h"


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
	       value_type& p_L);

#endif /* _WENO_H_ */

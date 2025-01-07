#ifndef _EIGENVECTORS_ROE_H_
#define _EIGENVECTORS_ROE_H_

#include "constants.h"


/*
 * @brief Calculate the eigenvectors of the flux Jacobian at the interface.
 *
 * This function computes the eigenvectors of the flux Jacobian matrix at the interface
 * using the Roe average of the left and right states. The inputs represent the 
 * thermodynamic properties (density, velocity components, and pressure) on the left 
 * and right sides of the interface, along with the normal vector components of the 
 * interface. The computed eigenvectors are stored in the matrices R and L, which 
 * are used for the diagonalization of the flux Jacobian.
 *
 * @param gamma: The specific heat ratio.
 * @param rho_L: The density on the left side of the interface.
 * @param u_L: The x-component of velocity on the left side.
 * @param v_L: The y-component of velocity on the left side.
 * @param w_L: The z-component of velocity on the left side.
 * @param p_L: The pressure on the left side of the interface.
 * @param rho_R: The density on the right side of the interface.
 * @param u_R The x-component of velocity on the right side.
 * @param v_R The y-component of velocity on the right side.
 * @param w_R The z-component of velocity on the right side.
 * @param p_R The pressure on the right side of the interface.
 * @param n_x The x-component of the unit normal vector at the interface.
 * @param n_y The y-component of the unit normal vector at the interface.
 * @param n_z The z-component of the unit normal vector at the interface.
 * @param R Output matrix for the right eigenvectors.
 * @param L Output matrix for the left eigenvectors.
 */
void calc_eigenvectors(const value_type gamma,
		       const value_type rho_L,
		       const value_type u_L,
		       const value_type v_L,
		       const value_type w_L,
		       const value_type p_L,
		       const value_type rho_R,
		       const value_type u_R,
		       const value_type v_R,
		       const value_type w_R,
		       const value_type p_R,
		       value_type n_x, value_type n_y, value_type n_z, 
		       value_type R[constant::NEQ][constant::NEQ],
		       value_type L[constant::NEQ][constant::NEQ]
		       );

#endif /* _EIGENVECTORS_ROE_H_ */

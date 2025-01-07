#ifndef _LF_FLUX_H_
#define _LF_FLUX_H_

#include "data_type.h"


/**
 * @brief Computes the Lax-Friedrichs flux as an approximate Riemann solver.
 *
 * This function calculates the numerical flux at interface using the Lax-Friedrichs
 * method, which is an approximate Riemann solver commonly used in numerical methods
 * for solving hyperbolic conservation laws. The flux is computed based on the state
 * variables on the left (L) and right (R) sides of the interface, as well as the
 * face-normal direction specified by the components (nx, ny, nz).
 *
 * @param gamma   Ratio of specific heats.
 * @param rho_R   Density on the right side of the interface.
 * @param u_R     Velocity in the x-direction on the right side.
 * @param v_R     Velocity in the y-direction on the right side.
 * @param w_R     Velocity in the z-direction on the right side.
 * @param p_R     Pressure on the right side of the interface.
 * @param rho_L   Density on the left side of the interface.
 * @param u_L     Velocity in the x-direction on the left side.
 * @param v_L     Velocity in the y-direction on the left side.
 * @param w_L     Velocity in the z-direction on the left side.
 * @param p_L     Pressure on the left side of the interface.
 * @param nx      x-component of the face-normal vector.
 * @param ny      y-component of the face-normal vector.
 * @param nz      z-component of the face-normal vector.
 * @param flux    Pointer to the array where the computed flux components will be stored.
 */
void lf_flux(const value_type gamma,
	     value_type rho_R,
	     value_type u_R, value_type v_R, value_type w_R,
	     value_type p_R,
	     value_type rho_L,
	     value_type u_L, value_type v_L, value_type w_L,
	     value_type p_L,
	     value_type nx, value_type ny, value_type nz,
	     value_type *flux);

#endif /* _LF_FLUX_H_ */

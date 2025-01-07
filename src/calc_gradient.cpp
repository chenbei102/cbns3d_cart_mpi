#include "Block3d.h"


/*
 * Compute the gradient of field variables within a 3D block.
 *
 * This function calculates the spatial derivatives of the input field variables 
 * using a fourth-order central differencing scheme.
 *
 * Parameters:
 * - blk: Pointer to the Block3d object representing the 3D computational block.
 * - f: Pointer to the input field variables array.
 * - f_x: Pointer to the output array for the x-direction gradient components.
 * - f_y: Pointer to the output array for the y-direction gradient components.
 * - f_z: Pointer to the output array for the z-direction gradient components.
 */
void calc_gradient(Block3d* blk,
		   const value_type* f,
		   value_type* f_x,
		   value_type* f_y,
		   value_type* f_z) {

  static const index_type NG = blk->NGV;

  static const index_type IM = blk->IM;
  static const index_type JM = blk->JM;
  static const index_type KM = blk->KM;

  for (index_type k = -NG; k < KM+NG; k++) {
    for (index_type j = -NG; j < JM+NG; j++) {
      for (index_type i = -NG; i < IM+NG; i++) {

	int idx0 = blk->get_idx_ux(i, j, k);
 
	int idx_1 = blk->get_idx(i - 1, j, k);
	int idx_2 = blk->get_idx(i - 2, j, k);
	int idx1  = blk->get_idx(i + 1, j, k);
	int idx2  = blk->get_idx(i + 2, j, k);

	f_x[idx0] = (1.0/12.0)*f[idx_2] - 2.0/3.0*f[idx_1] + (2.0/3.0)*f[idx1] - 1.0/12.0*f[idx2];


	idx_1 = blk->get_idx(i, j - 1, k);
	idx_2 = blk->get_idx(i, j - 2, k);
	idx1  = blk->get_idx(i, j + 1, k);
	idx2  = blk->get_idx(i, j + 2, k);

	f_y[idx0] = (1.0/12.0)*f[idx_2] - 2.0/3.0*f[idx_1] + (2.0/3.0)*f[idx1] - 1.0/12.0*f[idx2];


	idx_1 = blk->get_idx(i, j, k - 1);
	idx_2 = blk->get_idx(i, j, k - 2);
	idx1  = blk->get_idx(i, j, k + 1);
	idx2  = blk->get_idx(i, j, k + 2);

	f_z[idx0] = (1.0/12.0)*f[idx_2] - 2.0/3.0*f[idx_1] + (2.0/3.0)*f[idx1] - 1.0/12.0*f[idx2];

	f_x[idx0] /= blk->dx;
	f_y[idx0] /= blk->dy;
	f_z[idx0] /= blk->dz;

      }
    }
  }

}

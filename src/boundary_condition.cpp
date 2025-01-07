#include "Block3d.h"


void boundary_condition(Block3d* blk, value_type* Q) {

  // apply boundary conditions to the 3D block 

  static const size_type NEQ = blk->NEQ;

  static const index_type NG = blk->NG;

  static const index_type IM = blk->IM;
  static const index_type JM = blk->JM;
  static const index_type KM = blk->KM;

  static const index_type IM_G = blk->IM_G;
  static const index_type JM_G = blk->JM_G;

  // --------------------------------------------------------------------------
  // x-direction
  for (size_type k = 0; k < KM; k++) {
    for (size_type j = 0; j < JM; j++) {

      for (size_type i = 0; i < NG; i++) {

	size_type idx1 = NEQ * (i + NG * (j + JM * k));
	size_type idx2 = blk->get_idx_Q(0, i, j, k);

	blk->sbuf_left[idx1  ] = Q[idx2  ];
	blk->sbuf_left[idx1+1] = Q[idx2+1];
	blk->sbuf_left[idx1+2] = Q[idx2+2];
	blk->sbuf_left[idx1+3] = Q[idx2+3];
	blk->sbuf_left[idx1+4] = Q[idx2+4];

	idx2 = blk->get_idx_Q(0, IM-NG+i, j, k);

	blk->sbuf_right[idx1  ] = Q[idx2  ];
	blk->sbuf_right[idx1+1] = Q[idx2+1];
	blk->sbuf_right[idx1+2] = Q[idx2+2];
	blk->sbuf_right[idx1+3] = Q[idx2+3];
	blk->sbuf_right[idx1+4] = Q[idx2+4];

      }

    }
  }

  MPI_Request requests[4];

  static const int buf_size_x = NEQ * NG * JM * KM * sizeof(value_type);

  MPI_Isend(blk->sbuf_left, buf_size_x, MPI_BYTE, blk->proc_info->left, 0,
	    blk->proc_info->comm, &requests[0]);
  MPI_Irecv(blk->rbuf_right, buf_size_x, MPI_BYTE, blk->proc_info->right, 0,
	    blk->proc_info->comm, &requests[1]);

  MPI_Isend(blk->sbuf_right, buf_size_x, MPI_BYTE, blk->proc_info->right, 1,
	    blk->proc_info->comm, &requests[2]);
  MPI_Irecv(blk->rbuf_left, buf_size_x, MPI_BYTE, blk->proc_info->left, 1,
	    blk->proc_info->comm, &requests[3]);

  MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

  for (size_type k = 0; k < KM; k++) {
    for (size_type j = 0; j < JM; j++) {

      for (index_type i = 0; i < NG; i++) {

	size_type idx1 = NEQ * (i + NG * (j + JM * k));
	size_type idx2 = blk->get_idx_Q(0, i-NG, j, k);

	Q[idx2  ] = blk->rbuf_left[idx1  ];
	Q[idx2+1] = blk->rbuf_left[idx1+1];
	Q[idx2+2] = blk->rbuf_left[idx1+2];
	Q[idx2+3] = blk->rbuf_left[idx1+3];
	Q[idx2+4] = blk->rbuf_left[idx1+4];

	idx2 = blk->get_idx_Q(0, IM+i, j, k);

	Q[idx2  ] = blk->rbuf_right[idx1  ];
	Q[idx2+1] = blk->rbuf_right[idx1+1];
	Q[idx2+2] = blk->rbuf_right[idx1+2];
	Q[idx2+3] = blk->rbuf_right[idx1+3];
	Q[idx2+4] = blk->rbuf_right[idx1+4];

      }

    }
  }

  // --------------------------------------------------------------------------
  // y-direction
  for(size_type k = 0; k < KM; k++) {
    for(index_type i = -NG; i < IM+NG; i++) {

      for(size_type j = 0; j < NG; j++) {

	size_type idx1 = NEQ * ((NG+i) + IM_G * (j + NG * k));
	size_type idx2 = blk->get_idx_Q(0, i, j, k);

	blk->sbuf_down[idx1  ] = Q[idx2  ];
	blk->sbuf_down[idx1+1] = Q[idx2+1];
	blk->sbuf_down[idx1+2] = Q[idx2+2];
	blk->sbuf_down[idx1+3] = Q[idx2+3];
	blk->sbuf_down[idx1+4] = Q[idx2+4];

	idx2 = blk->get_idx_Q(0, i, JM-NG+j, k);

	blk->sbuf_up[idx1  ] = Q[idx2  ];
	blk->sbuf_up[idx1+1] = Q[idx2+1];
	blk->sbuf_up[idx1+2] = Q[idx2+2];
	blk->sbuf_up[idx1+3] = Q[idx2+3];
	blk->sbuf_up[idx1+4] = Q[idx2+4];

      }

    }
  }

  static const int buf_size_y = NEQ * IM_G * NG * KM * sizeof(value_type);

  MPI_Isend(blk->sbuf_down, buf_size_y, MPI_BYTE, blk->proc_info->down, 0,
	    blk->proc_info->comm, &requests[0]);
  MPI_Irecv(blk->rbuf_up, buf_size_y, MPI_BYTE, blk->proc_info->up, 0,
	    blk->proc_info->comm, &requests[1]);

  MPI_Isend(blk->sbuf_up, buf_size_y, MPI_BYTE, blk->proc_info->up, 1,
	    blk->proc_info->comm, &requests[2]);
  MPI_Irecv(blk->rbuf_down, buf_size_y, MPI_BYTE, blk->proc_info->down, 1,
	    blk->proc_info->comm, &requests[3]);

  MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

  for(size_type k = 0; k < KM; k++) {
    for(index_type i = -NG; i < IM+NG; i++) {

      for(index_type j = 0; j < NG; j++) {

	size_type idx1 = NEQ * ((NG+i) + IM_G * (j + NG * k));
	size_type idx2 = blk->get_idx_Q(0, i, j-NG, k);

	Q[idx2  ] = blk->rbuf_down[idx1  ];
	Q[idx2+1] = blk->rbuf_down[idx1+1];
	Q[idx2+2] = blk->rbuf_down[idx1+2];
	Q[idx2+3] = blk->rbuf_down[idx1+3];
	Q[idx2+4] = blk->rbuf_down[idx1+4];

	idx2 = blk->get_idx_Q(0, i, JM+j, k);

	Q[idx2  ] = blk->rbuf_up[idx1  ];
	Q[idx2+1] = blk->rbuf_up[idx1+1];
	Q[idx2+2] = blk->rbuf_up[idx1+2];
	Q[idx2+3] = blk->rbuf_up[idx1+3];
	Q[idx2+4] = blk->rbuf_up[idx1+4];

      }

    }
  }

  // --------------------------------------------------------------------------
  // z-direction
  for(index_type j = -NG; j < JM+NG; j++) {
    for(index_type i = -NG; i < IM+NG; i++) {

      for(size_type k = 0; k < NG; k++) {

	size_type idx1 = NEQ * ((NG+i) + IM_G * ((NG+j) + JM_G * k));
	size_type idx2 = blk->get_idx_Q(0, i, j, k);

	blk->sbuf_back[idx1  ] = Q[idx2  ];
	blk->sbuf_back[idx1+1] = Q[idx2+1];
	blk->sbuf_back[idx1+2] = Q[idx2+2];
	blk->sbuf_back[idx1+3] = Q[idx2+3];
	blk->sbuf_back[idx1+4] = Q[idx2+4];

	idx2 = blk->get_idx_Q(0, i, j, KM-NG+k);

	blk->sbuf_front[idx1  ] = Q[idx2  ];
	blk->sbuf_front[idx1+1] = Q[idx2+1];
	blk->sbuf_front[idx1+2] = Q[idx2+2];
	blk->sbuf_front[idx1+3] = Q[idx2+3];
	blk->sbuf_front[idx1+4] = Q[idx2+4];

      }

    }
  }

  static const int buf_size_z = NEQ * IM_G * JM_G * NG * sizeof(value_type);

  MPI_Isend(blk->sbuf_back, buf_size_z, MPI_BYTE, blk->proc_info->back, 0,
	    blk->proc_info->comm, &requests[0]);
  MPI_Irecv(blk->rbuf_front, buf_size_z, MPI_BYTE, blk->proc_info->front, 0,
	    blk->proc_info->comm, &requests[1]);

  MPI_Isend(blk->sbuf_front, buf_size_z, MPI_BYTE, blk->proc_info->front, 1,
	    blk->proc_info->comm, &requests[2]);
  MPI_Irecv(blk->rbuf_back, buf_size_z, MPI_BYTE, blk->proc_info->back, 1,
	    blk->proc_info->comm, &requests[3]);

  MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

  for(index_type j = -NG; j < JM+NG; j++) {
    for(index_type i = -NG; i < IM+NG; i++) {

      for(index_type k = 0; k < NG; k++) {

	size_type idx1 = NEQ * ((NG+i) + IM_G * ((NG+j) + JM_G * k));
	size_type idx2 = blk->get_idx_Q(0, i, j, k-NG);

	Q[idx2  ] = blk->rbuf_back[idx1  ];
	Q[idx2+1] = blk->rbuf_back[idx1+1];
	Q[idx2+2] = blk->rbuf_back[idx1+2];
	Q[idx2+3] = blk->rbuf_back[idx1+3];
	Q[idx2+4] = blk->rbuf_back[idx1+4];

	idx2 = blk->get_idx_Q(0, i, j, KM+k);

	Q[idx2  ] = blk->rbuf_front[idx1  ];
	Q[idx2+1] = blk->rbuf_front[idx1+1];
	Q[idx2+2] = blk->rbuf_front[idx1+2];
	Q[idx2+3] = blk->rbuf_front[idx1+3];
	Q[idx2+4] = blk->rbuf_front[idx1+4];

      }

    }
  }

}

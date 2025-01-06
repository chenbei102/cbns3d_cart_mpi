# cbns3d_cart_mpi

This repository hosts a **High-Order 3D Compressible Navier-Stokes Solver**. The solver utilizes high-order finite difference schemes to discretize and solve the 3D compressible Navier-Stokes equations on Cartesian grids. Leveraging MPI, the code enables parallel execution across multiple cores or processors, enhancing computational efficiency.

![taylor_green_vortex](fig/taylor_green_vortex.png)

<p align="center"><strong>Taylor-Green Vortex Test Case</strong></p>

## Features:
- **High-Order Accuracy:** Employs a 5th-order WENO scheme for inviscid fluxes, coupled with a 4th-order central finite difference scheme for viscous fluxes. Time integration utilizes the robust and efficient 3rd-order TVD Runge-Kutta method.
- **Parallel Performance:** Leverages MPI with 3D Cartesian topology for efficient parallel execution across multiple cores or processors, enhancing computational performance.
- **CMake Build System:** Utilizes CMake for streamlined build processes and seamless cross-platform compatibility.
- **VTK Output**: Simulation results are saved in the VTK XML format, enabling efficient parallel output and facilitating multi-block visualization.

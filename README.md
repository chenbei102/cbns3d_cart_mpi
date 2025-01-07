# cbns3d_cart_mpi

This repository hosts a **High-Order 3D Compressible Navier-Stokes Solver**. The solver utilizes high-order finite difference schemes to discretize and solve the 3D compressible Navier-Stokes equations on Cartesian grids. Leveraging MPI, the code enables parallel execution across multiple cores or processors, enhancing computational efficiency.

![taylor_green_vortex](fig/taylor_green_vortex.png)

<p align="center"><strong>Taylor-Green Vortex Test Case</strong></p>

## Features:
- **High-Order Accuracy:** Employs a 5th-order WENO scheme for inviscid fluxes, coupled with a 4th-order central finite difference scheme for viscous fluxes. Time integration utilizes the robust and efficient 3rd-order TVD Runge-Kutta method.
- **Parallel Performance:** Leverages MPI with 3D Cartesian topology for efficient parallel execution across multiple cores or processors, enhancing computational performance.
- **CMake Build System:** Utilizes CMake for streamlined build processes and seamless cross-platform compatibility.
- **VTK Output**: Simulation results are saved in the VTK XML format, enabling efficient parallel output and facilitating multi-block visualization.

## Getting Started

### Prerequisites

- **C++17** or higher
- **CMake 3.15** or higher
- **C++ compiler** (e.g., `g++`)
- **MPI library** (e.g., OpenMPI)
- **Python** (for running scripts)

### Clone the Repository

Clone the repository to your local machine:
```sh
git clone https://github.com/chenbei102/cbns3d_cart_mpi.git
```

### Build the Project

1. **Create and navigate to the build directory:**
    ```sh
    mkdir build
    cd build
    ```

2. **Generate build files with CMake:**
    ```sh
    cmake <path/to/source/code>
    ```

3. **Compile the project:**
    ```sh
    make
    ```

### Run the Solver

1. **Edit the configuration file:**
   Customize the `.ini` file to define simulation parameters. You can use `python/taylor_green.ini` as an example.

2. **Parse the configuration file:**
   Use the `parser_ini.py` script to parser the `.ini` file:
    ```sh
    python <path/to/source/code>/python/parser_ini.py <ini_file_name>
    ```

3. **Run the solver:**
   Execute the solver using `mpirun`. For example, 
    ```sh
    mpirun -np <numProc> ./cbns3d_cart_mpi <dimX> <dimY> <dimZ>
    ```
   Replace `<numProc>`, `<dimX>`, `<dimY>`, and `<dimZ>` with your specific values.
   
   **Note:** Ensure `<numProc>` equals `<dimX> * <dimY> * <dimZ>`.

### Visualization and Postprocessing

The solver outputs flow field data in the VTK XML format, organized as follows:  
- **MultiBlock file**: `output.vtm`  
- **Individual block files**: `output???.vtr` (one for each block)  

These files are stored in directories named `output???`, corresponding to different time series outputs.

To visualize and analyze the results:
1. Open an `output.vtm` file using **ParaView**.
2. Perform visualization and postprocessing as required.

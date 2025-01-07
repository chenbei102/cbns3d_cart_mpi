#ifndef _READ_INPUT_H_
#define _READ_INPUT_H_

#include "simulation_params.h"


/**
 * Reads the input file and initializes the simulation parameters.
 *
 * @param params A reference to a SimulationParams instance where the
 *               input parameters will be stored.
 * @return An integer indicating the success or error code of the operation.
 */
int read_input(SimulationParams& params);

#endif /* _READ_INPUT_H_ */

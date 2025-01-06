#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

// This header file defines global constants to:
//  - Avoid hardcoding values directly in the code.
//  - Facilitate future modifications and updates.

#include "data_type.h"


namespace constant {

  const size_type NEQ = 5;

  const size_type NGI = 3;
  const size_type NGV = 2;
  const size_type NG = (2*NGV > NGI ? 2*NGV : NGI);

  const value_type PI = 3.14159265358979323846;

}

#endif /* _CONSTANTS_H_ */

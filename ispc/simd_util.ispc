#include "real.h"

export void dummy (uniform int n)
{
  REAL x = 0.0, sum;

  foreach (i = 0 ... n)
  {
    sum += x;
  }
}

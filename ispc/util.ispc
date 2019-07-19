#include "real.h"
#include "alg.h"

export void translate (uniform REAL vector[], uniform REAL x[], uniform REAL y[], uniform REAL z[], uniform size_t size)
{
  foreach (i = 0 ... size)
  {
    #pragma ignore warning
    x[i] += vector[0];
    #pragma ignore warning
    y[i] += vector[1];
    #pragma ignore warning
    z[i] += vector[2];
  }
}

export void axis_rotate (uniform REAL point[], uniform REAL rotation[], uniform REAL x[], uniform REAL y[], uniform REAL z[], uniform size_t size)
{
  REAL v[3], w[3];

  foreach (i = 0 ... size)
  {
    #pragma ignore warning
    v[0] = x[i] - point[0];
    #pragma ignore warning
    v[1] = y[i] - point[1];
    #pragma ignore warning
    v[2] = z[i] - point[2];
    NVMUL (rotation, v, w);
    SCC (w, v);
    #pragma ignore warning
    x[i] += w[0];
    #pragma ignore warning
    y[i] += w[1];
    #pragma ignore warning
    z[i] += w[2];
  }
}

export void matrix_transform (uniform REAL transform[], uniform REAL x[], uniform REAL y[], uniform REAL z[], uniform size_t size)
{
  REAL v[3], w[3];

  foreach (i = 0 ... size)
  {
    #pragma ignore warning
    v[0] = x[i];
    #pragma ignore warning
    v[1] = y[i];
    #pragma ignore warning
    v[2] = z[i];
    NVMUL (transform, v, w);
    #pragma ignore warning
    x[i] = w[0];
    #pragma ignore warning
    y[i] = w[1];
    #pragma ignore warning
    z[i] = w[2];
  }
}

export void scale (uniform REAL point[], uniform REAL coef, uniform REAL x[], uniform REAL y[], uniform REAL z[], uniform size_t size)
{
  REAL v[3];

  foreach (i = 0 ... size)
  {
    #pragma ignore warning
    v[0] = x[i] - point[0];
    #pragma ignore warning
    v[1] = y[i] - point[1];
    #pragma ignore warning
    v[2] = z[i] - point[2];
    SCALE (v, coef);
    #pragma ignore warning
    x[i] = point[0] + v[0];
    #pragma ignore warning
    y[i] = point[1] + v[1];
    #pragma ignore warning
    z[i] = point[2] + v[2];
  }
}

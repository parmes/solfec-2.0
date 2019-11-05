/*
The MIT License (MIT)

Copyright (c) 2019 EDF Energy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* Contributors: Tomasz Koziara */

#include "real.h"
#include "alg.h"
#include "bulkmat.h"
#include "integration.h"

/* linear hexahedron shape functions */
inline static void hex1_shapes (REAL point[], REAL shapes[])
{
  shapes [0] = 0.125 * (1.0 - point [0]) * (1.0 - point [1]) * (1.0 - point [2]);
  shapes [1] = 0.125 * (1.0 + point [0]) * (1.0 - point [1]) * (1.0 - point [2]);
  shapes [2] = 0.125 * (1.0 + point [0]) * (1.0 + point [1]) * (1.0 - point [2]);
  shapes [3] = 0.125 * (1.0 - point [0]) * (1.0 + point [1]) * (1.0 - point [2]);
  shapes [4] = 0.125 * (1.0 - point [0]) * (1.0 - point [1]) * (1.0 + point [2]);
  shapes [5] = 0.125 * (1.0 + point [0]) * (1.0 - point [1]) * (1.0 + point [2]);
  shapes [6] = 0.125 * (1.0 + point [0]) * (1.0 + point [1]) * (1.0 + point [2]);
  shapes [7] = 0.125 * (1.0 - point [0]) * (1.0 + point [1]) * (1.0 + point [2]);
}

/* linear hexahedron shape derivatives */
inline static void hex1_derivs (REAL point[], REAL derivs[])
{
  derivs[0] = -0.125 * (1 - point[1]) * (1 - point[2]);
  derivs[1] = -0.125 * (1 - point[0]) * (1 - point[2]);
  derivs[2] = -0.125 * (1 - point[0]) * (1 - point[1]);

  derivs[3] =  0.125 * (1 - point[1]) * (1 - point[2]);
  derivs[4] = -0.125 * (1 + point[0]) * (1 - point[2]);
  derivs[5] = -0.125 * (1 + point[0]) * (1 - point[1]);

  derivs[6] =  0.125 * (1 + point[1]) * (1 - point[2]);
  derivs[7] =  0.125 * (1 + point[0]) * (1 - point[2]);
  derivs[8] = -0.125 * (1 + point[0]) * (1 + point[1]);

  derivs[9]  = -0.125 * (1 + point[1]) * (1 - point[2]);
  derivs[10] =  0.125 * (1 - point[0]) * (1 - point[2]);
  derivs[11] = -0.125 * (1 - point[0]) * (1 + point[1]);

  derivs[12] = -0.125 * (1 - point[1]) * (1 + point[2]);
  derivs[13] = -0.125 * (1 - point[0]) * (1 + point[2]);
  derivs[14] =  0.125 * (1 - point[0]) * (1 - point[1]);

  derivs[15] =  0.125 * (1 - point[1]) * (1 + point[2]);
  derivs[16] = -0.125 * (1 + point[0]) * (1 + point[2]);
  derivs[17] =  0.125 * (1 + point[0]) * (1 - point[1]);

  derivs[18] =  0.125 * (1 + point[1]) * (1 + point[2]);
  derivs[19] =  0.125 * (1 + point[0]) * (1 + point[2]);
  derivs[20] =  0.125 * (1 + point[0]) * (1 + point[1]);

  derivs[21] = -0.125 * (1 + point[1]) * (1 + point[2]);
  derivs[22] =  0.125 * (1 - point[0]) * (1 + point[2]);
  derivs[23] =  0.125 * (1 - point[0]) * (1 + point[1]);  
}

/* coordinates transformation determinant at local point */
inline static REAL hex1_det (REAL nodes[8][3], REAL point[], REAL F[])
{
  REAL derivs [24];
  int i, j, k;

  hex1_derivs (point, derivs);

  SET9 (F, 0.0);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 3; j ++)
      for (k = 0; k < 8; k ++)
       {
#pragma ignore warning
         F [3*j+i] += nodes[k][i] * derivs [3*k+j];
       }

  return DET (F);
}

/* element deformation determinant at local point */
inline static void hex1_gradient (REAL q[8][3], REAL point[], REAL F0[], REAL derivs[], REAL F[])
{
  REAL local_derivs[24], IF0[9], det;
  varying REAL * varying l;
  varying REAL * varying d;
  int i, j, k;

  hex1_derivs (point, local_derivs);

  INVERT (F0, IF0, det);
  assert (det > 0.0);
  for (k = 0, l = local_derivs, d = derivs; k < 8; k ++, l += 3, d += 3)
  { 
#pragma ignore warning
    TVMUL (IF0, l, d);
  }

  IDENTITY (F);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 3; j ++)
      for (k = 0; k < 8; k ++)
      {
#pragma ignore warning
        F [3*j+i] += q[k][i] * derivs [3*k+j];
      }
}

/* calculate element volume */
export void hex1_volumes (uniform int n, /* number of elements */
  uniform REAL * uniform x[8],
  uniform REAL * uniform y[8],
  uniform REAL * uniform z[8],
  uniform REAL volume[]) /* output volumes per element */
{
  foreach (i = 0 ... n)
  {
#pragma ignore warning
    REAL nodes [8][3] = {{x[0][i],y[0][i],z[0][i]}, {x[1][i],y[1][i],z[1][i]}, {x[2][i],y[2][i],z[2][i]}, {x[3][i],y[3][i],z[3][i]},
#pragma ignore warning
                         {x[4][i],y[4][i],z[4][i]}, {x[5][i],y[5][i],z[5][i]}, {x[6][i],y[6][i],z[6][i]}, {x[7][i],y[7][i],z[7][i]}};

    REAL vol = 0.0;

    INTEGRAL3D_BEGIN (8, 4)

      REAL F[9];
      REAL J = hex1_det (nodes, point, F);
      vol += J * weight;

    INTEGRAL3D_END()

#pragma ignore warning
    volume[i] = vol;
  }
}

/* calculate element lumped mass contribution */
export void hex1_lump_mass (uniform int n, /* number of elements */
  uniform REAL * uniform x[8],
  uniform REAL * uniform y[8],
  uniform REAL * uniform z[8],
  uniform REAL density[],
  uniform REAL * uniform lump_mass[8]) /* output nodal mass -- same along x, y, z -- per element */

{
  foreach (i = 0 ... n)
  {
#pragma ignore warning
    REAL nodes [8][3] = {{x[0][i],y[0][i],z[0][i]}, {x[1][i],y[1][i],z[1][i]}, {x[2][i],y[2][i],z[2][i]}, {x[3][i],y[3][i],z[3][i]},
#pragma ignore warning
                         {x[4][i],y[4][i],z[4][i]}, {x[5][i],y[5][i],z[5][i]}, {x[6][i],y[6][i],z[6][i]}, {x[7][i],y[7][i],z[7][i]}};

    REAL integral[8] = {0., 0., 0., 0., 0., 0., 0., 0.};

    REAL J, coef, rho, shapes [8], F[9];
    int j, k;

#pragma ignore warning
    rho = density[i];

    INTEGRAL3D_BEGIN (8, 4)

      hex1_shapes (point, shapes);
      J = hex1_det (nodes, point, F);
      coef = rho * J * weight;

      for (j = 0; j < 8; j ++)
      {
        for (k = 0; k < 8; k ++)
        {
#pragma ignore warning
          integral[j] += coef * shapes [j] * shapes [k];
        }
      }

    INTEGRAL3D_END() 

    for (uniform int l = 0; l < 8; l ++)
    {
#pragma ignore warning
      lump_mass[l][i] = integral[l];
    }
  }
}

/* compute element body force */
export void hex1_body_force (uniform int n, /* number of elements */
  uniform REAL * uniform x[8],
  uniform REAL * uniform y[8],
  uniform REAL * uniform z[8],
  uniform REAL density[],
  uniform REAL * uniform body_force[8]) /* output nodal body force -- to be scaled by gravity along x, y, z -- per element */
{
  foreach (i = 0 ... n)
  {
#pragma ignore warning
    REAL nodes [8][3] = {{x[0][i],y[0][i],z[0][i]}, {x[1][i],y[1][i],z[1][i]}, {x[2][i],y[2][i],z[2][i]}, {x[3][i],y[3][i],z[3][i]},
#pragma ignore warning
                         {x[4][i],y[4][i],z[4][i]}, {x[5][i],y[5][i],z[5][i]}, {x[6][i],y[6][i],z[6][i]}, {x[7][i],y[7][i],z[7][i]}};

    REAL integral[8] = {0., 0., 0., 0., 0., 0., 0., 0.};

    REAL J, coef, rho, shapes [8], F[9];
    int j, k;

#pragma ignore warning
    rho = density[i];

    INTEGRAL3D_BEGIN (8, 2)

      hex1_shapes (point, shapes);
      J = hex1_det (nodes, point, F);
      coef = rho * J * weight;

      for (j = 0; j < 8; j ++)
      {
#pragma ignore warning
        integral[j] += coef * shapes [j];
      }

    INTEGRAL3D_END() 

    for (uniform int l = 0; l < 8; l ++)
    {
#pragma ignore warning
      body_force[l][i] = integral[l];
    }
  }
}

/* compute element internal force */
export void hex1_internal_force (uniform int n, /* number of elements */
  uniform REAL * uniform x[8],
  uniform REAL * uniform y[8],
  uniform REAL * uniform z[8],
  uniform REAL * uniform dx[8],
  uniform REAL * uniform dy[8],
  uniform REAL * uniform dz[8],
  uniform REAL young[],
  uniform REAL poisson[],
  uniform REAL * uniform internal_force[24]) /* output internal force per element */
{
  foreach (i = 0 ... n)
  {
#pragma ignore warning
    REAL nodes [8][3] = {{x[0][i],y[0][i],z[0][i]}, {x[1][i],y[1][i],z[1][i]}, {x[2][i],y[2][i],z[2][i]}, {x[3][i],y[3][i],z[3][i]},
#pragma ignore warning
                         {x[4][i],y[4][i],z[4][i]}, {x[5][i],y[5][i],z[5][i]}, {x[6][i],y[6][i],z[6][i]}, {x[7][i],y[7][i],z[7][i]}};
#pragma ignore warning
    REAL q [8][3] = {{dx[0][i],dy[0][i],dz[0][i]}, {dx[1][i],dy[1][i],dz[1][i]}, {dx[2][i],dy[2][i],dz[2][i]}, {dx[3][i],dy[3][i],dz[3][i]},
#pragma ignore warning
                     {dx[4][i],dy[4][i],dz[4][i]}, {dx[5][i],dy[5][i],dz[5][i]}, {dx[6][i],dy[6][i],dz[6][i]}, {dx[7][i],dy[7][i],dz[7][i]}};

    REAL integral[24] = {0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0.};
    REAL J, coef, E, mu, derivs [24], F0[9], F[9], P[9];
    varying REAL * varying B;
    varying REAL * varying p;
    int j;

#pragma ignore warning
    E = young[i], mu = poisson[i];

    INTEGRAL3D_BEGIN (8, 2)

      J = hex1_det (nodes, point, F0);
      hex1_gradient (q, point, F0, derivs, F);
      coef = J * weight;

      BULK_MATERIAL_P (E, mu, F, coef, P);

      for (j = 0, B = derivs, p = integral; j < 8; j ++, B += 3, p += 3)
      { 
#pragma ignore warning
        NVADDMUL (p, P, B, p);
      }

    INTEGRAL3D_END() 

    for (uniform int l = 0; l < 24; l ++)
    {
#pragma ignore warning
      internal_force[l][i] = integral[l];
    }
  }
}

/* compute element stiffness matrix */
export void hex1_stiffness_matrix (uniform int n, /* number of elements */
  uniform REAL * uniform x[8],
  uniform REAL * uniform y[8],
  uniform REAL * uniform z[8],
  uniform REAL * uniform dx[8],
  uniform REAL * uniform dy[8],
  uniform REAL * uniform dz[8],
  uniform REAL young[],
  uniform REAL poisson[],
  uniform REAL * uniform stiffness_matrix[576]) /* output stiffness matrix per element */
{
  foreach (i = 0 ... n)
  {
#pragma ignore warning
    REAL nodes [8][3] = {{x[0][i],y[0][i],z[0][i]}, {x[1][i],y[1][i],z[1][i]}, {x[2][i],y[2][i],z[2][i]}, {x[3][i],y[3][i],z[3][i]},
#pragma ignore warning
                         {x[4][i],y[4][i],z[4][i]}, {x[5][i],y[5][i],z[5][i]}, {x[6][i],y[6][i],z[6][i]}, {x[7][i],y[7][i],z[7][i]}};
#pragma ignore warning
    REAL q [8][3] = {{dx[0][i],dy[0][i],dz[0][i]}, {dx[1][i],dy[1][i],dz[1][i]}, {dx[2][i],dy[2][i],dz[2][i]}, {dx[3][i],dy[3][i],dz[3][i]},
#pragma ignore warning
                     {dx[4][i],dy[4][i],dz[4][i]}, {dx[5][i],dy[5][i],dz[5][i]}, {dx[6][i],dy[6][i],dz[6][i]}, {dx[7][i],dy[7][i],dz[7][i]}};

    REAL integral[576] = {0.};
    REAL J, coef, E, mu, derivs [24], F0[9], F[9], K[81], KB[9];
    varying REAL * varying B;
    varying REAL * varying p;
    int j, k;

#pragma ignore warning
    E = young[i], mu = poisson[i];

    INTEGRAL3D_BEGIN (8, 2)

      J = hex1_det (nodes, point, F0);
      hex1_gradient (q, point, F0, derivs, F);
      coef = J * weight;

      BULK_MATERIAL_K (E, mu, F, coef, K);

      for (j = 0; j < 24; j ++) /* see doc/notes.lyx for details */
      {
	SET9 (KB, 0);
	for (k = 0; k < 3; k ++)
	{
          int l = j-((j/3)*3); /* j%3 since ISPC warns: "Modulus operator with varying types is very inefficient" */
	  p = &K [9*(l + (3*k))];
#pragma ignore warning
	  coef = derivs [3*(j/3)+k];
#pragma ignore warning
	  NNADDMUL (KB, coef, p, KB);
	}

	for (k = 0, B = derivs, p = &integral[24*j]; k < 8; k ++, B += 3, p += 3)
        { 
#pragma ignore warning
          NVADDMUL (p, KB, B, p);
        }
      }

    INTEGRAL3D_END() 

    for (uniform int l = 0; l < 576; l ++)
    {
#pragma ignore warning
      stiffness_matrix[l][i] = integral[l];
    }
  }
}

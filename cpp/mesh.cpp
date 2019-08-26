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

#include <metis/include/metis.h>
#include <algorithm>
#include <vector>
#include <tuple>
#include <map>
#include "real.h"
#include "alg.h"
#include "mesh.hpp"
#include "solfec.hpp"
#include "simplex.h"

static int tet [][4] =  /* 1-based indexing as the element lists start with the number of nodes */
  {{1,3,2,0},
   {1,2,4,0},
   {2,3,4,0},
   {3,1,4,0}}, pyr [][4] = 
  {{1,4,3,2},
   {1,2,5,0},
   {2,3,5,0},
   {3,4,5,0},
   {4,1,5,0}}, wed [][4] =
  {{1,3,2,0},
   {4,5,6,0},
   {1,2,5,4},
   {2,3,6,5},
   {3,1,4,6}}, hex [][4] =
  {{1,4,3,2},
   {1,2,6,5},
   {2,3,7,6},
   {3,4,8,7},
   {1,5,8,4},
   {5,6,7,8}};

/* create METIS-style all faces pointer and index vectors, and colors,
 * from elements definitions, a global color, and colors of selected faces;
 * returns [number of faces, face pointers, face indices, face colors] */
void mesh_create_metis_faces (const std::vector<uint64_t> &elements, uint64_t gcolor, const std::vector<uint64_t> &colors, /* input */
                              int64_t &nf, std::vector<int64_t> &fptr, std::vector<int64_t> &find, std::vector<uint64_t> &color) /* output */
{
  std::array<int(*)[4],9> ef = {NULL, NULL, NULL, NULL, tet, pyr, wed, NULL, hex};
  std::map<std::array<uint64_t,4>,int64_t> colored_faces;
  std::map<std::array<uint64_t,4>,int64_t> all_faces;
  int nef[] = {0, 0, 0, 0, 4, 5, 5, 0, 6};

  for(auto e = elements.begin(); e != elements.end(); e += e[0] + 2)
  {
    for (int i = 0; i < nef[e[0]]; i ++)
    {
      std::array<uint64_t,4> f = {e[ef[e[0]][i][0]]+1, e[ef[e[0]][i][1]]+1, e[ef[e[0]][i][2]]+1, e[ef[e[0]][i][3]] > 0 ? e[ef[e[0]][i][3]]+1 : 0};
      std::sort(f.begin(), f.end());
      all_faces[f] ++; /* see zero value-initialization: https://en.cppreference.com/w/cpp/language/zero_initialization */
    }
  }

  for(auto e = colors.begin(); e != colors.end(); e += e[0] + 2)
  {
    std::array<uint64_t,4> f = {e[1]+1, e[2]+1, e[3]+1, e[0] == 4 ? e[4]+1 : 0}; /* 1-based indexing */
    std::sort(f.begin(), f.end());
    colored_faces[f] = e[e[0]+1];
  }

  nf = 0;
  fptr.push_back(0);
  for (auto& [f, count] : all_faces)
  {
    if (count == 1) /* surface face */
    {
      if (f[0]) find.push_back(f[0]-1); /* if a quad */
      find.push_back(f[1]-1); /* 0-based indexing */
      find.push_back(f[2]-1);
      find.push_back(f[3]-1);
      fptr.push_back(find.size());
      color.push_back(colored_faces.count(f) ? colored_faces[f]: gcolor);
      nf ++;
    }
  }
}

/* sum up element charactersitics */
static inline void element_char_add (uint64_t matnum, int type, int64_t *enod,
  std::array<std::vector<REAL>,3> &nodes, REAL *me, REAL *sx, REAL *sy, REAL *sz, REAL *euler)
{
  REAL rho = solfec::materials[matnum].density;
  REAL J, *zero, *a, *b, *c;
  int (*ver) [4], nv[8], nf, i, j;

  switch (type)
  {
  case 4:
    nf = 4;
    ver = tet;
    nv[0] = nv[1] = nv[2] = nv[3] = 4;
  break;
  case 5:
    nf = 5;
    ver = pyr;
    nv[0] = 4; nv[1] = nv[2] = nv[3] = nv[4] = 3;
  break;
  case 6:
    nf = 5;
    ver = wed;
    nv[0] = nv[1] = 3; nv[2] = nv[3] = nv[4] = 4;
  break;
  case 8:
    nf = 6;
    ver = hex;
    nv[0] = nv[1] = nv[2] = nv[3] = 
    nv[4] = nv[5] = nv[6] = nv[7] = 4;
  break;
  }

  zero = &nodes[0][0];

  for (i = 0; i < nf; i++)
  {
    a = &nodes[enod[ver[i][0]-1]][0];

    for (j = 1; j < nv[i]-1; j ++)
    {
      b = &nodes[enod[ver[i][j]-1]][0];
      c = &nodes[enod[ver[i][j+1]-1]][0];

      J = rho * simplex_J (zero, a, b, c);
      *me += simplex_1 (J, zero, a, b, c);
      *sx += simplex_x (J, zero, a, b, c);
      *sy += simplex_y (J, zero, a, b, c);
      *sz += simplex_z (J, zero, a, b, c);
      euler [0] += simplex_xx (J, zero, a, b, c);
      euler [3] += simplex_xy (J, zero, a, b, c);
      euler [4] += simplex_yy (J, zero, a, b, c);
      euler [6] += simplex_xz (J, zero, a, b, c);
      euler [7] += simplex_yz (J, zero, a, b, c);
      euler [8] += simplex_zz (J, zero, a, b, c);
    }
  }
}

/* calculate mass characteristics: scalar mass, mass center, euler tensor */
void mesh_char (uint64_t bodnum, std::vector<uint64_t> &material, std::vector<int64_t> eptr, std::vector<int64_t> eind,
                REAL *mass, REAL *center, REAL *euler, REAL *inertia)
{
  struct mesh &mesh = solfec::meshes[bodnum];
  REAL me, sx, sy, sz, center0[3], euler0[9];

  me = sx = sy = sz = 0.0;
  SET9 (euler, 0.0);

  for (auto e = eptr.begin(); e != eptr.end()-1; e++)
  {
    auto i = *e, j = *(e+1);
    element_char_add (material[i], j-i, &eind[i], mesh.nodes, &me, &sx, &sy, &sz, euler0);
  }

  if (mass) mass[0] = me;

  center0 [0] = sx / me;
  center0 [1] = sy / me;
  center0 [2] = sz / me;

  if (center) COPY (center0, center);

  if (euler || inertia)
  {
    euler0 [0] -= (2*sx - center0[0]*me)*center0[0];
    euler0 [4] -= (2*sy - center0[1]*me)*center0[1];
    euler0 [8] -= (2*sz - center0[2]*me)*center0[2];
    euler0 [3] -= center0[0]*sy + center0[1]*sx - center0[0]*center0[1]*me;
    euler0 [6] -= center0[0]*sz + center0[2]*sx - center0[0]*center0[2]*me;
    euler0 [7] -= center0[1]*sz + center0[2]*sy - center0[1]*center0[2]*me;
    euler0 [1] = euler0[3];
    euler0 [2] = euler0[6];
    euler0 [5] = euler0[7];
  }

  if (euler) NNCOPY (euler0, euler);

  if (inertia)
  {
    /* convert Euler tensor to the inertia tensor */
    REAL trace = TRACE (euler0);
    IDENTITY (inertia);
    SCALE9 (inertia, trace);
    NNSUB (inertia, euler0, inertia); /* inertia = tr(euler)*one - euler */
  }
}

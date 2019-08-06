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
#include "mesh.hpp"

namespace mesh
{
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
};

/* create METIS-style all faces pointer and index vectors, and colors,
 * from elements definitions, a global color, and colors of selected faces;
 * returns [number of faces, face pointers, face indices, face colors] */
void mesh_create_metis_faces (const std::vector<uint64_t> &elements, uint64_t gcolor, const std::vector<uint64_t> &colors, /* input */
                              idx_t &nf, std::vector<idx_t> &fptr, std::vector<idx_t> &find, std::vector<uint64_t> &color) /* output */
{
  std::array<int(*)[4],9> ef = {NULL, NULL, NULL, NULL, mesh::tet, mesh::pyr, mesh::wed, NULL, mesh::hex};
  std::map<std::array<uint64_t,4>,int64_t> colored_faces;
  std::map<std::array<uint64_t,4>,int64_t> all_faces;
  int nef[] = {0, 0, 0, 0, 4, 5, 6, 0, 8};

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

  find.push_back(0);
  for (auto& [f, count] : all_faces)
  {
    if (count == 1) /* surface face */
    {
      if (f[0]) fptr.push_back(f[0]-1); /* if a quad */
      fptr.push_back(f[1]-1); /* 0-based indexing */
      fptr.push_back(f[2]-1);
      fptr.push_back(f[3]-1);
      find.push_back(fptr.size());
      color.push_back(colored_faces.count(f) ? colored_faces[f]: gcolor);
    }
  }
}

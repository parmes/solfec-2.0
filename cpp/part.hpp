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

#include <vector>
#include <tuple>
#include <map>

#ifndef __part__
#define __part__

/* mesh partitioning */
struct part
{
  int64_t nn, ne, nf, neparts, nfparts; /* number of: nodes, elements, faces, element parts, face parts */
  std::vector<int64_t> eptr, eind, epart, npart; /* element pointers, element indices, element partitioning, node partitioning */
  std::vector<int64_t> fptr, find, fpart; /* face pointers, face indices, face partitioning */
  std::vector<uint64_t> material; /* element materials */
  std::vector<uint64_t> color; /* face colors */
};

/* mesh mapping */
struct mapping
{
  std::vector<int> nrank; /* node MPI rank mapping */
  std::vector<uint64_t> nindex; /* node index within MPI rank */
  std::vector<int> erank; /* element MPI rank mapping */
  std::vector<int> frank; /* face MPI rank mapping */
};

/* partition input meshes and turn parts data */
std::map<uint64_t, part> partition_meshes(const std::set<uint64_t> &bodnum_subset);

/* map mesh partitioning to MPI ranks */
std::map<uint64_t, mapping> map_parts(const std::map<uint64_t, part> &parts);

/* return [maxnodes, maxeles, maxfaces] */
std::tuple<uint64_t, uint64_t, uint64_t> max_per_rank (const std::map<uint64_t, mapping> &maps);

#endif

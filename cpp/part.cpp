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

#include <taskflow/taskflow.hpp>
#include <metis/include/metis.h>
#include <mpi.h>
#include <map>
#include <set>
#include "part.hpp"
#include "real.h"
#include "mesh.hpp"
#include "solfec.hpp"
#include "compute.hpp"

/* Contributors: Tomasz Koziara */

/* partition input meshes and turn parts data */
std::map<uint64_t, part> partition_meshes(const std::set<uint64_t> &bodnum_subset, int elements_perpart, int faces_perpart)
{
  tf::Executor executor;
  tf::Taskflow taskflow;

  std::map<uint64_t, part> parts;

  for (auto& bodnum : bodnum_subset)
  {
    parts[bodnum]; /* populate map so that tasks only READ it */
  }

  taskflow.parallel_for(bodnum_subset.begin(), bodnum_subset.end(), [&] (const uint64_t bodnum)
  { 
    const struct mesh &mesh = solfec::meshes[bodnum];
    auto &part = parts[bodnum];
    std::vector<int64_t> temp;
    int64_t ncommon, objval;

    part.nn = mesh.nodes[0].size();
    part.ne = mesh.nhex+mesh.nwed+mesh.npyr+mesh.ntet;
    part.eptr.push_back(0);
    for(auto e = mesh.elements.begin(); e != mesh.elements.end();)
    {
      part.material.push_back(e[e[0]+1]);
      for (auto i = e+1; i != e+1+e[0]; i++)
      {
	part.eind.push_back(*i);
      }
      part.eptr.push_back(part.eind.size());
      e += e[0]+2;
    }

    /* partition elements into <= elements_perpart sized sets */
    part.epart.resize(part.ne);
    part.npart.resize(part.nn);
    ncommon = 3;
    part.neparts = 1+part.ne/elements_perpart;
    if (part.neparts > 1)
    {
      METIS_PartMeshDual(&part.ne, &part.nn, &part.eptr[0], &part.eind[0], NULL, NULL, &ncommon,
			 &part.neparts, NULL, NULL, &objval, &part.epart[0], &part.npart[0]);
    }
    else
    {
      std::fill(part.epart.begin(), part.epart.end(), 0);
      std::fill(part.npart.begin(), part.npart.end(), 0);
    }

    /* using mesh.elements, mesh.gcolor and mesh.colors create surface faces and colors */
    mesh_create_metis_faces (mesh.elements, mesh.gcolor, mesh.colors,
                             part.nf, part.fptr, part.find, part.color);

    /* partition faces into <= faces_perpart sized sets */
    part.fpart.resize(part.nf);
    temp.resize(part.nn);
    ncommon = 2;
    part.nfparts = 1+part.nf/faces_perpart;
    if (part.nfparts > 1)
    {
      METIS_PartMeshDual(&part.nf, &part.nn, &part.fptr[0], &part.find[0], NULL, NULL, &ncommon,
			 &part.nfparts, NULL, NULL, &objval, &part.fpart[0], &temp[0]);
    }
    else
    {
      std::fill(part.fpart.begin(), part.fpart.end(), 0);
    }
  });

  executor.run(taskflow).get();

  return parts;
}

/* map mesh partitioning to MPI ranks */
std::map<uint64_t, rankmap> map_parts(const std::map<uint64_t, part> &parts)
{
  int size;

  MPI_Comm_size (MPI_COMM_WORLD, &size);

  uint64_t etot = 0, ftot = 0;


  for (auto& [bodnum, part] : parts)
  {
    etot += part.neparts; /* total number of element parititons */
    ftot += part.nfparts; /* total number of face partitions */
  }

  uint64_t esplit = 1 + etot / size, fsplit = 1 + ftot / size;

  uint64_t esum = 0, fsum = 0;

  std::map<uint64_t, rankmap> output;

  for (auto& [bodnum, part] : parts)
  {
    struct rankmap &rankmap =  output[bodnum];

    for (auto& e : part.epart)
    {
      rankmap.erank.push_back ((esum + e)/esplit);
    }

    for (auto& f : part.fpart)
    {
      rankmap.frank.push_back ((fsum + f)/fsplit);
    }

    for (auto& n : part.npart)
    {
      int nrank = (esum + n)/esplit; /* |element partitions| = |node partitions| */
      rankmap.nrank.push_back (nrank);
    }

    esum += part.neparts;
    fsum += part.nfparts;
  }

  return output;
}

/* return [maxnodes, maxeles, maxfaces] */
std::tuple<uint64_t, uint64_t, uint64_t> max_per_rank (const std::map<uint64_t, rankmap> &maps)
{
  if (maps.empty()) return {0, 0, 0};

  std::map<int,uint64_t> nrank, erank, frank;

  for (auto& [bodnum, rankmap]: maps)
  {
    for (auto& r : rankmap.nrank)
    {
      nrank[r] ++;
    }

    for (auto& r : rankmap.erank)
    {
      erank[r] ++;
    }

    for (auto& r : rankmap.frank)
    {
      frank[r] ++;
    }
  }

  if (compute::debug_print)
  {
    std::cout << "DBG: In mesh partitioning:" << std::endl;

    for (auto& [r, n] : nrank)
    {
      std::cout << "DBG:   nodes on rank " << r << " = " << n << std::endl;
    }
    for (auto& [r, n] : erank)
    {
      std::cout << "DBG:   elements on rank " << r << " = " << n << std::endl;
    }
    for (auto& [r, n] : frank)
    {
      std::cout << "DBG:   faces on rank " << r << " = " << n << std::endl;
    }
  }

  using pair_type = decltype(nrank)::value_type;

  auto pair_cmp = [] (const pair_type &a, const pair_type &b) { return a.second < b.second; };

  auto n = std::max_element (std::begin(nrank), std::end(nrank), pair_cmp);
  auto e = std::max_element (std::begin(erank), std::end(erank), pair_cmp);
  auto f = std::max_element (std::begin(frank), std::end(frank), pair_cmp);

  return {n->second, e->second, f->second};
}

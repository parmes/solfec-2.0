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

#include <taskflow/taskflow.hpp>
#include <metis/include/metis.h>
#include <mpi.h>
#include <map>
#include <set>
#include "real.h"
#include "mesh.hpp"
#include "solfec.hpp"
#include "compute.hpp"

struct nodes_window /* nodes in mpi_nodes_window */
{
  REAL *refpos[3];
  REAL *curpos[3];
  REAL *velo[3];
  REAL *disp[3];
  size_t size;
};

#ifndef ELEMENTS_BUNCH
#define ELEMENTS_BUNCH 16
#endif

struct elements_bunch
{
  size_t node[2][ELEMENTS_BUNCH]; /* node[0] - MPI rank, node[1] - window index, of element nodes */

  REAL refpos[3][ELEMENTS_BUNCH]; /* local copy of nodes refpos */
  REAL curpos[3][ELEMENTS_BUNCH]; /* local copy of nodes curpos */
  REAL velo[3][ELEMENTS_BUNCH]; /* local copy of nodes velo */
  REAL disp[3][ELEMENTS_BUNCH]; /* local copy of nodes disp */

  short type[ELEMENTS_BUNCH]; /* 4, 5, 6, 8 => tetrahedron, pyramid, wedge, hexahedron */
  short lnod[8][ELEMENTS_BUNCH]; /* local node indices */

  size_t matnum[ELEMENTS_BUNCH]; /* material number */
  REAL density[ELEMENTS_BUNCH]; /* local material data */
  REAL young[ELEMENTS_BUNCH];
  REAL poisson[ELEMENTS_BUNCH];
  REAL viscosity[ELEMENTS_BUNCH];

  /* TODO: restrains and prescribes */

  size_t bodnum; /* body number */

  short size; /* bunch's size */
};

struct elements_window /* element bunches in mpi_elements_window */
{
  elements_bunch *bunch;
  size_t size;
};

#ifndef FACES_BUNCH
#define FACES_BUNCH 16
#endif

struct faces_bunch
{
  size_t node[2][FACES_BUNCH]; /* node[0] - MPI rank, node[1] - window index, of face nodes */

  REAL curpos[3][FACES_BUNCH]; /* local copy of nodes curpos */
  
  short type[FACES_BUNCH]; /* 3, 4 => triangle, quadrilateral */
  short lnod[4][FACES_BUNCH]; /* local node indices */
  short color[FACES_BUNCH];
  REAL normal[3][FACES_BUNCH];

  size_t elem[3][FACES_BUNCH]; /* elem[0] - MPI rank, elem[1] - window index, elem[2] - bunch index, of face element */

  size_t bodnum; /* body number */

  short size; /* bunch's size */
};

struct faces_window /* face bunches in mpi_faces_window */
{
  faces_bunch *bunch;
  size_t size;
};

/* compute gobal variables */
namespace compute
{
/* rank 0 --- */
bool partitioned = false; /* initially paritioned */

std::set<size_t> inserted_meshes;
std::set<size_t> deleted_meshes;
std::set<size_t> inserted_ellips;
std::set<size_t> deleted_ellips;
std::set<size_t> inserted_restrains;
std::set<size_t> deleted_restrains;
std::set<size_t> inserted_prescribes;
std::set<size_t> deleted_prescribes;
/* --- rank 0 */

/* all ranks --- */
std::map<size_t,material> materials; /* local MPI rank copy of all nominal materials */

MPI_Win mpi_nodes_window; /* storing nodes and dgreees of freedom data */
MPI_Win mpi_elements_window; /* storing elements data */
MPI_Win mpi_faces_window; /* storing faces data */
/* --- all ranks */
};

/* mesh partitioning */
struct part
{
  /* partition_meshes --- */
  idx_t nn, ne, nf, neparts, nfparts; /* number of: nodes, elements, faces, element parts, face parts */
  std::vector<idx_t> eptr, eind, epart, npart; /* element pointers, element indices, element partitioning, node partitioning */
  std::vector<idx_t> fptr, find, fpart; /* face pointers, face indices, face partitioning */
  std::vector<size_t> material; /* element materials */
  std::vector<size_t> color; /* face colors */
  /* --- partition_meshes */

  /* rank and bunch assignments --- */
  std::vector<int> erank, frank; /* element and face MPI rank assignment */
  std::vector<size_t> ebunch, fbunch; /* element and face data bunch assignment */
  /* --- rank and bunch assignments */
};

/* partition input meshes and turn parts data */
static std::map<size_t, part> partition_meshes(const std::set<size_t> &bodnum_subset)
{
  tf::Executor executor;
  tf::Taskflow taskflow;

  std::map<size_t, part> parts;

  for (auto& bodnum : bodnum_subset)
  {
    parts[bodnum]; /* populate map so that tasks only READ it */
  }

  taskflow.parallel_for(bodnum_subset.begin(), bodnum_subset.end(), [&] (const size_t bodnum)
  { 
    const struct mesh &mesh = solfec::meshes[bodnum];
    auto &part = parts[bodnum];
    std::vector<idx_t> temp;
    idx_t ncommon, objval;

    part.nn = mesh.nodes.size();
    part.ne = mesh.nhex+mesh.nwed+mesh.npyr+mesh.ntet;
    part.eptr.push_back(0);
    for(auto e = mesh.elements.begin(); e != mesh.elements.end();)
    {
      part.material.push_back(e[e[0]+1]);
      for (auto i = e+1; i != e+1+e[0]; i++)
      {
	part.eind.push_back(*i);
      }
      e += e[0]+2;
      part.eptr.push_back(e-mesh.elements.begin());
    }

    /* partition elements into <= ELEMENTS_BUNCH sized sets */
    part.epart.resize(part.ne);
    part.npart.resize(part.nn);
    ncommon = 3;
    part.neparts = 1+part.ne/ELEMENTS_BUNCH;
    METIS_PartMeshDual(&part.ne, &part.nn, &part.eptr[0], &part.eind[0], NULL, NULL, &ncommon,
                       &part.neparts, NULL, NULL, &objval, &part.epart[0], &part.npart[0]);

    /* using mesh.elements, mesh.gcolor and mesh.colors create surface faces and colors */
    mesh_create_metis_faces (mesh.elements, mesh.gcolor, mesh.colors,
                             part.nf, part.fptr, part.find, part.color);

    /* partition faces into <= FACES_BUNCH sized sets */
    part.fpart.resize(part.ne);
    temp.resize(part.nn);
    ncommon = 2;
    part.nfparts = 1+part.nf/ELEMENTS_BUNCH;
    METIS_PartMeshDual(&part.nf, &part.nn, &part.fptr[0], &part.find[0], NULL, NULL, &ncommon,
                       &part.nfparts, NULL, NULL, &objval, &part.fpart[0], &temp[0]);
  });

  executor.run(taskflow).get();

  return parts;
}

/* insert solfec::meshes[bodnum] into computation */
void compute_insert_mesh(size_t bodnum)
{
  compute::inserted_meshes.insert(bodnum);
}

/* delete mesh from computation */
void compute_delete_mesh(size_t bodnum)
{
  if (compute::inserted_meshes.count(bodnum))
    compute::inserted_meshes.erase(bodnum);
  else compute::deleted_meshes.insert(bodnum);
}

/* insert solfec::ellips[bodnum] into computation */
void compute_insert_ellip(size_t bodnum)
{
  compute::inserted_ellips.insert(bodnum);
}

/* delete ellip from computation */
void compute_delete_ellip(size_t bodnum)
{
  if (compute::inserted_ellips.count(bodnum))
    compute::inserted_ellips.erase(bodnum);
  else compute::deleted_ellips.insert(bodnum);
}

/* insert solfec::restrains[resnum] into computation */
void compute_insert_restrain(size_t resnum)
{
  compute::inserted_restrains.insert(resnum);
}

/* delete restrain from computation */
void compute_delete_restrain(size_t resnum)
{
  if (compute::inserted_restrains.count(resnum))
    compute::inserted_restrains.erase(resnum);
  else compute::deleted_restrains.insert(resnum);
}

/* insert solfec::prescribes[prenum] into computation */
void compute_insert_prescribe(size_t prenum)
{
  compute::inserted_prescribes.insert(prenum);
}

/* delete prescribe from computation */
void compute_delete_prescribe(size_t prenum)
{
  if (compute::inserted_prescribes.count(prenum))
    compute::inserted_prescribes.erase(prenum);
  else compute::deleted_prescribes.insert(prenum);
}

/* join compute main loop */
void compute_main_loop()
{
  using namespace compute;
  int rank, size;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    if (partitioned && (!inserted_meshes.empty() || !deleted_meshes.empty() ||
                        !inserted_ellips.empty() || !deleted_ellips.empty() ||
                        !inserted_restrains.empty() || !deleted_restrains.empty() ||
                        !inserted_prescribes.empty() || !deleted_prescribes.empty()))
    {
      if (!deleted_meshes.empty())
      {
       /*
	  TODO:
	  delete data bunches from individual MPI windows

	  rebalances bunches in windows for uniform load
        */
      }
      if (!deleted_ellips.empty())
      {
	/* TODO */
      }
      if (!deleted_restrains.empty())
      {
	/* TODO */
      }
      if (!deleted_prescribes.empty())
      {
	/* TODO */
      }

      if (!inserted_meshes.empty())
      {
       /*
          TODO:
          parition inserted meshes
	 
	  if (not enough space) resize windows

	  uniformly migrate inserted bunches into MPI windows with most free space
       */
      }
      if (!inserted_ellips.empty())
      {
	/* TODO */
      }
      if (!inserted_restrains.empty())
      {
	/* TODO */
      }
      if (!inserted_prescribes.empty())
      {
	/* TODO */
      }

      
    }
    else if (!partitioned)
    {
      std::map<size_t, part> parts = partition_meshes(inserted_meshes);

      size_t etot = 0, ftot = 0, ntot = 0;

      for (auto& [bodnum, part] : parts)
      {
        etot += part.neparts; /* total number of element parititons */
	ftot += part.nfparts; /* total number of face partitions */
	ntot += part.nn; /* total numbre of nodes */
      }

      size_t esplit = etot / size, fsplit = ftot / size;

      int currank = 0;

      for (auto& [bodnum, part] : parts)
      {
	/* TODO: map local to global data */
      }

      /* TODO: use one MPI_Scatterv for uint64_t data and one for REAL */

      /* TODO: use MPI_Scatterv to distribute mesh partitioning information to all ranks and then construct data bunches locally */

      /* TODO: use MPI_Scatterv to distribute ellip partitioning ... */

      /* TODO: use MPI to distribute restrains and prescribes ... */

      partitioned = true;
    }
  }
}

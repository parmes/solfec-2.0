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
#include "real.h"
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
std::map<size_t,material> materials; /* local MPI rank copy of all nominal materials */
MPI_Win mpi_nodes_window; /* storing nodes and dgreees of freedom data */
MPI_Win mpi_elements_window; /* storing elements data */
MPI_Win mpi_faces_window; /* storing faces data */
};

/* repartition all solfec::meshes from MPI rank 0 process and distribute them into
 * compute::mpi_nodes_window, compute::mpi_elements_window, compute::mpi_faces_window
 structures */
static void repartition_meshes()
{
  tf::Executor executor;
  tf::Taskflow taskflow;

  taskflow.parallel_for(solfec::meshes.begin(), solfec::meshes.end(), [&] (const std::pair<size_t,mesh> &it)
  { 
    std::vector<idx_t> eptr, eind, epart, npart;
    idx_t nn, ne, ncommon, nparts, objval;
    const struct mesh &mesh = it.second;

    nn = mesh.nodes.size();
    ne = mesh.nhex+mesh.nwed+mesh.npyr+mesh.ntet;
    ncommon = 3;
    nparts = 1+ne/ELEMENTS_BUNCH;
    eptr.push_back(0);
    for(auto e = mesh.elements.begin(); e != mesh.elements.end();)
    {
      for (auto i = e+1; i != e+1+e[0]; i++)
      {
	eind.push_back(*i);
      }
      e += e[0]+2;
      eptr.push_back(e-mesh.elements.begin());
    }
    epart.resize(ne);
    npart.resize(nn);

    /* partition elements into <= ELEMENTS_BUNCH sized sets */
    METIS_PartMeshDual(&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &ncommon,
                       &nparts, NULL, NULL, &objval, &epart[0], &npart[0]);

    /* TODO: using mesh.gcolor and mesh.colors create vector of all mesh faces */

    /* TODO: partition faces into <= FACES_BUNCH sized sets */

    /* TODO: create element and face bunch structures */
  });

  executor.run(taskflow).get();
}

/* insert solfec::meshes[bodnum] into computation */
void compute_insert_mesh(size_t bodnum)
{
}

/* delete mesh from computation */
void compute_delete_mesh(size_t bodnum)
{
}

/* insert solfec::ellips[bodnum] into computation */
void compute_insert_ellip(size_t bodnum)
{
}

/* delete ellip from computation */
void compute_delete_ellip(size_t bodnum)
{
}

/* insert solfec::restrains[resnum] into computation */
void compute_insert_restrain(size_t resnum)
{
}

/* delete restrain from computation */
void compute_delete_restrain(size_t resnum)
{
}

/* insert solfec::prescribes[prenum] into computation */
void compute_insert_prescribe(size_t prenum)
{
}

/* delete prescribe from computation */
void compute_delete_prescribe(size_t prenum)
{
}

/* join compute main loop */
void compute_main_loop()
{
  /* TODO rank 0:

    if (partitioned and inserted/deleted bodies)

      if (something to delete)

	delete bunches from individual MPI windows
	
	rebalances bunches in windows for uniform load

      if (inserted meshes) parition them

      if (not enough space)
	 
	 resize windows

      uniformly migrate inserted bunches into MPI windows with most free space

    else if (not partitioned)
       
       partition meshes
       
       create node sets, element bunches and face bunches and distribute them into MPI windows,
       in a way that balances parallel processing and memory usage
  */
}

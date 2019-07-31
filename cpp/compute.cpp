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
static void compute_repartition_meshes()
{
}

/* insert solfec::meshes[bodnum] into the existing window partitioning */
void compute_insert_mesh(size_t bodnum)
{
}

/* delete solfec::meshes[bodnum] from the existing window partitioning */
void compute_delete_mesh(size_t bodnum)
{
}

/* join compute main loop */
void compute_main_loop()
{
}

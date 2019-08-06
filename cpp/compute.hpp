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

#ifndef __compute__
#define __compute__

struct nodes_window /* nodes in mpi_nodes_window */
{
  REAL *refpos[3];
  REAL *curpos[3];
  REAL *velo[3];
  REAL *disp[3];
  uint64_t size;
};

#ifndef ELEMENTS_BUNCH
#define ELEMENTS_BUNCH 16
#endif

struct elements_bunch
{
  uint64_t node[2][ELEMENTS_BUNCH]; /* node[0] - MPI rank, node[1] - window index, of element nodes */

  REAL refpos[3][ELEMENTS_BUNCH]; /* local copy of nodes refpos */
  REAL curpos[3][ELEMENTS_BUNCH]; /* local copy of nodes curpos */
  REAL velo[3][ELEMENTS_BUNCH]; /* local copy of nodes velo */
  REAL disp[3][ELEMENTS_BUNCH]; /* local copy of nodes disp */

  short type[ELEMENTS_BUNCH]; /* 4, 5, 6, 8 => tetrahedron, pyramid, wedge, hexahedron */
  short lnod[8][ELEMENTS_BUNCH]; /* local node indices */

  uint64_t matnum[ELEMENTS_BUNCH]; /* material number */
  REAL density[ELEMENTS_BUNCH]; /* local material data */
  REAL young[ELEMENTS_BUNCH];
  REAL poisson[ELEMENTS_BUNCH];
  REAL viscosity[ELEMENTS_BUNCH];

  /* TODO: restrains and prescribes */

  uint64_t bodnum; /* body number */

  short size; /* bunch's size */
};

struct elements_window /* element bunches in mpi_elements_window */
{
  elements_bunch *bunch;
  uint64_t size;
};

#ifndef FACES_BUNCH
#define FACES_BUNCH 16
#endif

struct faces_bunch
{
  uint64_t node[2][FACES_BUNCH]; /* node[0] - MPI rank, node[1] - window index, of face nodes */

  REAL curpos[3][FACES_BUNCH]; /* local copy of nodes curpos */
  
  short type[FACES_BUNCH]; /* 3, 4 => triangle, quadrilateral */
  short lnod[4][FACES_BUNCH]; /* local node indices */
  short color[FACES_BUNCH];
  REAL normal[3][FACES_BUNCH];

  uint64_t elem[3][FACES_BUNCH]; /* elem[0] - MPI rank, elem[1] - window index, elem[2] - bunch index, of face element */

  uint64_t bodnum; /* body number */

  short size; /* bunch's size */
};

struct faces_window /* face bunches in mpi_faces_window */
{
  faces_bunch *bunch;
  uint64_t size;
};

/* insert solfec::meshes[bodnum] into computation */
void compute_insert_mesh(uint64_t bodnum);

/* delete mesh from computation */
void compute_delete_mesh(uint64_t bodnum);

/* insert solfec::ellips[bodnum] into computation */
void compute_insert_ellip(uint64_t bodnum);

/* delete ellip from computation */
void compute_delete_ellip(uint64_t bodnum);

/* insert solfec::restrains[resnum] into computation */
void compute_insert_restrain(uint64_t resnum);

/* delete restrain from computation */
void compute_delete_restrain(uint64_t resnum);

/* insert solfec::prescribes[prenum] into computation */
void compute_insert_prescribe(uint64_t prenum);

/* delete prescribe from computation */
void compute_delete_prescribe(uint64_t prenum);

/* join compute main loop */
void compute_main_loop();

#endif

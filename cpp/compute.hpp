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

#include "ga.hpp"

#ifndef __compute__
#define __compute__

struct mapping /* global array per rank entity ranges mapping */
{
  std::vector<std::array<uint64_t,3>> ga_nranges; /* externally populated global array ranges of same rank nodes */
  std::vector<std::array<uint64_t,3>> ga_eranges; /* externally populated global array ranges of same rank elements */
  std::vector<std::array<uint64_t,3>> ga_franges; /* externally populated global array ranges of same rank faces */
};

namespace compute
{
/* all ranks --- */
extern int ELEMENTS_BUNCH; /* elements SIMD bunch size */

extern int FACES_BUNCH; /* faces SIMD bunch size */

extern bool debug_print; /* enable debug printing */

extern GA *ga_counters; /* global array of MPI_UINT64_T counters; per rank:
			[count of materials
			 size of materials,
			 new size of material,
			 count of nodes,
			 size of nodes,
			 new size of nodes,
			 count of elements,
			 size of elements,
			 new size of elements,
			 count of faces,
			 size of faces,
			 new size of faces,
			 count of ellips,
			 size of ellips,
			 new size of ellips] */
enum {cn_materials, sz_materials, sz_materials_new,
      cn_nodes, sz_nodes, sz_nodes_new,
      cn_elements, sz_elements, sz_elements_new,
      cn_faces, sz_faces, sz_faces_new,
      cn_ellips, sz_ellips, sz_ellips_new,
      cn_last};

extern GA *ga_materials; /* global array of materials */
enum {mt_density,
      mt_young,
      mt_poisson,
      mt_viscosity,
      mt_last};

extern GA *ga_nodes; /* global array of nodal data */
enum {nd_vx, nd_vy, nd_vz,   /* linear velocity */
      nd_x, nd_y, nd_z,      /* current position */
      nd_X, nd_Y, nd_Z,      /* reference position */
      nd_last};

extern GA *ga_elements; /* global array of element data */
enum {el_bodnum,
      el_matnum,
      el_type, /* 4, 5, 6, 8 => tetrahedron, pyramid, wedge, hexahedron */
      el_nd0_rnk, el_nd0_idx, el_nd1_rnk, el_nd1_idx, el_nd2_rnk, el_nd2_idx, el_nd3_rnk, el_nd3_idx,
      el_nd4_rnk, el_nd4_idx, el_nd5_rnk, el_nd5_idx, el_nd6_rnk, el_nd6_idx, el_nd7_rnk, el_nd7_idx,
      el_last};

extern GA *ga_faces; /* global array of face data */
enum {fa_bodnum,
      fa_color,
      fa_type, /* 3, 4 => triangle, quadrilateral */
      fa_nd0_rnk, fa_nd0_idx, fa_nd1_rnk, fa_nd1_idx, fa_nd2_rnk, fa_nd2_idx, fa_nd3_rnk, fa_nd3_idx,
      fa_last};

extern GA *ga_elldata; /* global array of ellipsoids data */
enum {ll_vF0, ll_vF1, ll_vF2, ll_vF3, ll_vF4, ll_vF5, ll_vF6, ll_vF7, ll_vF8, /* deformation gradient velocity */
      ll_vx, ll_vy, ll_vz,                                                    /* linear velocity */
      ll_F0, ll_F1, ll_F2, ll_F3, ll_F4, ll_F5, ll_F6, ll_F7, ll_F8,          /* deformation gradient */
      ll_x, ll_y, ll_z,                                                       /* current position */
      ll_a, ll_b, ll_c,                                                       /* current radii */
      ll_r0, ll_r1, ll_r2, ll_r3, ll_r4, ll_r5, ll_r6, ll_r7, ll_r8,          /* current rotation */
      ll_X, ll_Y, ll_Z,                                                       /* reference position */
      ll_A, ll_B, ll_C,                                                       /* reference radii */
      ll_R0, ll_R1, ll_R2, ll_R3, ll_R4, ll_R5, ll_R6, ll_R7, ll_R8,          /* reference rotation */
      ll_last0};

extern GA *ga_ellips; /* global array of ellipsoids */
enum {ll_bodnum,
      ll_matnum,
      ll_color,
      ll_rnk, ll_idx, /* ellipsoid data rank and index */
      ll_last1};

extern GA *ga_contact; /* global array of contact REAL entities */
enum {co_rank, /* origin rank */
      co_index, /* origin rankn index */ /* XXX: uint64_t to REAL */
      co_type, /* 1.0 --> ellipsoid, 3.0 --> triangle, 4.0 --> quad */
      co_en0, co_en1, co_en2, /* (ll_x, ll_y, ll_z) or face node (x0, y0, z0) */
      co_en3, co_en4, co_en5, /* (ll_a, ll_b, ll_c) or face node (x1, y1, z1) */
      co_en6, co_en7, co_en8, /* (ll_r0, ll_r1, ll_r2) or face node (x2, y2, z2) */
      co_en9, co_en10, co_en11, /* (ll_r3, ll_r4, ll_r5) or unused */
      co_en12, co_en13, co_en14, /* (ll_r6, ll_r7, ll_r8) or unused */
      co_last};
/* --- all ranks */
};

/* insert solfec::materials[matnum] into computation */
void compute_insert_material(uint64_t matnum);

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
void compute_main_loop(REAL duration, REAL step);

/* finalize compute memory */
void compute_finalize();

#endif

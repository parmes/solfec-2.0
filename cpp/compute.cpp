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

#include <algorithm>
#include <mpi.h>
#include <map>
#include <set>
#include "real.h"
#include "err.h"
#include "ga.hpp"
#include "part.hpp"
#include "solfec.hpp"
#include "compute.hpp"

/* compute gobal variables */
namespace compute
{
/* rank 0 --- */
std::set<uint64_t> inserted_materials;
std::set<uint64_t> inserted_meshes;
std::set<uint64_t> deleted_meshes;
std::set<uint64_t> inserted_ellips;
std::set<uint64_t> deleted_ellips;
std::set<uint64_t> inserted_restrains;
std::set<uint64_t> deleted_restrains;
std::set<uint64_t> inserted_prescribes;
std::set<uint64_t> deleted_prescribes;
/* --- rank 0 */

/* all ranks --- */
int ELEMENTS_BUNCH = 16; /* elements SIMD bunch size */

int FACES_BUNCH = 16; /* faces SIMD bunch size */

bool partitioned = false; /* initially paritioned */

GA *ga_counters; /* global array of MPI_UINT64_T counters; per rank:
		    [count of materials
		     size of materials,
		     count of nodes,
		     size of nodes,
		     count of elements,
		     size of elements,
		     count of faces,
		     size of faces] */
enum {cn_materials, sz_materials,
      cn_nodes, sz_nodes,
      cn_elements, sz_elements,
      cn_faces, sz_faces, cn_last};

GA *ga_materials; /* global array of materials */
enum {mt_density,
      mt_young,
      mt_poisson,
      mt_viscosity,
      mt_last};

GA *ga_nodes; /* global array of nodal data */
enum {nd_X, nd_Y, nd_Z,
      nd_x, nd_y, nd_z,
      nd_dx, nd_dy, nd_dz,
      nd_vx, nd_vy, nd_vz,
      nd_last};

GA *ga_elements; /* global array of element data */
enum {el_bodnum,
      el_matnum,
      el_type, /* 4, 5, 6, 8 => tetrahedron, pyramid, wedge, hexahedron */
      el_nd0_rnk, el_nd0_idx, el_nd1_rnk, el_nd1_idx, el_nd2_rnk, el_nd2_idx, el_nd3_rnk, el_nd3_idx,
      el_nd4_rnk, el_nd4_idx, el_nd5_rnk, el_nd5_idx, el_nd6_rnk, el_nd6_idx, el_nd7_rnk, el_nd7_idx,
      el_last};

GA *ga_faces; /* global array of face data */
enum {fa_type, /* 3, 4 => triangle, quadrilateral */
      fa_nd0_rnk, fa_nd0_idx, fa_nd1_rnk, fa_nd1_idx, fa_nd2_rnk, fa_nd2_idx, fa_nd3_rnk, fa_nd3_idx,
      fa_color, fa_bodnum,
      fa_last};

/* --- all ranks */
};

/* insert solfec::materials[matnum] into computation */
void compute_insert_material(uint64_t matnum)
{
  compute::inserted_materials.insert(matnum);
}

/* insert solfec::meshes[bodnum] into computation */
void compute_insert_mesh(uint64_t bodnum)
{
  compute::inserted_meshes.insert(bodnum);
}

/* delete mesh from computation */
void compute_delete_mesh(uint64_t bodnum)
{
  if (compute::inserted_meshes.count(bodnum))
    compute::inserted_meshes.erase(bodnum);
  else compute::deleted_meshes.insert(bodnum);
}

/* insert solfec::ellips[bodnum] into computation */
void compute_insert_ellip(uint64_t bodnum)
{
  compute::inserted_ellips.insert(bodnum);
}

/* delete ellip from computation */
void compute_delete_ellip(uint64_t bodnum)
{
  if (compute::inserted_ellips.count(bodnum))
    compute::inserted_ellips.erase(bodnum);
  else compute::deleted_ellips.insert(bodnum);
}

/* insert solfec::restrains[resnum] into computation */
void compute_insert_restrain(uint64_t resnum)
{
  compute::inserted_restrains.insert(resnum);
}

/* delete restrain from computation */
void compute_delete_restrain(uint64_t resnum)
{
  if (compute::inserted_restrains.count(resnum))
    compute::inserted_restrains.erase(resnum);
  else compute::deleted_restrains.insert(resnum);
}

/* insert solfec::prescribes[prenum] into computation */
void compute_insert_prescribe(uint64_t prenum)
{
  compute::inserted_prescribes.insert(prenum);
}

/* delete prescribe from computation */
void compute_delete_prescribe(uint64_t prenum)
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

  if (!partitioned)
  {
    ga_counters = new GA(MPI_COMM_WORLD, 8, size, MPI_UINT64_T);

    ERRMEM (ga_counters);

    auto GA_ALL_CREATE = [](auto rank, auto size)
    {
      ga_counters->fence(); /* synchronize rank 0 puts */

      /* allocate global arrays */ 

      uint64_t counts[cn_last];

      ga_counters->get(0, cn_last, rank, rank+1, counts);

      ga_materials = new GA(MPI_COMM_WORLD, counts[sz_materials], mt_last*size, MPI_REAL);
      ga_nodes = new GA(MPI_COMM_WORLD, counts[sz_nodes], nd_last*size, MPI_REAL);
      ga_elements = new GA(MPI_COMM_WORLD, counts[sz_elements], el_last*size, MPI_UINT64_T);
      ga_faces = new GA(MPI_COMM_WORLD, counts[sz_faces], fa_last*size, MPI_UINT64_T);

      ERRMEM (ga_materials);
      ERRMEM (ga_nodes);
      ERRMEM (ga_elements);
      ERRMEM (ga_faces);
    };

    if (rank == 0) /* partition meshes and estimate array sizes */
    {
      std::map<uint64_t, part> parts = partition_meshes(inserted_meshes, ELEMENTS_BUNCH, FACES_BUNCH);

      std::map<uint64_t, mapping> maps = map_parts (parts);

      auto [maxnodes, maxeles, maxfaces] = max_per_rank (maps);

      for (int r = 0; r < size; r ++) /* initialize counters */
      {
	uint64_t counts[] = {0,
	                     inserted_materials.size() * 2,
                             0,
		             maxnodes * 2,
		             0,
		             maxeles * 2,
		             0,
		             maxfaces * 2};

	ga_counters->put (0, 8, r, r+1, counts);
      }

      GA_ALL_CREATE (rank, size); /* create global arrays */

      /* write materials */

      uint64_t matsize = inserted_materials.size();
      REAL *matdata = new REAL [matsize * mt_last];
      uint64_t matidx = 0;
      ERRMEM (matdata);

      for (auto& matnum : inserted_materials)
      {
	struct material &material = solfec::materials[matnum];
	matdata[mt_density*matsize + matidx] = material.density;
	matdata[mt_young*matsize + matidx] = material.young;
	matdata[mt_poisson*matsize + matidx] = material.poisson;
	matdata[mt_viscosity*matsize + matidx] = material.viscosity;
	matidx ++;
      }

      for (int r = 0; r < size; r ++)
      {
	ga_materials->put(0, matidx, r*mt_last, (r+1)*mt_last, matdata);
	ga_counters->acc(cn_materials, cn_materials+1, r, r+1, &matidx);
      }

      delete[] matdata;

      for (auto& [bodnum, map] : maps)
      {
         /* wrie nodes */

	for (auto r = map.nrank.begin(); r != map.nrank.end(); )
	{
	  auto eqr = std::equal_range(r, map.nrank.end(), *r);

	  uint64_t nodsize = eqr.second - eqr.first;
	  REAL *noddata = new REAL [nodsize * nd_last];
	  uint64_t nodidx = 0;
	  ERRMEM (noddata);

	  for (auto k = eqr.first; k != eqr.second; k++)
	  {
	    auto i = k - map.nrank.begin();

	    struct mesh &mesh = solfec::meshes[bodnum];
	    REAL x = mesh.nodes[i][0],
	         y = mesh.nodes[i][1],
	         z = mesh.nodes[i][2];
	    noddata[nd_X*nodsize + nodidx] = x;
	    noddata[nd_Y*nodsize + nodidx] = y;
	    noddata[nd_Z*nodsize + nodidx] = z;
	    noddata[nd_x*nodsize + nodidx] = x;
	    noddata[nd_y*nodsize + nodidx] = y;
	    noddata[nd_z*nodsize + nodidx] = z;
	    noddata[nd_dx*nodsize + nodidx] = 0.;
	    noddata[nd_dy*nodsize + nodidx] = 0.;
	    noddata[nd_dz*nodsize + nodidx] = 0.;
	    noddata[nd_vx*nodsize + nodidx] = 0.;
	    noddata[nd_vy*nodsize + nodidx] = 0.;
	    noddata[nd_vz*nodsize + nodidx] = 0.;
	    nodidx ++;
	  }

	  uint64_t count;
	  ga_counters->get(cn_nodes, cn_nodes+1, *r, (*r)+1, &count);
	  ga_nodes->put(count, count+nodidx, (*r)*nd_last, ((*r)+1)*nd_last, noddata);
	  ga_counters->acc(cn_nodes, cn_nodes+1, *r, (*r)+1, &nodidx);

	  delete[] noddata;

	  r = eqr.second;
	}

	/* write elements */

	struct part &part = parts[bodnum];

        for (auto r = map.erank.begin(); r != map.erank.end(); )
	{
	  auto eqr = std::equal_range(r, map.erank.end(), *r);

	  uint64_t elesize = eqr.second - eqr.first;
	  uint64_t *eledata = new uint64_t [elesize * el_last];
	  uint64_t eleidx = 0;
	  ERRMEM (eledata);

	  for (auto k = eqr.first; k != eqr.second; k++)
	  {
	    auto i = k - map.erank.begin();
	    auto j = part.eptr[i];
	    auto eltype = part.eptr[i+1]-j;

	    eledata[el_bodnum*elesize + eleidx] = bodnum;
	    eledata[el_matnum*elesize + eleidx] = part.material[i];
	    eledata[el_type*elesize + eleidx] = eltype;

	    eledata[el_nd0_rnk*elesize + eleidx] = map.nrank[j];
	    eledata[el_nd0_idx*elesize + eleidx] = map.nindex[j];
	    eledata[el_nd1_rnk*elesize + eleidx] = map.nrank[j+1];
	    eledata[el_nd1_idx*elesize + eleidx] = map.nindex[j+1];
	    eledata[el_nd2_rnk*elesize + eleidx] = map.nrank[j+2];
	    eledata[el_nd2_idx*elesize + eleidx] = map.nindex[j+2];
	    eledata[el_nd3_rnk*elesize + eleidx] = map.nrank[j+3];
	    eledata[el_nd3_idx*elesize + eleidx] = map.nindex[j+3];
	    if (eltype > 4)
	    { eledata[el_nd4_rnk*elesize + eleidx] = map.nrank[j+4];
	      eledata[el_nd4_idx*elesize + eleidx] = map.nindex[j+4]; }
	    if (eltype > 5)
	    { eledata[el_nd5_rnk*elesize + eleidx] = map.nrank[j+5];
	      eledata[el_nd5_idx*elesize + eleidx] = map.nindex[j+5]; }
	    if (eltype > 7)
	    { eledata[el_nd6_rnk*elesize + eleidx] = map.nrank[j+6];
	      eledata[el_nd6_idx*elesize + eleidx] = map.nindex[j+6];
	      eledata[el_nd7_rnk*elesize + eleidx] = map.nrank[j+7];
	      eledata[el_nd7_idx*elesize + eleidx] = map.nindex[j+7]; }

	    eleidx ++;
	  }

	  uint64_t count;
	  ga_counters->get(cn_elements, cn_elements+1, *r, (*r)+1, &count);
	  ga_nodes->put(count, count+eleidx, (*r)*el_last, ((*r)+1)*el_last, eledata);
	  ga_counters->acc(cn_elements, cn_elements+1, *r, (*r)+1, &eleidx);

	  delete[] eledata;

	  r = eqr.second;
	}

        /* write faces */

        /* TODO */
      }
    }
    else /* create global arrays */
    {
      GA_ALL_CREATE (rank, size);
    }

    ga_materials->fence();
    ga_nodes->fence();
    ga_elements->fence();
    ga_faces->fence();
    partitioned = true;
  }
  else
  {
    if (rank == 0)
    {
      if (!deleted_meshes.empty())
      {
       /* TODO */
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

      if (!inserted_materials.empty())
      {
       /* TODO */
      }
      if (!inserted_meshes.empty())
      {
       /* TODO */
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
  }

  /* TODO: compute */
}

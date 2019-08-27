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

#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <limits>
#include <mpi.h>
#include <map>
#include <set>
#include "real.h"
#include "err.h"
#include "alg.h"
#include "part.hpp"
#include "mesh.hpp"
#include "dynlb.hpp"
#include "debug.hpp"
#include "solfec.hpp"
#include "compute.hpp"
#include "alloc_ispc.h"

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

std::map<uint64_t, mapping> mesh_mapping; /* body number to mesh rank data mapping */
std::map<int, std::vector<std::pair<uint64_t, uint64_t>>>  deleted_nodes; /* rank mapping of deleted node ranges */
std::vector<uint64_t> deleted_nodes_count; /* deleted nodes counts per rank */
std::vector<uint64_t> deleted_node_index; /* deleted node index per rank */
std::vector<uint64_t> free_node_index; /* free node index per rank */
std::unordered_map<uint64_t, int>  ellip_mapping; /* ellipsoid body number to rank mapping */
/* --- rank 0 */

/* all ranks --- */
int ELEMENTS_BUNCH = 16; /* elements SIMD bunch size */

int FACES_BUNCH = 16; /* faces SIMD bunch size */

bool partitioned = false; /* initially paritioned */

bool debug_print = false; /* enable debug printing */

GA *ga_counters; /* global array of counters */
GA *ga_materials; /* global array of materials */
GA *ga_nodes; /* global array of nodal data */
GA *ga_elements; /* global array of element data */
GA *ga_faces; /* global array of face data */
GA *ga_elldata; /* global array of ellipsoids data */
GA *ga_ellips; /* global array of ellipsoids */
GA *ga_contact; /* global array of contact REAL entities */
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
void compute_main_loop(REAL duration, REAL step)
{
  while (true)
  {
    /* sync duration and step */

    REAL msg[2] = {duration, step};
    MPI_Bcast (msg, 2, MPI_REAL, 0, MPI_COMM_WORLD);
    duration = msg[0]; step = msg[1];

    if (duration == 0.0) break; /* rank 0 process communicated termination */

    using namespace compute;
    int rank, size;

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    auto GA_INSERT_MATERIALS = [](auto size)
    {
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
	uint64_t count;
	ga_counters->get(r, cn_materials, cn_materials+1, 0, 1, &count);
	ga_materials->put(r, count, count+matidx, 0, mt_last, matdata);
	ga_counters->acc(r, cn_materials, cn_materials+1, 0, 1, &matidx);
      }

      delete[] matdata;

      inserted_materials.clear();
    };

    auto GA_INSERT_MESHES = [](auto parts, auto maps)
    {
      for (auto& [bodnum, map] : maps)
      {
	struct part &part = parts[bodnum];

	/* wrie nodes */

	REAL center[3] = {0., 0., 0.};
	REAL linear[3] = {0., 0., 0.};
	REAL angular[3] = {0., 0., 0.};
	auto vi = solfec::velocities.find(bodnum);
	if (vi != solfec::velocities.end())
	{
	  mesh_char (bodnum, part.material, part.eptr, part.eind, NULL, center, NULL, NULL);

	  for (auto& it : (*vi).second)
	  {
	    ACC (it.linear_values, linear);
	    ACC (it.angular_values, angular);
	  }
	}

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

	    if (deleted_nodes_count[*r] > 0)
	    {
	      auto &stack = deleted_nodes[*r];
	      auto &back = stack.back();

	      if (deleted_node_index[*r] == std::numeric_limits<uint64_t>::max())
	      {
		deleted_node_index[*r] = back.first;
	      }
	      else if (deleted_node_index[*r]+1 < back.second)
	      {
		deleted_node_index[*r] ++;
	      }
	      else
	      {
		stack.pop_back();
		deleted_node_index[*r] = stack.back().first;
	      }

	      map.nindex[i] = deleted_node_index[*r];
	      deleted_nodes_count[*r] --;
	    }
	    else
	    {
	      map.nindex[i] = free_node_index[*r];
	      free_node_index[*r] ++;
	    }

	    struct mesh &mesh = solfec::meshes[bodnum];
	    REAL x = mesh.nodes[0][i],
		 y = mesh.nodes[1][i],
		 z = mesh.nodes[2][i];
	    REAL a[3] = {x - center[0],
			 y - center[1],
			 z - center[2]};
	    REAL v[3] = {linear[0],
			 linear[1],
			 linear[2]};
	    PRODUCTADD (v, a, angular);
	    noddata[nd_vx*nodsize + nodidx] = v[0];
	    noddata[nd_vy*nodsize + nodidx] = v[1];
	    noddata[nd_vz*nodsize + nodidx] = v[2];
	    noddata[nd_x*nodsize + nodidx] = x;
	    noddata[nd_y*nodsize + nodidx] = y;
	    noddata[nd_z*nodsize + nodidx] = z;
	    noddata[nd_X*nodsize + nodidx] = x;
	    noddata[nd_Y*nodsize + nodidx] = y;
	    noddata[nd_Z*nodsize + nodidx] = z;
	    nodidx ++;
	  }

	  uint64_t count;
	  ga_counters->get(*r, cn_nodes, cn_nodes+1, 0, 1, &count);
	  ga_nodes->put(*r, count, count+nodidx, 0, nd_last, noddata);
	  ga_counters->acc(*r, cn_nodes, cn_nodes+1, 0, 1, &nodidx);

	  std::array<uint64_t,3> rng = {(unsigned)*r, count, count+nodidx};
	  map.ga_nranges.push_back(rng); /* handy in deletion code */

	  delete[] noddata;

	  r = eqr.second;
	}

	/* write elements */

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
	    eledata[el_flags*elesize + eleidx] = 0;
	    eledata[el_type*elesize + eleidx] = eltype;

	    eledata[el_nd0_rnk*elesize + eleidx] = map.nrank[part.eind[j]];
	    eledata[el_nd0_idx*elesize + eleidx] = map.nindex[part.eind[j]];
	    eledata[el_nd1_rnk*elesize + eleidx] = map.nrank[part.eind[j+1]];
	    eledata[el_nd1_idx*elesize + eleidx] = map.nindex[part.eind[j+1]];
	    eledata[el_nd2_rnk*elesize + eleidx] = map.nrank[part.eind[j+2]];
	    eledata[el_nd2_idx*elesize + eleidx] = map.nindex[part.eind[j+2]];
	    eledata[el_nd3_rnk*elesize + eleidx] = map.nrank[part.eind[j+3]];
	    eledata[el_nd3_idx*elesize + eleidx] = map.nindex[part.eind[j+3]];
	    if (eltype > 4)
	    { eledata[el_nd4_rnk*elesize + eleidx] = map.nrank[part.eind[j+4]];
	      eledata[el_nd4_idx*elesize + eleidx] = map.nindex[part.eind[j+4]]; }
	    else
	    { eledata[el_nd4_rnk*elesize + eleidx] = 0;
	      eledata[el_nd4_idx*elesize + eleidx] = 0; }
	    if (eltype > 5)
	    { eledata[el_nd5_rnk*elesize + eleidx] = map.nrank[part.eind[j+5]];
	      eledata[el_nd5_idx*elesize + eleidx] = map.nindex[part.eind[j+5]]; }
	    { eledata[el_nd5_rnk*elesize + eleidx] = 0;
	      eledata[el_nd5_idx*elesize + eleidx] = 0; }
	    if (eltype > 7)
	    { eledata[el_nd6_rnk*elesize + eleidx] = map.nrank[part.eind[j+6]];
	      eledata[el_nd6_idx*elesize + eleidx] = map.nindex[part.eind[j+6]];
	      eledata[el_nd7_rnk*elesize + eleidx] = map.nrank[part.eind[j+7]];
	      eledata[el_nd7_idx*elesize + eleidx] = map.nindex[part.eind[j+7]]; }
	    else
	    { eledata[el_nd6_rnk*elesize + eleidx] = 0;
	      eledata[el_nd6_idx*elesize + eleidx] = 0;
	      eledata[el_nd7_rnk*elesize + eleidx] = 0;
	      eledata[el_nd7_idx*elesize + eleidx] = 0; }

	    eleidx ++;
	  }

	  uint64_t count;
	  ga_counters->get(*r, cn_elements, cn_elements+1, 0, 1, &count);
	  ga_elements->put(*r, count, count+eleidx, 0, el_last, eledata);
	  ga_counters->acc(*r, cn_elements, cn_elements+1, 0, 1, &eleidx);

	  std::array<uint64_t,3> rng = {(unsigned)*r, count, count+eleidx};
	  map.ga_eranges.push_back(rng); /* handy in deletion code */

	  delete[] eledata;

	  r = eqr.second;
	}

	/* write faces */

	for (auto r = map.frank.begin(); r != map.frank.end(); )
	{
	  auto eqr = std::equal_range(r, map.frank.end(), *r);

	  uint64_t facsize = eqr.second - eqr.first;
	  uint64_t *facdata = new uint64_t [facsize * fa_last];
	  uint64_t facidx = 0;
	  ERRMEM (facdata);

	  for (auto k = eqr.first; k != eqr.second; k++)
	  {
	    auto i = k - map.frank.begin();
	    auto j = part.fptr[i];
	    auto factype = part.fptr[i+1]-j;

	    facdata[fa_bodnum*facsize + facidx] = bodnum;
	    facdata[fa_color*facsize + facidx] = part.color[i];
	    facdata[fa_flags*facsize + facidx] = 0;
	    facdata[fa_type*facsize + facidx] = factype;

	    facdata[fa_nd0_rnk*facsize + facidx] = map.nrank[j];
	    facdata[fa_nd0_idx*facsize + facidx] = map.nindex[j];
	    facdata[fa_nd1_rnk*facsize + facidx] = map.nrank[j+1];
	    facdata[fa_nd1_idx*facsize + facidx] = map.nindex[j+1];
	    facdata[fa_nd2_rnk*facsize + facidx] = map.nrank[j+2];
	    facdata[fa_nd2_idx*facsize + facidx] = map.nindex[j+2];
	    if (factype > 3)
	    { facdata[fa_nd3_rnk*facsize + facidx] = map.nrank[j+3];
	      facdata[fa_nd3_idx*facsize + facidx] = map.nindex[j+3]; }
	    else
	    { facdata[fa_nd3_rnk*facsize + facidx] = 0;
	      facdata[fa_nd3_idx*facsize + facidx] = 0; }

	    facidx ++;
	  }

	  uint64_t count;
	  ga_counters->get(*r, cn_faces, cn_faces+1, 0, 1, &count);
	  ga_faces->put(*r, count, count+facidx, 0, fa_last, facdata);
	  ga_counters->acc(*r, cn_faces, cn_faces+1, 0, 1, &facidx);

	  std::array<uint64_t,3> rng = {(unsigned)*r, count, count+facidx};
	  map.ga_franges.push_back(rng); /* handy in deletion code */

	  delete[] facdata;

	  r = eqr.second;
	}
      }

      /* update global mesh mapping */

      mesh_mapping.merge (maps); /* move temporary mappings into the global container */

      inserted_meshes.clear();
    };

    auto GA_INSERT_ELLIPS = [](auto size)
    {
      uint64_t ellsize = inserted_ellips.size();
      REAL *elldata = new REAL [ellsize * ll_last0];
      uint64_t *ellips = new uint64_t [ellsize * ll_last1];
      uint64_t ellidx = 0, i;
      using namespace ispc;
      ERRMEM (elldata);
      ERRMEM (ellips);
      REAL *point[3] = {aligned_real_alloc(ellsize),
			aligned_real_alloc(ellsize),
			aligned_real_alloc(ellsize)};
      dynlb lb;
     
      i = 0;
      for (auto& bodnum : inserted_ellips)
      {
	struct ellip &ellip = solfec::ellips[bodnum];
	point[0][i] = ellip.center[0];
	point[1][i] = ellip.center[1];
	point[2][i] = ellip.center[2];
	i ++;
      }

      lb.local_create (ellidx, point); /* create local load balancer */
      aligned_real_free(point[0]);
      aligned_real_free(point[1]);
      aligned_real_free(point[2]);

      std::vector<int> ellip_rank(ellsize);
      i = 0;
      for (auto& bodnum : inserted_ellips)
      {
	struct ellip &ellip = solfec::ellips[bodnum];
	ellip_rank[i] = lb.point_assign(ellip.center); /* assign ranks */
	ellip_mapping[bodnum] = ellip_rank[i];
	i ++;
      }

      std::vector<uint64_t> ellinrank(size, 0);
      std::vector<uint64_t> ellip_index(ellsize);
      for (i = 0; i < ellsize; i ++)
      {
	ellip_index[i] = ellinrank[ellip_rank[i]];
	ellinrank[ellip_rank[i]] ++;
      }

      for (auto& bodnum : inserted_ellips)
      {
	struct ellip &ellip = solfec::ellips[bodnum];

	REAL linear[3] = {0., 0., 0.};
	REAL angular[3] = {0., 0., 0.};
	auto vi = solfec::velocities.find(bodnum);
	if (vi != solfec::velocities.end())
	{
	  for (auto& it : (*vi).second)
	  {
	    ACC (it.linear_values, linear);
	    ACC (it.angular_values, angular);
	  }
	}
	REAL L[9], F[9], vF[9];
	VECSKEW(angular, L);
	IDENTITY(F);
	NNMUL (L, F, vF);

	elldata[ll_vF0*ellsize + ellidx] = vF[0];
	elldata[ll_vF1*ellsize + ellidx] = vF[1];
	elldata[ll_vF2*ellsize + ellidx] = vF[2];
	elldata[ll_vF3*ellsize + ellidx] = vF[3];
	elldata[ll_vF4*ellsize + ellidx] = vF[4];
	elldata[ll_vF5*ellsize + ellidx] = vF[5];
	elldata[ll_vF6*ellsize + ellidx] = vF[6];
	elldata[ll_vF7*ellsize + ellidx] = vF[7];
	elldata[ll_vF8*ellsize + ellidx] = vF[8];
	elldata[ll_vx*ellsize + ellidx] = linear[0];
	elldata[ll_vy*ellsize + ellidx] = linear[1];
	elldata[ll_vz*ellsize + ellidx] = linear[2];

	elldata[ll_F0*ellsize + ellidx] = 1.;
	elldata[ll_F1*ellsize + ellidx] = 0.;
	elldata[ll_F2*ellsize + ellidx] = 0.;
	elldata[ll_F3*ellsize + ellidx] = 0.;
	elldata[ll_F4*ellsize + ellidx] = 1.;
	elldata[ll_F5*ellsize + ellidx] = 0.;
	elldata[ll_F6*ellsize + ellidx] = 0.;
	elldata[ll_F7*ellsize + ellidx] = 0.;
	elldata[ll_F8*ellsize + ellidx] = 1.;

	elldata[ll_x*ellsize + ellidx] = ellip.center[0];
	elldata[ll_y*ellsize + ellidx] = ellip.center[1];
	elldata[ll_z*ellsize + ellidx] = ellip.center[2];
	elldata[ll_a*ellsize + ellidx] = ellip.radius[0];
	elldata[ll_b*ellsize + ellidx] = ellip.radius[1];
	elldata[ll_c*ellsize + ellidx] = ellip.radius[2];
	elldata[ll_r0*ellsize + ellidx] = ellip.rotation[0];
	elldata[ll_r1*ellsize + ellidx] = ellip.rotation[1];
	elldata[ll_r2*ellsize + ellidx] = ellip.rotation[2];
	elldata[ll_r3*ellsize + ellidx] = ellip.rotation[3];
	elldata[ll_r4*ellsize + ellidx] = ellip.rotation[4];
	elldata[ll_r5*ellsize + ellidx] = ellip.rotation[5];
	elldata[ll_r6*ellsize + ellidx] = ellip.rotation[6];
	elldata[ll_r7*ellsize + ellidx] = ellip.rotation[7];
	elldata[ll_r8*ellsize + ellidx] = ellip.rotation[8];

	elldata[ll_X*ellsize + ellidx] = ellip.center[0];
	elldata[ll_Y*ellsize + ellidx] = ellip.center[1];
	elldata[ll_Z*ellsize + ellidx] = ellip.center[2];
	elldata[ll_A*ellsize + ellidx] = ellip.radius[0];
	elldata[ll_B*ellsize + ellidx] = ellip.radius[1];
	elldata[ll_C*ellsize + ellidx] = ellip.radius[2];
	elldata[ll_R0*ellsize + ellidx] = ellip.rotation[0];
	elldata[ll_R1*ellsize + ellidx] = ellip.rotation[1];
	elldata[ll_R2*ellsize + ellidx] = ellip.rotation[2];
	elldata[ll_R3*ellsize + ellidx] = ellip.rotation[3];
	elldata[ll_R4*ellsize + ellidx] = ellip.rotation[4];
	elldata[ll_R5*ellsize + ellidx] = ellip.rotation[5];
	elldata[ll_R6*ellsize + ellidx] = ellip.rotation[6];
	elldata[ll_R7*ellsize + ellidx] = ellip.rotation[7];
	elldata[ll_R8*ellsize + ellidx] = ellip.rotation[8];

	ellips[ll_bodnum*ellsize + ellidx] = bodnum;
	ellips[ll_matnum*ellsize + ellidx] = ellip.matnum;
	ellips[ll_color*ellsize + ellidx] = ellip.gcolor;
	ellips[ll_rnk*ellsize + ellidx] = ellip_rank[ellidx];
	ellips[ll_idx*ellsize + ellidx] = ellip_index[ellidx];

	ellidx ++;
      }

      for (int r = 0; r < size; r ++)
      {
	ga_elldata->put(r, 0, ellidx, 0, ll_last0, elldata);
	ga_ellips->put(r, 0, ellidx, 0, ll_last1, ellips);
	ga_counters->acc(r, cn_ellips, cn_ellips+1, 0, 1, &ellidx);
      }

      delete[] elldata;
      delete[] ellips;

      inserted_ellips.clear();
    };

    if (!partitioned)
    {
      deleted_nodes_count.resize(size);
      deleted_node_index.resize(size);
      free_node_index.resize(size);

      std::fill(deleted_nodes_count.begin(), deleted_nodes_count.end(), 0);
      std::fill(deleted_node_index.begin(), deleted_node_index.end(), std::numeric_limits<uint64_t>::max());
      std::fill(free_node_index.begin(), free_node_index.end(), 0);

      ga_counters = new GA(MPI_COMM_WORLD, cn_last, 1, MPI_UINT64_T);

      ERRMEM (ga_counters);

      auto GA_ALL_CREATE = [](auto rank, auto size)
      {
	ga_counters->fence(); /* synchronize rank 0 puts */

	/* allocate global arrays */ 

	uint64_t counts[cn_last];

	ga_counters->get(rank, 0, cn_last, 0, 1, counts);

	ga_materials = new GA(MPI_COMM_WORLD, counts[sz_materials], mt_last, MPI_REAL);
	ga_nodes = new GA(MPI_COMM_WORLD, counts[sz_nodes], nd_last, MPI_REAL);
	ga_elements = new GA(MPI_COMM_WORLD, counts[sz_elements], el_last, MPI_UINT64_T);
	ga_faces = new GA(MPI_COMM_WORLD, counts[sz_faces], fa_last, MPI_UINT64_T);
	ga_elldata = new GA(MPI_COMM_WORLD, counts[sz_ellips], ll_last0, MPI_REAL);
	ga_ellips = new GA(MPI_COMM_WORLD, counts[sz_ellips], ll_last1, MPI_UINT64_T);

	ERRMEM (ga_materials);
	ERRMEM (ga_nodes);
	ERRMEM (ga_elements);
	ERRMEM (ga_faces);
	ERRMEM (ga_elldata);
	ERRMEM (ga_ellips);
      };

      if (rank == 0) /* partition meshes and estimate array sizes */
      {
	std::map<uint64_t, part> parts = partition_meshes(inserted_meshes, ELEMENTS_BUNCH, FACES_BUNCH);

	std::map<uint64_t, mapping> maps = map_parts (parts);

	auto [maxnodes, maxeles, maxfaces] = max_per_rank (maps);

	for (int r = 0; r < size; r ++) /* initialize counters */
	{
	  uint64_t counts[] = {0, inserted_materials.size() * 2, 0,
			       0, maxnodes * 2, 0,
			       0, maxeles * 2, 0,
			       0, maxfaces * 2, 0,
			       0, inserted_ellips.size() * 2, 0};

	  ga_counters->put (r, 0, cn_last, 0, 1, counts);
	}

	GA_ALL_CREATE (rank, size); /* create global arrays */

	if (!inserted_materials.empty())
	{
	  GA_INSERT_MATERIALS (size);
	}

	if (!inserted_meshes.empty())
	{
	  GA_INSERT_MESHES (parts, maps);
	}

	if (!inserted_ellips.empty())
	{
	  GA_INSERT_ELLIPS(size);
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
      else /* create global arrays */
      {
	GA_ALL_CREATE (rank, size);
      }

      partitioned = true; /* initially partitioned */
    }
    else /* already initially partitioned */
    {
      auto GA_RESIZE_UP = [](auto rank)
      {
	ga_counters->fence();

	uint64_t counts[cn_last];

	ga_counters->get(rank, 0, cn_last, 0, 1, counts);

	if (counts[sz_materials_new] > counts[sz_materials])
	{
	  GA *ga = new GA(MPI_COMM_WORLD, counts[sz_materials_new], mt_last, MPI_REAL);
	  REAL *data = new REAL [counts[cn_materials] * mt_last];
	  ERRMEM (data);
	  ga_materials->get(rank, 0, counts[cn_materials], 0, mt_last, data);
	  ga->put(rank, 0, counts[cn_materials], 0, mt_last, data);
	  GA *tmp = ga_materials;
	  ga_materials = ga;
	  delete tmp;
	}

        if (counts[sz_nodes_new] > counts[sz_nodes])
	{
	  GA *ga = new GA(MPI_COMM_WORLD, counts[sz_nodes_new], nd_last, MPI_REAL);
	  REAL *data = new REAL [counts[cn_nodes] * nd_last];
	  ERRMEM (data);
	  ga_nodes->get(rank, 0, counts[cn_nodes], 0, nd_last, data);
	  ga->put(rank, 0, counts[cn_nodes], 0, nd_last, data);
	  GA *tmp = ga_nodes;
	  ga_nodes = ga;
	  delete tmp;
	}

	if (counts[sz_elements_new] > counts[sz_elements])
	{
	  GA *ga = new GA(MPI_COMM_WORLD, counts[sz_elements_new], el_last, MPI_UINT64_T);
	  uint64_t *data = new uint64_t [counts[cn_elements] * el_last];
	  ERRMEM (data);
	  ga_elements->get(rank, 0, counts[cn_elements], 0, el_last, data);
	  ga->put(rank, 0, counts[cn_elements], 0, el_last, data);
	  GA *tmp = ga_elements;
	  ga_elements = ga;
	  delete tmp;
	}

	if (counts[sz_faces_new] > counts[sz_faces])
	{
	  GA *ga = new GA(MPI_COMM_WORLD, counts[sz_faces_new], fa_last, MPI_UINT64_T);
	  uint64_t *data = new uint64_t [counts[cn_faces] * fa_last];
	  ERRMEM (data);
	  ga_faces->get(rank, 0, counts[cn_faces], 0, fa_last, data);
	  ga->put(rank, 0, counts[cn_faces], 0, fa_last, data);
	  GA *tmp = ga_faces;
	  ga_faces = ga;
	  delete tmp;
	}

	if (counts[sz_ellips_new] > counts[sz_ellips])
	{
	  GA *ga0 = new GA(MPI_COMM_WORLD, counts[sz_ellips_new], ll_last0, MPI_REAL);
	  GA *ga1 = new GA(MPI_COMM_WORLD, counts[sz_ellips_new], ll_last1, MPI_UINT64_T);
	  REAL *data0 = new REAL [counts[cn_ellips] * ll_last0];
	  uint64_t *data1 = new uint64_t [counts[cn_ellips] * ll_last1];
	  ERRMEM (data0);
	  ERRMEM (data1);
	  ga_elldata->get(rank, 0, counts[cn_ellips], 0, ll_last0, data0);
	  ga_ellips->get(rank, 0, counts[cn_ellips], 0, ll_last1, data1);
	  ga0->put(rank, 0, counts[cn_ellips], 0, ll_last0, data0);
	  ga1->put(rank, 0, counts[cn_ellips], 0, ll_last1, data1);
	  GA *tmp = ga_elldata;
	  ga_elldata = ga0;
	  delete tmp;
	  tmp = ga_ellips;
	  ga_ellips = ga1;
	  delete tmp;
	}

	counts[sz_materials] = counts[sz_materials_new];
	counts[sz_nodes] = counts[sz_nodes_new];
	counts[sz_elements] = counts[sz_elements_new];
	counts[sz_faces] = counts[sz_faces_new];
	counts[sz_ellips] = counts[sz_ellips_new];

	ga_counters->put(rank, 0, cn_last, 0, 1, counts);
      };
	
      if (rank == 0)
      {
	auto GA_RESIZE_DOWN = [](auto cn_name, auto ga_name, auto xx_last,
	  auto xx_flags, auto xx_flags_deleted, auto rank_deleted_ranges)
	{
	  for (auto& [r, vec] : rank_deleted_ranges)
	  {
	    uint64_t count0, count1;
	    ga_counters->get(r, cn_name, cn_name+1, 0, 1, &count0);

	    uint64_t *data = new uint64_t [count0 * xx_last];
	    ga_name->get(r, 0, count0, 0, xx_last, data);

	    count1 = count0;
	    for (auto& rng : vec)
	    {
	      for (auto i = rng.first; i < rng.second; i ++)
	      {
		uint64_t *item = &data[i*xx_last];
		item[xx_flags] |= xx_flags_deleted;
		count1 --;
	      }
	    }

	    for (uint64_t i = 0, j = 0, *item = data; i < count0; i ++, item += xx_last)
	    {
	      if (item[xx_flags] & xx_flags_deleted) continue;

	      if (j < i)
	      {
		std::memmove (&data[j*xx_last], item, sizeof(uint64_t [xx_last]));
	      }

	      j ++;
	    }

	    ga_name->put(r, 0, count1, 0, xx_last, data);
	    ga_counters->put(r, cn_name, cn_name+1, 0, 1, &count1);

	    delete [] data;
	  }
	};

	if (!deleted_meshes.empty())
	{
	  std::map<uint64_t, mapping>  deleted_mesh_mapping;
	  std::map<int, std::vector<std::pair<uint64_t,uint64_t>>> deleted_elements;
	  std::map<int, std::vector<std::pair<uint64_t,uint64_t>>> deleted_faces;

	  for (auto& bodnum : deleted_meshes)
	  {
	    auto nh = mesh_mapping.extract(bodnum);
	    deleted_mesh_mapping.insert(std::move(nh));
	    mapping &map = deleted_mesh_mapping[bodnum];

	    for (auto& r : map.ga_nranges)
	    {
	      deleted_nodes[r[0]].push_back(std::make_pair(r[1],r[2])); /* used during insertion of new meshes */
	      deleted_nodes_count[r[0]] += r[2]-r[1];
	    }

	    for (auto& r : map.ga_eranges)
	    {
	      deleted_elements[r[0]].push_back(std::make_pair(r[1],r[2]));
	    }

	    for (auto& r : map.ga_franges)
	    {
	      deleted_faces[r[0]].push_back(std::make_pair(r[1],r[2]));
	    }
	  }

	  /* resize element and face arrays */
	  GA_RESIZE_DOWN (cn_elements, ga_elements, el_last, el_flags, el_flags_deleted, deleted_elements);
	  GA_RESIZE_DOWN (cn_faces, ga_faces, fa_last, fa_flags, fa_flags_deleted, deleted_faces);
	}

	if (!deleted_ellips.empty())
	{
	  std::unordered_map<int, std::set<uint64_t>> deleted_ellips_onrank;

	  for (auto& bodnum : deleted_ellips)
	  {
	    deleted_ellips_onrank[ellip_mapping[bodnum]].insert(bodnum);
	  }

	  for (auto& [r, set] : deleted_ellips_onrank)
	  {
	    uint64_t count0, count1;
	    ga_counters->get(r, cn_ellips, cn_ellips+1, 0, 1, &count0);

	    REAL *data0 = new REAL [count0 * ll_last0];
	    uint64_t *data1 = new uint64_t [count0 * ll_last1];
	    ga_elldata->get(r, 0, count0, 0, ll_last0, data0);
	    ga_ellips->get(r, 0, count0, 0, ll_last1, data1);

	    count1 = count0;
	    REAL *item = data0;
	    uint64_t i = 0, j = 0, *jtem = data1;

	    for (; i < count0; i ++, item += ll_last0, jtem += ll_last1)
	    {
	      if (set.find(jtem[ll_bodnum]) != set.end())
	      {
		ellip_mapping.erase(jtem[ll_bodnum]);
		count1 --;
	      }

	      if (j < i)
	      {
		std::memmove (&data0[j*ll_last0], item, sizeof(REAL [ll_last0]));
		std::memmove (&data1[j*ll_last1], jtem, sizeof(uint64_t [ll_last1]));
	      }

	      j ++;
	    }

	    ga_elldata->put(r, 0, count1, 0, ll_last0, data0);
	    ga_ellips->put(r, 0, count1, 0, ll_last1, data1);
	    ga_counters->put(r, cn_ellips, cn_ellips+1, 0, 1, &count1);

	    delete [] data0;
	    delete [] data1;
	  }
	}

	if (!deleted_restrains.empty())
	{
	  /* TODO */
	}

	if (!deleted_prescribes.empty())
	{
	  /* TODO */
	}

        /* partition new meshes */

	std::map<uint64_t, part> parts = partition_meshes(inserted_meshes, ELEMENTS_BUNCH, FACES_BUNCH);

	std::map<uint64_t, mapping> maps = map_parts (parts);

	auto [maxnodes, maxeles, maxfaces] = max_per_rank (maps);

	/* first resize arrays if needed */

	uint64_t counts0[cn_last];

	ga_counters->get (0, 0, cn_last, 0, 1, counts0);

	counts0[sz_materials_new] = counts0[sz_materials];
	counts0[sz_nodes_new] = counts0[sz_nodes];
	counts0[sz_elements_new] = counts0[sz_elements];
	counts0[sz_faces_new] = counts0[sz_faces];
	counts0[sz_ellips_new] = counts0[sz_ellips];

	for (int r = 0; r < size; r ++) /* update counters */
	{
	  uint64_t counts[cn_last];

	  ga_counters->get (r, 0, cn_last, 0, 1, counts);

	  counts[sz_materials_new] = counts[cn_materials] + inserted_materials.size();
	  counts[sz_nodes_new] = counts[cn_nodes] - deleted_nodes_count[r] + maxnodes;
	  counts[sz_elements_new] = counts[cn_elements] + maxeles;
	  counts[sz_faces_new] = counts[cn_faces] + maxfaces,
	  counts[sz_ellips_new] = counts[cn_ellips] + inserted_ellips.size()/size + 1;

	  if (counts[sz_materials_new] > counts0[sz_materials_new]) counts0[sz_materials_new] = counts[sz_materials_new] * 2;
	  if (counts[sz_nodes_new] > counts0[sz_nodes_new]) counts0[sz_nodes_new] = counts[sz_nodes_new] * 2;
	  if (counts[sz_elements_new] > counts0[sz_elements_new]) counts0[sz_elements_new] = counts[sz_elements_new] * 2;
	  if (counts[sz_faces_new] > counts0[sz_faces_new]) counts0[sz_faces_new] = counts[sz_faces_new] * 2;
	  if (counts[sz_ellips_new] > counts0[sz_ellips_new]) counts0[sz_ellips_new] = counts[sz_ellips_new] * 2;
	}

	for (int r = 0; r < size; r ++) /* update counters */
	{
	  uint64_t counts[cn_last];

	  ga_counters->get (r, 0, cn_last, 0, 1, counts);

	  counts[sz_materials_new] = counts0[sz_materials_new];
	  counts[sz_nodes_new] = counts0[sz_nodes_new];
	  counts[sz_elements_new] = counts0[sz_elements_new];
	  counts[sz_faces_new] = counts0[sz_faces_new];
	  counts[sz_ellips_new] = counts0[sz_ellips_new];

	  ga_counters->put (r, 0, cn_last, 0, 1, counts);

	  free_node_index[r] = counts[cn_nodes];
	}

	GA_RESIZE_UP(rank);

	/* insert new itmes into arrays */

	if (!inserted_materials.empty())
	{
	  GA_INSERT_MATERIALS (size);
	}

	if (!inserted_meshes.empty())
	{
	  GA_INSERT_MESHES (parts, maps);
	}

	if (!inserted_ellips.empty())
	{
	  GA_INSERT_ELLIPS(size);
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
      else /* resize on ranks > 0 */
      {
	GA_RESIZE_UP(rank);
      }
    }

    /* sync all RMA ops on arrays */

    ga_counters->fence();
    ga_materials->fence();
    ga_nodes->fence();
    ga_elements->fence();
    ga_faces->fence();
    ga_elldata->fence();
    ga_ellips->fence();

    if (debug_print)
    {
      debug_print_compute_data();
    }

    /* TODO: create compute task graph including: */

    /* { */

      /* compute step */

      /* output solfec::histories at solfec::interval */

      /* output solfec::outputs at solfec::interval */

    /* } */

    for (REAL end = solfec::simulation_time + duration;
	 solfec::simulation_time < end; solfec::simulation_time += step)
    {
      /* TODO: exectue compute task graph */
    }

    if (rank == 0) break; /* exit and rejoin upon next RUN(), while others wait at MPI_Bcast */
  }
}

/* finalize compute memory */
void compute_finalize()
{
  using namespace compute;
  delete ga_counters;
  delete ga_nodes;
  delete ga_elements;
  delete ga_faces;
  delete ga_ellips;
}

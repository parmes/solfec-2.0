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
#include <iostream>
#include <cstring>
#include <numeric>
#include <limits>
#include <mpi.h>
#include <list>
#include <map>
#include <set>
#include "real.h"
#include "err.h"
#include "alg.h"
#include "fmt.hpp"
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

std::vector<uint64_t> nodes_per_rank; /* nodes per rank counts */
std::vector<uint64_t> elements_per_rank; /* element per rank counts */
std::vector<uint64_t> faces_per_rank; /* faces per rank counts */
std::map<uint64_t, mapping> mesh_mapping; /* body number to mesh rank data range mapping */
std::map<int, std::list<std::pair<uint64_t, uint64_t>>>  deleted_nodes; /* rank mapping of deleted node ranges */
std::vector<uint64_t> deleted_nodes_count; /* deleted nodes counts per rank */
std::unordered_map<uint64_t, int>  ellip_mapping; /* ellipsoid body number to rank mapping */
std::vector<uint64_t> ellips_per_rank; /* ellipsoids per rank counts */
/* --- rank 0 */

/* all ranks --- */
int ELEMENTS_BUNCH = 16; /* elements SIMD bunch size */

int FACES_BUNCH = 16; /* faces SIMD bunch size */

REAL IMBALANCE_TOLERANCE = 1.1; /* data imbalance tolerance */

REAL ELEMIGRATIO_THRESHOLD = 0.05; /* element migration ratio below which rank-0 GA read-write access is used */

bool partitioned = false; /* initially paritioned */

bool debug_print = false; /* enable debug printing */

bool debug_files = false; /* enable debug files output */

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

/* RANK 0: calculate mesh migration mapping;
 *** in/out: ***
 * r_b_rng - MPI_COMM_WORLD size vector of maps of bodnums to vectors of element/face ranges per rank;
 * modified internally and undefined upon return;
 * rng_type - range type as either 'elements' or 'faces';
 * append - append sendbuf, sendcounts, displs vectors flag;
 *** out: ***
 * (sendbuf, sendcounts, displs) as appropriate for MPI_Scatterv, where
 * sendbuf stores continguous tuples (bodnum, rank, type, org0, org1, dst0, dst1)
 * bodnums, their destination ranks, migrating entity types - 0/node, 1/element,
 * 2/face, origin global array ranges org0, org1, and destination global array
 * ranges dst0, dst1;
 *** modified globals: ***
 * nodes_per_rank, elements_per_rank, faces_per_rank, mesh_mapping, deleted_nodes,
 * deleted_nodes_counts */
static void mesh_migration_map (std::vector<std::map<uint64_t,std::vector<std::array<uint64_t,3>*>>> r_b_rng, /* (r)ank-(b)ody-(r)a(ng)e */
  std::string rng_type, std::vector<uint64_t> &sendbuf, std::vector<int> &sendcounts, std::vector<int> &displs, bool append=false)
{
  /* XXX 1: note that this is an element-balance driven migration scheme, where nodes and faces are migrated with elements;
            a more refined approach can supersede it, where nodes, elements, and faces are migrated independently */

  /* XXX 2: also note that currenlty, pre-existing per-rank continguous ranges of entities (elements) are used as data
            transfer increments, even though these originally coalesced from the smaller BUNCHes; this constraints the
            achievable degree of balance; a refinement is possible, where these coalesced ranges are again dismemberd into
            BUNCHes for the purpose of a finer balancing */

  /* XXX 3: finally note, that currently node movement is not implemented */

  std::vector<std::vector<std::pair<uint64_t,std::array<uint64_t,3>*>>> send (r_b_rng.size()); /* post-migrtion new rank to (body, origin ranges) map */
  bool rng_ele = rng_type == "elements" ? true : false;
  std::set<std::pair<uint64_t,int>> sorted; /* updatable set of pairs of (element count, rank) */
  using namespace compute;

  for (auto it = r_b_rng.begin(); it != r_b_rng.end(); it ++)
  {
    int rank = (int)(it - r_b_rng.begin());
    uint64_t totalcount = 0;
    for (auto& [bodnum, vector] : *it)
    {
      for (auto rng : vector)
      {
        totalcount += (*rng)[2]-(*rng)[1];
      }
    }
    sorted.insert(std::make_pair(totalcount, rank));
  }

  REAL im0 = (REAL)(--sorted.end())->first/(REAL)(sorted.begin()->first+1);

  if (debug_print) std::cout << "DBG: " << rng_type << " migration algo im0 = " << im0 << std::endl;

  while (true)
  {
    auto big = sorted.extract(--sorted.end()); /* largest element count rank */

    auto& val0 = big.value(); /* element count, rank */

    auto rnk = val0.second;

    std::array<uint64_t,3>* min_rng = nullptr;
    typedef std::vector<decltype(min_rng)>::iterator it_type;
    uint64_t min_bodnum;
    it_type min_it;

    for (auto& [bodnum, vector] : r_b_rng[rnk])
    {
      it_type it = min_element(vector.begin(), vector.end(), /* range of minimum element count */
                   [](auto l, auto r) { return ((*l)[2]-(*l)[1]) < ((*r)[2] - (*r)[1]); });

      if (it == vector.end()) break; /* empty range set case */

      if (min_rng == nullptr || ((**it)[2]-(**it)[1]) <= ((*min_rng)[2]-(*min_rng)[1]))
      {
        min_bodnum = bodnum;
        min_rng = *it;
        min_it = it;
      }
    }

    auto small = sorted.extract(sorted.begin()); /* smallest element count rank */

    auto& val1 = small.value();

    if (min_rng == nullptr) break; /* empty rank range set case */

    uint64_t count = (*min_rng)[2]-(*min_rng)[1];

    if (debug_print)
    {
      std::cout << "DBG: Largest per rank " << rng_type << " count: " << val0.first << ", smallest: " << val1.first;
      std::cout << ", transferred count: " << count << std::endl;
    }

    val0.first -= count; /* subtract from largest */
    val1.first += count; /* add to smallest */

    /* INFO: this approach is inherently safe in terms of the global array storage space, since the largest possible
       space was already uniformly allocated on all ranks, and the algorithm only makes the largest used space smaller */

    REAL im1 = (REAL)std::max(val0.first,val1.first)/(REAL)(std::min(val0.first,val1.first)+1);

    if (debug_print) std::cout << "DBG: migration algo im1 = " << im1 << std::endl;

    if (im1 >= im0) break; /* no more refinement possible */

    im0 = im1;
 
    r_b_rng[val0.second][min_bodnum].erase (min_it);
    send[val0.second].push_back(std::make_pair(min_bodnum,min_rng));
    /* RANK 0 --- */
    (*min_rng)[0] = val1.second; /* overwrite source rank with desitnation rank (via pointer with mesh_mapping[bodnum]) */
    /* --- RANK 0 */

    sorted.insert(std::move(val0)); /* reinsert modified pairs */
    sorted.insert(std::move(val1));
  }

  std::vector<uint64_t> in_sendbuf = sendbuf;
  std::vector<int> in_sendcounts = sendcounts;
  std::vector<int> in_displs = displs;

  sendbuf.clear();
  sendcounts.clear();
  displs.resize(1,0);

  for (auto it = send.begin(); it != send.end(); it ++)
  {
    auto rank0 = it-send.begin(); /* source range */

    size_t in_count = 0;

    if (append) /* append existing current rank data */
    {
      in_count = in_sendcounts[rank0];

      for (auto x = &in_sendbuf[in_displs[rank0]], y = &in_sendbuf[in_displs[rank0+1]]; x != y; x++)
      {
        sendbuf.push_back(*x);
      }
    }

    for (auto jt : *it)
    {
      auto bodnum = jt.first;
      auto rng = jt.second;
      auto rank1 = (*rng)[0]; /* destination range */

      sendbuf.push_back(bodnum);
      sendbuf.push_back(rank1);
      if (rng_ele) sendbuf.push_back(1); /* entity type: element */
      else sendbuf.push_back(2); /* entity type: element */

      sendbuf.push_back((*rng)[1]); /* source range */
      sendbuf.push_back((*rng)[2]);

      auto sz = (*rng)[2]-(*rng)[1];

      uint64_t dst0 = rng_ele ? elements_per_rank[rank1] : faces_per_rank[rank1];
      uint64_t dst1 = dst0 + sz;

      sendbuf.push_back(dst0); /* destination range */
      sendbuf.push_back(dst1);

      /* RANK 0 --- */
      (*rng)[1] = dst0; /* overwrite source range with destination range (via pointer with mesh_mapping[bodnum]) */
      (*rng)[2] = dst1;
      if (rng_ele)
      {
        elements_per_rank[rank0] -= sz;
        elements_per_rank[rank1] += sz;
      }
      else
      {
        faces_per_rank[rank0] -= sz;
        faces_per_rank[rank1] += sz;
      }
      /* --- RANK 0 */
    }

    sendcounts.push_back(7*it->size()+in_count);
    displs.push_back(displs.back()+7*it->size()+in_count);
  }
}

/* RANK 0: resize down compute::ga_elements or compute::ga_faces arrays
 *** in: ***
 * cn_name - cn_elements or cn_faces
 * ga_name - ga_elments or ga_faces
 * xx_last - el_last or fa_last
 * rank_deleted_ranges - std::vector<std::pair<uint64_t,uint64_t>> vector of entity ranges;
 * xx_bodnum - el_bodnum or fa_bodnum
 * xx_remap - mesh_mapping[bodnum] element/face ranges remapping flag */
auto GA_RESIZE_DOWN_RANK_MESH_ELE_FAC = [](int rank, auto cn_name, auto ga_name, auto xx_last, auto rank_deleted_ranges, auto xx_bodnum, bool xx_remap = true)
{
  using namespace compute;

  uint64_t count0, count1;
  ga_counters->get(rank, cn_name, cn_name+1, 0, 1, &count0);

  uint64_t *data = new uint64_t [count0 * xx_last];
  ga_name->get(rank, 0, count0, 0, xx_last, data);

  std::set<std::pair<uint64_t,uint64_t>> deleted_ranges;
  count1 = count0;
  for (auto& rng : rank_deleted_ranges)
  {
    count1 -= rng.second-rng.first;
    deleted_ranges.insert(rng);
  }

  std::set<std::pair<uint64_t,uint64_t>> left_ranges;
  uint64_t j = 0;
  for (auto& rng : deleted_ranges)
  {
    if (j < rng.first)
    {
      left_ranges.insert(std::pair<uint64_t,uint64_t>(j, rng.first));
    }

    j = rng.second;
  }
  if (j < count0)
  {
    left_ranges.insert(std::pair<uint64_t,uint64_t>(j, count0));
  }

  j = 0;
  for (uint64_t i = 0, *item = data; i < xx_last; i ++, item += count0)
  {
    for (auto & rng : left_ranges)
    {
      std::memmove (&data[j], &item[rng.first], sizeof(uint64_t [rng.second-rng.first]));
      j += rng.second-rng.first;
    }
  }

  ga_name->put(rank, 0, count1, 0, xx_last, data);
  ga_counters->put(rank, cn_name, cn_name+1, 0, 1, &count1);

  if (xx_remap) /* only used during rank-0-driven element/face deletion */
  {
    for (uint64_t *start = &data[count1*xx_bodnum], *item = start,
         *jtem = item, *end = start+count1; item < end; item = jtem)
    {
      uint64_t bodnum = *jtem;
      while (bodnum == *jtem && jtem < end) jtem ++;
      auto& xx_ranges = (ga_name == ga_elements ? mesh_mapping[bodnum].ga_eranges : mesh_mapping[bodnum].ga_franges);
      auto new_end = std::remove_if(xx_ranges.begin(), xx_ranges.end(), [rank](auto& rng) { return rng[0] == rank; });
      xx_ranges.erase(new_end, xx_ranges.end()); /* erase current rank ranges */
    }

    for (uint64_t *start = &data[count1*xx_bodnum], *item = start,
         *jtem = item, *end = start+count1; item < end; item = jtem)
    {
      uint64_t bodnum = *jtem;
      while (bodnum == *jtem && jtem < end) jtem ++;
      auto& xx_ranges = (ga_name == ga_elements ? mesh_mapping[bodnum].ga_eranges : mesh_mapping[bodnum].ga_franges);
      std::array<uint64_t,3> rng = {(unsigned)rank, (unsigned)(item-start), (unsigned)(jtem-start)}; /* remapped range */
      xx_ranges.push_back(rng); /* insert remapped range */
    }
  }

  delete [] data;
};

/* RANK 0: migrate mesh data by modifying global arrays */
static void migrate_meshes (std::vector<uint64_t> &sendbuf, std::vector<int> &sendcounts, std::vector<int> &displs)
{
  using namespace compute;

  for (auto jt = sendcounts.begin(); jt != sendcounts.end(); jt ++)
  {
    auto rank0 = jt - sendcounts.begin();

    std::vector<std::pair<uint64_t,uint64_t>> deleted_elements, deleted_faces;

    ASSERT_DEBUG ((displs[rank0+1]-displs[rank0]) % 7 == 0, "Mesh migratio data size is corrupted");
    for (auto it = &sendbuf[displs[rank0]]; it < &sendbuf[displs[rank0+1]]; it += 7)
    {
      uint64_t bodnum = it[0],
               rank1  = it[1],
               type   = it[2],
               org0   = it[3],
               org1   = it[4],
               dst0   = it[5],
               dst1   = it[6];
      uint64_t count = org1-org0;

      switch (type)
      {
      case 1:
      {
        std::vector<uint64_t> data(count*el_last);
        ga_elements->get(rank0, org0, org1, 0, el_last, &data[0]);
        ga_elements->put(rank1, dst0, dst1, 0, el_last, &data[0]);
        ga_counters->acc(rank1, cn_elements, cn_elements+1, 0, 1, &count);
        deleted_elements.push_back(std::make_pair(org0, org1));
      }
      break;
      case 2:
      {
        std::vector<uint64_t> data(count*fa_last);
        ga_faces->get(rank0, org0, org1, 0, fa_last, &data[0]);
        ga_faces->put(rank1, dst0, dst1, 0, fa_last, &data[0]);
        ga_counters->acc(rank1, cn_faces, cn_faces+1, 0, 1, &count);
        deleted_faces.push_back(std::make_pair(org0, org1));
      }
      break;
      }
    }

    GA_RESIZE_DOWN_RANK_MESH_ELE_FAC (rank0, cn_elements, ga_elements, el_last, deleted_elements, el_bodnum);
    GA_RESIZE_DOWN_RANK_MESH_ELE_FAC (rank0, cn_faces, ga_faces, fa_last, deleted_faces, fa_bodnum);
    /* TODO: nodes */
  }
}

/* COLLECTIVE: migrate mesh data by modifying global arrays
 *** in: ***
 * rank0 - current rank
 * data - array of 7-tuples of (bodnum, rank, type, org0, org1, dst0, dst1)
 *        as defined in mesh_migration_map() */
static void ALL_migrate_meshes (int rank0, const std::vector<uint64_t> &data)
{
  ASSERT_DEBUG (data.size() % 7 == 0, "Mesh migratio data size is corrupted");
  using namespace compute;

  std::vector<std::pair<uint64_t,uint64_t>> deleted_elements, deleted_faces;

  for (auto it = data.begin(); it != data.end(); it += 7)
  {
    uint64_t bodnum = it[0],
             rank1  = it[1],
             type   = it[2],
             org0   = it[3],
             org1   = it[4],
             dst0   = it[5],
             dst1   = it[6];
    uint64_t count = org1-org0;

    switch (type)
    {
    case 1:
    {
      std::vector<uint64_t> data(count*el_last);
      ga_elements->get(rank0, org0, org1, 0, el_last, &data[0]);
      ga_elements->put(rank1, dst0, dst1, 0, el_last, &data[0]);
      ga_counters->acc(rank1, cn_elements, cn_elements+1, 0, 1, &count);
      deleted_elements.push_back(std::make_pair(org0, org1));
    }
    break;
    case 2:
    {
      std::vector<uint64_t> data(count*fa_last);
      ga_faces->get(rank0, org0, org1, 0, fa_last, &data[0]);
      ga_faces->put(rank1, dst0, dst1, 0, fa_last, &data[0]);
      ga_counters->acc(rank1, cn_faces, cn_faces+1, 0, 1, &count);
      deleted_faces.push_back(std::make_pair(org0, org1));
    }
    break;
    }
  }

  ga_counters->fence();
  ga_elements->fence();
  ga_faces->fence();
  ga_nodes->fence();

  GA_RESIZE_DOWN_RANK_MESH_ELE_FAC (rank0, cn_elements, ga_elements, el_last, deleted_elements, el_bodnum, false);
  GA_RESIZE_DOWN_RANK_MESH_ELE_FAC (rank0, cn_faces, ga_faces, fa_last, deleted_faces, fa_bodnum, false);
  /* TODO: nodes */

  ga_counters->fence();
  ga_elements->fence();
  ga_faces->fence();
  ga_nodes->fence();
}

/* RANK 0: remap element/face ranges to reflect compactification undergone in GA_RESIZE_DOWN_RANK_MESH_ELE_FAC */
static void post_ALL_migrate_meshes (int size)
{
  using namespace compute;

  std::vector<std::set<std::array<uint64_t,3>>> rank_to_erange(size); /* rank-transposed mesh_mapping[*].ga_eranges */

  for (auto& [bodnum, mapping] : mesh_mapping)
  {
    for (auto& rng : mapping.ga_eranges)
    {
      std::array<uint64_t,3> r = {rng[1], rng[2], bodnum};
      rank_to_erange[rng[0]].insert(r);
    }

    mapping.ga_eranges.clear();
  }

  for (auto ir = rank_to_erange.begin(); ir != rank_to_erange.end(); ir ++)
  {
    std::vector<std::array<uint64_t,3>> ranges;

    auto rank = ir-rank_to_erange.begin();

    uint64_t pos = 0;

    if (debug_print) std::cout << "DBG: Unmodified ranges of rank " << rank << " are:";
    for (auto& rng : *ir)
    {
      if (debug_print) std::cout << " (" << rng[0] << "," << rng[1] << ")";
      std::array<uint64_t,3> r = {rng[2], pos, pos+(rng[1]-rng[0])};
      ranges.push_back(r);
      pos = r[2];
    }
    if (debug_print) std::cout << std::endl;

    if (debug_print) std::cout << "DBG: Modified ranges of rank " << rank << " are:";
    for (auto& rng : ranges)
    {
      if (debug_print) std::cout << " (" << rng[1] << "," << rng[2] << ")";
      std::array<uint64_t,3> r = {(unsigned)rank, rng[1], rng[2]};
      mesh_mapping[rng[0]].ga_eranges.push_back(r);
    }
    if (debug_print) std::cout << std::endl;
  }
}

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

    if (duration == -1.) break; /* rank 0 process communicated termination */

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
        std::vector<uint64_t> nindex(map.nrank.size(), 0); /* node index within MPI rank */

	struct part &part = parts[bodnum];

	struct mapping &mapping = mesh_mapping[bodnum]; /* create mesh range mapping */

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

	std::map<int,std::vector<std::pair<uint64_t, std::vector<uint64_t>>>> replace_nodes; /* replaced nodes on rank mapping;
                                      vectors of pairs of (replaced node range start index, vector of inserted node indices) */
	std::map<int,std::vector<uint64_t>> new_nodes; /* new nodes on rank mapping */

	for (auto r = map.nrank.begin(); r != map.nrank.end(); r++)
	{
          uint64_t i = r - map.nrank.begin();

          bool deleted_nodes_on_r = deleted_nodes.find(*r) != deleted_nodes.end();

          if (deleted_nodes_on_r && !deleted_nodes[*r].empty() > 0) /* use current range of deleted nodes */
          {
            auto &deleted_nodes_of_r = deleted_nodes[*r];

            if (auto it{ replace_nodes.find(*r) }; it != std::end(replace_nodes)) /* if mapped */
            {
              auto& vec = it->second;
              vec.back().second.push_back(i); /* add current node index */
            }
            else /* not yet mapped */
            {
              std::vector<uint64_t> vec(1, i);
              replace_nodes[*r].push_back(std::make_pair(
                 deleted_nodes_of_r.front().first, vec)); /* pair up initial resued node index
                                                      with a vector of inserted node indices */
            }

            deleted_nodes_of_r.front().first ++;

            if (deleted_nodes_of_r.front().first == deleted_nodes_of_r.front().second) /* no more space */
            {
              deleted_nodes_of_r.pop_front();

              if (!deleted_nodes_of_r.empty())
              {
                std::vector<uint64_t> vec; /* empty */
                replace_nodes[*r].push_back(std::make_pair(
                   deleted_nodes_of_r.front().first, vec)); /* pair up initial resued node index with a an empty vector */
              }
            }
          }
          else
          {
	    new_nodes[*r].push_back (i);
          }
	}

        for (auto& [r, vec0] : replace_nodes)
	{
          for (auto& [start, vec1] : vec0)
          {
            if (vec1.empty()) continue;
            struct mesh &mesh = solfec::meshes[bodnum];
            uint64_t nodsize = vec1.size();
            REAL *noddata = new REAL [nodsize * nd_last];
            uint64_t nodidx = 0;
            ERRMEM (noddata);

            for (auto& i  : vec1)
            {
              REAL x = mesh.nodes[0][i],
                   y = mesh.nodes[1][i],
                   z = mesh.nodes[2][i];
              REAL a[3] = {x - center[0],
                           y - center[1],
                           z - center[2]};
              REAL v[3] = {linear[0],
                           linear[1],
                           linear[2]};
              PRODUCTADD (angular, a, v);
              noddata[nd_vx*nodsize + nodidx] = v[0];
              noddata[nd_vy*nodsize + nodidx] = v[1];
              noddata[nd_vz*nodsize + nodidx] = v[2];
              noddata[nd_x*nodsize + nodidx] = x;
              noddata[nd_y*nodsize + nodidx] = y;
              noddata[nd_z*nodsize + nodidx] = z;
              noddata[nd_X*nodsize + nodidx] = x;
              noddata[nd_Y*nodsize + nodidx] = y;
              noddata[nd_Z*nodsize + nodidx] = z;
              noddata[nd_unused*nodsize + nodidx] = 0.;
              nindex[i] = start + nodidx;
              nodidx ++;
            }

            ga_nodes->put(r, start, start+nodidx, 0, nd_last, noddata);

            std::array<uint64_t,3> rng = {(unsigned)r, start, start+nodidx};
            mapping.ga_nranges.push_back(rng); /* handy in deletion code */
            nodes_per_rank[rng[0]] += rng[2]-rng[1];

	    deleted_nodes_count[r] -= nodidx; /* this number of deleted nodes is now reused */

            delete[] noddata;
          }
	}

	for (auto& [r, vec] : new_nodes)
	{
	  struct mesh &mesh = solfec::meshes[bodnum];
	  uint64_t nodsize = vec.size();
	  REAL *noddata = new REAL [nodsize * nd_last];
	  uint64_t nodidx = 0;
	  ERRMEM (noddata);

	  uint64_t count;
	  ga_counters->get(r, cn_nodes, cn_nodes+1, 0, 1, &count);

	  for (auto& i  : vec)
	  {
	    REAL x = mesh.nodes[0][i],
		 y = mesh.nodes[1][i],
		 z = mesh.nodes[2][i];
	    REAL a[3] = {x - center[0],
			 y - center[1],
			 z - center[2]};
	    REAL v[3] = {linear[0],
			 linear[1],
			 linear[2]};
	    PRODUCTADD (angular, a, v);
	    noddata[nd_vx*nodsize + nodidx] = v[0];
	    noddata[nd_vy*nodsize + nodidx] = v[1];
	    noddata[nd_vz*nodsize + nodidx] = v[2];
	    noddata[nd_x*nodsize + nodidx] = x;
	    noddata[nd_y*nodsize + nodidx] = y;
	    noddata[nd_z*nodsize + nodidx] = z;
	    noddata[nd_X*nodsize + nodidx] = x;
	    noddata[nd_Y*nodsize + nodidx] = y;
	    noddata[nd_Z*nodsize + nodidx] = z;
            noddata[nd_unused*nodsize + nodidx] = 0.;
	    nindex[i] = count + nodidx;
	    nodidx ++;
	  }

	  ga_nodes->put(r, count, count+nodidx, 0, nd_last, noddata);
	  ga_counters->acc(r, cn_nodes, cn_nodes+1, 0, 1, &nodidx);

	  std::array<uint64_t,3> rng = {(unsigned)r, count, count+nodidx};
	  mapping.ga_nranges.push_back(rng); /* handy in deletion code */
          nodes_per_rank[rng[0]] += rng[2]-rng[1];

	  delete[] noddata;
	}

	/* write elements */

	std::map<int,std::vector<uint64_t>> eonrank; /* elements on rank mapping */

	for (auto r = map.erank.begin(); r != map.erank.end(); r++)
	{
	  eonrank[*r].push_back (r - map.erank.begin());
	}

        for (auto& [r, vec] : eonrank)
	{
	  uint64_t elesize = vec.size();
	  uint64_t *eledata = new uint64_t [elesize * el_last];
	  uint64_t eleidx = 0;
	  ERRMEM (eledata);

	  for (auto& i : vec)
	  {
	    auto j = part.eptr[i];
	    auto eltype = part.eptr[i+1]-j;

	    eledata[el_bodnum*elesize + eleidx] = bodnum;
	    eledata[el_matnum*elesize + eleidx] = part.material[i];
	    eledata[el_type*elesize + eleidx] = eltype;

	    eledata[el_nd0_rnk*elesize + eleidx] = map.nrank[part.eind[j]];
	    eledata[el_nd0_idx*elesize + eleidx] = nindex[part.eind[j]];
	    eledata[el_nd1_rnk*elesize + eleidx] = map.nrank[part.eind[j+1]];
	    eledata[el_nd1_idx*elesize + eleidx] = nindex[part.eind[j+1]];
	    eledata[el_nd2_rnk*elesize + eleidx] = map.nrank[part.eind[j+2]];
	    eledata[el_nd2_idx*elesize + eleidx] = nindex[part.eind[j+2]];
	    eledata[el_nd3_rnk*elesize + eleidx] = map.nrank[part.eind[j+3]];
	    eledata[el_nd3_idx*elesize + eleidx] = nindex[part.eind[j+3]];
	    if (eltype > 4)
	    { eledata[el_nd4_rnk*elesize + eleidx] = map.nrank[part.eind[j+4]];
	      eledata[el_nd4_idx*elesize + eleidx] = nindex[part.eind[j+4]]; }
	    else
	    { eledata[el_nd4_rnk*elesize + eleidx] = 0;
	      eledata[el_nd4_idx*elesize + eleidx] = 0; }
	    if (eltype > 5)
	    { eledata[el_nd5_rnk*elesize + eleidx] = map.nrank[part.eind[j+5]];
	      eledata[el_nd5_idx*elesize + eleidx] = nindex[part.eind[j+5]]; }
	    else
	    { eledata[el_nd5_rnk*elesize + eleidx] = 0;
	      eledata[el_nd5_idx*elesize + eleidx] = 0; }
	    if (eltype > 7)
	    { eledata[el_nd6_rnk*elesize + eleidx] = map.nrank[part.eind[j+6]];
	      eledata[el_nd6_idx*elesize + eleidx] = nindex[part.eind[j+6]];
	      eledata[el_nd7_rnk*elesize + eleidx] = map.nrank[part.eind[j+7]];
	      eledata[el_nd7_idx*elesize + eleidx] = nindex[part.eind[j+7]]; }
	    else
	    { eledata[el_nd6_rnk*elesize + eleidx] = 0;
	      eledata[el_nd6_idx*elesize + eleidx] = 0;
	      eledata[el_nd7_rnk*elesize + eleidx] = 0;
	      eledata[el_nd7_idx*elesize + eleidx] = 0; }

	    eleidx ++;
	  }

	  uint64_t count;
	  ga_counters->get(r, cn_elements, cn_elements+1, 0, 1, &count);
	  ga_elements->put(r, count, count+eleidx, 0, el_last, eledata);
	  ga_counters->acc(r, cn_elements, cn_elements+1, 0, 1, &eleidx);

	  std::array<uint64_t,3> rng = {(unsigned)r, count, count+eleidx};
	  mapping.ga_eranges.push_back(rng); /* handy in deletion code */
          elements_per_rank[rng[0]] += rng[2]-rng[1];

	  delete[] eledata;
	}

	/* write faces */

	std::map<int,std::vector<uint64_t>> fonrank; /* faces on rank mapping */

	for (auto r = map.frank.begin(); r != map.frank.end(); r++)
	{
	  fonrank[*r].push_back (r - map.frank.begin());
	}

        for (auto& [r, vec] : fonrank)
	{
	  uint64_t facsize = vec.size();
	  uint64_t *facdata = new uint64_t [facsize * fa_last];
	  uint64_t facidx = 0;
	  ERRMEM (facdata);

	  for (auto& i : vec)
	  {
	    auto j = part.fptr[i];
	    auto factype = part.fptr[i+1]-j;

	    facdata[fa_bodnum*facsize + facidx] = bodnum;
	    facdata[fa_color*facsize + facidx] = part.color[i];
	    facdata[fa_type*facsize + facidx] = factype;

	    facdata[fa_nd0_rnk*facsize + facidx] = map.nrank[j];
	    facdata[fa_nd0_idx*facsize + facidx] = nindex[j];
	    facdata[fa_nd1_rnk*facsize + facidx] = map.nrank[j+1];
	    facdata[fa_nd1_idx*facsize + facidx] = nindex[j+1];
	    facdata[fa_nd2_rnk*facsize + facidx] = map.nrank[j+2];
	    facdata[fa_nd2_idx*facsize + facidx] = nindex[j+2];
	    if (factype > 3)
	    { facdata[fa_nd3_rnk*facsize + facidx] = map.nrank[j+3];
	      facdata[fa_nd3_idx*facsize + facidx] = nindex[j+3]; }
	    else
	    { facdata[fa_nd3_rnk*facsize + facidx] = 0;
	      facdata[fa_nd3_idx*facsize + facidx] = 0; }

	    facidx ++;
	  }

	  uint64_t count;
	  ga_counters->get(r, cn_faces, cn_faces+1, 0, 1, &count);
	  ga_faces->put(r, count, count+facidx, 0, fa_last, facdata);
	  ga_counters->acc(r, cn_faces, cn_faces+1, 0, 1, &facidx);

	  std::array<uint64_t,3> rng = {(unsigned)r, count, count+facidx};
	  mapping.ga_franges.push_back(rng); /* handy in deletion code */
          faces_per_rank[rng[0]] += rng[2]-rng[1];

	  delete[] facdata;
	}
      }

      inserted_meshes.clear();
    };

    auto GA_INSERT_ELLIPS = [](auto size)
    {
      std::set<std::pair<uint64_t, int>> equeue;
      for (auto ec = ellips_per_rank.begin(); ec != ellips_per_rank.end(); ec++)
      {
        equeue.insert(std::make_pair(*ec, ec - ellips_per_rank.begin()));
      }

      std::map<int,std::vector<uint64_t>> ellips_on_rank;
      for (auto& bodnum : inserted_ellips)
      {
        auto nh = equeue.extract(equeue.begin()); /* smallest size rank */
        int rank = nh.value().second; /* insert to this rank */
        nh.value().first ++; /* update rank size */
        equeue.insert(std::move(nh)); /* update queue */

	ellip_mapping[bodnum] = rank;
        ellips_on_rank[rank].push_back(bodnum);
        ellips_per_rank[rank] ++;
      }

      for (auto& [r, vec] : ellips_on_rank)
      {
        uint64_t ellsize = vec.size();
        REAL *elldata = new REAL [ellsize * ll_last0];
        uint64_t *ellips = new uint64_t [ellsize * ll_last1];
        uint64_t ellidx = 0;
        ERRMEM (elldata);
        ERRMEM (ellips);

        for (auto& bodnum : vec)
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

          ellidx ++;
        }

        uint64_t count;
        ga_counters->get(r, cn_ellips, cn_ellips+1, 0, 1, &count);
	ga_elldata->put(r, count, count+ellidx, 0, ll_last0, elldata);
	ga_ellips->put(r, count, count+ellidx, 0, ll_last1, ellips);
	ga_counters->acc(r, cn_ellips, cn_ellips+1, 0, 1, &ellidx);

        delete[] elldata;
        delete[] ellips;
      }

      inserted_ellips.clear();
    };

    if (!partitioned)
    {
      nodes_per_rank.resize(size);
      elements_per_rank.resize(size);
      faces_per_rank.resize(size);
      ellips_per_rank.resize(size);

      deleted_nodes_count.resize(size);

      std::fill(deleted_nodes_count.begin(), deleted_nodes_count.end(), 0);

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

	std::map<uint64_t, rankmap> maps = map_parts (parts, nodes_per_rank, elements_per_rank, faces_per_rank);

	auto [maxnodes, maxeles, maxfaces] = max_per_rank (maps);

	for (int r = 0; r < size; r ++) /* initialize counters */
	{
	  uint64_t counts[] = {0, inserted_materials.size()*2, 0,
			       0, maxnodes*2, 0,
			       0, maxeles*2, 0,
			       0, maxfaces*2, 0,
			       0, (inserted_ellips.size()/size+1)*2, 0};

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
          ERRMEM (ga);
	  ga_materials->get(rank, 0, counts[cn_materials], 0, mt_last, data);
	  ga->put(rank, 0, counts[cn_materials], 0, mt_last, data);
	  GA *tmp = ga_materials;
	  ga_materials = ga;
	  delete tmp;

          counts[sz_materials] = counts[sz_materials_new];
	}

        if (counts[sz_nodes_new] > counts[sz_nodes])
	{
	  GA *ga = new GA(MPI_COMM_WORLD, counts[sz_nodes_new], nd_last, MPI_REAL);
	  REAL *data = new REAL [counts[cn_nodes] * nd_last];
	  ERRMEM (data);
          ERRMEM (ga);
	  ga_nodes->get(rank, 0, counts[cn_nodes], 0, nd_last, data);
	  ga->put(rank, 0, counts[cn_nodes], 0, nd_last, data);
	  GA *tmp = ga_nodes;
	  ga_nodes = ga;
	  delete tmp;

          counts[sz_nodes] = counts[sz_nodes_new];
	}

	if (counts[sz_elements_new] > counts[sz_elements])
	{
	  GA *ga = new GA(MPI_COMM_WORLD, counts[sz_elements_new], el_last, MPI_UINT64_T);
	  uint64_t *data = new uint64_t [counts[cn_elements] * el_last];
	  ERRMEM (data);
          ERRMEM (ga);
	  ga_elements->get(rank, 0, counts[cn_elements], 0, el_last, data);
	  ga->put(rank, 0, counts[cn_elements], 0, el_last, data);
	  GA *tmp = ga_elements;
	  ga_elements = ga;
	  delete tmp;

          counts[sz_elements] = counts[sz_elements_new];
	}

	if (counts[sz_faces_new] > counts[sz_faces])
	{
	  GA *ga = new GA(MPI_COMM_WORLD, counts[sz_faces_new], fa_last, MPI_UINT64_T);
	  uint64_t *data = new uint64_t [counts[cn_faces] * fa_last];
	  ERRMEM (data);
          ERRMEM (ga);
	  ga_faces->get(rank, 0, counts[cn_faces], 0, fa_last, data);
	  ga->put(rank, 0, counts[cn_faces], 0, fa_last, data);
	  GA *tmp = ga_faces;
	  ga_faces = ga;
	  delete tmp;

	  counts[sz_faces] = counts[sz_faces_new];
	}

	if (counts[sz_ellips_new] > counts[sz_ellips])
	{
	  GA *ga0 = new GA(MPI_COMM_WORLD, counts[sz_ellips_new], ll_last0, MPI_REAL);
	  GA *ga1 = new GA(MPI_COMM_WORLD, counts[sz_ellips_new], ll_last1, MPI_UINT64_T);
	  REAL *data0 = new REAL [counts[cn_ellips] * ll_last0];
	  uint64_t *data1 = new uint64_t [counts[cn_ellips] * ll_last1];
	  ERRMEM (data0);
	  ERRMEM (data1);
          ERRMEM (ga0);
          ERRMEM (ga1);
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

	  counts[sz_ellips] = counts[sz_ellips_new];
	}

	ga_counters->put(rank, 0, cn_last, 0, 1, counts);

	ga_counters->fence();
      };
	
      if (rank == 0)
      {
        auto GA_RESIZE_DOWN = [](auto cn_name, auto ga_name, auto xx_last, auto rank_deleted_ranges, auto xx_bodnum)
        {
          for (auto& [rnk, vec] : rank_deleted_ranges)
          {
            GA_RESIZE_DOWN_RANK_MESH_ELE_FAC (rnk, cn_name, ga_name, xx_last, vec, xx_bodnum);
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
              nodes_per_rank[r[0]] -= r[2]-r[1];
              std::vector<REAL> unused(r[2]-r[1],1.);
              ga_nodes->put(r[0], r[1], r[2], nd_unused, nd_unused+1, &unused[0]);
	    }

	    for (auto& r : map.ga_eranges)
	    {
	      deleted_elements[r[0]].push_back(std::make_pair(r[1],r[2]));
              elements_per_rank[r[0]] -= r[2]-r[1];
	    }

	    for (auto& r : map.ga_franges)
	    {
	      deleted_faces[r[0]].push_back(std::make_pair(r[1],r[2]));
              faces_per_rank[r[0]] -= r[2]-r[1];
	    }
	  }

          for (auto& [r, list] : deleted_nodes)
          {
            list.sort(); /* sort deleted ranges list per rank */
          }

	  /* resize element and face arrays */
	  GA_RESIZE_DOWN (cn_elements, ga_elements, el_last, deleted_elements, el_bodnum);
	  GA_RESIZE_DOWN (cn_faces, ga_faces, fa_last, deleted_faces, fa_bodnum);

          deleted_meshes.clear();
	}

	if (!deleted_ellips.empty())
	{
	  std::unordered_map<int, std::set<uint64_t>> deleted_ellips_onrank;

	  for (auto& bodnum : deleted_ellips)
	  {
            int rank = ellip_mapping[bodnum];
	    deleted_ellips_onrank[rank].insert(bodnum);
            ellip_mapping.erase(bodnum);
            ellips_per_rank[rank] --;
	  }

	  for (auto& [r, set] : deleted_ellips_onrank)
	  {
	    uint64_t count0, count1;
	    ga_counters->get(r, cn_ellips, cn_ellips+1, 0, 1, &count0);

	    REAL *data0 = new REAL [count0 * ll_last0];
	    uint64_t *data1 = new uint64_t [count0 * ll_last1];
	    ga_elldata->get(r, 0, count0, 0, ll_last0, data0);
	    ga_ellips->get(r, 0, count0, 0, ll_last1, data1);

	    std::vector<uint64_t> left_indices;

            for (uint64_t i = 0; i < count0; i ++)
            {
              uint64_t bodnum = data1[count0*ll_bodnum+i];

              if (set.find(bodnum) == set.end())
              {
                left_indices.push_back(i);
              }
            }

            count1 = left_indices.size();

            uint64_t i = 0, j = 0;
	    for (REAL *item = data0; i < ll_last0; i ++, item += count0)
	    {
              for (auto & k : left_indices)
              {
                data0[j] = item[k];
                j ++;
              }
	    }

            i = 0, j = 0;
	    for (uint64_t *item = data1; i < ll_last1; i ++, item += count0)
	    {
              for (auto & k : left_indices)
              {
                data1[j] = item[k];
                j ++;
              }
	    }

	    ga_elldata->put(r, 0, count1, 0, ll_last0, data0);
	    ga_ellips->put(r, 0, count1, 0, ll_last1, data1);
	    ga_counters->put(r, cn_ellips, cn_ellips+1, 0, 1, &count1);

	    delete [] data0;
	    delete [] data1;
	  }

          deleted_ellips.clear();
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

	std::map<uint64_t, rankmap> maps = map_parts (parts, nodes_per_rank, elements_per_rank, faces_per_rank);

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
	  counts[sz_ellips_new] = counts[cn_ellips] + (inserted_ellips.size()/size+1);

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

    /* rebalance data if required */

    REAL imb[2]; /* elements and ellipsoids imbalances */

    if (rank == 0)
    {
      auto x = std::minmax_element (elements_per_rank.begin(), elements_per_rank.end());
      auto y = std::minmax_element (ellips_per_rank.begin(), ellips_per_rank.end());
      imb[0] = *x.first == 0 ? std::numeric_limits<REAL>::infinity() : (REAL)*x.second / (REAL)(*x.first);
      imb[1] = *y.first == 0 ? std::numeric_limits<REAL>::infinity() : (REAL)*y.second / (REAL)(*y.first);

      if (debug_print)
      {
        for (auto it = elements_per_rank.begin(); it != elements_per_rank.end(); it ++)
        {
          std::cout << "DBG: pre-migration --> elements on rank " << it-elements_per_rank.begin() << " = " << *it << std::endl;
        }
        for (auto it = ellips_per_rank.begin(); it != ellips_per_rank.end(); it ++)
        {
          std::cout << "DBG: pre-migration --> ellipsoids on rank " << it-ellips_per_rank.begin() << " = " << *it << std::endl;
        }

        std::cout << "DBG: pre-migration --> ELEMENTs imbalance: " << FMT("%.2f") << imb[0] << std::endl;
        std::cout << "DBG: pre-migration --> ELLIPs imbalance: " << FMT("%.2f") << imb[1] << std::endl;
      }
    }

    MPI_Bcast (imb, 2, MPI_REAL, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
      if (imb[0] > IMBALANCE_TOLERANCE)
      {
        std::vector<std::map<uint64_t,std::vector<std::array<uint64_t,3>*>>> eranges_on_ranks(size),
                                                                             franges_on_ranks(size);

        for (auto& [bodnum, mapping] : mesh_mapping)
        {
          for (auto& rng : mapping.ga_eranges)
          {
            eranges_on_ranks[rng[0]][bodnum].push_back(&rng); /* collect element ranges per mesh per rank */
          }

          for (auto& rng : mapping.ga_franges)
          {
            franges_on_ranks[rng[0]][bodnum].push_back(&rng); /* collect face ranges per mesh per rank */
          }
        }

        /* pre-migration bodnums */
        if (debug_print)
        {
          std::cout << "DBG: pre-migration element range bodnums:" << std::endl;
          for (auto it = eranges_on_ranks.begin(); it != eranges_on_ranks.end(); it ++)
          {
            std::cout << "DBG:  on rank " << (it-eranges_on_ranks.begin()) << ":";
            for (auto& [bodnum, v] : *it)
            {
              std::cout << " " << bodnum;
            }
            std::cout << std::endl;
          }

          std::cout << "DBG: pre-migration face range bodnums:" << std::endl;
          for (auto it = franges_on_ranks.begin(); it != franges_on_ranks.end(); it ++)
          {
            std::cout << "DBG:  on rank " << (it-franges_on_ranks.begin()) << ":";
            for (auto& [bodnum, v] : *it)
            {
              std::cout << " " << bodnum;
            }
            std::cout << std::endl;
          }
        }

        /* calculate migration mapping based on eranges_on_ranks, franges_on_ranks */
        std::vector<uint64_t> sendbuf;
        std::vector<int> sendcounts;
        std::vector<int> displs;
        mesh_migration_map (eranges_on_ranks, "elements", sendbuf, sendcounts, displs);
        mesh_migration_map (franges_on_ranks, "faces", sendbuf, sendcounts, displs, true);

        if (debug_print)
        {
          for (int i = 0; i < size; i ++)
          {
            std::cout << "DBG: Migrating " << sendcounts[i]/7 << " meshes from rank " << i << std::endl;
            for (int j = displs[i]; j < displs[i+1]; j += 7)
            {
              std::cout << "DBG:    bodnum " << sendbuf[j] << " to rank " << sendbuf[j+1] << std::endl;
            }
          }
        }

        uint64_t etot = std::accumulate(elements_per_rank.begin(), elements_per_rank.end(), 0);
        uint64_t esent = 0;
        for (auto it = &sendbuf[0], end = &sendbuf[displs[size]]; it != end; it += 7) { if (it[2] == 1) esent += it[4]-it[3]; }

        REAL elemigratio = (REAL)esent/(REAL)(etot+1);
        MPI_Bcast (&elemigratio, 1, MPI_REAL, 0, MPI_COMM_WORLD);

        if (debug_print) std::cout << "DBG: migration (rank-0/ranks-all) approach elemigratio: " << elemigratio << "; and threshold: " << ELEMIGRATIO_THRESHOLD << std::endl;

        if (elemigratio > ELEMIGRATIO_THRESHOLD) /* all ranks GA read-write access */
        {
          /* broadcast migration mapping */
          int recvcount;
          ASSERT_DEBUG (sendcounts.size() == size, "Inconsistent sendcounts size");
          MPI_Scatter (&sendcounts[0], 1, MPI_INT, &recvcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
          std::vector<uint64_t> recvbuf(recvcount);
          MPI_Scatterv (&sendbuf[0], &sendcounts[0], &displs[0], MPI_UINT64_T, &recvbuf[0], recvcount, MPI_UINT64_T, 0, MPI_COMM_WORLD);

          /* migrate meshes */
          ALL_migrate_meshes (rank, recvbuf);

          post_ALL_migrate_meshes (size);
        }
        else /* rank-0 GA read-write access */
        {
          migrate_meshes (sendbuf, sendcounts, displs);

          ga_counters->fence();
          ga_elements->fence();
          ga_faces->fence();
          ga_nodes->fence();
        }

        /* post-migration bodnums */
        if (debug_print)
        {
          std::vector<std::set<uint64_t>> bodnums_on_ranks(size);

          for (auto& [bodnum, mapping] : mesh_mapping)
          {
            for (auto& rng : mapping.ga_eranges)
            {
              bodnums_on_ranks[rng[0]].insert(bodnum);
            }
          }

          std::cout << "DBG: post-migration bodnums:" << std::endl;
          for (auto it = bodnums_on_ranks.begin(); it != bodnums_on_ranks.end(); it ++)
          {
            std::cout << "DBG:  on rank " << (it-bodnums_on_ranks.begin()) << ":";
            for (auto bodnum : *it)
            {
              std::cout << " " << bodnum;
            }
            std::cout << std::endl;
          }
        }
      }

      if (imb[1] > IMBALANCE_TOLERANCE)
      {
        auto sum = std::accumulate(ellips_per_rank.begin(), ellips_per_rank.end(), 0);
        auto avg = sum/size;

        std::vector<std::vector<uint64_t>> ellips_on_ranks(size);

        for (auto& [bodnum, rank] : ellip_mapping)
        {
          ellips_on_ranks[rank].push_back(bodnum);
        }

        /* TODO: calculate migration mapping based on (avg, ellips_on_ranks) */

        /* TODO: broadcast migration mapping */

        /* TODO: migrate ellipsoids */
      }
    }
    else /* rank > 0 */
    {
      if (imb[0] > IMBALANCE_TOLERANCE)
      {
        REAL elemigratio;
        MPI_Bcast (&elemigratio, 1, MPI_REAL, 0, MPI_COMM_WORLD);

        if (elemigratio > ELEMIGRATIO_THRESHOLD) /* all ranks GA read-write access */
        {
          /* broadcast migration mapping */
          int recvcount;
          MPI_Scatter (NULL, 0, MPI_INT, &recvcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
          std::vector<uint64_t> recvbuf(recvcount);
          MPI_Scatterv (NULL, NULL, NULL, MPI_UINT64_T, &recvbuf[0], recvcount, MPI_UINT64_T, 0, MPI_COMM_WORLD);

          /* migrate meshes */
          ALL_migrate_meshes (rank, recvbuf);
        }
        else /* rank-0 GA read-write access */
        {
          ga_counters->fence();
          ga_elements->fence();
          ga_faces->fence();
          ga_nodes->fence();
        }
      }

      if (imb[1] > IMBALANCE_TOLERANCE)
      {
        /* TODO: broadcast migration mapping */

        /* TODO: migrate ellipsoids */
      }
    }

    if (rank == 0 && debug_print)
    {
      for (auto it = nodes_per_rank.begin(); it != nodes_per_rank.end(); it ++)
      {
        std::cout << "DBG: nodes on rank " << it-nodes_per_rank.begin() << " = " << *it << std::endl;
      }
      for (auto it = elements_per_rank.begin(); it != elements_per_rank.end(); it ++)
      {
        std::cout << "DBG: elements on rank " << it-elements_per_rank.begin() << " = " << *it << std::endl;
      }
      for (auto it = faces_per_rank.begin(); it != faces_per_rank.end(); it ++)
      {
        std::cout << "DBG: faces on rank " << it-faces_per_rank.begin() << " = " << *it << std::endl;
      }
      for (auto it = ellips_per_rank.begin(); it != ellips_per_rank.end(); it ++)
      {
        std::cout << "DBG: ellipsoids on rank " << it-ellips_per_rank.begin() << " = " << *it << std::endl;
      }

      auto z = std::minmax_element (nodes_per_rank.begin(), nodes_per_rank.end());
      REAL im = *z.first == 0 ? std::numeric_limits<REAL>::infinity() : (REAL)*z.second / (REAL)(*z.first);
      std::cout << "DBG: NODES imbalance: " << FMT("%.2f") << im << std::endl;

      z = std::minmax_element (elements_per_rank.begin(), elements_per_rank.end());
      im = *z.first == 0 ? std::numeric_limits<REAL>::infinity() : (REAL)*z.second / (REAL)(*z.first);
      std::cout << "DBG: ELEMENTS imbalance: " << FMT("%.2f") << im << std::endl;

      z = std::minmax_element (faces_per_rank.begin(), faces_per_rank.end());
      im = *z.first == 0 ? std::numeric_limits<REAL>::infinity() : (REAL)*z.second / (REAL)(*z.first);
      std::cout << "DBG: FACES imbalance: " << FMT("%.2f") << im << std::endl;

      z = std::minmax_element (ellips_per_rank.begin(), ellips_per_rank.end());
      im = *z.first == 0 ? std::numeric_limits<REAL>::infinity() : (REAL)*z.second / (REAL)(*z.first);
      std::cout << "DBG: ELLIPSOIDS imbalance: " << FMT("%.2f") << im << std::endl;
    }

    if (debug_files)
    {
      debug_output_compute_data();
    }

    /* TODO: create compute task graph */

    if (solfec::simulation_time == 0. || duration == 0.) output_current_results(); /* all ranks */

    for (REAL end = solfec::simulation_time + duration;
	 solfec::simulation_time < end; solfec::simulation_time += step)
    {
      /* TODO: execute compute task graph */

      if (0 /* all ranks -> time to write output files or histories */)
      {
        ga_nodes->fence(); /* XXX: may not be needed (depending on compute graph) */
      }

      if (0 /* all ranks -> time to write output files */)
      {
        output_current_results();
      }

      if (0 /* time to write histories */)
      {
        /* TODO */
      }
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

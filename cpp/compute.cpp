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
#include <set>
#include "real.h"
#include "part.hpp"
#include "solfec.hpp"
#include "compute.hpp"

/* compute gobal variables */
namespace compute
{
/* rank 0 --- */
bool partitioned = false; /* initially paritioned */

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
std::map<uint64_t,material> materials; /* local MPI rank copy of all nominal materials */

MPI_Win mpi_nodes_window; /* storing nodes and dgreees of freedom data */
MPI_Win mpi_elements_window; /* storing elements data */
MPI_Win mpi_faces_window; /* storing faces data */
/* --- all ranks */
};

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
  int rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

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
      std::map<uint64_t, part> parts = partition_meshes(inserted_meshes);

      std::map<uint64_t, mapping> maps = map_parts (parts);

      auto [maxnodes, maxeles, maxfaces] = max_per_rank (maps);

      /* TODO: using max and bunch sizes create node, element, face windows */

      /* TODO: start putting node data to remote ranks */

      /* TODO: put element and  face data to ranks > 0, overallping
               communication and creation of bunch structures if possible */

      partitioned = true;
    }
  }
}

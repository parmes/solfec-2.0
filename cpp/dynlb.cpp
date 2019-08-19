/*
The MIT License (MIT)

Copyright (c) 2016 EDF Energy

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include "real.h"
#include "err.h"
#include "alg.h"
#include "dynlb.hpp"

/* create local load balancer */
struct dynlb* dynlb_local_create (uint64_t n, REAL *point[3], int64_t cutoff, REAL epsilon)
{
  using namespace ispc;
  int i, size, *rank_size;
  partitioning *ptree;
  struct dynlb *lb;
  REAL *gpoint[3];
  uint64_t leaf_count;

  MPI_Comm_size (MPI_COMM_WORLD, &size);

  ERRMEM (lb = (struct dynlb*)malloc (sizeof(struct dynlb)));
  lb->cutoff = cutoff;
  lb->epsilon = epsilon;

  if (cutoff <= 0)
  {
    cutoff = -size; /* as many leaves as ranks */
  }

  ptree = partitioning_create (0, n, point, cutoff, &lb->ptree_size, &leaf_count);

  partitioning_assign_ranks (ptree, leaf_count / size, leaf_count % size);

  partitioning_store (0, ptree, n, point);

  /* determine initial imbalance */

  ERRMEM (rank_size = (int*)calloc (size, sizeof (int)));

  for (i = 0; i < lb->ptree_size; i ++)
  {
    if (ptree[i].dimension < 0) /* leaf */
    {
      rank_size[ptree[i].rank] += ptree[i].size;
    }
  }

  int min_size = INT_MAX, max_size = 0;

  for (i = 0; i < size; i ++)
  {
    min_size = MIN (min_size, rank_size[i]);
    max_size = MAX (max_size, rank_size[i]);
  }

  lb->imbalance = (REAL)max_size/(REAL)min_size;

  if (isnan(lb->imbalance)) lb->imbalance = (REAL)1/(REAL)0; /* inf istead */

  lb->ptree = ptree;

  free (rank_size);

  return lb;
}

/* assign an MPI rank to a point; return this rank */
int dynlb_point_assign (struct dynlb *lb, REAL point[])
{
  return ispc::partitioning_point_assign (lb->ptree, 0, point);
}

/* assign MPI ranks to a box spanned between lo and hi points; return the number of ranks assigned */
int dynlb_box_assign (struct dynlb *lb, REAL lo[], REAL hi[], int ranks[])
{
  int count = 0;

  ispc::partitioning_box_assign (lb->ptree, 0, lo, hi, ranks, &count);

  return count;
}

/* update local load balancer */
void dynlb_local_update (struct dynlb *lb, uint64_t n, REAL *point[3])
{
  using namespace ispc;
  partitioning *ptree = lb->ptree;
  int i, size, *rank_size;

  MPI_Comm_size (MPI_COMM_WORLD, &size);

  partitioning_store (0, ptree, n, point);

  ERRMEM (rank_size = (int*)calloc (size, sizeof (int)));

  for (i = 0; i < lb->ptree_size; i ++)
  {
    if (ptree[i].dimension < 0) /* leaf */
    {
      rank_size[ptree[i].rank] += ptree[i].size;
    }
  }

  int min_size = INT_MAX, max_size = 0;

  for (i = 0; i < size; i ++)
  {
    min_size = MIN (min_size, rank_size[i]);
    max_size = MAX (max_size, rank_size[i]);
  }

  lb->imbalance = (REAL)max_size/(REAL)min_size;

  free (rank_size);

  if (isnan (lb->imbalance) || isinf(lb->imbalance) || lb->imbalance > 1.0 + lb->epsilon) /* update partitioning */
  {
    struct dynlb *dy = dynlb_local_create (n, point, lb->cutoff, lb->epsilon);

    partitioning_destroy (lb->ptree);
    lb->ptree = dy->ptree;
    lb->ptree_size = dy->ptree_size;

    lb->imbalance = dy->imbalance;

    free (dy);
  }
}

/* destroy load balancer */
void dynlb_destroy (struct dynlb *lb)
{
  ispc::partitioning_destroy (lb->ptree);
  free (lb);
}

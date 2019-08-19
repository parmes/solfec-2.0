/*
The MIT License (MIT)

Copyright (c) 2015 Tomasz Koziara

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

#include "real.h"
#include "rcb.h"

/* count leaves */
static void count_leaves (uniform struct partitioning ptree[], uniform int node, uniform unsigned int64 * uniform leaf_count)
{
  if (ptree[node].dimension >= 0) /* node */
  {
    count_leaves (ptree, ptree[node].left, leaf_count);
    count_leaves (ptree, ptree[node].right, leaf_count);
  }
  else /* leaf */
  {
    ptree[node].size = 0;
    (*leaf_count) ++;
  }
}

/* zero leaf data */
static void zero_leaves (uniform struct partitioning ptree[], uniform int node)
{
  if (ptree[node].dimension >= 0) /* node */
  {
    zero_leaves (ptree, ptree[node].left);
    zero_leaves (ptree, ptree[node].right);
  }
  else /* leaf */
  {
    ptree[node].size = 0;
  }
}

/* drop point down the partitioning tree */
static void drop_point (uniform partitioning ptree[], uniform int node, uniform unsigned int64 i, uniform REAL * uniform point[3])
{
  uniform int d = ptree[node].dimension;

  if (d >= 0) /* node */
  {
    if (point[d][i] < ptree[node].coord)
      drop_point (ptree, ptree[node].left, i, point);
    else drop_point (ptree, ptree[node].right, i, point);
  }
  else /* leaf */
  {
    atomic_add_global (&ptree[node].size, 1);
  }
}

/* store points at tree leaves */
task void store_points (uniform unsigned int64 span, uniform partitioning ptree[], uniform unsigned int64 n, uniform REAL * uniform point[3])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? n: start+span;

  for (uniform int i = start; i < end; i ++)
  {
    drop_point (ptree, 0, i, point);
  }
}

/* assign ranks to partitioning tree leaves */
static void assign_ranks (uniform partitioning * uniform ptree, uniform int node,
  uniform unsigned int64 leaves_per_rank, uniform unsigned int64 * uniform remainder,
  uniform unsigned int64 * uniform leaf, uniform int * uniform rank)
{
  if (ptree[node].dimension >= 0) /* node */
  {
    assign_ranks (ptree, ptree[node].left, leaves_per_rank, remainder, leaf, rank);
    assign_ranks (ptree, ptree[node].right, leaves_per_rank, remainder, leaf, rank);
  }
  else /* leaf */
  {
    ptree[node].rank =  (*rank);

    (*leaf) ++;

    uniform unsigned int64 m = leaves_per_rank + ((*remainder) ? 1 : 0); /* remainder is distributed into initial ranks */

    if ((*leaf) == m) /* finished with this rank */
    {
      (*leaf) = 0; /* zero leaf per rank counter */

      (*rank) ++; /* increment rank */

      if ((*remainder)) (*remainder) --; /* if remainder left then decrement it */
    }
  }
}

/* create partitioning tree based on rcb tree;
 * uniformly bisect untill leaf size <= cutoff; or if cutoff < 0 then create -cutoff equal size leaves */
export uniform partitioning * uniform partitioning_create (uniform int ntasks, uniform unsigned int64 n, uniform REAL * uniform point[3],
  uniform int64 cutoff, uniform int * uniform tree_size, uniform unsigned int64 * uniform leaf_count)
{
  uniform partitioning * uniform ptree = rcb_tree_create (ntasks, n, point, cutoff, tree_size);

  *leaf_count = 0;

  count_leaves (ptree, 0, leaf_count);

  return ptree;
}

/* allocate partitioning tree memory */
export uniform partitioning * uniform partitioning_alloc (uniform int tree_size)
{
  uniform partitioning * uniform ptree = uniform new uniform partitioning[tree_size];

  return ptree;
}

/* assign ranks to partitioning tree leaves */
export void partitioning_assign_ranks (uniform partitioning * uniform ptree, uniform unsigned int64 leaves_per_rank, uniform unsigned int64 remainder)
{
  uniform unsigned int64 leaf = 0;

  uniform int rank = 0;

  assign_ranks (ptree, 0, leaves_per_rank, &remainder, &leaf, &rank);
}

/* store points in the partitioning tree leaves */
export void partitioning_store (uniform int ntasks, uniform partitioning * uniform ptree, uniform unsigned int64 n, uniform REAL * uniform point[3])
{
  uniform int num = ntasks < 1 ? num_cores () : ntasks;

  zero_leaves (ptree, 0);

  launch [num] store_points (n/num, ptree, n, point);
}

/* assign leaf rank to a point */
export uniform int partitioning_point_assign (uniform partitioning ptree[], uniform int node, uniform REAL point[])
{
  uniform int d = ptree[node].dimension;

  if (d >= 0) /* node */
  {
    if (point[d] < ptree[node].coord)
      return partitioning_point_assign (ptree, ptree[node].left, point);
    else return partitioning_point_assign (ptree, ptree[node].right, point);
  }
  else /* leaf */
  {
    return ptree[node].rank;
  }
}

/* assign leaf ranks to a box */
export void partitioning_box_assign (uniform partitioning ptree[], uniform int node,
  uniform REAL lo[], uniform REAL hi[], uniform int ranks[], uniform int * uniform rank_count)
{
  uniform int d = ptree[node].dimension;

  if (d >= 0) /* node */
  {
    if (lo[d] < ptree[node].coord)
      partitioning_box_assign (ptree, ptree[node].left, lo, hi, ranks, rank_count);
    if (hi[d] > ptree[node].coord)
      partitioning_box_assign (ptree, ptree[node].right, lo, hi, ranks, rank_count);
  }
  else /* leaf */
  {
    uniform int i, r = ptree[node].rank;
    
    for (i = 0; i < (*rank_count); i ++)
    {
      if (ranks[i] == r) break;
    }

    if (i == (*rank_count))
    {
      ranks[i] = r;
      (*rank_count) ++;
    }
  }
}

/* destroy partitioning tree */
export void partitioning_destroy (uniform partitioning * uniform ptree)
{
  delete ptree;
}

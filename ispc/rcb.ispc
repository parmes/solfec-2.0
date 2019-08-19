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

#include "real.h"
#include "alg.h"
#include "rcb.h"

/* calculate extrema of x, y, z */
task void extrema (uniform unsigned int64 span, uniform unsigned int64 n, uniform REAL x[], uniform REAL y[], uniform REAL z[], uniform REAL extents[])
{
  uniform unsigned int64 start = taskIndex*span;
  uniform unsigned int64 end = taskIndex == taskCount-1 ? n : start+span;

  REAL e[6] = {REAL_MAX,REAL_MAX,REAL_MAX,-REAL_MAX,-REAL_MAX,-REAL_MAX};

  foreach (i = start ... end)
  {
    #pragma ignore warning
    if (x[i] < e[0]) e[0] = x[i];
    #pragma ignore warning
    if (y[i] < e[1]) e[1] = y[i];
    #pragma ignore warning
    if (z[i] < e[2]) e[2] = z[i];
    #pragma ignore warning
    if (x[i] > e[3]) e[3] = x[i];
    #pragma ignore warning
    if (y[i] > e[4]) e[4] = y[i];
    #pragma ignore warning
    if (z[i] > e[5]) e[5] = z[i];
  }

  uniform REAL * uniform out = &extents [6*taskIndex];

  out[0] = reduce_min (e[0]);
  out[1] = reduce_min (e[1]);
  out[2] = reduce_min (e[2]);
  out[3] = reduce_max (e[3]);
  out[4] = reduce_max (e[4]);
  out[5] = reduce_max (e[5]);
}

/* task based and vectorized extents of points */
void extents_of_points (uniform int ntasks, uniform unsigned int64 n, uniform REAL * uniform point[3], uniform REAL extents[])
{
  uniform int num = ntasks < 1 ? num_cores () : ntasks;
  uniform unsigned int64 span = n / num;

  uniform REAL * uniform task_extents = uniform new uniform REAL [6*num];

  launch[num] extrema (span, n, point[0], point[1], point[2], task_extents);
  sync;

  extents[0] = task_extents[0];
  extents[1] = task_extents[1];
  extents[2] = task_extents[2];
  extents[3] = task_extents[3];
  extents[4] = task_extents[4];
  extents[5] = task_extents[5];

  for (uniform int i = 1; i < num; i ++)
  {
    uniform REAL * uniform e = &task_extents [6*i];

    if (e[0] < extents[0]) extents[0] = e[0];
    if (e[1] < extents[1]) extents[1] = e[1];
    if (e[2] < extents[2]) extents[2] = e[2];
    if (e[3] > extents[3]) extents[3] = e[3];
    if (e[4] > extents[4]) extents[4] = e[4];
    if (e[5] > extents[5]) extents[5] = e[5];
  }

  delete task_extents;
}

static void rcb_tree_size (uniform unsigned int64 n, uniform int64 cutoff, uniform int * uniform tree_size)
{
  if (n > cutoff)
  {
    (*tree_size) += 2; /* two potential siblings */

    rcb_tree_size (n/2, cutoff, tree_size);
    rcb_tree_size (n-n/2, cutoff, tree_size);
  }
}

static void rcb_tree_init (uniform unsigned int64 n, uniform int64 cutoff, uniform partitioning tree[], uniform int node, uniform int * uniform i)
{
  if (n > cutoff) /* node */
  {
    tree[node].dimension = 3; /* mark as node node; actual dimension 0,1 or 2 will be determined in rcb_tree_task */
    tree[node].left = ++(*i);
    tree[node].right = ++(*i);

    rcb_tree_init (n/2, cutoff, tree, tree[node].left, i);
    rcb_tree_init (n-n/2, cutoff, tree, tree[node].right, i);
  }
  else /* leaf */
  {
    tree[node].dimension = -1; /* mark as leaf */
    tree[node].left = tree[node].right = -1;
    tree[node].rank = -1;
    tree[node].size = 0;
  }
}

static void leaf_count (uniform partitioning tree[], uniform int node, uniform int * uniform i)
{
  if (tree[node].dimension >= 0)
  {
    leaf_count (tree, tree[node].left, i);
    leaf_count (tree, tree[node].right, i);
  }
  else (*i) ++;
}

/* O(n) split of point[] such that point[d][i<k] <= point[d][i>=k]; other dimensions are copied accordingly */
static uniform REAL quick_split (uniform unsigned int64 n, uniform REAL * uniform point[3], uniform int d, uniform unsigned int64 k)
{
  uniform REAL pivot = point[d][n/2]; /* random index might be more effective; XXX */

  uniform int d1, d2;

  switch (d)
  {
  case 0:
    d1 = 1;
    d2 = 2;
    break;
  case 1:
    d1 = 0;
    d2 = 2;
    break;
  case 2:
    d1 = 0;
    d2 = 1;
    break;
  }

  uniform unsigned int64 i = 0, j = n-1;

  while (i < j)
  {
    while (i < j && point[d][i] <= pivot) i++;
    while (i < j && point[d][j] > pivot) j--;

    if (i < j)
    {
      uniform REAL temp = point[d][i];
      point[d][i] = point[d][j];
      point[d][j] = temp;

      temp = point[d1][i];
      point[d1][i] = point[d1][j];
      point[d1][j] = temp;

      temp = point[d2][i];
      point[d2][i] = point[d2][j];
      point[d2][j] = temp;
    }
  }

  if (i == k) return pivot;
  else if (i > k) return quick_split (i, point, d, k);
  else
  {
    uniform REAL * uniform rpoint[3] = {point[0]+i, point[1]+i, point[2]+i};
    return quick_split (n-i, rpoint, d, k-i);
  }
}

task void rcb_tree_task (uniform int ntasks, uniform unsigned int64 n, uniform REAL * uniform point[3], uniform partitioning tree[], uniform int node)
{
  if (tree[node].dimension >= 0)
  {
    uniform REAL extents[6];

    extents_of_points (ntasks, n, point, extents);

    extents[0] = extents[3] - extents[0];
    extents[1] = extents[4] - extents[1];
    extents[2] = extents[5] - extents[2];

    uniform int dimension = 0;

    if (extents[1] > extents[0]) dimension = 1;
    if (extents[2] > extents[dimension]) dimension = 2;

    tree[node].dimension = dimension;

    uniform int left_count = 0, right_count = 0;

    leaf_count (tree, tree[node].left, &left_count);
    leaf_count (tree, tree[node].right, &right_count);

    uniform unsigned int64 k = (REAL) n * (REAL) left_count / (REAL) (left_count + right_count);

    tree[node].coord = quick_split (n, point, dimension, k);

    uniform REAL * uniform rpoint[3] = {point[0]+k, point[1]+k, point[2]+k};

    launch rcb_tree_task (ntasks, k, point, tree, tree[node].left);
    launch rcb_tree_task (ntasks, n-k, rpoint, tree, tree[node].right);
  }
}

/* create rcb tree; uniformly bisect untill leaf size <= cutoff; or if cutoff < 0 then create -cutoff equal size leaves */
uniform partitioning * uniform rcb_tree_create (uniform int ntasks, uniform unsigned int64 n,
  uniform REAL * uniform point[3], uniform int64 cutoff, uniform int * uniform tree_size)
{
  *tree_size = 1;
  
  if (cutoff < 0)
  {
    rcb_tree_size (-cutoff, 1, tree_size);
  }
  else
  {
    cutoff = MAX (cutoff, 1); /* 0 would cause trouble */

    rcb_tree_size (n, cutoff, tree_size);
  }

  uniform partitioning * uniform tree = uniform new uniform partitioning [*tree_size];

  uniform int i = 0;

  if (cutoff < 0)
  {
    rcb_tree_init (-cutoff, 1, tree, 0, &i);
  }
  else
  {
    rcb_tree_init (n, cutoff, tree, 0, &i);
  }

  launch rcb_tree_task (ntasks, n, point, tree, 0);
  sync;

  return tree;
}

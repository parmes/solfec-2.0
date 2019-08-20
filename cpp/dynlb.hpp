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

#include "part_ispc.h"

#ifndef __dynlb__
#define __dynlb__

struct dynlb /* load balancer interface */
{
  int64_t cutoff; /* partitioning tree cutoff; 0 means use default selection */
  REAL epsilon; /* imbalance epsilon; rebalance when imbalance > 1.0 + epsilon */
  ispc::partitioning *ptree; /* partitioning tree; used internally */
  int ptree_size; /* partitioning tree size; used internally */
  REAL imbalance; /* current imbalance */

  /* default constructor */
  dynlb(): cutoff(0), epsilon(0.), ptree(NULL), ptree_size(0), imbalance(1.) { }

  /* create local load balancer */
  void local_create (uint64_t n, REAL *point[3], int64_t cutoff = 0, REAL epsilon = 0.1);

  /* assign an MPI rank to a point; return this rank */
  int point_assign (REAL point[]);

  /* assign MPI ranks to a box spanned between lo and hi points; return the number of ranks assigned */
  int box_assign (REAL lo[], REAL hi[], int ranks[]);

  /* update local load balancer */
  void local_update (uint64_t n, REAL *point[3]);

  /* destroy load balancer */
  ~dynlb();
};

#endif

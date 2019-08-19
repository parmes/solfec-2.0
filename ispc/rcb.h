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

#ifndef __rcb__
#define __rcb__

struct partitioning /* binary space partitioning recurisve bisection tree */
{
  uniform REAL coord;
  uniform int dimension;
  uniform int left;
  uniform int right;

  uniform int rank; /* > 0 for leaves */
  uniform unsigned int64 size; /* leaf size */
};

/* create rcb tree; uniformly bisect untill leaf size <= cutoff; or if cutoff < 0 then create -cutoff equal size leaves */
uniform partitioning * uniform rcb_tree_create (uniform int ntasks, uniform unsigned int64 n,
  uniform REAL * uniform point[3], uniform int64 cutoff, uniform int * uniform tree_size);

#endif


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

/* aligned real allocator */
export uniform REAL * uniform  aligned_real_alloc (uniform unsigned int64 n)
{
  return uniform new uniform REAL [n];
}

export void aligned_real_free (uniform REAL * uniform ptr)
{
  delete ptr;
}

/* aligned int allocator */
export uniform int * uniform  aligned_int_alloc (uniform unsigned int64 n)
{
  return uniform new uniform int [n];
}

export void aligned_int_free (uniform int * uniform ptr)
{
  delete ptr;
}

/* aligned unsigned int allocator */
export uniform unsigned int * uniform  aligned_uint_alloc (uniform unsigned int64 n)
{
  typedef unsigned int uint;

  return uniform new uniform uint [n];
}

export void aligned_uint_free (uniform unsigned int * uniform ptr)
{
  delete ptr;
}

/* aligned int64 allocator */
export uniform int64 * uniform  aligned_int64_alloc (uniform unsigned int64 n)
{
  return uniform new uniform int64 [n];
}

export void aligned_int64_free (uniform int64 * uniform ptr)
{
  delete ptr;
}

/* aligned unsigned int64 allocator */
export uniform unsigned int64 * uniform  aligned_uint64_alloc (uniform unsigned int64 n)
{
  typedef unsigned int64 uint64;

  return uniform new uniform uint64 [n];
}

export void aligned_uint64_free (uniform unsigned int64 * uniform ptr)
{
  delete ptr;
}

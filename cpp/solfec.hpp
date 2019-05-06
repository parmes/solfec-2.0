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

#ifndef __solfec__
#define __solfec__

/* compilation-level constants */
#define BLOCK_SIZE 64
#define REAL float
/* --------------------------- */

enum fe_type {TET, PYR, WED, HEX};

struct fe_block
{
  uint32_t bid; /* body identifier */

  fe_type type; /* fixed type blocks */

  REAL data[BLOCK_SIZE]; /* static per-component data blocks for ISPC routines */

  uint8_t size;

  /* XXX: one map<string,block> instead ? */

  /* XXX: a global mempool for all element data blocks */

  /* XXX: per block algebra ? (additive-Schwarz-overlapped) */

  REAL extents[6]; /* for contact detection and MPI rank migration */
};

struct solfec
{
  REAL *nodal_data; /* global nodal data via MPI-3.0 shared memory blocks */

  /* XXX: use one map<string,vector> instead ? */

  std::list<fe_block> fe_blocks;
};

#endif

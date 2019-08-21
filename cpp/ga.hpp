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

/* Global Array implementation based on the source
 * code include with the book "Using Advanced MPI":
 * [1] https://mitpress.mit.edu/books/using-advanced-mpi
 * [2] https://books.google.pl/books?id=GY5IBQAAQBAJ
 * [3] http://wgropp.cs.illinois.edu/usingmpiweb/index.html
 */

#include <mpi.h>
#include "mutex.hpp"

#ifndef __ga__
#define __ga__

struct GA
{
  MPI_Win ga_win;
  MPE_MUTEX mutex;
  MPI_Win lock_win;
  /* datatype and size */
  MPI_Datatype dtype;
  int dtype_size;
  /* local sizes of the global array */
  uint64_t dim1, dim2;
  GA(MPI_Comm comm, uint64_t dim1, uint64_t dim2, MPI_Datatype dtype);
  ~GA();
  /* 0-based indexing [ilo, ihigh), [jlo, jhigh) */
  void acc(int rank, uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf);
  void get(int rank, uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf);
  void put(int rank, uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf);
  /* fence synchronization on all ranks */
  void fence();
};

#endif

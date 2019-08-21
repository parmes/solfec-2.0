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

/* Contributors: Tomasz Koziara */

#include <climits>
#include "ga.hpp"

GA::GA(MPI_Comm comm, uint64_t dim1, uint64_t dim2, MPI_Datatype dtype): dim1(dim1), dim2(dim2)
{
  MPI_Aint local_size;
  void *ga_win_ptr;
  MPI_Info info;
  int size;

  /* Determine sizes */
  MPI_Comm_size(comm, &size);
  MPI_Type_size(dtype, &dtype_size);
  local_size = dim1 * dim2 * dtype_size;

  /* Specify ordering of accumulate operations (this is the default behavior in MPI-3) */
  MPI_Info_create(&info);
  MPI_Info_set(info,"accumulate_ordering", "rar,raw,war,waw");

  /* Allocate memory and create window */
  MPI_Win_allocate(local_size, dtype_size, info, comm, &ga_win_ptr, &ga_win);
  MPI_Info_free(&info);

  /* Create critical section window */
  mutex.create(comm, size, &lock_win);
}

GA::~GA()
{
  void *ga_win_ptr;
  int flag;

  MPI_Win_get_attr(ga_win, MPI_WIN_BASE, &ga_win_ptr, &flag);

  MPI_Win_free(&ga_win);

  if (flag && ga_win_ptr) MPI_Free_mem(ga_win_ptr);

  mutex.free(&lock_win);
}

void GA::acc(int rank, uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf)
{
  MPI_Datatype vtype;
  MPI_Aint disp;
  uint64_t j;

  /* In order to ensure that the entire update is atomic,
     we must first mutex-lock the window that we will access */
  mutex.acquire(lock_win, rank);

  MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK, ga_win);

  if (dim1 < INT_MAX)
  {
    MPI_Type_vector(jhigh-jlo, ihigh-ilo, (int)dim1, dtype, &vtype);

    MPI_Type_commit(&vtype);

    disp = jlo * dim1 + ilo;

    MPI_Accumulate(buf, (ihigh - ilo)*(jhigh - jlo), dtype, rank, disp, 1, vtype, MPI_SUM, ga_win);

    MPI_Type_free(&vtype);
  }
  else for (j=jlo; j<jhigh; j++)
  {
    disp = j * dim1 + ilo;

    MPI_Accumulate(buf, ihigh - ilo, dtype, rank, disp, ihigh - ilo, dtype, MPI_SUM, ga_win);

    buf = (void *)(((char *)buf) + (ihigh - ilo) *  dtype_size);
  }

  MPI_Win_unlock(rank, ga_win);

  mutex.release(lock_win, rank);
}

void GA::get(int rank, uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf)
{
  MPI_Datatype vtype;
  MPI_Aint disp;
  uint64_t j;

  mutex.acquire(lock_win,rank);

  /* Using lock_shared allows get accesses to proceed */
  MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK, ga_win);

  if (dim1 < INT_MAX)
  {
    MPI_Type_vector(jhigh-jlo, ihigh-ilo, (int)dim1, dtype, &vtype);

    MPI_Type_commit(&vtype);

    disp = jlo * dim1 + ilo;

    MPI_Get(buf, (ihigh - ilo)*(jhigh - jlo), dtype, rank, disp, 1, vtype, ga_win);

    MPI_Type_free(&vtype);
  }
  else for (j=jlo; j<jhigh; j++)
  {
    disp = j * dim1 + ilo;

    MPI_Get(buf, ihigh - ilo, dtype, rank, disp, ihigh - ilo, dtype, ga_win);

    buf = (void *)(((char *)buf) + (ihigh - ilo) *  dtype_size);
  }

  MPI_Win_unlock(rank, ga_win);

  mutex.release(lock_win,rank);
}

void GA::put(int rank, uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf)
{
  MPI_Datatype vtype;
  MPI_Aint disp;
  uint64_t j;

  mutex.acquire(lock_win,rank);

  /* Using lock_shared allows get accesses to proceed */
  MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK, ga_win);

  if (dim1 < INT_MAX)
  {
    MPI_Type_vector(jhigh-jlo, ihigh-ilo, (int)dim1, dtype, &vtype);

    MPI_Type_commit(&vtype);

    disp = jlo * dim1 + ilo;

    MPI_Put(buf, (ihigh - ilo)*(jhigh - jlo), dtype, rank, disp, 1, vtype, ga_win);

    MPI_Type_free(&vtype);
  }
  else for (j=jlo; j<jhigh; j++)
  {
    disp = j * dim1 + ilo;

    MPI_Put(buf, ihigh - ilo, dtype, rank, disp, ihigh - ilo, dtype, ga_win);

    buf = (void *)(((char *)buf) + (ihigh - ilo) *  dtype_size);
  }

  MPI_Win_unlock(rank, ga_win);

  mutex.release(lock_win,rank);
}

void GA::fence()
{
  MPI_Win_fence (0, ga_win);
}

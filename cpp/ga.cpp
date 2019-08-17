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

#include "ga.hpp"

GA::GA(MPI_Comm comm, uint64_t dim1, uint64_t dim2, MPI_Datatype dtype): dim1(dim1), dim2(dim2)
{
  MPI_Aint local_size;
  void *ga_win_ptr;
  MPI_Info info;
  int size;

  /* Determine size of GA memory */
  MPI_Comm_size(comm, &size);
  chunk2 = dim2 / size;
  /* Require size to exactly divide dim2 */
  if ((dim2 % size) != 0) MPI_Abort(comm, 1);
  MPI_Type_size(dtype, &dtype_size);
  local_size = dim1 * chunk2 * dtype_size;

  /* Specify ordering of accumulate operations (this is the
     default behavior in MPI-3) */
  MPI_Info_create(&info);
  MPI_Info_set(info,"accumulate_ordering", "rar,raw,war,waw");

  /* Allocate memory and create window */
  MPI_Win_allocate(local_size, dtype_size, info, comm,
		   &ga_win_ptr, &ga_win);
  MPI_Info_free(&info);

  /* Create critical section window */
  mutex.create(comm, size, &lock_win);
}

GA::~GA()
{
  int flag;
  void *ga_win_ptr;

  MPI_Win_get_attr(ga_win, MPI_WIN_BASE, &ga_win_ptr, &flag);
  MPI_Win_free(&ga_win);
  if (flag && ga_win_ptr)
      MPI_Free_mem(ga_win_ptr);
  mutex.free(&lock_win);
}

int GA::acc(uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf)
{
  uint64_t jcur, jfirst, jlast, j;
  int rank, rank_first, rank_last;
  MPI_Aint disp;

  /* In order to ensure that the entire update is atomic, we must
     first mutex-lock all of the windows that we will access */
  rank_first = jlo / chunk2;
  rank_last = jhigh / chunk2;
  for (rank = rank_first; rank < rank_last; rank++) {
      mutex.acquire(lock_win, rank);
  }

  jcur = jlo;
  while (jcur < jhigh) {
    rank   = (jcur - 1) /chunk2;
    jfirst = rank * chunk2 + 1;
    jlast  = (rank + 1) * chunk2;
    if (jlast >= jhigh) jlast = jhigh;

    MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK,
                 ga_win);
    for (j=jcur; j<=jlast; j++) {
      disp = (j - jfirst) * dim1 + ilo;
      MPI_Accumulate(buf, ihigh - ilo, dtype,
                     rank, disp, ihigh - ilo, dtype,
                     MPI_SUM, ga_win);
      buf = (void *)( ((char *)buf) +
                      (ihigh - ilo) *  dtype_size);
    }
    MPI_Win_unlock(rank, ga_win);

    mutex.release(lock_win, rank);
    jcur = jlast + 1;
  }
  return 0;
}

int GA::get(uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf)
{
  uint64_t jcur, jfirst, jlast, j;
  MPI_Aint disp;
  int rank;

  jcur = jlo;
  while (jcur < jhigh) {
      rank   = (jcur - 1) /chunk2;
      jfirst = rank * chunk2 + 1;
      jlast  = (rank + 1) * chunk2;
      if (jlast >= jhigh) jlast = jhigh;

      mutex.acquire(lock_win,rank);

      /* Using lock_shared allows get accesses to proceed */
      MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK,
		   ga_win);
      for (j=jcur; j<=jlast; j++) {
	  disp = (j - jfirst) * dim1 + ilo;
	  MPI_Get(buf, ihigh - ilo, dtype,
		  rank, disp, ihigh - ilo, dtype,
		  ga_win);
	  buf = (void *)( ((char *)buf) +
			  (ihigh - ilo) *  dtype_size );
      }
      MPI_Win_unlock(rank, ga_win);

      mutex.release(lock_win,rank);
      jcur = jlast + 1;
  }
  return 0;
}

int GA::put(uint64_t ilo, uint64_t ihigh, uint64_t jlo, uint64_t jhigh, void *buf)
{
  uint64_t jcur, jfirst, jlast, j;
  MPI_Aint disp;
  int rank;

  jcur = jlo;
  while (jcur < jhigh) {
      rank   = (jcur - 1) /chunk2;
      jfirst = rank * chunk2 + 1;
      jlast  = (rank + 1) * chunk2;
      if (jlast >= jhigh) jlast = jhigh;

      mutex.acquire(lock_win,rank);

      /* Using lock_shared allows get accesses to proceed */
      MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK,
		   ga_win);
      for (j=jcur; j<=jlast; j++) {
	  disp = (j - jfirst) * dim1 + ilo;
	  MPI_Put(buf, ihigh - ilo, dtype,
		  rank, disp, ihigh - ilo, dtype,
		  ga_win);
	  buf = (void *)( ((char *)buf) +
			  (ihigh - ilo) *  dtype_size );
      }
      MPI_Win_unlock(rank, ga_win);

      mutex.release(lock_win,rank);
      jcur = jlast + 1;
  }
  return 0;
}

void GA::fence()
{
  MPI_Win_fence (0, ga_win);
}
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

/* MPI mutex implementation based on the source
 * code include with the book "Using Advanced MPI":
 * [1] https://mitpress.mit.edu/books/using-advanced-mpi
 * [2] https://books.google.pl/books?id=GY5IBQAAQBAJ
 * [3] http://wgropp.cs.illinois.edu/usingmpiweb/index.html
 */

/* Contributors: Tomasz Koziara */

#include "mutex.hpp"

int MPE_MUTEX::create(MPI_Comm comm, int num, MPI_Win *mutex_win)
{
  int rank, size, *counterMem = 0;
  int lnum, lleft, i;
  MPI_Aint counterSize = 0;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  lnum  = num / size;
  lleft = num % size;
  if (rank < lleft) lnum++;
  counterSize = lnum * sizeof(int);
  if (counterSize > 0) {
    MPI_Alloc_mem(counterSize, MPI_INFO_NULL, &counterMem);
    for (i=0; i<lnum; i++) counterMem[i] = 0;
  }
  MPI_Win_create(counterMem, counterSize, sizeof(int),
		 MPI_INFO_NULL, MPI_COMM_WORLD, mutex_win);
  if (MPE_MUTEX_KEYVAL == MPI_KEYVAL_INVALID) {
    MPI_Win_create_keyval(MPI_WIN_NULL_COPY_FN, MPI_WIN_NULL_DELETE_FN,
			  &MPE_MUTEX_KEYVAL, (void*)0);
  }
  MPI_Win_set_attr(*mutex_win, MPE_MUTEX_KEYVAL, (void*)(MPI_Aint)(size));
  return 0;
}

int MPE_MUTEX::free(MPI_Win *mutex_win)
{
  int flag, *counterMem;

  MPI_Win_get_attr(*mutex_win, MPI_WIN_BASE, &counterMem, &flag);
  if (flag && counterMem != 0) {
    MPI_Free_mem(counterMem);
  }
  MPI_Win_free(mutex_win);
  return 0;
}

int MPE_MUTEX::acquire(MPI_Win mutex_win, int num)
{
  int mone = -1, one=1, oldval;
  int lrank, flag, size, *attrval;
  MPI_Aint lidx;

  /* Compute the location of the counter */
  MPI_Win_get_attr(mutex_win, MPE_MUTEX_KEYVAL, &attrval, &flag);
  if (!flag) return -1;  /* Error: counterWin not setup */
  size = (int)(MPI_Aint)attrval;  /* We stored the integer as a
                                     pointer */
  lrank = num % size; lidx  = num / size;

  MPI_Win_lock(MPI_LOCK_SHARED, lrank, 0, mutex_win);
  do {
    MPI_Fetch_and_op(&one, &oldval, MPI_INT,
                     lrank, lidx, MPI_SUM, mutex_win);
    MPI_Win_flush(lrank, mutex_win);
    if (oldval == 0) break;
    MPI_Accumulate(&mone, 1, MPI_INT, lrank, lidx, 1, MPI_INT,
                   MPI_SUM, mutex_win);
    MPI_Win_flush(lrank, mutex_win);
    /* We could wait a little bit, depending on oldval */
  } while (1);
  MPI_Win_unlock(lrank, mutex_win);
  return 0;
}

int MPE_MUTEX::release(MPI_Win mutex_win, int num)
{
  int mone = -1;
  int lrank, flag, size, *attrval;
  MPI_Aint lidx;

  /* Compute the location of the counter */
  MPI_Win_get_attr(mutex_win, MPE_MUTEX_KEYVAL, &attrval, &flag);
  if (!flag) return -1;           /* Error: counterWin setup */
  size = (int)(MPI_Aint)attrval;  /* We stored the integer as a
                                     pointer */
  lrank = num % size; lidx  = num / size;

  MPI_Win_lock(MPI_LOCK_SHARED, lrank, 0, mutex_win);
  MPI_Accumulate(&mone, 1, MPI_INT, lrank, lidx, 1, MPI_INT,
                 MPI_SUM, mutex_win);
  MPI_Win_unlock(lrank, mutex_win);
  return 0;
}

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

/* MSC lock implementation based on the source
 * code include with the book "Using Advanced MPI":
 * [1] https://mitpress.mit.edu/books/using-advanced-mpi
 * [2] https://books.google.pl/books?id=GY5IBQAAQBAJ
 * [3] http://wgropp.cs.illinois.edu/usingmpiweb/index.html
 */

/* Contributors: Tomasz Koziara */

#include "lock.hpp"

void MSC_LOCK::init(MPI_Comm comm, MPI_Win *lock_win)
{
  int *lmem, rank;
  MPI_Aint winsize;
  MPI_Comm_rank (comm, &rank);

  if (MSC_LOCKRANK == MPI_KEYVAL_INVALID)
    MPI_Win_create_keyval (MPI_WIN_NULL_COPY_FN,
                           MPI_WIN_NULL_DELETE_FN,
			   &MSC_LOCKRANK, (void*)0);

  winsize = 2 * sizeof(int);
  if (rank == 0) winsize += sizeof(int);
  MPI_Win_allocate(winsize, sizeof(int), MPI_INFO_NULL, comm,
                   &lmem, lock_win);
  lmem[nextRank] = -1;
  lmem[blocked] = 0;
  if (rank == 0) {
    lmem[lockTail] = -1;
  }
  MPI_Win_set_attr(*lock_win, MSC_LOCKRANK, (void*)(MPI_Aint)rank);
  MPI_Barrier(comm);
}

void MSC_LOCK::acquire(MPI_Win lock_win)
{
  int flag, myrank, predecessor, *lmem;
  void *attrval;

  MPI_Win_get_attr(lock_win, MSC_LOCKRANK, &attrval, &flag);
  myrank = (int)(MPI_Aint) attrval;
  MPI_Win_get_attr(lock_win, MPI_WIN_BASE, &lmem, &flag);
  lmem[blocked] = 1; /* in case we are blocked */
  MPI_Win_lock_all(0, lock_win);
  MPI_Fetch_and_op(&myrank, &predecessor, MPI_INT,
                   0, lockTail, MPI_REPLACE, lock_win);
  MPI_Win_flush (0, lock_win);
  if (predecessor != -1) {
    /* We didn't get the lock. Add us to the tail of the list */
    MPI_Accumulate(&myrank, 1, MPI_INT, predecessor,
                   nextRank, 1, MPI_INT, MPI_REPLACE, lock_win);
    /* Now spin on our local value "blocked" until we are
       given the lock */
    do {
      MPI_Win_sync(lock_win);
    } while(lmem[blocked] == 1);
  }
  // else we have the lock
  MPI_Win_unlock_all(lock_win);
}

void MSC_LOCK::release(MPI_Win lock_win)
{
  int nullrank=-1, zero=0, myrank, curtail, flag, *lmem;
  void *attrval;

  MPI_Win_get_attr(lock_win, MSC_LOCKRANK, &attrval, &flag);
  myrank = (int)(MPI_Aint)attrval;
  MPI_Win_get_attr(lock_win, MPI_WIN_BASE, &lmem, &flag);
  MPI_Win_lock_all(0, lock_win);
  if (lmem[nextRank] == -1) {
    /* See if we're waiting for th enext to notify us */
    MPI_Compare_and_swap (&nullrank, &myrank, &curtail, MPI_INT,
                          0, lockTail, lock_win);
    if (curtail == myrank) {
      /* We are the only process in the list */
      MPI_Win_unlock_all(lock_win);
      return;
    }
    /* Otherwise, someone else has added themeselves to the list. */
    do {
      MPI_Win_sync(lock_win);
    } while (lmem[nextRank] == -1);
  }
  /* Now we can notify them. Use accumulate with replace instead
     of put since we want an atomic update of the location */
  MPI_Accumulate (&zero, 1, MPI_INT, lmem[nextRank], blocked,
                  1, MPI_INT, MPI_REPLACE, lock_win);
  MPI_Win_unlock_all(lock_win);
}

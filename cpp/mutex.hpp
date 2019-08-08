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

#include <mpi.h>

#ifndef __mutex__
#define __mutex__

struct MPE_MUTEX
{
  int MPE_MUTEX_KEYVAL;
  MPE_MUTEX(): MPE_MUTEX_KEYVAL(MPI_KEYVAL_INVALID) {};
  int create(MPI_Comm comm, int num, MPI_Win *mutex_win);
  int free(MPI_Win *mutex_win);
  int acquire(MPI_Win mutex_win, int num);
  int release(MPI_Win mutex_win, int num);
};

#endif

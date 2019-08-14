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

#ifndef __rmadds_vect__
#define __rmadds_vect__

/* FIXME: unedited below */

typedef struct _GA {
  MPI_Win ga_win;
  MPI_Win lock_win;
  /* Datatype and size */
  MPI_Datatype dtype;
  int dtype_size;
  /* sizes of the global array */
  int dim1, dim2, chunk2;
} *GA;

int ga_create(MPI_Comm comm, int dim1, int dim2,
              MPI_Datatype dtype, GA *ga)
{
    GA       new_ga;
    int      size, chunk2, sizeoftype;
    MPI_Aint local_size;
    MPI_Info info;
    void     *ga_win_ptr;

    /* Get a new structure */
    new_ga = (GA)malloc(sizeof(struct _GA));
    if (!new_ga) return 0;
    /* Determine size of GA memory */
    MPI_Comm_size(comm, &size);
    chunk2 = dim2 / size;
    /* Require size to exactly divide dim2 */
    if ((dim2 % size) != 0) MPI_Abort(comm, 1);
    MPI_Type_size(dtype, &sizeoftype);
    local_size = dim1 * chunk2 * sizeoftype;

    /* Specify ordering of accumulate operations (this is the
       default behavior in MPI-3) */
    MPI_Info_create(&info);
    MPI_Info_set(info,"accumulate_ordering", "rar,raw,war,waw");

    /* Allocate memory and create window */
    MPI_Win_allocate(local_size, sizeoftype, info, comm,
		     &ga_win_ptr, &new_ga->ga_win);
    MPI_Info_free(&info);

    /* Create critical section window */
    MPE_Mutex_create(comm, size, &new_ga->lock_win);

    /* Save other data and return */
    new_ga->dtype      = dtype;   new_ga->dtype_size = sizeoftype;
    new_ga->dim1       = dim1;    new_ga->dim2       = dim2;
    new_ga->chunk2     = chunk2;
    *ga                = new_ga;
    return 1;
}

int ga_free(GA ga)
{
    int flag;
    void *ga_win_ptr;

    MPI_Win_get_attr(ga->ga_win, MPI_WIN_BASE, &ga_win_ptr,
		     &flag);
    MPI_Win_free(&ga->ga_win);
    if (flag && ga_win_ptr)
        MPI_Free_mem(ga_win_ptr);
    MPE_Mutex_free(&ga->lock_win);

    free(ga);
    return 0;
}

int ga_acc(GA ga, int ilo, int ihigh, int jlo, int jhigh,
           void *buf)
{
 int      jcur, jfirst, jlast, j, rank, rank_first, rank_last;
 MPI_Aint disp;

 /* In order to ensure that the entire update is atomic, we must
    first mutex-lock all of the windows that we will access */
 rank_first = (jlo - 1) / ga->chunk2;
 rank_last  = (jhigh - 1) / ga->chunk2;
 for (rank = rank_first; rank <= rank_last; rank++) {
     MPE_Mutex_acquire(ga->lock_win, rank);
 }

 jcur = jlo;
 while (jcur <= jhigh) {
   rank   = (jcur - 1) /ga->chunk2;
   jfirst = rank * ga->chunk2 + 1;
   jlast  = (rank + 1) * ga->chunk2;
   if (jlast > jhigh) jlast = jhigh;

   MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK,
                ga->ga_win);
   for (j=jcur; j<=jlast; j++) {
       disp = (j - jfirst) * ga->dim1 + (ilo - 1);
       MPI_Accumulate(buf, ihigh - ilo + 1, ga->dtype,
                      rank, disp, ihigh - ilo + 1, ga->dtype,
                      MPI_SUM, ga->ga_win);
       buf = (void *)( ((char *)buf) +
                      (ihigh - ilo + 1) *  ga->dtype_size);
   }
   MPI_Win_unlock(rank, ga->ga_win);

   MPE_Mutex_release(ga->lock_win, rank);
   jcur = jlast + 1;
 }
 return 0;
}

int ga_put(GA ga, int ilo, int ihigh, int jlo, int jhigh,
	   void *buf)
{
    int      jcur, jfirst, jlast, j, rank;
    MPI_Aint disp;

    jcur = jlo;
    while (jcur <= jhigh) {
	rank   = (jcur - 1) /ga->chunk2;
	jfirst = rank * ga->chunk2 + 1;
	jlast  = (rank + 1) * ga->chunk2;
	if (jlast > jhigh) jlast = jhigh;

	MPE_Mutex_acquire(ga->lock_win,rank);

	/* Using lock_shared allows get accesses to proceed */
	MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK,
		     ga->ga_win);
	for (j=jcur; j<=jlast; j++) {
	    disp = (j - jfirst) * ga->dim1 + (ilo - 1);
	    MPI_Get(buf, ihigh - ilo + 1, ga->dtype,
		    rank, disp, ihigh - ilo + 1, ga->dtype,
		    ga->ga_win);
	    buf = (void *)( ((char *)buf) +
			    (ihigh - ilo + 1) *  ga->dtype_size );
	}
	MPI_Win_unlock(rank, ga->ga_win);

	MPE_Mutex_release(ga->lock_win,rank);
	jcur = jlast + 1;
    }
    return 0;
}

int ga_put(GA ga, int ilo, int ihigh, int jlo, int jhigh,
	   void *buf)
{
    int      jcur, jfirst, jlast, j, rank;
    MPI_Aint disp;

    jcur = jlo;
    while (jcur <= jhigh) {
	rank   = (jcur - 1) /ga->chunk2;
	jfirst = rank * ga->chunk2 + 1;
	jlast  = (rank + 1) * ga->chunk2;
	if (jlast > jhigh) jlast = jhigh;

	MPE_Mutex_acquire(ga->lock_win,rank);

	/* Using lock_shared allows get accesses to proceed */
	MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK,
		     ga->ga_win);
	for (j=jcur; j<=jlast; j++) {
	    disp = (j - jfirst) * ga->dim1 + (ilo - 1);
	    MPI_Put(buf, ihigh - ilo + 1, ga->dtype,
		    rank, disp, ihigh - ilo + 1, ga->dtype,
		    ga->ga_win);
	    buf = (void *)( ((char *)buf) +
			    (ihigh - ilo + 1) *  ga->dtype_size );
	}
	MPI_Win_unlock(rank, ga->ga_win);

	MPE_Mutex_release(ga->lock_win,rank);
	jcur = jlast + 1;
    }
    return 0;
}

#endif

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

/* Distributed list implementation based on the source
 * code include with the book "Using Advanced MPI":
 * [1] https://mitpress.mit.edu/books/using-advanced-mpi
 * [2] https://books.google.pl/books?id=GY5IBQAAQBAJ
 * [3] http://wgropp.cs.illinois.edu/usingmpiweb/index.html
 */

/* Contributors: Tomasz Koziara */

#include <cstring>
#include <mpi.h>
#include "lock.h"

#ifndef __dlist__
#define __dlist__

struct RemotePointer
{
  MPI_Aint disp; /* Displacement in window */
  int rank; /* Rank (process) */
  void *local_pointer; /* Local address of data pointed at (if data local) */
};

template<typename key_type, typename value_type>
struct DListElm
{
  RemotePointer next;
  key_type key;
  value_type value;
};

#define DispInListElm( _dptr, _field ) \
        (MPI_Aint)&(((DListElm<key_type,value_type> *)((_dptr).disp))->_field)

template <typename key_type, typename value_type>
struct DList
{
  int wrank;
  MPI_Win listwin;
  MPI_Win lock_win;
  MSC_LOCK mcslock;
  RemotePointer headDptr;
  RemotePointer nullDptr;
  MPI_Datatype listelmType, dptrType;
  int MPE_LISTWIN_KEY_RANK;

  DList(): headDptr{0,-1,0}, nullDptr{0,-1,0}, MPE_LISTWIN_KEY_RANK(MPI_KEYVAL_INVALID) { };

  ~DList() { /* TODO */ };

  int Init() /* Initialize the mutex, datatype, window, and list head  */
  {
    DListElm<key_type,value_type> sampleElm;
    MPI_Datatype dtypes[3];
    MPI_Aint displ[3];
    int blens[3];

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    mcslock.init(MPI_COMM_WORLD, &lock_win);

    blens[0] = 1; blens[1] = 1;
    MPI_Get_address( &nullDptr.disp, &displ[0] );
    MPI_Get_address( &nullDptr.rank, &displ[1] );
    displ[1] = displ[1] - displ[0];
    displ[0] = 0;
    dtypes[0] = MPI_AINT;
    dtypes[1] = MPI_INT;
    MPI_Type_create_struct(2, blens, displ, dtypes, &dptrType);
    MPI_Type_commit(&dptrType);

    blens[0] = 1;                  dtypes[0] = dptrType;
    blens[1] = sizeof(key_type);   dtypes[1] = MPI_CHAR;
    blens[2] = sizeof(value_type); dtypes[2] = MPI_CHAR;
    MPI_Get_address(&sampleElm.next,&displ[0]);
    MPI_Get_address(&sampleElm.key[0],&displ[1]);
    MPI_Get_address(&sampleElm.value[0],&displ[2]);
    displ[2] -= displ[0];
    displ[1] -= displ[0];
    displ[0] = 0;
    MPI_Type_create_struct(3, blens, displ, dtypes, &listelmType);
    MPI_Type_commit(&listelmType);

    MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &listwin);

    MPI_Win_create_keyval(MPI_WIN_NULL_COPY_FN, MPI_WIN_NULL_DELETE_FN,
			  &MPE_LISTWIN_KEY_RANK, (void*)0);
    MPI_Win_set_attr(listwin, MPE_LISTWIN_KEY_RANK, (void*)(MPI_Aint)wrank);

    headDptr.rank = 0;
    if (wrank == 0) {
      DListElm<key_type,value_type> *headLptr;
      MPI_Alloc_mem(sizeof(DListElm), MPI_INFO_NULL, &headLptr);
      MPI_Get_address(headLptr,&headDptr.disp);
      headLptr->next.rank = -1;
      headLptr->next.disp = (MPI_Aint)MPI_BOTTOM;
      MPI_Win_attach(listwin, headLptr, sizeof(DListElm));
    }
    MPI_Bcast(&headDptr.disp, 1, MPI_AINT, 0, MPI_COMM_WORLD);
    return 0;
  }

  int Insert (const key_type &key, const value_type &value) /* duplicate keys are allowed */
  {
    RemotePointer dptr, last_dptr, new_dptr;
    DListElm elmOfptr, *new_lptr;
    int *attrval, myrank, flag;

    MPI_Win_get_attr(listwin, MPE_LISTWIN_KEY_RANK, &attrval, &flag);

    if (!flag) {
      /* Listwin not properly initialized */
      return 1;
    }
    myrank = (int)(MPI_Aint)attrval;  /* We store the rank in the
					 attrval, which is an
					 address-sized value */
    mcslock.acquire(lock_win);

    last_dptr = headDptr;
    MPI_Win_lock_all(0, listwin);
    MPI_Get(&dptr, 1, dptrType, last_dptr.rank,
      DispInListElm(last_dptr,next), 1, dptrType, listwin);
    MPI_Win_flush(last_dptr.rank, listwin);

    while (dptr.rank != nullDptr.rank) {
      MPI_Get(&elmOfptr, 1, listelmType, dptr.rank, dptr.disp,
	      1, listelmType, listwin);
      MPI_Win_flush(dptr.rank, listwin);
      /* elm is what ptr points to (i.e., *ptr) */
      #if 0
      if (elmOfptr.key == key)
      {
	/* Duplicate key.  Ignore */
	MPI_Win_unlock_all(listwin);
	return 0;
      }
      #endif
      if (key <= elmOfptr.key) break; /* Insert in front of this */
      last_dptr = dptr;
      dptr = elmOfptr.next; /* i.e., ptr->next */
    }

    /* Create new element */
    MPI_Alloc_mem(sizeof(DListElm), MPI_INFO_NULL, &new_lptr);
    std::memcpy(new_lptr->key, key, sizeof(key));
    std::memcpy(new_lptr->value, value, sizeof(value));
    new_lptr->next = dptr;
    MPI_Win_attach(listwin, new_lptr, sizeof(DListElm));

    new_dptr.rank = myrank;
    MPI_Get_address(new_lptr,&new_dptr.disp);
    MPI_Put(&new_dptr, 1, dptrType, last_dptr.rank, 
	DispInListElm(last_dptr,next), 1, dptrType, listwin);
    MPI_Win_unlock_all(listwin);
    mcslock.release(lock_win);
    return 0;
  }

  int Delete (const key_type &key) /* delete all items with key */
  {
    /* TODO */
  }

  /* find next value with key, given previously found value */
  std::pair<RemotePointer,const value_type*> Find (const key_type &key, RemoteValue ptr=headDptr)
  {
    static DListElm local_copy;
    int myrank, *attrval, flag;
    DListElm *local_copy_ptr;
    MPI_Group win_group;

    MPI_Win_get_attr(listwin, MPE_LISTWIN_KEY_RANK, &attrval, &flag);
    if (!flag) {
      /* listwin not properly initialized */
      return {nullDptr, NULL};
    }
    myrank = (int)(MPI_Aint)attrval;  /* We store the rank in the attrval,
                                         which is an address-sized value */
    MPI_Win_lock_all(0, listwin);
    while (ptr.rank >= 0) {
      /* Make sure we have the data */
      if (ptr.rank != myrank) {
        MPI_Get(&local_copy, 1, listelmType,
                ptr.rank, ptr.disp, 1, listelmType, listwin);
	MPI_Win_flush(ptr.rank, listwin);
	local_copy_ptr = &local_copy;
      } else
	local_copy_ptr = (DListElm *)(ptr.local_pointer); /* FIXME: local_pointer is never set */

      if (key == local_copy_ptr->key) {
        MPI_Win_unlock_all(listwin);
	return {local_copy_ptr->next, &local_copy_ptr->value}; /* FIXME: is this OK? does ->value memory persist? */
      }
      ptr = local_copy_ptr->next;
    }
    MPI_Win_unlock_all(listwin);
    return {nullDptr, NULL};  /* Did not find key */
  }

  /* get next list element */
  RemotePointer Next(DListElem &nextElm, RemotePointer ptr=headDptr)
  {
    RemotePointer nextDptr;

    mcslock.acquire(lock_win); /* FIXME: is this locking necessary? */
    MPI_Win_lock_all(0, listwin); /* XXX: also here? */
    if (ptr.rank != nullDptr.rank) {
      /* Get the list element pointed at by ptr */
      MPI_Get(&nextElm, 1, listelmType, ptr.rank, ptr.disp,
	      1, listelmType, listwin);
      MPI_Win_flush(ptr.rank, listwin);
      nextDptr = nextElm.next;
    }
    MPI_Win_unlock_all(listwin);
    mcslock.release(lock_win);
    return nextDptr;
  }
};

#endif

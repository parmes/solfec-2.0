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

#ifndef __compute__
#define __compute__

/* insert solfec::materials[matnum] into computation */
void compute_insert_material(uint64_t matnum);

/* insert solfec::meshes[bodnum] into computation */
void compute_insert_mesh(uint64_t bodnum);

/* delete mesh from computation */
void compute_delete_mesh(uint64_t bodnum);

/* insert solfec::ellips[bodnum] into computation */
void compute_insert_ellip(uint64_t bodnum);

/* delete ellip from computation */
void compute_delete_ellip(uint64_t bodnum);

/* insert solfec::restrains[resnum] into computation */
void compute_insert_restrain(uint64_t resnum);

/* delete restrain from computation */
void compute_delete_restrain(uint64_t resnum);

/* insert solfec::prescribes[prenum] into computation */
void compute_insert_prescribe(uint64_t prenum);

/* delete prescribe from computation */
void compute_delete_prescribe(uint64_t prenum);

/* join compute main loop */
void compute_main_loop(REAL duration, REAL step);

/* finalize compute memory */
void compute_finalize();

#endif

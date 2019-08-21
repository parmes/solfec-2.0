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

#ifndef __mesh__
#define __mesh__

/* create METIS-style all faces pointer and index vectors, and colors,
 * from elements definitions, a global color, and colors of selected faces;
 * returns [number of faces, face pointers, face indices, face colors] */
void mesh_create_metis_faces (const std::vector<uint64_t> &elements, uint64_t gcolor, const std::vector<uint64_t> &colors, /* input */
                              int64_t &nf, std::vector<int64_t> &fptr, std::vector<int64_t> &find, std::vector<uint64_t> &color); /* output */

/* calculate mass characteristics: scalar mass, mass center, euler tensor and inertia tensor */
void mesh_char (uint64_t bodnum, std::vector<uint64_t> &material, std::vector<int64_t> eptr, std::vector<int64_t> eind,
                REAL *mass, REAL *center, REAL *euler, REAL *inertia);

#endif

/*
The MIT License (MIT)

Copyright (c) 2016 Tomasz Koziara

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

#ifndef __simplex__
#define __simplex__

/* simplex integration of few
 * usefull global coordinate functions;
 * self-explanatory
 */

inline static REAL simplex_J (REAL *a, REAL *b, REAL *c, REAL *d)
{
  REAL q [9], J;

  q [0] = b [0] - a [0];
  q [1] = c [0] - a [0];
  q [2] = d [0] - a [0];
  q [3] = b [1] - a [1];
  q [4] = c [1] - a [1];
  q [5] = d [1] - a [1];
  q [6] = b [2] - a [2];
  q [7] = c [2] - a [2];
  q [8] = d [2] - a [2];

  J = q [0]*q [4]*q [8] + q [3]*q [7]*q [2] + q [6]*q [1]*q [5] -
      q [6]*q [4]*q [2] - q [0]*q [7]*q [5] - q [3]*q [1]*q [8];

  return J;
}

#define simplex_1(J, a, b, c, d) ((J)/6.)

#define simplex_x(J, a, b, c, d) (((J)/24.)*((a) [0] + (b) [0] + (c) [0] + (d) [0]))

#define simplex_y(J, a, b, c, d) (((J)/24.)*((a) [1] + (b) [1] + (c) [1] + (d) [1]))

#define simplex_z(J, a, b, c, d) (((J)/24.)*((a) [2] + (b) [2] + (c) [2] + (d) [2]))

#define simplex_xx(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [0]*(a) [0] + (a) [0]*(b) [0] + (a) [0]*(c) [0] + (a) [0]*(d) [0] +\
(b) [0]*(a) [0] + 2.*(b) [0]*(b) [0] + (b) [0]*(c) [0] + (b) [0]*(d) [0] +\
(c) [0]*(a) [0] + (c) [0]*(b) [0] + 2.*(c) [0]*(c) [0] + (c) [0]*(d) [0] +\
(d) [0]*(a) [0] + (d) [0]*(b) [0] + (d) [0]*(c) [0] + 2.*(d) [0]*(d) [0] \
))

#define simplex_xy(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [0]*(a) [1] + (a) [0]*(b) [1] + (a) [0]*(c) [1] + (a) [0]*(d) [1] +\
(b) [0]*(a) [1] + 2.*(b) [0]*(b) [1] + (b) [0]*(c) [1] + (b) [0]*(d) [1] +\
(c) [0]*(a) [1] + (c) [0]*(b) [1] + 2.*(c) [0]*(c) [1] + (c) [0]*(d) [1] +\
(d) [0]*(a) [1] + (d) [0]*(b) [1] + (d) [0]*(c) [1] + 2.*(d) [0]*(d) [1] \
))

#define simplex_xz(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [0]*(a) [2] + (a) [0]*(b) [2] + (a) [0]*(c) [2] + (a) [0]*(d) [2] +\
(b) [0]*(a) [2] + 2.*(b) [0]*(b) [2] + (b) [0]*(c) [2] + (b) [0]*(d) [2] +\
(c) [0]*(a) [2] + (c) [0]*(b) [2] + 2.*(c) [0]*(c) [2] + (c) [0]*(d) [2] +\
(d) [0]*(a) [2] + (d) [0]*(b) [2] + (d) [0]*(c) [2] + 2.*(d) [0]*(d) [2] \
))

#define simplex_yy(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [1]*(a) [1] + (a) [1]*(b) [1] + (a) [1]*(c) [1] + (a) [1]*(d) [1] +\
(b) [1]*(a) [1] + 2.*(b) [1]*(b) [1] + (b) [1]*(c) [1] + (b) [1]*(d) [1] +\
(c) [1]*(a) [1] + (c) [1]*(b) [1] + 2.*(c) [1]*(c) [1] + (c) [1]*(d) [1] +\
(d) [1]*(a) [1] + (d) [1]*(b) [1] + (d) [1]*(c) [1] + 2.*(d) [1]*(d) [1] \
))

#define simplex_yz(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [1]*(a) [2] + (a) [1]*(b) [2] + (a) [1]*(c) [2] + (a) [1]*(d) [2] +\
(b) [1]*(a) [2] + 2.*(b) [1]*(b) [2] + (b) [1]*(c) [2] + (b) [1]*(d) [2] +\
(c) [1]*(a) [2] + (c) [1]*(b) [2] + 2.*(c) [1]*(c) [2] + (c) [1]*(d) [2] +\
(d) [1]*(a) [2] + (d) [1]*(b) [2] + (d) [1]*(c) [2] + 2.*(d) [1]*(d) [2] \
))

#define simplex_zz(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [2]*(a) [2] + (a) [2]*(b) [2] + (a) [2]*(c) [2] + (a) [2]*(d) [2] +\
(b) [2]*(a) [2] + 2.*(b) [2]*(b) [2] + (b) [2]*(c) [2] + (b) [2]*(d) [2] +\
(c) [2]*(a) [2] + (c) [2]*(b) [2] + 2.*(c) [2]*(c) [2] + (c) [2]*(d) [2] +\
(d) [2]*(a) [2] + (d) [2]*(b) [2] + (d) [2]*(c) [2] + 2.*(d) [2]*(d) [2] \
))

#endif

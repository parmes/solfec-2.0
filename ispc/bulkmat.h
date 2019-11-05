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

/* Contributors: Tomasz Koziara */

#ifndef __bulkmat__
#define __bulkmat__

#define SQ(X) ((X)*(X))
#define F(I, J) F[(I-1) + 3 * (J-1)]
#define K(I, J) K[(I-1) + 9 * (J-1)]
/* Saint Venant - Kirchhoff material tangent assuming column-wise F;
 * below code were generated using Maxima's 'itensor' package; */
inline static void SVK_Tangent_C (REAL lambda, REAL mi, REAL coef, REAL F[], REAL K[])
{
  K (1, 1) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+3*SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,3))+2*SQ(F(1,2))+6*SQ(F(1,1))-2)*mi)*0.5;
  K (2, 1) = F(1,1)*F(2,1)*lambda+(F(1,3)*F(2,3)+F(1,2)*F(2,2)+2*F(1,1)*F(2,1))*mi;
  K (3, 1) = F(1,1)*F(3,1)*lambda+(F(1,3)*F(3,3)+F(1,2)*F(3,2)+2*F(1,1)*F(3,1))*mi;
  K (4, 1) = F(1,1)*F(1,2)*lambda+(F(3,1)*F(3,2)+F(2,1)*F(2,2)+2*F(1,1)*F(1,2))*mi;
  K (5, 1) = F(1,1)*F(2,2)*lambda+F(1,2)*F(2,1)*mi;
  K (6, 1) = F(1,1)*F(3,2)*lambda+F(1,2)*F(3,1)*mi;
  K (7, 1) = F(1,1)*F(1,3)*lambda+(F(3,1)*F(3,3)+F(2,1)*F(2,3)+2*F(1,1)*F(1,3))*mi;
  K (8, 1) = F(1,1)*F(2,3)*lambda+F(1,3)*F(2,1)*mi;
  K (9, 1) = F(1,1)*F(3,3)*lambda+F(1,3)*F(3,1)*mi;
  K (1, 2) = K(2,1);
  K (2, 2) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+3*SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(2,2))+6*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (3, 2) = F(2,1)*F(3,1)*lambda+(F(2,3)*F(3,3)+F(2,2)*F(3,2)+2*F(2,1)*F(3,1))*mi;
  K (4, 2) = F(1,2)*F(2,1)*lambda+F(1,1)*F(2,2)*mi;
  K (5, 2) = F(2,1)*F(2,2)*lambda+(F(3,1)*F(3,2)+2*F(2,1)*F(2,2)+F(1,1)*F(1,2))*mi;
  K (6, 2) = F(2,1)*F(3,2)*lambda+F(2,2)*F(3,1)*mi;
  K (7, 2) = F(1,3)*F(2,1)*lambda+F(1,1)*F(2,3)*mi;
  K (8, 2) = F(2,1)*F(2,3)*lambda+(F(3,1)*F(3,3)+2*F(2,1)*F(2,3)+F(1,1)*F(1,3))*mi;
  K (9, 2) = F(2,1)*F(3,3)*lambda+F(2,3)*F(3,1)*mi;
  K (1, 3) = K(3,1);
  K (2, 3) = K(3,2);
  K (3, 3) = ((SQ(F(3,3))+SQ(F(3,2))+3*SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,3))+2*SQ(F(3,2))+6*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (4, 3) = F(1,2)*F(3,1)*lambda+F(1,1)*F(3,2)*mi;
  K (5, 3) = F(2,2)*F(3,1)*lambda+F(2,1)*F(3,2)*mi;
  K (6, 3) = F(3,1)*F(3,2)*lambda+(2*F(3,1)*F(3,2)+F(2,1)*F(2,2)+F(1,1)*F(1,2))*mi;
  K (7, 3) = F(1,3)*F(3,1)*lambda+F(1,1)*F(3,3)*mi;
  K (8, 3) = F(2,3)*F(3,1)*lambda+F(2,1)*F(3,3)*mi;
  K (9, 3) = F(3,1)*F(3,3)*lambda+(2*F(3,1)*F(3,3)+F(2,1)*F(2,3)+F(1,1)*F(1,3))*mi;
  K (1, 4) = K(4,1);
  K (2, 4) = K(4,2);
  K (3, 4) = K(4,3);
  K (4, 4) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+3*SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,2))+2*SQ(F(1,3))+6*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (5, 4) = F(1,2)*F(2,2)*lambda+(F(1,3)*F(2,3)+2*F(1,2)*F(2,2)+F(1,1)*F(2,1))*mi;
  K (6, 4) = F(1,2)*F(3,2)*lambda+(F(1,3)*F(3,3)+2*F(1,2)*F(3,2)+F(1,1)*F(3,1))*mi;
  K (7, 4) = F(1,2)*F(1,3)*lambda+(F(3,2)*F(3,3)+F(2,2)*F(2,3)+2*F(1,2)*F(1,3))*mi;
  K (8, 4) = F(1,2)*F(2,3)*lambda+F(1,3)*F(2,2)*mi;
  K (9, 4) = F(1,2)*F(3,3)*lambda+F(1,3)*F(3,2)*mi;
  K (1, 5) = K(5,1);
  K (2, 5) = K(5,2);
  K (3, 5) = K(5,3);
  K (4, 5) = K(5,4);
  K (5, 5) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+3*SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,3))+6*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,2))-2)*mi)*0.5;
  K (6, 5) = F(2,2)*F(3,2)*lambda+(F(2,3)*F(3,3)+2*F(2,2)*F(3,2)+F(2,1)*F(3,1))*mi;
  K (7, 5) = F(1,3)*F(2,2)*lambda+F(1,2)*F(2,3)*mi;
  K (8, 5) = F(2,2)*F(2,3)*lambda+(F(3,2)*F(3,3)+2*F(2,2)*F(2,3)+F(1,2)*F(1,3))*mi;
  K (9, 5) = F(2,2)*F(3,3)*lambda+F(2,3)*F(3,2)*mi;
  K (1, 6) = K(6,1);
  K (2, 6) = K(6,2);
  K (3, 6) = K(6,3);
  K (4, 6) = K(6,4);
  K (5, 6) = K(6,5);
  K (6, 6) = ((SQ(F(3,3))+3*SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,3))+6*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,2))+2*SQ(F(1,2))-2)*mi)*0.5;
  K (7, 6) = F(1,3)*F(3,2)*lambda+F(1,2)*F(3,3)*mi;
  K (8, 6) = F(2,3)*F(3,2)*lambda+F(2,2)*F(3,3)*mi;
  K (9, 6) = F(3,2)*F(3,3)*lambda+(2*F(3,2)*F(3,3)+F(2,2)*F(2,3)+F(1,2)*F(1,3))*mi;
  K (1, 7) = K(7,1);
  K (2, 7) = K(7,2);
  K (3, 7) = K(7,3);
  K (4, 7) = K(7,4);
  K (5, 7) = K(7,5);
  K (6, 7) = K(7,6);
  K (7, 7) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+3*SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +  
             (2*SQ(F(3,3))+2*SQ(F(2,3))+6*SQ(F(1,3))+2*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (8, 7) = F(1,3)*F(2,3)*lambda+(2*F(1,3)*F(2,3)+F(1,2)*F(2,2)+F(1,1)*F(2,1))*mi;
  K (9, 7) = F(1,3)*F(3,3)*lambda+(2*F(1,3)*F(3,3)+F(1,2)*F(3,2)+F(1,1)*F(3,1))*mi;
  K (1, 8) = K(8,1);
  K (2, 8) = K(8,2);
  K (3, 8) = K(8,3);
  K (4, 8) = K(8,4);
  K (5, 8) = K(8,5);
  K (6, 8) = K(8,6);
  K (7, 8) = K(8,7);
  K (8, 8) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+3*SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +   
             (2*SQ(F(3,3))+6*SQ(F(2,3))+2*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,3))-2)*mi)*0.5;
  K (9, 8) = F(2,3)*F(3,3)*lambda+(2*F(2,3)*F(3,3)+F(2,2)*F(3,2)+F(2,1)*F(3,1))*mi;
  K (1, 9) = K(9,1);
  K (2, 9) = K(9,2);
  K (3, 9) = K(9,3);
  K (4, 9) = K(9,4);
  K (5, 9) = K(9,5);
  K (6, 9) = K(9,6);
  K (7, 9) = K(9,7);
  K (8, 9) = K(9,8);
  K (9, 9) = ((3*SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +    
             (6*SQ(F(3,3))+2*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(1,3))-2)*mi)*0.5;

#pragma ignore warning
  for (int i = 0; i < 81; i ++) K[i] *= coef;
}

/* Saint Venant - Kirchhoff material model; assuming column-wise F;
 * given Lame coefficients (lambda, mi) and deformation gradient 'F'
 * return det (F) end output first Piola-Kirchhoff stress 'P' (scaled by the 'coef'); */
inline static void SVK_Stress_C (REAL lambda, REAL mi, REAL coef, REAL F[], REAL P[])
{
  REAL J, E [9], S [9], trace;

  /* J = det (F) */
  J = F [0]*F [4]*F [8] + F [3]*F [7]*F [2] + F [6]*F [1]*F [5] -
      F [6]*F [4]*F [2] - F [0]*F [7]*F [5] - F [3]*F [1]*F [8];

  /* calculate Green tensor: E = (FF - 1) / 2 */
  E [0] = .5 * (F [0]*F [0] + F [1]*F [1] + F [2]*F [2] - 1.);
  E [1] = .5 * (F [3]*F [0] + F [4]*F [1] + F [5]*F [2]);
  E [2] = .5 * (F [6]*F [0] + F [7]*F [1] + F [8]*F [2]);
  E [3] = .5 * (F [0]*F [3] + F [1]*F [4] + F [2]*F [5]);
  E [4] = .5 * (F [3]*F [3] + F [4]*F [4] + F [5]*F [5] - 1.);
  E [5] = .5 * (F [6]*F [3] + F [7]*F [4] + F [8]*F [5]);
  E [6] = .5 * (F [0]*F [6] + F [1]*F [7] + F [2]*F [8]);
  E [7] = .5 * (F [3]*F [6] + F [4]*F [7] + F [5]*F [8]);
  E [8] = .5 * (F [6]*F [6] + F [7]*F [7] + F [8]*F [8] - 1.);

  /* obtain the second PK tensor
   * trough the Saint Venant - Kirchhoff law */
  trace = 2. * mi;
  S [0] = trace * E [0];
  S [1] = trace * E [1];
  S [2] = trace * E [2];
  S [3] = trace * E [3];
  S [4] = trace * E [4];
  S [5] = trace * E [5];
  S [6] = trace * E [6];
  S [7] = trace * E [7];
  S [8] = trace * E [8];

  trace = E [0] + E [4] + E [8];
  S [0] += lambda * trace;
  S [4] += lambda * trace;
  S [8] += lambda * trace;

  /* now conver S - which is the second PK tensor
   * into the first PK tensor: P = F S */
  P [0] = coef * (F [0]*S [0] + F [3]*S [1] + F [6]*S [2]);
  P [1] = coef * (F [1]*S [0] + F [4]*S [1] + F [7]*S [2]);
  P [2] = coef * (F [2]*S [0] + F [5]*S [1] + F [8]*S [2]);
  P [3] = coef * (F [0]*S [3] + F [3]*S [4] + F [6]*S [5]);
  P [4] = coef * (F [1]*S [3] + F [4]*S [4] + F [7]*S [5]);
  P [5] = coef * (F [2]*S [3] + F [5]*S [4] + F [8]*S [5]);
  P [6] = coef * (F [0]*S [6] + F [3]*S [7] + F [6]*S [8]);
  P [7] = coef * (F [1]*S [6] + F [4]*S [7] + F [7]*S [8]);
  P [8] = coef * (F [2]*S [6] + F [5]*S [7] + F [8]*S [8]);
}

inline static void BULK_MATERIAL_P (REAL young, REAL poisson, REAL F[], REAL coef, REAL P[])
{
  REAL lambda = young*poisson / ((1.0 + poisson)*(1.0 - 2.0*poisson));
  REAL mi = young / (2.0*(1.0 + poisson));

  SVK_Stress_C (lambda, mi, coef, F, P);
}

inline static void BULK_MATERIAL_K (REAL young, REAL poisson, REAL F[], REAL coef, REAL K[])
{
  REAL lambda = young*poisson / ((1.0 + poisson)*(1.0 - 2.0*poisson));
  REAL mi = young / (2.0*(1.0 + poisson));

  SVK_Tangent_C (lambda, mi, coef, F, K);
}

#endif

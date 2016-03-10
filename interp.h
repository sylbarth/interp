/*
Interpolation Library - Version 1.0 - April 1998

Copyright (c) 1998 Sylvain BARTHELEMY

Contact: sylbarth@gmail.com, www.sylbarth.com

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _INTERP_H_
#define _INTERP_H_


// --- cubic splines end condition
// iend = 1 : linear end condition, s(1)=s(n)=0
// iend = 2 : parabolic ends: s(1)=s(2),s(n-1)=s(n)
// iend = 3 : cubic end condition
// iend = 4 : first derivative is known at x(1)
//            and stored in a[1] and the first
//            derivative is know at x(n) and stored
//            in a[2].

const int CSPLINE_LINEAR    = 1;
const int CSPLINE_PARABOLIC = 2;
const int CSPLINE_CUBIC     = 3;
const int CSPLINE_FDEV      = 4;


// --- interpolations class
class interp {

  // --- cubic splines end condition
  int     m_iend;

  // --- divided differences polynome coef
  double* m_dd;

  // --- cspline polynome coef
  double* m_a;
  double* m_b;
  double* m_c;

  // --- orig data
  double* m_x;
  double* m_y;
  int     m_n;

  // --- booleans
  bool m_cspline;
  bool m_divdiff;

  void divdiff_coef();
  void cspline_coef();

public:
  interp(const double* x, const double* y, const int& n);
  ~interp();

  double lagrange ( const double& u );
  double divdiff  ( const double& u );
  double cspline  ( const double& u, int* store);

};

#endif

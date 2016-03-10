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

#include "stdio.h"
#include "interp.h"

interp::interp(const double* x, const double* y, const int& n) :
  m_n      (n),
  m_iend   (CSPLINE_LINEAR),
  m_cspline(false),
  m_divdiff(false)
{
  m_x = new double[n];
  m_y = new double[n];
  memcpy(m_x,x,sizeof(double)*n);
  memcpy(m_y,y,sizeof(double)*n);
  m_a  = new double[n];
  m_b  = new double[n];
  m_c  = new double[n];
  m_dd = new double[n];
}

interp::~interp()
{
  if(m_x)  delete[] m_x;
  if(m_y)  delete[] m_y;
  if(m_a)  delete[] m_a;
  if(m_b)  delete[] m_b;
  if(m_c)  delete[] m_c;
  if(m_dd) delete[] m_dd;
}

// -------------------------------------------------------------------
// LAGRANGE
// -------------------------------------------------------------------
double interp::lagrange(const double& u)
{
  double interp = 0.0;
  
  for(int i=0; i<m_n; i++) {
    double p = 1.0;
    for(int j=0; j<m_n; j++) {
      if(i!=j) p *= (u-m_x[j])/(m_x[i]-m_x[j]);
    }
    interp += p * m_y[i];
  }

  return interp;
}

// -------------------------------------------------------------------
// DIVIDED DIFFRERENCES
// -------------------------------------------------------------------
void interp::divdiff_coef()
{
  double temp1, temp2;

  for ( int i=0; i<m_n; i++ )
    m_dd[i] = m_y[i];

  for ( int j=1; j<m_n; j++ ) { 
    temp1 = m_dd[j-1]; 
    for ( int k=j; k<m_n; k++ ) {
      temp2 = m_dd[k];
      m_dd[k] = (m_dd[k] - temp1)/(m_x[k] - m_x[k-j]);
      temp1 = temp2;   
    }
  }
}

double interp::divdiff(const double& u)
{
  double interp = 0.0;

  if(m_divdiff == false) {
    divdiff_coef();
    m_divdiff = true;
  }

   for ( int i=(m_n-1); i>=1; i-- )
      interp = (interp + m_dd[i]) * (u - m_x[i-1]);
   interp += m_dd[0];

  return interp;
}

// -------------------------------------------------------------------
// CUBIC SPLINE
// -------------------------------------------------------------------
void interp::cspline_coef()
{
  int n    = m_n;
  int iend = m_iend;

  // --- x[0] -> x[1]
  const double* x = m_x;
  const double* y = m_y;
  --x;
  --y;

  double* a = m_a;
  double* b = m_b;
  double* c = m_c;
  --a;
  --b;
  --c;

  // --- alloc des struct
  double*  s       = new double [n+1];
  double*  h       = new double [n+1];
  double** smatrix = new double*[n+1];
  for(int k=0;k<n;k++)
    smatrix[k] = new double[5];

  // --- init des locales
  double dx1, dy1, dx2, dy2, dxn1, dxn2;
  int    nm1, nm2, i, j, first, last;

  dx1  = .0; dy1  = .0;
  dx2  = .0; dy2  = .0;
  dxn1 = .0; dxn2 = .0;

  // --- COMPUTE FOR N-2 ROWS
  nm2 = n - 2; 
  nm1 = n - 1;
  dx1 = x[2] - x[1]; 
  dy1 = (y[2]-y[1])/dx1 * 6.0;
  for ( i = 1; i <= nm2; ++i ) { 
      dx2 = x[i+2] - x[i+1];
      dy2 = (y[i+2] - y[i+1])/dx2 * 6.0;
      smatrix[i][1] = dx1;
      smatrix[i][2] = 2.0*(dx1 + dx2);
      smatrix[i][3] = dx2;
      smatrix[i][4] = dy2 - dy1;
      dx1 = dx2; 
      dy1 = dy2;
    }     /* end of for i loop */
  first = 2; 
  last  = nm2;

  // --- ADJUST FIRST AND LAST ROWS ACCORDING TO IEND CONDITIONS
  switch (iend)
    {
    case 1 : /* natural end conditions */       ; /* no change needed */
      break;
    case 2 : /* parabolic end conditions: s(1)=s(2) and s(n)=s(n-1) */
      smatrix[  1][2] = smatrix[1][2]   + x[2]-x[  1];
      smatrix[nm2][2] = smatrix[nm2][2] + x[n]-x[nm1];
      break;
    case 3 : /* cubic ends--s[1], s[n] are extrapolated */
      dx1 = x[2] - x[1]; 
      dx2 = x[3] - x[2];
      smatrix[1][2] = (dx1+dx2)*(dx1+2.0*dx2)/dx2;
      smatrix[1][3] = (dx2*dx2-dx1*dx1)/dx2;
      dxn2 = x[nm1]-x[nm2]; 
      dxn1 = x[  n]-x[nm1];
      smatrix[nm2][1] = (dxn2*dxn2-dxn1*dxn1)/dxn2;
      smatrix[nm2][2] = (dxn1+dxn2)*(dxn1+2.0*dxn2)/dxn2;
      break;
    case 4 : /* Derivative end conditions */
      dx1 = x[2]-x[1]; 
      dy1 = (y[2]-y[1])/dx1;
      s[1] = a[1]; 
      s[n] = a[2];
      smatrix[0][1] = 0.0;
      smatrix[0][2] = 2.0*dx1; 
      smatrix[0][3] = dx1;
      smatrix[0][4] = (dy1-s[1])*6;
      dx1 = x[n]-x[n-1]; 
      dy1 = (y[n]-y[n-1])/dx1;
      smatrix[nm1][1] = dx1; 
      smatrix[nm1][2] = 2.0*dx1;
      smatrix[nm1][3] = 0.0;
      smatrix[nm1][4] = (s[n]-dy1)*6.0;
      first = 1; last = n-1;  
      break;
    }

  // --- NOW WE SOLVE THE TRIDIAGONAL SYSTEM
  for ( i = first; i <= last; ++i ) {
    smatrix[i][1] = smatrix[i][1]/smatrix[i-1][2];
    smatrix[i][2] = smatrix[i][2]-smatrix[i][1]*smatrix[i-1][3];
    smatrix[i][4] = smatrix[i][4]-smatrix[i][1]*smatrix[i-1][4];
  }

  // --- NOW WE BACK SUBSTITUTE
  smatrix[last][4] = smatrix[last][4]/smatrix[last][2];
  for ( j = last-1; j >= first-1; --j )
    smatrix[j][4] = (smatrix[j][4]-smatrix[j][3]*smatrix[j+1][4])/
      smatrix[j][2];

  // --- NOW PUT THE VALUES INTO THE S VECTOR
  for ( i = first-1; i <= last; ++i )
    s[i+1] = smatrix[i][4];

  // --- TAKE CARE OF SPECIAL END CONDITIONS  FOR S VECTOR
  switch (iend) {
  case 1 :  
    s[1] = 0.0; 
    s[n] = 0.0;
    break;
  case 2 :  
    s[1] = s[2]; 
    s[n] = s[n-1];
    break;
  case 3 :  
    s[1] = ((dx1+dx2)*s[2]-dx1*s[3])/dx2;
    s[n] = ((dxn2+dxn1)*s[nm1]-dxn1*s[nm2])/dxn2;
    break;
  case 4 :  /* already taken care of */;
    break;
  }

  // --- GENERATE THE COEFFICIENTS OF THE POLYNOMIALS
  //  printf("\n\n");
  //  for ( i = 1; i <= n; ++i )
  //    printf("%.3f ",s[i]);
  for ( i = 1; i <= n-1; ++i ) {
    h[i] = x[i+1]-x[i];
    a[i] = (s[i+1]-s[i])/(6.0*h[i]);
    b[i] = s[i]/2.0;
    c[i] = ((y[i+1]-y[i])/h[i]) - ((2.0*h[i]*s[i] + h[i]*s[i+1])/6.0);
  }

  // --- free des structs
  if(s) delete[] s;
  if(h) delete[] h;
  for(int k=0;k<n;k++)
    if(smatrix[k]) delete[] smatrix[k];
  if(smatrix) delete[] smatrix;
}



// --- THIS FUNCTION EVALUATES THE SPLINE AT A GIVEN VALUE OF U
//     eval = y(i) + c(i)*(u-x(i)) + b(i)(u-x(i))**2 + a(i)(u-x(i))**3
//            on the interval [x(i), x(i+1)].
double interp::cspline(const double& u, int* store)
{
  if(m_cspline == false) {
    cspline_coef();
    m_cspline = true;
  }

  int n = m_n;
  // --- x[0] -> x[1]
  const double* x = m_x;
  const double* y = m_y;
  --x;
  --y;

  double* a = m_a;
  double* b = m_b;
  double* c = m_c;
  --a;
  --b;
  --c;

  int    j, k;
  double dx;

  if (*store > n) *store = 1;

  if (u >= x[n])  
    {
      dx = u-x[n-1];
      return(y[n-1]+dx*(c[n-1]+dx*(b[n-1]+dx*a[n-1])));
    }
  else if ((u>= x[*store]) && (u<=x[*store+1])) 
    {
      dx = u - x[*store];
      return(y[*store]+dx*(c[*store]+dx*(b[*store]+dx*a[*store])));
    }
  else
    {
      // --- WE DO A BINARY SEARCH TO FIND THE APPROPRIATE INTERVAL
      *store = 1; j = n+1;
      while (j > *store+1) 
	{
	  k = (*store+j)/2;
	  if (u<x[k]) 
	    j = k;
	  else
	    *store = k;
	}
      dx = u - x[*store];
      return(y[*store]+dx*(c[*store]+dx*(b[*store]+dx*a[*store])));
    }
}


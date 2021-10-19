/* *****************************************************************************
 *	HF_EWO_fft.c                                                                
 *   Coded by    S. Matsuda  2016/3/26                                           
 *   Modofied by F. Tsuchiya 2021/10/15 
 *                replace double with float
 *                replace app_math functions with functions defined by math.h
 *                remove Japanese charactors
 *                change function name : app_EWO_* -> HF_EWO_* or HF_app_EWO_*
 ****************************************************************************** */
#include <math.h>
#include "HF_EWO_fft.h"

// -----------------------------------------------------------------------------------------------------
// --- MACRO: local ------------------------------------------------------------------------------------
void HF_app_EWO_makewt 	(int nw, int *ip, float *w);
void HF_app_EWO_makect 	(int nc, int *ip, float *c);
void HF_app_EWO_bitrv2 	(int n, int *ip, float *a);
void HF_app_EWO_bitrv2conj	(int n, int *ip, float *a);
void HF_app_EWO_cftfsub	(int n, float *a, float *w);
void HF_app_EWO_cftbsub	(int n, float *a, float *w);
void HF_app_EWO_rftfsub	(int n, float *a, int nc, float *c);
void HF_app_EWO_rftbsub	(int n, float *a, int nc, float *c);
void HF_app_EWO_cft1st		(int n, float *a, float *w);
void HF_app_EWO_cftmdl		(int n, int l, float *a, float *w);
// -----------------------------------------------------------------------------------------------------

/* ************************************* */

// ===========================================================================
// void HF_EWO_cdft(int n, int isgn, float *a, int *ip, float *w);
//  int     n       number of input data
//  int     isgn    direction of FFT (1: forward, -1 inverse)
//  float   *a      input/output data
//  int     *ip     work area
//  float   *w      sin table
// ===========================================================================
void HF_EWO_cdft(int n, int isgn, float *a, int *ip, float *w)
{

    if (n > (ip[0] << 2)) {
        HF_app_EWO_makewt(n >> 2, ip, w);
    }
    if (n > 4) {
        if (isgn >= 0) {
            HF_app_EWO_bitrv2(n, ip + 2, a);
            HF_app_EWO_cftfsub(n, a, w);
        } else {
            HF_app_EWO_bitrv2conj(n, ip + 2, a);
            HF_app_EWO_cftbsub(n, a, w);
        }
    } else if (n == 4) {
        HF_app_EWO_cftfsub(n, a, w);
    }
}

/* ************************************* */

// ===========================================================================
// void HF_EWO_rdft(int n, int isgn, float *a, int *ip, float *w);
//  int     n       number of input data
//  int     isgn    direction of FFT (1: forward, -1 inverse)
//  float   *a      input/output data
//  int     *ip     work area
//  float   *w      sin table
// ===========================================================================
void 	HF_EWO_rdft( int n, int isgn, float *a, int *ip, float *w ){
  	int 	nw, nc;
  	float 	xi;
  
  	nw = ip[0];
  	if (n > (nw << 2)) {
    	nw = n >> 2;
    	HF_app_EWO_makewt(nw, ip, w);
  	}
  	nc = ip[1];
  	if (n > (nc << 2)) {
    	nc = n >> 2;
    	HF_app_EWO_makect(nc, ip, w + nw);
  	}
  	if (isgn >= 0) {
    	if (n > 4) {
    		HF_app_EWO_bitrv2(n, ip + 2, a);
      		HF_app_EWO_cftfsub(n, a, w);
      		HF_app_EWO_rftfsub(n, a, nc, w + nw);
    	} else if (n == 4) {
      		HF_app_EWO_cftfsub(n, a, w);
    	}
 	xi = a[0] - a[1];
    a[0] += a[1];
    a[1] = xi;
  } else {
    a[1] = 0.5 * (a[0] - a[1]);
    a[0] -= a[1];
    if (n > 4) {
      HF_app_EWO_rftbsub(n, a, nc, w + nw);
      HF_app_EWO_bitrv2(n, ip + 2, a);
      HF_app_EWO_cftbsub(n, a, w);
    } else if (n == 4) {
      HF_app_EWO_cftfsub(n, a, w);
    }
  }
}

/* -------- initializing routines -------- */
void	HF_app_EWO_makewt(int nw, int *ip, float *w)
{
  void HF_app_EWO_bitrv2(int n, int *ip, float *a);
  int j, nwh;
  float delta, x, y;
  
  ip[0] = nw;
  ip[1] = 1;
  if (nw > 2) {
    nwh = nw >> 1;
    delta = atan(1.0) / nwh;
    w[0] = 1;
    w[1] = 0;
    w[nwh] = cos(delta * nwh);
    w[nwh + 1] = w[nwh];
    if (nwh > 2) {
      for (j = 2; j < nwh; j += 2) {
	      x = cos(delta * j);
	      y = sin(delta * j);
	      w[j] = x;
	      w[j + 1] = y;
	      w[nw - j] = y;
	      w[nw - j + 1] = x;
      }
      HF_app_EWO_bitrv2(nw, ip + 2, w);
    }
  }
}

void	HF_app_EWO_makect(int nc, int *ip, float *c)
{
  int j, nch;
  float delta;
  
  ip[1] = nc;
  if (nc > 1) {
    nch = nc >> 1;
    delta = atan(1.0) / nch;
    c[0] = cos(delta * nch);
    c[nch] = 0.5 * c[0];
    for (j = 1; j < nch; j++) {
      c[j] = 0.5 * cos(delta * j);
      c[nc - j] = 0.5 * sin(delta * j);
    }
  }
}


/* -------- child routines -------- */


void	HF_app_EWO_bitrv2(int n, int *ip, float *a)
{
  int j, j1, k, k1, l, m, m2;
  float xr, xi, yr, yi;
  
  ip[0] = 0;
  l = n;
  m = 1;
  while ((m << 3) < l) {
    l >>= 1;
    for (j = 0; j < m; j++) {
      ip[m + j] = ip[j] + l;
    }
    m <<= 1;
  }
  m2 = 2 * m;
  if ((m << 3) == l) {
    for (k = 0; k < m; k++) {
      for (j = 0; j < k; j++) {
	j1 = 2 * j + ip[k];
	k1 = 2 * k + ip[j];
	xr = a[j1];
	xi = a[j1 + 1];
	yr = a[k1];
	yi = a[k1 + 1];
	a[j1] = yr;
	a[j1 + 1] = yi;
	a[k1] = xr;
	a[k1 + 1] = xi;
	j1 += m2;
	k1 += 2 * m2;
	xr = a[j1];
	xi = a[j1 + 1];
	yr = a[k1];
	yi = a[k1 + 1];
	a[j1] = yr;
	a[j1 + 1] = yi;
	a[k1] = xr;
	a[k1 + 1] = xi;
	j1 += m2;
	k1 -= m2;
	xr = a[j1];
	xi = a[j1 + 1];
	yr = a[k1];
	yi = a[k1 + 1];
	a[j1] = yr;
	a[j1 + 1] = yi;
	a[k1] = xr;
	a[k1 + 1] = xi;
	j1 += m2;
	k1 += 2 * m2;
	xr = a[j1];
	xi = a[j1 + 1];
	yr = a[k1];
	yi = a[k1 + 1];
	a[j1] = yr;
	a[j1 + 1] = yi;
	a[k1] = xr;
	a[k1 + 1] = xi;
      }
      j1 = 2 * k + m2 + ip[k];
      k1 = j1 + m2;
      xr = a[j1];
      xi = a[j1 + 1];
      yr = a[k1];
      yi = a[k1 + 1];
      a[j1] = yr;
      a[j1 + 1] = yi;
      a[k1] = xr;
      a[k1 + 1] = xi;
    }
  } else {
    for (k = 1; k < m; k++) {
      for (j = 0; j < k; j++) {
	j1 = 2 * j + ip[k];
	k1 = 2 * k + ip[j];
	xr = a[j1];
	xi = a[j1 + 1];
	yr = a[k1];
	yi = a[k1 + 1];
	a[j1] = yr;
	a[j1 + 1] = yi;
	a[k1] = xr;
	a[k1 + 1] = xi;
	j1 += m2;
	k1 += m2;
	xr = a[j1];
	xi = a[j1 + 1];
	yr = a[k1];
	yi = a[k1 + 1];
	a[j1] = yr;
	a[j1 + 1] = yi;
	a[k1] = xr;
	a[k1 + 1] = xi;
      }
    }
  }
}

void HF_app_EWO_bitrv2conj(int n, int *ip, float *a)
{
    int j, j1, k, k1, l, m, m2;
    float xr, xi, yr, yi;
    
    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 3) < l) {
        l >>= 1;
        for (j = 0; j < m; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    m2 = 2 * m;
    if ((m << 3) == l) {
        for (k = 0; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 -= m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
            k1 = 2 * k + ip[k];
            a[k1 + 1] = -a[k1 + 1];
            j1 = k1 + m2;
            k1 = j1 + m2;
            xr = a[j1];
            xi = -a[j1 + 1];
            yr = a[k1];
            yi = -a[k1 + 1];
            a[j1] = yr;
            a[j1 + 1] = yi;
            a[k1] = xr;
            a[k1 + 1] = xi;
            k1 += m2;
            a[k1 + 1] = -a[k1 + 1];
        }
    } else {
        a[1] = -a[1];
        a[m2 + 1] = -a[m2 + 1];
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
            k1 = 2 * k + ip[k];
            a[k1 + 1] = -a[k1 + 1];
            a[k1 + m2 + 1] = -a[k1 + m2 + 1];
        }
    }
}

void 	HF_app_EWO_cftfsub(int n, float *a, float *w)
{
  void HF_app_EWO_cft1st(int n, float *a, float *w);
  void HF_app_EWO_cftmdl(int n, int l, float *a, float *w);
  int j, j1, j2, j3, l;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
  
  l = 2;
  if (n > 8) {
    HF_app_EWO_cft1st(n, a, w);
    l = 8;
    while ((l << 2) < n) {
      HF_app_EWO_cftmdl(n, l, a, w);
      l <<= 2;
    }
  }
  if ((l << 2) == n) {
    for (j = 0; j < l; j += 2) {
      j1 = j + l;
      j2 = j1 + l;
      j3 = j2 + l;
      x0r = a[j] + a[j1];
      x0i = a[j + 1] + a[j1 + 1];
      x1r = a[j] - a[j1];
      x1i = a[j + 1] - a[j1 + 1];
      x2r = a[j2] + a[j3];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2] - a[j3];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j] = x0r + x2r;
      a[j + 1] = x0i + x2i;
      a[j2] = x0r - x2r;
      a[j2 + 1] = x0i - x2i;
      a[j1] = x1r - x3i;
      a[j1 + 1] = x1i + x3r;
      a[j3] = x1r + x3i;
      a[j3 + 1] = x1i - x3r;
    }
  } else {
    for (j = 0; j < l; j += 2) {
      j1 = j + l;
      x0r = a[j] - a[j1];
      x0i = a[j + 1] - a[j1 + 1];
      a[j] += a[j1];
      a[j + 1] += a[j1 + 1];
      a[j1] = x0r;
      a[j1 + 1] = x0i;
    }
  }
}


void	HF_app_EWO_cftbsub(int n, float *a, float *w)
{
  int j, j1, j2, j3, l;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
  
  l = 2;
  if (n > 8) {
    HF_app_EWO_cft1st(n, a, w);
    l = 8;
    while ((l << 2) < n) {
      HF_app_EWO_cftmdl(n, l, a, w);
      l <<= 2;
    }
  }
  if ((l << 2) == n) {
    for (j = 0; j < l; j += 2) {
      j1 = j + l;
      j2 = j1 + l;
      j3 = j2 + l;
      x0r = a[j] + a[j1];
      x0i = -a[j + 1] - a[j1 + 1];
      x1r = a[j] - a[j1];
      x1i = -a[j + 1] + a[j1 + 1];
      x2r = a[j2] + a[j3];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2] - a[j3];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j] = x0r + x2r;
      a[j + 1] = x0i - x2i;
      a[j2] = x0r - x2r;
      a[j2 + 1] = x0i + x2i;
      a[j1] = x1r - x3i;
      a[j1 + 1] = x1i - x3r;
      a[j3] = x1r + x3i;
      a[j3 + 1] = x1i + x3r;
    }
  } else {
    for (j = 0; j < l; j += 2) {
      j1 = j + l;
      x0r = a[j] - a[j1];
      x0i = -a[j + 1] + a[j1 + 1];
      a[j] += a[j1];
      a[j + 1] = -a[j + 1] - a[j1 + 1];
      a[j1] = x0r;
      a[j1 + 1] = x0i;
    }
  }
}


void
HF_app_EWO_cft1st(int n, float *a, float *w)
{
  int j, k1, k2;
  float wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
  
  x0r = a[0] + a[2];
  x0i = a[1] + a[3];
  x1r = a[0] - a[2];
  x1i = a[1] - a[3];
  x2r = a[4] + a[6];
  x2i = a[5] + a[7];
  x3r = a[4] - a[6];
  x3i = a[5] - a[7];
  a[0] = x0r + x2r;
  a[1] = x0i + x2i;
  a[4] = x0r - x2r;
  a[5] = x0i - x2i;
  a[2] = x1r - x3i;
  a[3] = x1i + x3r;
  a[6] = x1r + x3i;
  a[7] = x1i - x3r;
  wk1r = w[2];
  x0r = a[8] + a[10];
  x0i = a[9] + a[11];
  x1r = a[8] - a[10];
  x1i = a[9] - a[11];
  x2r = a[12] + a[14];
  x2i = a[13] + a[15];
  x3r = a[12] - a[14];
  x3i = a[13] - a[15];
  a[8] = x0r + x2r;
  a[9] = x0i + x2i;
  a[12] = x2i - x0i;
  a[13] = x0r - x2r;
  x0r = x1r - x3i;
  x0i = x1i + x3r;
  a[10] = wk1r * (x0r - x0i);
  a[11] = wk1r * (x0r + x0i);
  x0r = x3i + x1r;
  x0i = x3r - x1i;
  a[14] = wk1r * (x0i - x0r);
  a[15] = wk1r * (x0i + x0r);
  k1 = 0;
  for (j = 16; j < n; j += 16) {
    k1 += 2;
    k2 = 2 * k1;
    wk2r = w[k1];
    wk2i = w[k1 + 1];
    wk1r = w[k2];
    wk1i = w[k2 + 1];
    wk3r = wk1r - 2 * wk2i * wk1i;
    wk3i = 2 * wk2i * wk1r - wk1i;
    x0r = a[j] + a[j + 2];
    x0i = a[j + 1] + a[j + 3];
    x1r = a[j] - a[j + 2];
    x1i = a[j + 1] - a[j + 3];
    x2r = a[j + 4] + a[j + 6];
    x2i = a[j + 5] + a[j + 7];
    x3r = a[j + 4] - a[j + 6];
    x3i = a[j + 5] - a[j + 7];
    a[j] = x0r + x2r;
    a[j + 1] = x0i + x2i;
    x0r -= x2r;
    x0i -= x2i;
    a[j + 4] = wk2r * x0r - wk2i * x0i;
    a[j + 5] = wk2r * x0i + wk2i * x0r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j + 2] = wk1r * x0r - wk1i * x0i;
    a[j + 3] = wk1r * x0i + wk1i * x0r;
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j + 6] = wk3r * x0r - wk3i * x0i;
    a[j + 7] = wk3r * x0i + wk3i * x0r;
    wk1r = w[k2 + 2];
    wk1i = w[k2 + 3];
    wk3r = wk1r - 2 * wk2r * wk1i;
    wk3i = 2 * wk2r * wk1r - wk1i;
    x0r = a[j + 8] + a[j + 10];
    x0i = a[j + 9] + a[j + 11];
    x1r = a[j + 8] - a[j + 10];
    x1i = a[j + 9] - a[j + 11];
    x2r = a[j + 12] + a[j + 14];
    x2i = a[j + 13] + a[j + 15];
    x3r = a[j + 12] - a[j + 14];
    x3i = a[j + 13] - a[j + 15];
    a[j + 8] = x0r + x2r;
    a[j + 9] = x0i + x2i;
    x0r -= x2r;
    x0i -= x2i;
    a[j + 12] = -wk2i * x0r - wk2r * x0i;
    a[j + 13] = -wk2i * x0i + wk2r * x0r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j + 10] = wk1r * x0r - wk1i * x0i;
    a[j + 11] = wk1r * x0i + wk1i * x0r;
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j + 14] = wk3r * x0r - wk3i * x0i;
    a[j + 15] = wk3r * x0i + wk3i * x0r;
  }
}


void 
HF_app_EWO_cftmdl(int n, int l, float *a, float *w)
{
  int j, j1, j2, j3, k, k1, k2, m, m2;
  float wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
  
  m = l << 2;
  for (j = 0; j < l; j += 2) {
    j1 = j + l;
    j2 = j1 + l;
    j3 = j2 + l;
    x0r = a[j] + a[j1];
    x0i = a[j + 1] + a[j1 + 1];
    x1r = a[j] - a[j1];
    x1i = a[j + 1] - a[j1 + 1];
    x2r = a[j2] + a[j3];
    x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2] - a[j3];
    x3i = a[j2 + 1] - a[j3 + 1];
    a[j] = x0r + x2r;
    a[j + 1] = x0i + x2i;
    a[j2] = x0r - x2r;
    a[j2 + 1] = x0i - x2i;
    a[j1] = x1r - x3i;
    a[j1 + 1] = x1i + x3r;
    a[j3] = x1r + x3i;
    a[j3 + 1] = x1i - x3r;
  }
  wk1r = w[2];
  for (j = m; j < l + m; j += 2) {
    j1 = j + l;
    j2 = j1 + l;
    j3 = j2 + l;
    x0r = a[j] + a[j1];
    x0i = a[j + 1] + a[j1 + 1];
    x1r = a[j] - a[j1];
    x1i = a[j + 1] - a[j1 + 1];
    x2r = a[j2] + a[j3];
    x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2] - a[j3];
    x3i = a[j2 + 1] - a[j3 + 1];
    a[j] = x0r + x2r;
    a[j + 1] = x0i + x2i;
    a[j2] = x2i - x0i;
    a[j2 + 1] = x0r - x2r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j1] = wk1r * (x0r - x0i);
    a[j1 + 1] = wk1r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    a[j3] = wk1r * (x0i - x0r);
    a[j3 + 1] = wk1r * (x0i + x0r);
  }
  k1 = 0;
  m2 = 2 * m;
  for (k = m2; k < n; k += m2) {
    k1 += 2;
    k2 = 2 * k1;
    wk2r = w[k1];
    wk2i = w[k1 + 1];
    wk1r = w[k2];
    wk1i = w[k2 + 1];
    wk3r = wk1r - 2 * wk2i * wk1i;
    wk3i = 2 * wk2i * wk1r - wk1i;
    for (j = k; j < l + k; j += 2) {
      j1 = j + l;
      j2 = j1 + l;
      j3 = j2 + l;
      x0r = a[j] + a[j1];
      x0i = a[j + 1] + a[j1 + 1];
      x1r = a[j] - a[j1];
      x1i = a[j + 1] - a[j1 + 1];
      x2r = a[j2] + a[j3];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2] - a[j3];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j] = x0r + x2r;
      a[j + 1] = x0i + x2i;
      x0r -= x2r;
      x0i -= x2i;
      a[j2] = wk2r * x0r - wk2i * x0i;
      a[j2 + 1] = wk2r * x0i + wk2i * x0r;
      x0r = x1r - x3i;
      x0i = x1i + x3r;
      a[j1] = wk1r * x0r - wk1i * x0i;
      a[j1 + 1] = wk1r * x0i + wk1i * x0r;
      x0r = x1r + x3i;
      x0i = x1i - x3r;
      a[j3] = wk3r * x0r - wk3i * x0i;
      a[j3 + 1] = wk3r * x0i + wk3i * x0r;
    }
    wk1r = w[k2 + 2];
    wk1i = w[k2 + 3];
    wk3r = wk1r - 2 * wk2r * wk1i;
    wk3i = 2 * wk2r * wk1r - wk1i;
    for (j = k + m; j < l + (k + m); j += 2) {
      j1 = j + l;
      j2 = j1 + l;
      j3 = j2 + l;
      x0r = a[j] + a[j1];
      x0i = a[j + 1] + a[j1 + 1];
      x1r = a[j] - a[j1];
      x1i = a[j + 1] - a[j1 + 1];
      x2r = a[j2] + a[j3];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2] - a[j3];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j] = x0r + x2r;
      a[j + 1] = x0i + x2i;
      x0r -= x2r;
      x0i -= x2i;
      a[j2] = -wk2i * x0r - wk2r * x0i;
      a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
      x0r = x1r - x3i;
      x0i = x1i + x3r;
      a[j1] = wk1r * x0r - wk1i * x0i;
      a[j1 + 1] = wk1r * x0i + wk1i * x0r;
      x0r = x1r + x3i;
      x0i = x1i - x3r;
      a[j3] = wk3r * x0r - wk3i * x0i;//no 
      a[j3 + 1] = wk3r * x0i + wk3i * x0r;//no
    }
  }
}


void 
HF_app_EWO_rftfsub(int n, float *a, int nc, float *c)
{
  int j, k, kk, ks, m;
  float wkr, wki, xr, xi, yr, yi;
  
  m = n >> 1;
  ks = 2 * nc / m;
  kk = 0;
  for (j = 2; j < m; j += 2) {
    k = n - j;
    kk += ks;
    wkr = 0.5 - c[nc - kk];
    wki = c[kk];
    xr = a[j] - a[k];
    xi = a[j + 1] + a[k + 1];
    yr = wkr * xr - wki * xi;
    yi = wkr * xi + wki * xr;
    a[j] -= yr;
    a[j + 1] -= yi;
    a[k] += yr;
    a[k + 1] -= yi;
  }
}


void	HF_app_EWO_rftbsub(int n, float *a, int nc, float *c)
{
  int j, k, kk, ks, m;
  float wkr, wki, xr, xi, yr, yi;
  
  a[1] = -a[1];
  m = n >> 1;
  ks = 2 * nc / m;
  kk = 0;
  for (j = 2; j < m; j += 2) {
    k = n - j;
    kk += ks;
    wkr = 0.5 - c[nc - kk];
    wki = c[kk];
    xr = a[j] - a[k];
    xi = a[j + 1] + a[k + 1];
    yr = wkr * xr + wki * xi;
    yi = wkr * xi - wki * xr;
    a[j] -= yr;
    a[j + 1] = yi - a[j + 1];
    a[k] += yr;
    a[k + 1] = yi - a[k + 1];
  }
  a[m + 1] = -a[m + 1];
}


void 
HF_app_EWO_dctsub(int n, float *a, int nc, float *c)
{
  int j, k, kk, ks, m;
  float wkr, wki, xr;
  
  m = n >> 1;
  ks = nc / n;
  kk = 0;
  for (j = 1; j < m; j++) {
    k = n - j;
    kk += ks;
    wkr = c[kk] - c[nc - kk];
    wki = c[kk] + c[nc - kk];
    xr = wki * a[j] - wkr * a[k];
    a[j] = wkr * a[j] + wki * a[k];
    a[k] = xr;
  }
  a[m] *= c[0];
}


// ===========================================================================
// HF_EWO_Hannig		**** Hanning function: 0.5-0.5 cos(2 pi x) [x=0-1] ****
//		[INPUT]		float					f_data[]			Original data
//				 	float 					f_work[]			Return   data
//					int 	  				i_num				Num of data
// ===========================================================================
void	HF_EWO_Hanning( 	float	f_data[],
						 	float 	f_work[],
							int 	i_num )
{
	int		i;
	float	f_ave, f_tmp;
	
  // remove bias
	f_ave = 0;
  for (i = 0; i < i_num; i++)
    f_ave += f_data[i];
  f_ave /= i_num;
	
	for (i=0; i < i_num; i++) {
		f_tmp     =  cos( 6.2831852 * (float)i/ ((float)i_num-1.0) );
		f_tmp     =  0.5 - (0.5 * f_tmp);
		f_work[i] = (f_data[i] - f_ave) * f_tmp;
	}
}

// ===========================================================================
// HF_EWO_Hammig		**** Hamming function: 0.54-0.46 cos(2 pi x) [x=0-1] ****
//		[INPUT]		float					f_data[]			Original data
//				 	float 					f_work[]			Return   data
//					int 	  				i_num				Num of data
// ===========================================================================
void	HF_EWO_Hamming( 	float	f_data[],
						 	float 	f_work[],
							int 	i_num )
{
	int		i;
	float	f_ave, f_tmp;
	
  // remove bias
	f_ave = 0;
  for (i = 0; i < i_num; i++)
    f_ave += f_data[i];
  f_ave /= (float)i_num;
	
	for (i=0; i < i_num; i++) {
		f_tmp     =  cos( 6.2831852 * (float)i/ ((float)i_num-1.0) );
		f_tmp     =  0.54 - (0.46 * f_tmp);
		f_work[i] = (f_data[i] - f_ave) * f_tmp;
	}
}

// ===========================================================================
// HF_EWO_Blackman		**** Blackman function: 0.42 - 0.5 cos(2 pi x) + 0.08 cos(4 pi x) ****
//		[INPUT]		float					f_data[]			Original & Return data
//					int 	  				i_num				Num of data
// ===========================================================================
void	HF_EWO_Blackman( 	float	f_data[],
							int 	i_num )
{
	int		i;
	float  f_tmp, f_tmp2;
	
	for (i=0; i < i_num; i++) {
		f_tmp     =  cos( 6.2831852  * (float)i/ ((float)i_num-1.0) );
		f_tmp2    =  cos( 12.5663704 * (float)i/ ((float)i_num-1.0) );
		f_tmp     =  0.42 - 0.50 * f_tmp + 0.08 * f_tmp2;
		f_data[i] = f_data[i] * f_tmp;
	}
}


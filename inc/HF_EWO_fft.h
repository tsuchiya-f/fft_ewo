/*******************************************************************************
 *	HF_EWO_fft.h
 *   Coded by    S. Matsuda  2016/3/26
 *   Modofied by F. Tsuchiya 2021/10/15
 *******************************************************************************/
#ifndef _HF_EWO_FFT_H_
#define _HF_EWO_FFT_H_

// ===========================================================================
// void HF_EWO_rdft(int n, int isgn, float *a, int *ip, float *w);
//  int     n       number of input data
//  int     isgn    direction of FFT (1: forward, -1 inverse)
//  float   *a      input/output data
//  int     *ip     work area
//  float   *w      sin table
/*
-------- Real DFT / Inverse of Real DFT --------
    [definition]
        <case1> RDFT
            R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
            I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
        <case2> IRDFT (excluding scale)
            a[k] = (R[0] + R[n/2]*cos(pi*k))/2 +
                   sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) +
                   sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            rdft(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            rdft(n, -1, a, ip, w);
    [parameters]
        n              :data length (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (double *)
                        <case1>
                            output data
                                a[2*k] = R[k], 0<=k<n/2
                                a[2*k+1] = I[k], 0<k<n/2
                                a[1] = R[n/2]
                        <case2>
                            input data
                                a[2*j] = R[j], 0<=j<n/2
                                a[2*j+1] = I[j], 0<j<n/2
                                a[1] = R[n/2]
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/2)
                        strictly,
                        length of ip >=
                            2+(1<<(int)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of
            rdft(n, 1, a, ip, w);
        is
            rdft(n, -1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .
*/
void HF_EWO_rdft(int n, int isgn, float *a, int *ip, float *w);

// ===========================================================================

// ===========================================================================
// void HF_EWO_cdft(int n, int isgn, float *a, int *ip, float *w);
//  int     n       number of input data
//  int     isgn    direction of FFT (1: forward, -1 inverse)
//  float   *a      input/output data
//  int     *ip     work area
//  float   *w      sin table
/*
-------- Complex DFT (Discrete Fourier Transform) --------
    [definition]
        <case1>
            X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
        <case2>
            X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n
        (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
    [usage]
        <case1>
            ip[0] = 0; // first time only
            cdft(2*n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            cdft(2*n, -1, a, ip, w);
    [parameters]
        2*n            :data length (int)
                        n >= 1, n = power of 2
        a[0...2*n-1]   :input/output data (double *)
                        input data
                            a[2*j] = Re(x[j]),
                            a[2*j+1] = Im(x[j]), 0<=j<n
                        output data
                            a[2*k] = Re(X[k]),
                            a[2*k+1] = Im(X[k]), 0<=k<n
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n)
                        strictly,
                        length of ip >=
                            2+(1<<(int)(log(n+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (double *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of
            cdft(2*n, -1, a, ip, w);
        is
            cdft(2*n, 1, a, ip, w);
            for (j = 0; j <= 2 * n - 1; j++) {
                a[j] *= 1.0 / n;
            }
        .
*/
// ===========================================================================
void HF_EWO_cdft(int n, int isgn, float *a, int *ip, float *w);

// ===========================================================================
// HF_EWO_Hannig		**** Hanning function: 0.5-0.5 cos(2 pi x) [x=0-1] ****
//		[INPUT]		float					f_data[]			Original data
//				 	float 					f_work[]			Return   data
//					int 	  				i_num				Num of data
// ===========================================================================
void HF_EWO_Hanning(float[], float[], int);

// ===========================================================================
// HF_EWO_Hammig		**** Hamming function: 0.54-0.46 cos(2 pi x) [x=0-1] ****
//		[INPUT]		float					s_data[]			Original data
//				 	float 					f_work[]			Return   data
//					int 	  				i_num				Num of data
// ===========================================================================
void HF_EWO_Hamming(float[], float[], int);

// ===========================================================================
// HF_EWO_Blackman		**** Blackman function: 0.42 - 0.5 cos(2 pi x) + 0.08 cos(4 pi x) ****
//		[INPUT]		float					s_data[]			Original & Output data
//					int 	  				i_num				Num of data
// ===========================================================================
void HF_EWO_Blackman(float[], int);

#endif /* _HF_EWO_FFT_H_ */

/*

calculation of covariance and correlation for pairs of vectors in N-D matrices

FORMAT:       [cv, r] = cov_nd(X, Y [, lag])

Input fields:

      X, Y        N-D numeric matrices of same size
      lag         optional lag, useful for autocorrelation cov_nd(d, d, 1)

Output fields:

      cv          covariance matrix
      r           correlation coefficients matrix

Note: cv and r are computed over the LAST dimension, so if X and Y
are 10-by-30-by-50 matrices, cv and r are 10-by-30 matrices.

% Version:  v0.9d
% Build:    14062016
% Date:     Jun-20 2014, 4:12 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, Jochen Weber
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "mex.h"
#include "math.h"

#define COV_1D_LOOP(P1, P2)							\
	if (!docr) {								\
		for(c = 0; c < ts; c += rs) {					\
			m1 = 0.0;						\
			m2 = 0.0;						\
			mp = 0.0;						\
			for(ic = c, ec = c + nl; ic < ec; ++ic) {		\
				m1 += (double) P1[ic];				\
				m2 += (double) P2[ic];				\
				mp += ((double) P1[ic]) * ((double) P2[ic]);	\
			}							\
			m1 /= dnl;						\
			m2 /= dnl;						\
			mp /= dnl;						\
			*t1++ = (double) ((mp - m1 * m2) * (dnl / dnx));	\
		}								\
	} else {								\
		plhs[1] = mxCreateNumericArray(rd, dr, mxDOUBLE_CLASS, mxREAL);	\
		t2 = (double *) mxGetPr(plhs[1]);				\
		for(c = 0;c < ts; c += rs) {					\
			m1 = 0.0;						\
			m2 = 0.0;						\
			x1 = 0.0;						\
			x2 = 0.0;						\
			mp = 0.0;						\
			for(ic = c, ec = c + nl; ic < ec; ++ic) {		\
				m1 += (v1 = (double) P1[ic]);			\
				m2 += (v2 = (double) P2[ic]);			\
				x1 += v1 * v1;					\
				x2 += v2 * v2;					\
				mp += v1 * v2;					\
			}							\
			m1 /= dnl;						\
			m2 /= dnl;						\
			x1 /= dnl;						\
			x2 /= dnl;						\
			mp /= dnl;						\
			cv = (mp - m1 * m2) * (dnl / dnx);			\
			*t1++ = (double) cv;					\
			cv /= sqrt((x1 - m1 * m1) * (x2 - m2 * m2)) * (dnl / dnx);	\
			if (cv < -1.0)						\
				cv = -1.0;					\
			else if (cv > 1.0)					\
				cv = 1.0;					\
			*t2++ = (double) cv;					\
		}								\
	}

#define COV_1D_IILOOP(P1, P2, C2)						\
	P2 = (const C2 *) mxGetPr(prhs[1]);					\
	if (lag > 0)								\
		P2 = &P2[lag];							\
	COV_1D_LOOP(P1, P2)							\
	break;

#define COV_1D_ILOOP(P1, C1)							\
	P1 = (const C1 *) mxGetPr(*prhs);					\
	switch (mxGetClassID(prhs[1])) {					\
		case mxDOUBLE_CLASS: COV_1D_IILOOP(P1, idbl2, double)		\
		case mxSINGLE_CLASS: COV_1D_IILOOP(P1, isng2, float)		\
		case  mxINT32_CLASS: COV_1D_IILOOP(P1, ii322, signed int)	\
		case mxUINT32_CLASS: COV_1D_IILOOP(P1, iu322, unsigned int)	\
		case  mxINT16_CLASS: COV_1D_IILOOP(P1, ii162, signed short)	\
		case mxUINT16_CLASS: COV_1D_IILOOP(P1, iu162, unsigned short)	\
		case   mxINT8_CLASS: COV_1D_IILOOP(P1, ii082, signed char)	\
		case  mxUINT8_CLASS:						\
		case mxLOGICAL_CLASS:COV_1D_IILOOP(P1, iu082, unsigned char)	\
		default: mexErrMsgTxt("Datatype of second argument not supported.");	\
			 break;							\
	}									\
	break;

#define COV_ND_LOOP(P1, P2)							\
	if (!docr) {								\
		for(c = 0; c < rs; ++c) {					\
			m1 = 0.0;						\
			m2 = 0.0;						\
			mp = 0.0;						\
			for(ic = c; ic < ts; ic += rs) {			\
				m1 += (double) P1[ic];				\
				m2 += (double) P2[ic];				\
				mp += ((double) P1[ic]) * ((double) P2[ic]);	\
			}							\
			m1 /= dnl;						\
			m2 /= dnl;						\
			mp /= dnl;						\
			*t1++ = (double) ((mp - m1 * m2) * (dnl / dnx));	\
		}								\
	} else {								\
		plhs[1] = mxCreateNumericArray(rd, dr, mxDOUBLE_CLASS, mxREAL);	\
		t2 = (double *) mxGetPr(plhs[1]);				\
		for(c = 0;c < rs; ++c) {					\
			m1 = 0.0;						\
			m2 = 0.0;						\
			x1 = 0.0;						\
			x2 = 0.0;						\
			mp = 0.0;						\
			for(ic = c; ic < ts; ic += rs) {			\
				m1 += (v1 = (double) P1[ic]);			\
				m2 += (v2 = (double) P2[ic]);			\
				x1 += v1 * v1;					\
				x2 += v2 * v2;					\
				mp += v1 * v2;					\
			}							\
			m1 /= dnl;						\
			m2 /= dnl;						\
			x1 /= dnl;						\
			x2 /= dnl;						\
			mp /= dnl;						\
			cv = (mp - m1 * m2) * (dnl / dnx);			\
			*t1++ = (double) cv;					\
			cv /= sqrt((x1 - m1 * m1) * (x2 - m2 * m2)) * (dnl / dnx);	\
			if (cv < -1.0)						\
				cv = -1.0;					\
			else if (cv > 1.0)					\
				cv = 1.0;					\
			*t2++ = (double) cv;					\
		}								\
	}

#define COV_ND_IILOOP(P1, P2, C2)						\
	P2 = (const C2 *) mxGetPr(prhs[1]);					\
	if (lag > 0)								\
		P2 = &P2[rs * lag];						\
	COV_ND_LOOP(P1, P2)							\
	break;

#define COV_ND_ILOOP(P1, C1)							\
	P1 = (const C1 *) mxGetPr(*prhs);					\
	switch (mxGetClassID(prhs[1])) {					\
		case mxDOUBLE_CLASS: COV_ND_IILOOP(P1, idbl2, double)		\
		case mxSINGLE_CLASS: COV_ND_IILOOP(P1, isng2, float)		\
		case  mxINT32_CLASS: COV_ND_IILOOP(P1, ii322, signed int)	\
		case mxUINT32_CLASS: COV_ND_IILOOP(P1, iu322, unsigned int)	\
		case  mxINT16_CLASS: COV_ND_IILOOP(P1, ii162, signed short)	\
		case mxUINT16_CLASS: COV_ND_IILOOP(P1, iu162, unsigned short)   \
		case   mxINT8_CLASS: COV_ND_IILOOP(P1, ii082, signed char)	\
		case  mxUINT8_CLASS:						\
		case mxLOGICAL_CLASS:COV_ND_IILOOP(P1, iu082, unsigned char)	\
		default: mexErrMsgTxt("Datatype of second argument not supported.");	\
		break;								\
	}									\
	break;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int docr, c, ec, ic, nd, nl, dr[64], rd, rs, ts, lag;
	const int *d1, *d2;
	double dnl, dnx, *t1, *t2;
	const double *idbl1, *idbl2;
	const float *isng1, *isng2;
	const signed int *ii321, *ii322;
	const unsigned int *iu321, *iu322;
	const signed short *ii161, *ii162;
	const unsigned short *iu161, *iu162;
	const signed char *ii081, *ii082;
	const unsigned char *iu081, *iu082;
	long double m1, m2, x1, x2, v1, v2, cv;
	double mp;

	if (nrhs < 2 || nlhs > 2)
		mexErrMsgTxt("Bad number of input/output arguments.");

	if (nlhs > 1)
		docr = 1;
	else
		docr = 0;

	nd = mxGetNumberOfDimensions(prhs[0]);
	if (nd != mxGetNumberOfDimensions(prhs[1]) || nd > 64)
		mexErrMsgTxt("Matrices must have the same size.");

	d1 = mxGetDimensions(prhs[0]);
	d2 = mxGetDimensions(prhs[1]);
	for(c = 0; c < nd; ++c)
		if (d1[c] != d2[c] || d1[c] == 0)
			mexErrMsgTxt("Matrices must have the same size and not be empty.");

	for(c = 0; c < 2; ++c)
		if (!mxIsNumeric(prhs[c]) ||
		     mxIsComplex(prhs[c]) ||
		     mxIsSparse(prhs[c]))
			mexErrMsgTxt("Matrices must be numeric, real, full and single/double.");

	rd = nd - 1;
	lag = 0;
	if (nrhs > 2) {
		if (mxIsDouble(prhs[2]) &&
		    (mxGetNumberOfElements(prhs[2]) == 1)) {
			idbl1 = (const double *) mxGetPr(prhs[2]);
			if (!mxIsInf(*idbl1) &&
			    !mxIsNaN(*idbl1) &&
			    *idbl1 >= (2 - d1[0]) &&
			    *idbl1 <= (d1[rd] - 2)) {
				lag = (int) *idbl1;
				if (lag <= 0) {
					lag = -lag;
					if (d1[0] <= (1 + lag))
						mexErrMsgTxt("The first matrix dimensions must not be singleton.");
					rs = d1[0];
					ts = rs;
					for(c = 0; c < rd; ++c) {
						dr[c] = d1[c+1];
						ts *= dr[c];
					}
					nl = rs - lag;
					dnl = (double) nl;
					dnx = dnl - 1.0;

					plhs[0] = mxCreateNumericArray(rd, dr, mxDOUBLE_CLASS, mxREAL);
					t1 = (double *) mxGetPr(plhs[0]);

					switch (mxGetClassID(*prhs)) {
						case mxDOUBLE_CLASS: COV_1D_ILOOP(idbl1, double)
						case mxSINGLE_CLASS: COV_1D_ILOOP(isng1, float)
						case  mxINT32_CLASS: COV_1D_ILOOP(ii321, signed int)
						case mxUINT32_CLASS: COV_1D_ILOOP(iu321, unsigned int)
						case  mxINT16_CLASS: COV_1D_ILOOP(ii161, signed short)
						case mxUINT16_CLASS: COV_1D_ILOOP(iu161, unsigned short)
						case   mxINT8_CLASS: COV_1D_ILOOP(ii081, signed char)
						case  mxUINT8_CLASS:
						case mxLOGICAL_CLASS:COV_1D_ILOOP(iu081, unsigned char)
						default: mexErrMsgTxt("Datatype of first argument not supported.");
							 break;
					}
					return;
				}
			} else
				mexErrMsgTxt("The lag value must be between 0 and size - 2.");
		}
	}

	if (d1[rd] <= (1 + lag))
		mexErrMsgTxt("The last matrix dimensions must not be singleton.");

	rs = 1;
	for(c = 0; c < rd; ++c) {
		dr[c] = d1[c];
		rs *= dr[c];
	}
	nl = d1[rd];
	if (lag > 0)
		nl = nl - lag;

	ts = rs * nl;
	dnl = (double) nl;
	dnx = dnl - 1.0;

	plhs[0] = mxCreateNumericArray(rd, dr, mxDOUBLE_CLASS, mxREAL);
	t1 = (double *) mxGetPr(plhs[0]);

	switch (mxGetClassID(*prhs)) {
		case mxDOUBLE_CLASS: COV_ND_ILOOP(idbl1, double)
		case mxSINGLE_CLASS: COV_ND_ILOOP(isng1, float)
		case  mxINT32_CLASS: COV_ND_ILOOP(ii321, signed int)
		case mxUINT32_CLASS: COV_ND_ILOOP(iu321, unsigned int)
		case  mxINT16_CLASS: COV_ND_ILOOP(ii161, signed short)
		case mxUINT16_CLASS: COV_ND_ILOOP(iu161, unsigned short)
		case   mxINT8_CLASS: COV_ND_ILOOP(ii081, signed char)
		case  mxUINT8_CLASS:
		case mxLOGICAL_CLASS:COV_ND_ILOOP(iu081, unsigned char)
		default: mexErrMsgTxt("Datatype of first argument not supported.");
			 break;
	}
}

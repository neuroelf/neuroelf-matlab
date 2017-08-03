/*

return min, max, and mean of array

FORMAT:         mmm = minmaxmean(v [, flags]);

Input fields:

      v             N-d data array
      flags         sum of optional flags (default 0):
                    1 - also compute unbiased variance
                    4 - "safe" mode (skip Inf/NaNs in v, only for floats)

Output fields:

      mmm           [min, max, mean, minpos, maxpos, var]

Note: all output is computed over all elements--as if used with v(:) !

% Version:  v0.9d
% Build:    14080113
% Date:     Aug-01 2014, 1:35 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2010, 2014, Jochen Weber
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Columbia University nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "mex.h"
#include "math.h"
#include "isinfnan.h"


/* diverse PREPROCESSOR patterns */
#define MIN_MAX_ADDTOMEAN tmean += tval;

#define MIN_MAX_ADDTOMEANVAR \
        ++tnum; \
        delta = tval - tmean; \
        tmean += delta / (double) tnum; \
        mvar += delta * (tval - tmean); \

#define MIN_MAX_LOOP_OPEN(VAR) \
    for (numc = ne - 1; numc >= 0; --numc) { \
        tval = (double) (*VAR++);

#define MIN_MAX_LOOP_CLOSE }

#define MIN_MAX_MAKEMEAN(NELEM) tmean /= (double) NELEM;

#define MIN_MAX_MAKEMEANVAR(NELEM) tvar = mvar / (double) (NELEM - 1);

#define MIN_MAX_SAFE \
        IF_IS_BAD_VAL( tval ) continue;

#define MIN_MAX_VALPOS \
        if (tval > tmax) { tmax = tval; maxpos = numc; } \
        if (tval < tmin) { tmin = tval; minpos = numc; }


/* PREPROCESSOR function-like calls */
#define MIN_MAX_MEAN_VALPOS(VAR) \
    MIN_MAX_LOOP_OPEN( VAR ) \
        MIN_MAX_ADDTOMEAN \
        MIN_MAX_VALPOS \
    MIN_MAX_LOOP_CLOSE \
    MIN_MAX_MAKEMEAN( ne )

#define MIN_MAX_MEANVAR_VALPOS(VAR) \
    MIN_MAX_LOOP_OPEN( VAR ) \
        MIN_MAX_ADDTOMEANVAR \
        MIN_MAX_VALPOS \
    MIN_MAX_LOOP_CLOSE \
    MIN_MAX_MAKEMEANVAR( tnum )

#define MIN_MAX_MEAN_VALPOS_SAFE(VAR) \
    MIN_MAX_LOOP_OPEN( VAR ) \
        MIN_MAX_SAFE \
        ++tnum; \
        MIN_MAX_ADDTOMEAN \
        MIN_MAX_VALPOS \
    MIN_MAX_LOOP_CLOSE \
    MIN_MAX_MAKEMEAN( tnum )

#define MIN_MAX_MEANVAR_VALPOS_SAFE(VAR) \
    MIN_MAX_LOOP_OPEN( VAR ) \
        MIN_MAX_SAFE \
        MIN_MAX_ADDTOMEANVAR \
        MIN_MAX_VALPOS \
    MIN_MAX_LOOP_CLOSE \
    MIN_MAX_MAKEMEANVAR( tnum )


/* main function called by Matlab */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxClassID cid;
    const void *dat;
    const signed char *csc;
    const unsigned char *cuc;
    const signed short *css;
    const unsigned short *cus;
    const signed int *csl;
    const unsigned int *cul;
    const float *csf;
    const double *clf;
    double dflags, delta, mvar = 0.0;

    unsigned int flags = 0;
    double *odbl;

    double tval, tmin, tmax, tmean, tvar;
    signed long ne = 0, numc = 0, tnum = 0, maxpos = -1, minpos = -1;
    
    VARS_FOR_ISINFNAN;

	if (nrhs < 1 || nrhs > 2 || nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");

    tval = 0.0;
    INIT_INF_NAN_BAD_VAL();

    cid = mxGetClassID(*prhs);
    ne = mxGetNumberOfElements(*prhs);
    dat = (const void *) mxGetData(*prhs);
    csc = (const signed char *) dat;
    cuc = (const unsigned char *) dat;
    css = (const signed short *) dat;
    cus = (const unsigned short *) dat;
    csl = (const signed int *) dat;
    cul = (const unsigned int *) dat;
    csf = (const float *) dat;
    clf = (const double *) dat;

    *plhs = mxCreateDoubleMatrix(1, 6, mxREAL);
    if ((nrhs > 1) && (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfElements(prhs[1]) == 1)) {
        dflags = *((const double *) mxGetData(prhs[1]));
        if (!mxIsInf(dflags) && !mxIsNaN(dflags) &&
            (dflags >= 0.0) && (dflags <= 7.0))
            flags = (unsigned int) dflags;
        if ((cid != mxDOUBLE_CLASS) && (cid != mxSINGLE_CLASS))
            flags &= 0x1;
        else
            flags &= 0x5;
    }

    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output argument.");
    odbl = (double *) mxGetData(*plhs);
    if (odbl == NULL)
        mexErrMsgTxt("Invalid data pointer in output argument.");

    /* start real code here */
    tmin = mxGetInf();
    tmax = -mxGetInf();
    tmean = tvar = 0.0;

    if (ne == 0) {
        *odbl++ = tmin;
        *odbl++ = tmax;
        *odbl = mxGetNaN();
        if (nrhs > 1)
            odbl[3] = *odbl;
        return;
    }

    if (cid == mxLOGICAL_CLASS)
        cid = mxUINT8_CLASS;
    switch (flags) {
        case 0:
            switch (cid) {
                case mxDOUBLE_CLASS:   MIN_MAX_MEAN_VALPOS( clf )  break;
                case mxSINGLE_CLASS:   MIN_MAX_MEAN_VALPOS( csf )  break;
                case mxINT32_CLASS:    MIN_MAX_MEAN_VALPOS( csl )  break;
                case mxUINT32_CLASS:   MIN_MAX_MEAN_VALPOS( cul )  break;
                case mxINT16_CLASS:    MIN_MAX_MEAN_VALPOS( css )  break;
                case mxUINT16_CLASS:   MIN_MAX_MEAN_VALPOS( cus )  break;
                case mxINT8_CLASS:     MIN_MAX_MEAN_VALPOS( csc )  break;
                case mxUINT8_CLASS:    MIN_MAX_MEAN_VALPOS( cuc )  break;
                default: mexErrMsgTxt("Invalid input class.");
            }
            break;
        case 1:
            switch (cid) {
                case mxDOUBLE_CLASS:   MIN_MAX_MEANVAR_VALPOS( clf )  break;
                case mxSINGLE_CLASS:   MIN_MAX_MEANVAR_VALPOS( csf )  break;
                case mxINT32_CLASS:    MIN_MAX_MEANVAR_VALPOS( csl )  break;
                case mxUINT32_CLASS:   MIN_MAX_MEANVAR_VALPOS( cul )  break;
                case mxINT16_CLASS:    MIN_MAX_MEANVAR_VALPOS( css )  break;
                case mxUINT16_CLASS:   MIN_MAX_MEANVAR_VALPOS( cus )  break;
                case mxINT8_CLASS:     MIN_MAX_MEANVAR_VALPOS( csc )  break;
                case mxUINT8_CLASS:    MIN_MAX_MEANVAR_VALPOS( cuc )  break;
                default: mexErrMsgTxt("Invalid input class.");
            }
            break;
        case 4:
            switch (cid) {
                case mxDOUBLE_CLASS:   MIN_MAX_MEAN_VALPOS_SAFE( clf )  break;
                case mxSINGLE_CLASS:   MIN_MAX_MEAN_VALPOS_SAFE( csf )  break;
                default: mexErrMsgTxt("Invalid input class.");
            }
            break;
        case 5:
            switch (cid) {
                case mxDOUBLE_CLASS:   MIN_MAX_MEANVAR_VALPOS_SAFE( clf )  break;
                case mxSINGLE_CLASS:   MIN_MAX_MEANVAR_VALPOS_SAFE( csf )  break;
                default: mexErrMsgTxt("Invalid input class.");
            }
            break;
    }
    *odbl++ = tmin;
    *odbl++ = tmax;
    *odbl++ = tmean;
    *odbl++ = (double) (ne - minpos);
    *odbl++ = (double) (ne - maxpos);
    *odbl = tvar;
}

/*
    normcdfc.c - compute Phi for standard normal distribution

  FORMAT:       phi = normcdfc(x)

  Input fields:

        x           x ordinate to get normcdf at

  Output fields:

        phi         normcdf at x

  Note: this function supports single and double input.

% Version:  v0.9b
% Build:    11050511
% Date:     Aug-29 2010, 1:47 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2010, Jochen Weber
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

#include "math.h"
#include "mex.h"
#include "isinfnan.h"

/*
   Based on:
   Marsaglia, George "Evaluating the Normal Distribution",
   Journal of Statistical Software 11, 4 (July 2004).
   http://www.jstatsoft.org/
 */

#define NORMCDF_PHI(INPUT)                                                  \
    for ( ; ne > 0; --ne) {                                                 \
        x = (double) *INPUT++;                                              \
        IF_IS_NAN(x)                                                        \
            *outp++ = x;                                                    \
        else if (x < -8.0)                                                  \
            *outp++ = 0.0;                                                  \
        else if (x > 8)                                                     \
            *outp++ = 1.0;                                                  \
        else {                                                              \
            b = s = x;                                                      \
            t = 0.0;                                                        \
            q = x * x;                                                      \
            e = 1.0;                                                        \
            while(s != t)                                                   \
                s = (t = s) + (b *= q / (e += 2.0));                        \
            *outp++ = .5 + s * exp(-.5 * q - .91893853320467274178L);       \
        }                                                                   \
	}

/* MEX function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const double *inpd;
    const float *inpf;
    double *outp;
    const int *idim;
    int ndim, ne;
    long double s, t, b, q, e;
    double x;
    /* char vstr[256]; */

    VARS_FOR_ISINFNAN

    /* argument check */
    if (nrhs != 1)
        mexErrMsgTxt("One input and one output required.");
    if ((mxGetClassID(*prhs) != mxDOUBLE_CLASS) &&
        (mxGetClassID(*prhs) != mxSINGLE_CLASS))
        mexErrMsgTxt("Input must be either double or single.");

    /* get dims, etc. */
    ndim = mxGetNumberOfDimensions(*prhs);
    idim = mxGetDimensions(*prhs);
    ne = mxGetNumberOfElements(*prhs);

    /* create double output */
    *plhs = mxCreateNumericArray(ndim, idim, mxDOUBLE_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output array.");
    outp = (double *) mxGetData(*plhs);
    if (outp == NULL)
        mexErrMsgTxt("Error getting output data pointer.");

    /* return if nothing to do */
    if (ne == 0)
        return;

    /* get data pointer */
    inpd = (const double*) mxGetData(*prhs);

    /* initialize IS_BAD_VAL values */
    INIT_INF_NAN_BAD_VAL()

    /* depending on class */
    if (mxGetClassID(*prhs) == mxDOUBLE_CLASS) {
        NORMCDF_PHI(inpd)
    } else {
        inpf = (const float*) inpd;
        NORMCDF_PHI(inpf)
    }
}

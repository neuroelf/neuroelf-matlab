/*

converting thresholded map to looked up colors

FORMAT:       l = threshlutc(map, lut)

Input fields:

      map         XxY double map with values
      lut         Cx3 lookup colors (C must be even!)

Output fields:

      l           looked-up color map (XxYx3 uint8)

% Version:  v0.9b
% Build:    11050512
% Date:     Sep-16 2010, 9:13 AM EST
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

#include "mex.h"
#include "math.h"
#include <stdio.h>
#include "isinfnan.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* dimensions */
    int nv = 0, vc = 0, nc = 0;
    int od[3] = {0, 0, 3};
    const int *id = NULL;
    short ci = 0, lc = 0;

    /* pointers */
    const double *map = NULL, *inval = NULL;
    unsigned char *lum = NULL, *lum2 = NULL, *lum3 = NULL,
        *lut1 = NULL, *lut2 = NULL, *lut3 = NULL,
        *nlut1 = NULL, *nlut2 = NULL, *nlut3 = NULL;

    /* thresholds and value */
    double dnc = 0.0, val = 0.0;

    /* isinfnan */
    VARS_FOR_ISINFNAN

    /* variable output string */
    /* char vstr[256]; */

    /* check number, type, fields of in/out arguments */
	if (nrhs != 2)
		mexErrMsgTxt("Bad number of input arguments.");
    if (!mxIsDouble(*prhs))
        mexErrMsgTxt("Input map must be of type double.");
    if (mxGetNumberOfDimensions(*prhs) != 2)
        mexErrMsgTxt("Input map must be 2D.");
    nv = mxGetNumberOfElements(*prhs);
    id = (const int *) mxGetDimensions(*prhs);
    *od = *id++;
    od[1] = *id;
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("LUT must be of type double.");
    if (mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("LUT must be 2D.");
    id = (const int *) mxGetDimensions(prhs[1]);
    if (id[1] != 3)
        mexErrMsgTxt("LUT must be of size Cx3.");
    nc = (int) (0.5 * ((double) (*id)));
    if ((2 * nc) != *id)
        mexErrMsgTxt("LUT must contain even number of colors.");
    if (nc > 32000)
        mexErrMsgTxt("Only up to 32000 colors (per tail) supported.");
    lc = (short) (nc - 1);
    dnc = (double) lc;

    /* get map and LUT */
    map = (const double *) mxGetData(*prhs);
    inval = (const double *) mxGetData(prhs[1]);

    /* create output */
    *plhs = mxCreateNumericArray(3, od, mxUINT8_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output.");
    lum = (unsigned char *) mxGetData(*plhs);
    if (lum == NULL)
        mexErrMsgTxt("Error allocating memory for output.");
    lum2 = &lum[nv];
    lum3 = &lum2[nv];

    /* create uint8 LUT copy */
    lut1 = mxCalloc(6 * nc, sizeof(unsigned char));
    if (lut1 == NULL) {
        mxDestroyArray(*plhs);
        mexErrMsgTxt("Error allocating uint8 LUT buffer.");
    }
    for (nlut1 = lut1, vc = 6 * nc; vc > 0; --vc)
        *nlut1++ = (unsigned char) *inval++;
    nlut1 = &lut1[nc];
    lut2 = &nlut1[nc];
    nlut2 = &lut2[nc];
    lut3 = &nlut2[nc];
    nlut3 = &lut3[nc];

    /* init isinfnan */
    INIT_INF_NAN_BAD_VAL()

    /* iterate over values */
    for (vc = nv; vc > 0; --vc) {

        /* get value */
        val = *map++;

        /* for zero values */
        if (val == 0.0) {
            ++lum;
            ++lum2;
            ++lum3;
            continue;
        }

        /* also for bad values */
        IF_IS_BAD_VAL(val) {
            ++lum;
            ++lum2;
            ++lum3;
            continue;
        }

        /* for positive */
        if (val > 0.0) {

            /* for values at the end of the range, set to last color */
            if (val >= 1.0) {
                *lum++ = lut1[lc];
                *lum2++ = lut2[lc];
                *lum3++ = lut3[lc];

            /* otherwise */
            } else {

                /* get index */
                ci = (short) (dnc * val);

                /* put color into output */
                *lum++ = lut1[ci];
                *lum2++ = lut2[ci];
                *lum3++ = lut3[ci];
            }

        /* for negative */
        } else {

            /* flip sign */
            val = -val;
            if (val >= 1.0) {
                *lum++ = nlut1[lc];
                *lum2++ = nlut2[lc];
                *lum3++ = nlut3[lc];
            } else {
                ci = (short) (dnc * val);
                *lum++ = nlut1[ci];
                *lum2++ = nlut2[ci];
                *lum3++ = nlut3[ci];
            }
        }
    }

    /* remove LUT buffer */
    mxFree(lut1);
}

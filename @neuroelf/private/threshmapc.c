/*

thresholding a double map to 0...1

FORMAT:       t = threshmapc(map, lt, ut, tails)

Input fields:

      map         XxY double map with values
      lt          lower threshold
      ut          upper threshold
      tails       tails flag (1, positive tail, 2 negative tail, 3 both)

Output fields:

      t           thresholded map (-1...1)

% Version:  v0.9a
% Build:    11050512
% Date:     May-17 2010, 10:48 AM EST
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
    int nd = 0, nv = 0, vc = 0;
    const int *id = NULL;

    /* pointers */
    const double *map = NULL, *inval = NULL;
    double *omap = NULL;

    /* thresholds and value */
    double lt = 0.0, ut = 0.0, lut = 0.0, nlt = 0.0, val = 0.0;

    /* boolean flags */
    bool ptail = 1, ntail = 1;

    /* isinfnan */
    VARS_FOR_ISINFNAN

    /* variable output string */
    /* char vstr[256]; */

    /* check number, type, fields of in/out arguments */
	if ((nrhs < 3) ||
        (nlhs > 1))
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (!mxIsDouble(*prhs))
        mexErrMsgTxt("Input map must be of type double.");
    nv = mxGetNumberOfElements(*prhs);
    nd = mxGetNumberOfDimensions(*prhs);
    if (nd > 3)
        mexErrMsgTxt("Input map must be 2D or 3D.");
    id = (const int *) mxGetDimensions(*prhs);
    map = (const double *) mxGetData(*prhs);
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("Input lt must be of type double.");
    if (mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgTxt("Input lt must be 1x1.");
    if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("Input ut must be of type double.");
    if (mxGetNumberOfElements(prhs[2]) != 1)
        mexErrMsgTxt("Input ut must be 1x1.");
    inval = (const double *) mxGetData(prhs[1]);
    lt = *inval;
    if (lt < 0.0)
        lt = -lt;
    else if (lt == 0.0)
        lt = 0.0000001;
    inval = (const double *) mxGetData(prhs[2]);
    ut = *inval;
    if (ut < 0)
        ut = -ut;
    else if (ut == 0.0)
        ut = 0.0000002;
    if (lt > ut) {
        val = ut;
        ut = lt;
        lt = val;
    } else if (lt == ut)
        ut = ut + 0.0000001;
    lut = ut - lt;
    lt -= 0.0000000000001;
    nlt = -lt;

    /* create outputs */
    *plhs = (mxArray *) mxCreateNumericArray(nd, id, mxDOUBLE_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating first output mxArray.");
    omap = (double *) mxGetData(*plhs);
    if (omap == NULL){
        mxDestroyArray(*plhs);
        mexErrMsgTxt("Error getting data pointers to output mxArray.");
    }

    /* tails */
    if (nrhs > 3) {
        if (mxIsDouble(prhs[3]) &&
            (mxGetNumberOfElements(prhs[3]) == 1)) {
            inval = (const double *) mxGetData(prhs[3]);
            if (*inval == 1.0)
                ntail = 0;
            else if (*inval == 2.0)
                ptail = 0;
            else if (*inval != 3.0)
                return;
        }
    }

    /* init isinfnan */
    INIT_INF_NAN_BAD_VAL()

    /* perform computation */
    if (!ntail) {
        for (vc = nv; vc > 0; --vc) {
            val = *map++;
            IF_IS_BAD_VAL(val) {
                *omap++ = 0.0;
                continue;
            }
            if (val < lt) {
                *omap++ = 0.0;
                continue;
            }
            val -= lt;
            val /= lut;
            *omap++ = (val > 1.0) ? 1.0 : val;
        }
    } else if (!ptail) {
        for (vc = nv; vc > 0; --vc) {
            val = *map++;
            IF_IS_BAD_VAL(val) {
                *omap++ = 0.0;
                continue;
            }
            if (val > nlt) {
                *omap++ = 0.0;
                continue;
            }
            val -= nlt;
            val /= lut;
            *omap++ = (val < -1.0) ? -1.0 : val;
        }
    } else {
        for (vc = nv; vc > 0; --vc) {
            val = *map++;
            IF_IS_BAD_VAL(val) {
                *omap++ = 0.0;
                continue;
            }
            if (val > 0.0) {
                if (val < lt) {
                    *omap++ = 0.0;
                    continue;
                }
                val -= lt;
                val /= lut;
                *omap++ = (val > 1.0) ? 1.0 : val;
            } else {
                if (val > nlt) {
                    *omap++ = 0.0;
                    continue;
                }
                val -= nlt;
                val /= lut;
                *omap++ = (val < -1.0) ? -1.0 : val;
            }
        }
    }
}

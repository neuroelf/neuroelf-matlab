/*

sprintfbx  - apply sprintf('%bx', ...) format

FORMAT:       s = sprintfbx(V)

Input fields:

      V           N-D double matrix

Output fields:

      s           output string

% Version:  v0.9c
% Build:    11050515
% Date:     May-05 2010, 3:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2011, Jochen Weber
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
#include "isinfnan.h"

const mxChar sbx_digits[16] = {
     48,  49,  50,  51,  52,  53,  54,  55,
     56,  57,  97,  98,  99, 100, 101, 102};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t c, ne;
    const double *idbl;
    const float *isng;
    unsigned int ii;
    const unsigned int *ii2;
    mxChar *out;

    VARS_FOR_ISINFNAN

    if (nrhs != 1 ||
        nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (mxIsComplex(*prhs) ||
        mxIsSparse(*prhs) ||
        !mxIsNumeric(*prhs) ||
        (mxGetClassID(*prhs) != mxDOUBLE_CLASS &&
         mxGetClassID(*prhs) != mxSINGLE_CLASS))
        mexErrMsgTxt("Input array must be real, full double array.");
    ne = mxGetNumberOfElements(*prhs);
    if (ne == 0) {
        *plhs = mxCreateNumericMatrix(0, 0, mxCHAR_CLASS, mxREAL);
        return;
    }

    INIT_INF_NAN_BAD_VAL()

    if (mxGetClassID(*prhs) == mxDOUBLE_CLASS) {
        idbl = (const double *) mxGetData(*prhs);
        if (idbl == NULL)
            mexErrMsgTxt("Error getting input array pointer.");
        *plhs = mxCreateNumericMatrix(1, 16 * (int) ne, mxCHAR_CLASS, mxREAL);
        if (*plhs == NULL)
            mexErrMsgTxt("Error creating output array.");
        out = (mxChar *) mxGetData(*plhs);
        if (out == NULL)
            mexErrMsgTxt("Error getting output array pointer.");
        ii2 = (const unsigned int*) &isbadval.ui[1];

        for (c = 0; c < ne; ++c) {
            isbadval.d = *idbl++;
            ii = *ii2;
            *out++ = sbx_digits[(ii & 0xF0000000) >> 28];
            *out++ = sbx_digits[(ii & 0x0F000000) >> 24];
            *out++ = sbx_digits[(ii & 0x00F00000) >> 20];
            *out++ = sbx_digits[(ii & 0x000F0000) >> 16];
            *out++ = sbx_digits[(ii & 0x0000F000) >> 12];
            *out++ = sbx_digits[(ii & 0x00000F00) >> 8];
            *out++ = sbx_digits[(ii & 0x000000F0) >> 4];
            *out++ = sbx_digits[(ii & 0x0000000F)];
            ii = *isbadval.ui;
            *out++ = sbx_digits[(ii & 0xF0000000) >> 28];
            *out++ = sbx_digits[(ii & 0x0F000000) >> 24];
            *out++ = sbx_digits[(ii & 0x00F00000) >> 20];
            *out++ = sbx_digits[(ii & 0x000F0000) >> 16];
            *out++ = sbx_digits[(ii & 0x0000F000) >> 12];
            *out++ = sbx_digits[(ii & 0x00000F00) >> 8];
            *out++ = sbx_digits[(ii & 0x000000F0) >> 4];
            *out++ = sbx_digits[(ii & 0x0000000F)];
        }
    } else {
        isng = (const float *) mxGetData(*prhs);
        if (isng == NULL)
            mexErrMsgTxt("Error getting input array pointer.");
        *plhs = mxCreateNumericMatrix(1, 8 * (int) ne, mxCHAR_CLASS, mxREAL);
        if (*plhs == NULL)
            mexErrMsgTxt("Error creating output array.");
        out = (mxChar *) mxGetData(*plhs);
        if (out == NULL)
            mexErrMsgTxt("Error getting output array pointer.");

        for (c = 0; c < ne; ++c) {
            *isbadval.f = *isng++;
            ii = *isbadval.ui;
            *out++ = sbx_digits[(ii & 0xF0000000) >> 28];
            *out++ = sbx_digits[(ii & 0x0F000000) >> 24];
            *out++ = sbx_digits[(ii & 0x00F00000) >> 20];
            *out++ = sbx_digits[(ii & 0x000F0000) >> 16];
            *out++ = sbx_digits[(ii & 0x0000F000) >> 12];
            *out++ = sbx_digits[(ii & 0x00000F00) >> 8];
            *out++ = sbx_digits[(ii & 0x000000F0) >> 4];
            *out++ = sbx_digits[(ii & 0x0000000F)];
        }
    }
}

/*

checking whether (any) element(s) are inf/nan

FORMAT:       o = isinfnan(V, fa, fi, fn)

Input fields:

      V           N-D single/double matrix (for other numeric arrays, false)
      fa          "any" flag (return scalar true/false; default: true)
      fi          "inf" flag (test for infs, default: true)
      fn          "nan" flag (test for nans, default: true)

Output fields:

      o           output

% Version:  v0.9a
% Build:    11050510
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
#include "isinfnan.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nd;
    size_t c, ne;
    const int *dim;
    const double *idbl;
    const float  *isng;
    const unsigned char *flag;
    unsigned char *out;
    bool fa = 0, fi = 1, fn = 1;

    VARS_FOR_ISINFNAN

    if (nrhs < 1 ||
        nrhs > 4 ||
        nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (mxIsComplex(*prhs) ||
        mxIsSparse(*prhs) ||
        !mxIsNumeric(*prhs))
        mexErrMsgTxt("Input array must be real, full and numeric.");
	nd = mxGetNumberOfDimensions(*prhs);
    ne = mxGetNumberOfElements(*prhs);
	dim = mxGetDimensions(*prhs);
    if ((nrhs > 1) &&
        (mxGetClassID(prhs[1]) == mxLOGICAL_CLASS) &&
        (mxGetNumberOfElements(prhs[1]) == 1)) {
        flag = (const unsigned char*) mxGetData(prhs[1]);
        fa = (*flag != 0);
    }
    if ((nrhs > 2) &&
        (mxGetClassID(prhs[2]) == mxLOGICAL_CLASS) &&
        (mxGetNumberOfElements(prhs[2]) == 1)) {
        flag = (const unsigned char*) mxGetData(prhs[2]);
        fi = (*flag != 0);
    }
    if ((nrhs > 3) &&
        (mxGetClassID(prhs[3]) == mxLOGICAL_CLASS) &&
        (mxGetNumberOfElements(prhs[3]) == 1)) {
        flag = (const unsigned char*) mxGetData(prhs[3]);
        fn = (*flag != 0);
    }
    if (fa)
        *plhs = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
    else
        *plhs = mxCreateNumericArray(nd, dim, mxLOGICAL_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output array.");
    out = (unsigned char *) mxGetData(*plhs);
    if (out == NULL)
        mexErrMsgTxt("Error getting output array pointer.");

    INIT_INF_NAN_BAD_VAL()

    if (fa) {
        if ((!mxIsDouble(*prhs)) &&
            (!mxIsSingle(*prhs))) {
            *out = 1;
            return;
        }
        if (mxIsDouble(*prhs)) {
            idbl = (const double*) mxGetData(*prhs);
            if (fi && fn) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_BAD_VAL(*idbl++) {
                        *out = 1;
                        return;
                    }
                }
            } else if (fi) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_INF(*idbl++) {
                        *out = 1;
                        return;
                    }
                }
            } else if (fn) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_NAN(*idbl++) {
                        *out = 1;
                        return;
                    }
                }
            }
        } else {
            isng = (const float*) mxGetData(*prhs);
            if (fi && fn) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_BAD_VAL(*isng++) {
                        *out = 1;
                        return;
                    }
                }
            } else if (fi) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_INF(*isng++) {
                        *out = 1;
                        return;
                    }
                }
            } else if (fn) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_NAN(*isng++) {
                        *out = 1;
                        return;
                    }
                }
            }
        }
    } else {
        if ((!mxIsDouble(*prhs)) &&
            (!mxIsSingle(*prhs))) {
            for (c = ne; c > 0; --c) {
                *out++ = 1;
                return;
            }
        }
        if (mxIsDouble(*prhs)) {
            idbl = (const double*) mxGetData(*prhs);
            if (fi && fn) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_BAD_VAL(*idbl++)
                        *out++ = 1;
                    else
                        ++out;
                }
            } else if (fi) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_INF(*idbl++)
                        *out++ = 1;
                    else
                        ++out;
                }
            } else if (fn) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_NAN(*idbl++)
                        *out++ = 1;
                    else
                        ++out;
                }
            }
        } else {
            isng = (const float*) mxGetData(*prhs);
            if (fi && fn) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_BAD_VAL(*isng++)
                        *out++ = 1;
                    else
                        ++out;
                }
            } else if (fi) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_INF(*isng++)
                        *out++ = 1;
                    else
                        ++out;
                }
            } else if (fn) {
                for (c = 0; c < ne; ++c) {
                    IF_IS_NAN(*isng++)
                        *out++ = 1;
                    else
                        ++out;
                }
            }
        }
    }
}

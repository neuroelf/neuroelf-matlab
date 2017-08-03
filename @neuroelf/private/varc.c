/*

calculation of variance for N-D matrices

FORMAT:       v = varc(X,d,s)

Input fields:

      X           N-D single/double matrix
      d           dim (otherwise first non-singleton)
      s           if any third argument is given, skip over Inf/NaN values

Output fields:

      v           variance

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
#include "isinfnan.h"

#define VARC_UNSAFE(P1,P2)                                              \
    if (i1 == 1) {                                                      \
        for (c3 = 0; c3 < ne; c3 += i2) {                               \
            for (P2 = P1, s = 0.0, c1 = 0; c1 < ti; ++c1)               \
                s += ((double) *P2++);                                  \
            s /= ddnum;                                                 \
            for (P2 = P1, vr = 0.0, c1 = 0; c1 < ti; ++c1) {            \
                v = ((double) *P2++) - s;                               \
                vr += v * v;                                            \
            }                                                           \
            *varo++ = vr / ddnom;                                       \
            P1 = &P1[ti];                                               \
        }                                                               \
    } else {                                                            \
        for (c3 = 0; c3 < ne; c3 += i2) {                               \
            for (c2 = 0; c2 < i1; ++c2) {                               \
                for (P2 = &P1[c2+c3], s = 0.0, c1 = 0; c1 < ti; ++c1) { \
                    s += ((double) *P2);                                \
                    P2 = &P2[i1];                                       \
                }                                                       \
                s /= ddnum;                                             \
                for (P2 = &P1[c2+c3], vr = 0.0, c1 = 0; c1 < ti; ++c1) {\
                    v = ((double) *P2) - s;                             \
                    P2 = &P2[i1];                                       \
                    vr += v * v;                                        \
                }                                                       \
                *varo++ = vr / ddnom;                                   \
            }                                                           \
        }                                                               \
    }

#define VARC_SAFE(P1,P2)                                                \
    if (i1 == 1) {                                                      \
        for (c3 = 0; c3 < ne; c3 += i2) {                               \
            for (P2 = P1, dnom = -1, s = 0.0, c1 = 0; c1 < ti; ++c1) {  \
                v = *P2++;                                              \
                IF_IS_BAD_VAL(v) continue;                              \
                ++dnom;                                                 \
                s += v;                                                 \
            }                                                           \
            s /= ((double) (dnom + 1));                                 \
            for (P2 = P1, vr = 0.0, c1 = 0; c1 < ti; ++c1) {            \
                v = *P2++ - s;                                          \
                IF_IS_BAD_VAL(v) continue;                              \
                vr += v * v;                                            \
            }                                                           \
            *varo++ = (dnom > 0) ? (vr / ((double) dnom)) : 0.0;        \
            P1 = &P1[ti];                                               \
        }                                                               \
    } else {                                                            \
        for (c3 = 0; c3 < ne; c3 += i2) {                               \
            for (c2 = 0; c2 < i1; ++c2) {                               \
                for (P2 = &P1[c2+c3], dnom = -1, s = 0.0, c1 = 0; c1 < ti; ++c1) { \
                    v = *P2;                                            \
                    P2 = &P2[i1];                                       \
                    IF_IS_BAD_VAL(v) continue;                          \
                    ++dnom;                                             \
                    s += v;                                             \
                }                                                       \
                s /= ((double) (dnom + 1));                             \
                for (P2 = &P1[c2+c3], vr = 0.0, c1 = 0; c1 < ti; ++c1) {\
                    v = *P2 - s;                                        \
                    P2 = &P2[i1];                                       \
                    IF_IS_BAD_VAL(v) continue;                          \
                    vr += v * v;                                        \
                }                                                       \
                *varo++ = (dnom > 0) ? (vr / ((double) dnom)) : 0.0;    \
            }                                                           \
        }                                                               \
    }


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nd, ne, td, ti, c1, c2, c3, i1, i2, dnom;
    int odim[64];
    const int *dim;
	double s, v, vr, ddnom, ddnum, *varo;
    const double *idbl, *idb2;
    const float  *isng, *isn2;
    const int    *iint, *iin2;
    const short  *ishr, *ish2;
    const char   *ichr, *ich2;
    const unsigned int   *iuin, *iui2;
    const unsigned short *iush, *ius2;
    const unsigned char  *iuch, *iuc2;
    /* char vstr[256]; */

    VARS_FOR_ISINFNAN

    if (nrhs < 1 ||
        nrhs > 3 ||
        nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (mxIsComplex(*prhs) ||
        mxIsSparse(*prhs) ||
        !mxIsNumeric(*prhs))
        mexErrMsgTxt("Input array must be real, full and single/double.");
	nd = mxGetNumberOfDimensions(*prhs);
    if (nd > 64)
        mexErrMsgTxt("Input array cannot be more than 64-D.");
    ne = mxGetNumberOfElements(*prhs);
	dim = mxGetDimensions(*prhs);
    td = 0;
    INIT_INF_NAN_BAD_VAL()
    if (nrhs > 1) {
        if (mxIsDouble(prhs[1]) &&
            (mxGetNumberOfElements(prhs[1]) == 1)) {
            idbl = (const double *) mxGetPr(prhs[1]);
            IF_IS_GOOD_VAL(*idbl) {
                if ((*idbl >= 1.0) &&
                    (*idbl <= ((double) nd)))
                    td = (int) *idbl;
            }
        }
    }
    if (td == 0) {
        for (c1 = 0; c1 < nd; ++c1)
            if (dim[c1] > 1) {
                td = c1 + 1;
                break;
            }
    }
    if ((td == 0) ||
        (dim[td-1] < 2))
        mexErrMsgTxt("Invalid dimension provided or array not usable for varc.");
    --td;
    ti = dim[td];
    ddnom = (double) (ti - 1);
    ddnum = (double) ti;
	for (c1 = 0; c1 < nd; ++c1)
        odim[c1] = *dim++;
    odim[td] = 1;
    *plhs = mxCreateNumericArray(nd, odim, mxDOUBLE_CLASS, mxREAL);
    if (ne == 0)
        return;
    idbl = (const double*) mxGetData(*prhs);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output array.");
    varo = (double *) mxGetData(*plhs);
    if (varo == NULL)
        mexErrMsgTxt("Error getting output pointer.");
    i1 = 1;
    for (c1 = 0; c1 < td; ++c1)
        i1 *= odim[c1];
    i2 = i1 * ti;

    switch (mxGetClassID(*prhs)) {
        case mxDOUBLE_CLASS:
            if (nrhs < 3)
                VARC_UNSAFE(idbl, idb2)
            else
                VARC_SAFE(idbl, idb2)
            break;
        case mxSINGLE_CLASS:
            isng = (const float *) idbl;
            if (nrhs < 3)
                VARC_UNSAFE(isng, isn2)
            else
                VARC_SAFE(isng, isn2)
            break;
        case mxINT32_CLASS:
            iint = (const int *) idbl;
            VARC_UNSAFE(iint, iin2)
            break;
        case mxINT16_CLASS:
            ishr = (const short *) idbl;
            VARC_UNSAFE(ishr, ish2)
            break;
        case mxINT8_CLASS:
            ichr = (const char *) idbl;
            VARC_UNSAFE(ichr, ich2)
            break;
        case mxUINT32_CLASS:
            iuin = (const unsigned int *) idbl;
            VARC_UNSAFE(iuin, iui2)
            break;
        case mxUINT16_CLASS:
            iush = (const unsigned short *) idbl;
            VARC_UNSAFE(iush, ius2)
            break;
        case mxUINT8_CLASS:
        case mxLOGICAL_CLASS:
            iuch = (const unsigned char *) idbl;
            VARC_UNSAFE(iuch, iuc2)
            break;
        default:
            mexErrMsgTxt("Unsupported input class.");
    }
}

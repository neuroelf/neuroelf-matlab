/*

limit data to a range

FORMAT:       V = limitrangec(V, rmin, rmax [, rinv])

Input fields:

      V           N-D numeric matrix
      rmin        minimum value
      rmax        maximum value
      rinv        if given, replace invalid values (Inf/NaN) with this

Output fields:

      V           limited output

% Version:  v0.9a
% Build:    11050511
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

const double singlerealmax = 3.4028235e+38;

#define LIMIT_RANGE(IN, OUT, INVAL, MINVAL, MAXVAL, REPVAL)             \
    if (dinv) {                                                         \
        if (MINVAL >= MAXVAL) {                                         \
            for (c = ne; c > 0; --c) {                                  \
                INVAL = *IN++;                                          \
                IF_IS_BAD_VAL(INVAL) {                                  \
                    *OUT++ = REPVAL;                                    \
                } else {                                                \
                    *OUT++ = MAXVAL;                                    \
                }                                                       \
            }                                                           \
        } else {                                                        \
            if (!dmin) {                                                \
                if (!dmax) {                                            \
                    for (c = ne; c > 0; --c) {                          \
                        INVAL = *IN++;                                  \
                        IF_IS_BAD_VAL(INVAL) {                          \
                            *OUT++ = REPVAL;                            \
                        } else {                                        \
                            *OUT++ = INVAL;                             \
                        }                                               \
                    }                                                   \
                } else {                                                \
                    for (c = ne; c > 0; --c) {                          \
                        INVAL = *IN++;                                  \
                        IF_IS_BAD_VAL(INVAL) {                          \
                            *OUT++ = REPVAL;                            \
                        } else {                                        \
                            *OUT++ = (INVAL < MAXVAL) ? INVAL : MAXVAL; \
                        }                                               \
                    }                                                   \
                }                                                       \
            } else if (!dmax) {                                         \
                for (c = ne; c > 0; --c) {                              \
                    INVAL = *IN++;                                      \
                    IF_IS_BAD_VAL(INVAL) {                              \
                        *OUT++ = REPVAL;                                \
                    } else {                                            \
                        *OUT++ = (INVAL > MINVAL) ? INVAL : MINVAL;     \
                    }                                                   \
                }                                                       \
            } else {                                                    \
                for (c = ne; c > 0; --c) {                              \
                    INVAL = *IN++;                                      \
                    IF_IS_BAD_VAL(INVAL) {                              \
                        *OUT++ = REPVAL;                                \
                    } else if (INVAL < MINVAL)                          \
                        *OUT++ = MINVAL;                                \
                    else                                                \
                        *OUT++ = (INVAL > MAXVAL) ? MAXVAL : INVAL;     \
                }                                                       \
            }                                                           \
        }                                                               \
    } else {                                                            \
        if (MINVAL >= MAXVAL) {                                         \
            for (c = ne; c > 0; --c) {                                  \
                *OUT++ = MAXVAL;                                        \
            }                                                           \
        } else {                                                        \
            if (!dmin) {                                                \
                if (!dmax) {                                            \
                    for (c = ne; c > 0; --c)                            \
                        *OUT++ = *IN++;                                 \
                } else {                                                \
                    for (c = ne; c > 0; --c) {                          \
                        INVAL = *IN++;                                  \
                        *OUT++ = (INVAL < MAXVAL) ? INVAL : MAXVAL;     \
                    }                                                   \
                }                                                       \
            } else if (!dmax) {                                         \
                for (c = ne; c > 0; --c) {                              \
                    INVAL = *IN++;                                      \
                    *OUT++ = (INVAL > MINVAL) ? INVAL : MINVAL;         \
                }                                                       \
            } else {                                                    \
                for (c = ne; c > 0; --c) {                              \
                    INVAL = *IN++;                                      \
                    if (INVAL < MINVAL)                                 \
                        *OUT++ = MINVAL;                                \
                    else                                                \
                        *OUT++ = (INVAL > MAXVAL) ? MAXVAL : INVAL;     \
                }                                                       \
            }                                                           \
        }                                                               \
    }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nd, ne, c;
    mxClassID cid;
    double rmin, rmax, rinv = 0.0, rval;
    float srmin = 0.0, srmax = 0.0, srinv = 0.0, srval;
    int imin = 0, imax = 0, iinv = 0, ival;
    short smin = 0, smax = 0, sinv = 0, sval;
    char cmin = 0, cmax = 0, cinv = 0, cval;
    unsigned int uimin = 0, uimax = 0, uiinv = 0, uival;
    unsigned short usmin = 0, usmax = 0, usinv = 0, usval;
    unsigned char ucmin = 0, ucmax = 0, ucinv = 0, ucval;
    const int *dim;
    const double *idbl;
    const float  *isng;
    const int    *iint;
    const short  *ishr;
    const char   *ichr;
    const unsigned int   *uint;
    const unsigned short *ushr;
    const unsigned char  *uchr;
    double *odbl;
    float  *osng;
    int    *oint;
    short  *oshr;
    char   *ochr;
    unsigned int   *ouint;
    unsigned short *oushr;
    unsigned char  *ouchr;
    bool dinv = 0, dmin = 1, dmax = 1;

    VARS_FOR_ISINFNAN


    if (nrhs < 3 ||
        nrhs > 4 ||
        nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (mxIsComplex(*prhs) ||
        mxIsSparse(*prhs) ||
        !mxIsNumeric(*prhs))
        mexErrMsgTxt("Input array must be real, full and numeric.");
    cid = mxGetClassID(*prhs);
    if ((cid < mxDOUBLE_CLASS) ||
        (cid > mxUINT32_CLASS))
        mexErrMsgTxt("Input must be one of the numeric classes.");
    if (!mxIsDouble(prhs[1]) ||
        (mxGetNumberOfElements(prhs[1]) != 1))
        mexErrMsgTxt("Input rmin must 1x1 double.");
    if (!mxIsDouble(prhs[2]) ||
        (mxGetNumberOfElements(prhs[2]) != 1))
        mexErrMsgTxt("Input rmax must 1x1 double.");

    INIT_INF_NAN_BAD_VAL()

    idbl = (const double*) mxGetData(prhs[1]);
    rmin = *idbl;
    IF_IS_NAN(rmin)
        mexErrMsgTxt("Input rmin must not be NaN.");
    IF_IS_POS_INF(rmin)
        mexErrMsgTxt("Input rmin must not be +Inf.");
    IF_IS_NEG_INF(rmin)
        dmin = 0;

    idbl = (const double*) mxGetData(prhs[2]);
    rmax = *idbl;
    IF_IS_NAN(rmax)
        mexErrMsgTxt("Input rmax must not be NaN.");
    IF_IS_NEG_INF(rmax)
        mexErrMsgTxt("Input rmax must not be -Inf.");
    IF_IS_POS_INF(rmax)
        dmax = 0;

    if (rmin > rmax)
        mexErrMsgTxt("Invalid range.");

    if (nrhs > 3) {
        if (!mxIsDouble(prhs[3]) ||
            (mxGetNumberOfElements(prhs[3]) != 1))
            mexErrMsgTxt("If given, rinv must be 1x1 double.");
        idbl = (const double*) mxGetData(prhs[3]);
        rinv = *idbl;
        IF_IS_BAD_VAL(rinv)
            mexErrMsgTxt("If given, rinv must not be Inf/NaN.");
        dinv = 1;
        IF_IS_BAD_VAL(rinv)
            mexErrMsgTxt("If given, rinv must not be Inf/NaN.");
    }

	nd = mxGetNumberOfDimensions(*prhs);
    ne = mxGetNumberOfElements(*prhs);
	dim = mxGetDimensions(*prhs);
    idbl = mxGetData(*prhs);
    *plhs = mxCreateNumericArray(nd, dim, cid, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output array.");
    odbl = (double *) mxGetData(*plhs);
    if (odbl == NULL)
        mexErrMsgTxt("Error getting output array pointer.");

    switch (cid) {
        case mxDOUBLE_CLASS:
            LIMIT_RANGE(idbl, odbl, rval, rmin, rmax, rinv)
            break;
        case mxSINGLE_CLASS:
            isng = (const float*) idbl;
            osng = (float*) odbl;
            if (rmin <= (-1.0 * singlerealmax))
                dmin = 0;
            else
                srmin = (float) rmin;
            if (rmax >= singlerealmax)
                dmax = 0;
            else
                srmax = (float) rmax;
            if (dinv)
                srinv = (float) rinv;
            LIMIT_RANGE(isng, osng, srval, srmin, srmax, srinv)
            break;
        case mxINT32_CLASS:
            iint = (const int*) idbl;
            oint = (int*) idbl;
            if (rmin < -2147483647.0)
                dmin = 0;
            else
                imin = (int) rmin;
            if (rmax >= 2147483647.0)
                dmax = 0;
            else
                imax = (int) rmax;
            if (dinv)
                iinv = (int) rinv;
            LIMIT_RANGE(iint, oint, ival, imin, imax, iinv)
            break;
        case mxUINT32_CLASS:
            uint = (const unsigned int*) idbl;
            ouint = (unsigned int*) idbl;
            if (rmin < 1.0)
                dmin = 0;
            else
                uimin = (unsigned int) rmin;
            if (rmax >= 4294967295.0)
                dmax = 0;
            else
                uimax = (unsigned int) rmax;
            if (dinv)
                uiinv = (unsigned int) rinv;
            LIMIT_RANGE(uint, ouint, uival, uimin, uimax, uiinv)
            break;
        case mxINT16_CLASS:
            ishr = (const short*) idbl;
            oshr = (short*) idbl;
            if (rmin < -32767.0)
                dmin = 0;
            else
                smin = (short) rmin;
            if (rmax >= 32767.0)
                dmax = 0;
            else
                smax = (short) rmax;
            if (dinv)
                sinv = (short) rinv;
            LIMIT_RANGE(ishr, oshr, sval, smin, smax, sinv)
            break;
        case mxUINT16_CLASS:
            ushr = (const unsigned short*) idbl;
            oushr = (unsigned short*) idbl;
            if (rmin < 1.0)
                dmin = 0;
            else
                usmin = (unsigned short) rmin;
            if (rmax >= 65535.0)
                dmax = 0;
            else
                usmax = (unsigned short) rmax;
            if (dinv)
                usinv = (unsigned short) rinv;
            LIMIT_RANGE(ushr, oushr, usval, usmin, usmax, usinv)
            break;
        case mxINT8_CLASS:
            ichr = (const char*) idbl;
            ochr = (char*) idbl;
            if (rmin < -127.0)
                dmin = 0;
            else
                cmin = (char) rmin;
            if (rmax >= 127.0)
                dmax = 0;
            else
                cmax = (char) rmax;
            if (dinv)
                cinv = (char) rinv;
            LIMIT_RANGE(ichr, ochr, cval, cmin, cmax, cinv)
            break;
        case mxUINT8_CLASS:
            uchr = (const unsigned char*) idbl;
            ouchr = (unsigned char*) idbl;
            if (rmin < 1.0)
                dmin = 0;
            else
                ucmin = (unsigned char) rmin;
            if (rmax >= 255.0)
                dmax = 0;
            else
                ucmax = (unsigned char) rmax;
            if (dinv)
                ucinv = (unsigned char) rinv;
            LIMIT_RANGE(uchr, ouchr, ucval, ucmin, ucmax, ucinv)
            break;
        default:
            mexErrMsgTxt("Unsupported class.");
            break;
    }
}

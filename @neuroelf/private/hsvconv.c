/*

convert from and to HSV/RGB schemes

FORMAT:       col = hsvconv(col, dir, scale)

Input fields:

      col         Xx3 or X*Y*3 or X*Y*Z*3 data (uint8, single, double)
      dir         either {1} (from HSV to RGB) or 2 (from RGB to HSV)
      scale       scaling for hue value, {1} ([0 .. 1]) or 2 ([0 .. 360))

Output fields:

      col         converted data (uint8->double for hsv, uint8 for rgb)

% Version:  v0.9b
% Build:    11050511
% Date:     Aug-11 2010, 5:28 PM EST
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
#include <stdio.h>

#define DOWNSCALE 0.00392156862745098 *
#define HSVHSCALE 0.00277777777777778 *
#define UPSCALE 255.0 *
#define SIXTH 0.16666666666666667
#define THIRD 0.33333333333333333
#define TWOTHIRDS 0.66666666666666667
#define FIVESIXTH 0.83333333333333333

#define CONVERT_HSV2RGB_U8(IP1, IP2, IP3, OUTTYPE, OP1, OP2, OP3)           \
    if (hsvscale1)                                                          \
        for (c = nv; c > 0; --c) {                                          \
            h = (double) *IP1++;                                            \
            s = (double) *IP2++;                                            \
            v = (double) *IP3++;                                            \
            x = 6.0 * h;                                                    \
            ndims = (int) x;                                                \
            x -= (double) ndims;                                            \
            if (!(ndims & 1))                                               \
                x = 1.0 - x;                                                \
            r = v * (1.0 - s);                                              \
            g = v * (1.0 - s * x);                                          \
            switch (ndims) {                                                \
                case 6:                                                     \
                case 0:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE v);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE g);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE r);                         \
                    break;                                                  \
                case 1:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE g);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE v);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE r);                         \
                    break;                                                  \
                case 2:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE r);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE v);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE g);                         \
                    break;                                                  \
                case 3:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE r);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE g);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE v);                         \
                    break;                                                  \
                case 4:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE g);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE r);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE v);                         \
                    break;                                                  \
                case 5:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE v);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE r);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE g);                         \
                    break;                                                  \
            }                                                               \
        }                                                                   \
    else                                                                    \
        for (c = nv; c > 0; --c) {                                          \
            h = HSVHSCALE ((double) *IP1++);                                \
            s = (double) *IP2++;                                            \
            v = (double) *IP3++;                                            \
            x = 6.0 * h;                                                    \
            ndims = (int) x;                                                \
            x -= (double) ndims;                                            \
            if (!(ndims & 1))                                               \
                x = 1.0 - x;                                                \
            r = v * (1.0 - s);                                              \
            g = v * (1.0 - s * x);                                          \
            switch (ndims) {                                                \
                case 6:                                                     \
                case 0:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE v);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE g);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE r);                         \
                    break;                                                  \
                case 1:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE g);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE v);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE r);                         \
                    break;                                                  \
                case 2:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE r);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE v);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE g);                         \
                    break;                                                  \
                case 3:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE r);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE g);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE v);                         \
                    break;                                                  \
                case 4:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE g);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE r);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE v);                         \
                    break;                                                  \
                case 5:                                                     \
                    *OP1++ = (OUTTYPE) (UPSCALE v);                         \
                    *OP2++ = (OUTTYPE) (UPSCALE r);                         \
                    *OP3++ = (OUTTYPE) (UPSCALE g);                         \
                    break;                                                  \
            }                                                               \
        }

#define CONVERT_RGB2HSV(IP1, IP2, IP3, OUTTYPE, OP1, OP2, OP3)              \
    if (hsvscale1)                                                          \
        for (c = nv; c > 0; --c) {                                          \
            r = (double) *IP1++;                                            \
            g = (double) *IP2++;                                            \
            b = (double) *IP3++;                                            \
            if (r > g) {                                                    \
                x = ((g < b) ? g : b);                                      \
                v = ((r > b) ? r : b);                                      \
            } else {                                                        \
                x = ((r < b) ? r : b);                                      \
                v = ((g > b) ? g : b);                                      \
            }                                                               \
            *OP3++ = (OUTTYPE) v;                                           \
            if (v == x) {                                                   \
                *OP1++ = (OUTTYPE) 0;                                       \
                *OP2++ = (OUTTYPE) 0;                                       \
                continue;                                                   \
            }                                                               \
            s = v - x;                                                      \
            *OP2++ = (OUTTYPE) (s / v);                                     \
            if (r == x) {                                                   \
                h = g - b;                                                  \
                x = 0.5;                                                    \
            } else if (b == x) {                                            \
                h = r - g;                                                  \
                x = SIXTH;                                                  \
            } else {                                                        \
                h = b - r;                                                  \
                x = FIVESIXTH;                                              \
            }                                                               \
            *OP1++ = (OUTTYPE) (x - h / (6.0 * s));                         \
        }                                                                   \
    else                                                                    \
        for (c = nv; c > 0; --c) {                                          \
            r = (double) *IP1++;                                            \
            g = (double) *IP2++;                                            \
            b = (double) *IP3++;                                            \
            if (r > g) {                                                    \
                x = ((g < b) ? g : b);                                      \
                v = ((r > b) ? r : b);                                      \
            } else {                                                        \
                x = ((r < b) ? r : b);                                      \
                v = ((g > b) ? g : b);                                      \
            }                                                               \
            *OP3++ = (OUTTYPE) v;                                           \
            if (v == x) {                                                   \
                *OP1++ = (OUTTYPE) 0;                                       \
                *OP2++ = (OUTTYPE) 0;                                       \
                continue;                                                   \
            }                                                               \
            s = v - x;                                                      \
            *OP2++ = (OUTTYPE) (s / v);                                     \
            if (r == x) {                                                   \
                h = g - b;                                                  \
                x = 0.5;                                                    \
            } else if (b == x) {                                            \
                h = r - g;                                                  \
                x = SIXTH;                                                  \
            } else {                                                        \
                h = b - r;                                                  \
                x = FIVESIXTH;                                              \
            }                                                               \
            *OP1++ = (OUTTYPE) 360.0 * (x - h / (6.0 * s));                 \
        }

#define CONVERT_RGB2HSV_SC(IP1, IP2, IP3, OUTTYPE, OP1, OP2, OP3)           \
    if (hsvscale1)                                                          \
        for (c = nv; c > 0; --c) {                                          \
            r = DOWNSCALE ((double) *IP1++);                                \
            g = DOWNSCALE ((double) *IP2++);                                \
            b = DOWNSCALE ((double) *IP3++);                                \
            if (r > g) {                                                    \
                x = ((g < b) ? g : b);                                      \
                v = ((r > b) ? r : b);                                      \
            } else {                                                        \
                x = ((r < b) ? r : b);                                      \
                v = ((g > b) ? g : b);                                      \
            }                                                               \
            *OP3++ = (OUTTYPE) v;                                           \
            if (v == x) {                                                   \
                *OP1++ = (OUTTYPE) 0;                                       \
                *OP2++ = (OUTTYPE) 0;                                       \
                continue;                                                   \
            }                                                               \
            s = v - x;                                                      \
            *OP2++ = (OUTTYPE) (s / v);                                     \
            if (r == x) {                                                   \
                h = g - b;                                                  \
                x = 0.5;                                                    \
            } else if (b == x) {                                            \
                h = r - g;                                                  \
                x = SIXTH;                                                  \
            } else {                                                        \
                h = b - r;                                                  \
                x = FIVESIXTH;                                              \
            }                                                               \
            *OP1++ = (OUTTYPE) (x - h / (6.0 * s));                         \
        }                                                                   \
    else                                                                    \
        for (c = nv; c > 0; --c) {                                          \
            r = DOWNSCALE ((double) *IP1++);                                \
            g = DOWNSCALE ((double) *IP2++);                                \
            b = DOWNSCALE ((double) *IP3++);                                \
            if (r > g) {                                                    \
                x = ((g < b) ? g : b);                                      \
                v = ((r > b) ? r : b);                                      \
            } else {                                                        \
                x = ((r < b) ? r : b);                                      \
                v = ((g > b) ? g : b);                                      \
            }                                                               \
            *OP3++ = (OUTTYPE) v;                                           \
            if (v == x) {                                                   \
                *OP1++ = (OUTTYPE) 0;                                       \
                *OP2++ = (OUTTYPE) 0;                                       \
                continue;                                                   \
            }                                                               \
            s = v - x;                                                      \
            *OP2++ = (OUTTYPE) (s / v);                                     \
            if (r == x) {                                                   \
                h = g - b;                                                  \
                x = 0.5;                                                    \
            } else if (b == x) {                                            \
                h = r - g;                                                  \
                x = SIXTH;                                                  \
            } else {                                                        \
                h = b - r;                                                  \
                x = FIVESIXTH;                                              \
            }                                                               \
            *OP1++ = (OUTTYPE) 360.0 * (x - h / (6.0 * s));                 \
        }


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* dimensions */
    const int *idim;
    int ndims, nv = 0, c;

    /* pointers */
    const double *dp = NULL, *dp2 = NULL, *dp3 = NULL;
    const unsigned char *pp = NULL, *pp2 = NULL, *pp3 = NULL;
    const float *fp = NULL, *fp2 = NULL, *fp3 = NULL;
    double *odp, *odp2, *odp3;
    unsigned char *opp, *opp2, *opp3;
    float *ofp, *ofp2, *ofp3;

    /* values */
    double r, g, b, h, s, v, x;

    /* direction and scaling flags */
    bool fromhsv = 1, hsvscale1 = 1;

    /* variable output string */
    /* char vstr[256]; */

    /* check number, type, fields of in/out arguments */
	if ((nrhs < 1) ||
        (nlhs > 1))
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (!mxIsDouble(*prhs) &&
        !mxIsSingle(*prhs) &&
        !mxIsUint8(*prhs))
        mexErrMsgTxt("First input must be of either double, single, or uint8 type.");
    if (mxGetNumberOfElements(*prhs) < 3)
        mexErrMsgTxt("Too few elements in input data.");
    ndims = mxGetNumberOfDimensions(*prhs);
    if (ndims > 4)
        mexErrMsgTxt("First input can only be up to 4-D");
    idim = mxGetDimensions(*prhs);
    if (idim[ndims-1] != 3)
        mexErrMsgTxt("First input must have the last non-singleton dimension = 3.");

    /* compute number of values */
    nv = 1;
    for (c = ndims - 2; c >= 0; --c)
        nv *= idim[c];

    /* get optional settings */
    if ((nrhs > 1) &&
        mxIsDouble(prhs[1]) &&
        (mxGetNumberOfElements(prhs[1]) == 1)) {
        dp = (const double *) mxGetData(prhs[1]);
        if (dp != NULL) {
            if (!mxIsInf(*dp) &&
                !mxIsNaN(*dp) &&
                (*dp == 2.0))
                    fromhsv = 0;
        }
    }
    if ((nrhs > 2) &&
        mxIsDouble(prhs[2]) &&
        (mxGetNumberOfElements(prhs[2]) == 1)) {
        dp = (const double *) mxGetData(prhs[2]);
        if (dp != NULL) {
            if (!mxIsInf(*dp) &&
                !mxIsNaN(*dp) &&
                (*dp == 2.0))
                    hsvscale1 = 0;
        }
    }

    /* create output */
    if (fromhsv)
        *plhs = mxCreateNumericArray(ndims, idim, mxUINT8_CLASS, mxREAL);
    else if (mxGetClassID(*prhs) != mxSINGLE_CLASS)
        *plhs = mxCreateNumericArray(ndims, idim, mxDOUBLE_CLASS, mxREAL);
    else
        *plhs = mxCreateNumericArray(ndims, idim, mxSINGLE_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output array.");

    /* switch over input class */
    switch (mxGetClassID(*prhs)) {
        case mxDOUBLE_CLASS:
            dp = (const double *) mxGetData(*prhs);
            if (dp == NULL)
                mexErrMsgTxt("Error getting input data pointer.");
            dp2 = &dp[nv];
            dp3 = &dp2[nv];
            if (fromhsv) {
                opp = (unsigned char*) mxGetData(*plhs);
                if (opp == NULL)
                    mexErrMsgTxt("Error getting output data pointer.");
                opp2 = &opp[nv];
                opp3 = &opp2[nv];
                CONVERT_HSV2RGB_U8(dp, dp2, dp3, unsigned char, opp, opp2, opp3);
            } else {
                odp = (double *) mxGetData(*plhs);
                if (odp == NULL)
                    mexErrMsgTxt("Error getting output data pointer.");
                odp2 = &odp[nv];
                odp3 = &odp2[nv];
                CONVERT_RGB2HSV(dp, dp2, dp3, double, odp, odp2, odp3);
            }
            break;
        case mxSINGLE_CLASS:
            fp = (const float *) mxGetData(*prhs);
            if (fp == NULL)
                mexErrMsgTxt("Error getting input data pointer.");
            fp2 = &fp[nv];
            fp3 = &fp2[nv];
            odp = (double *) mxGetData(*plhs);
            if (odp == NULL)
                mexErrMsgTxt("Error getting output data pointer.");
            if (fromhsv) {
                opp = (unsigned char*) mxGetData(*plhs);
                if (opp == NULL)
                    mexErrMsgTxt("Error getting output data pointer.");
                opp2 = &opp[nv];
                opp3 = &opp2[nv];
                CONVERT_HSV2RGB_U8(fp, fp2, fp3, unsigned char, opp, opp2, opp3);
            } else {
                ofp = (float *) mxGetData(*plhs);
                if (ofp == NULL)
                    mexErrMsgTxt("Error getting output data pointer.");
                ofp2 = &ofp[nv];
                ofp3 = &ofp2[nv];
                CONVERT_RGB2HSV(fp, fp2, fp3, float, ofp, ofp2, ofp3);
            }
            break;
        case mxUINT8_CLASS:
            pp = (const unsigned char *) mxGetData(*prhs);
            if (pp == NULL)
                mexErrMsgTxt("Error getting input data pointer.");
            pp2 = &pp[nv];
            pp3 = &pp2[nv];
            if (fromhsv)
                mexErrMsgTxt("Conversion from HSV to RGB with uint8 input unsupported");
            else {
                odp = (double *) mxGetData(*plhs);
                if (odp == NULL)
                    mexErrMsgTxt("Error getting output data pointer.");
                odp2 = &odp[nv];
                odp3 = &odp2[nv];
                CONVERT_RGB2HSV_SC(pp, pp2, pp3, double, odp, odp2, odp3);
            }
            break;
        default:
            mexErrMsgTxt("Unsupported class.");
            break;
    }
}

/*

return kendall tau ties-correction value

FORMAT:         ps = kendtaupairsign(m, d);

Input fields:

      m             matrix
      d             dim

Output fields:

      ps            sign of all pairs

Note: the order in which pairs are made is (1;2), (1;3), ..., (1;n),
      (2;3), ... (n-1;n)

% Version:  v0.9a
% Build:    11102110
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
#include "isinfnan.h"

#define OPEN_SIGN_LOOP(DTYPE, PTR1, PTR2)                               \
            PTR1 = (const DTYPE *) mxGetData(*prhs);                    \
            for (c1 = 0; c1 < nepre; ++c1) {                            \
                for (c2 = c1, tc2 = c1; c2 < ne; c2 += inc, tc2 += tinc) { \
                    PTR2 = &PTR1[c2];                                   \
                    otg = &odbl[tc2];                                   \
                    for (i1 = 0; i1 < len; i1 += nepre) {               \
                        for (i2 = i1 + nepre; i2 < len; i2 += nepre) {

#define OPEN_SIGN_LOOP_FL(DTYPE, PTR1, PTR2)                            \
    OPEN_SIGN_LOOP(DTYPE, PTR1, PTR2)                                   \
                            dv = ((double) PTR2[i2]) - ((double) PTR2[i1]); \
                            IF_IS_NAN(dv)                               \
                                dv = 0.0;                               \
                            else if (dv < 0.0)                          \
                                dv = -1.0;                              \
                            else if (dv > 0.0)                          \
                                dv = 1.0;                               \

#define OPEN_SIGN_LOOP_INT(DTYPE, PTR1, PTR2)                           \
    OPEN_SIGN_LOOP(DTYPE, PTR1, PTR2)                                   \
                            dv = ((double) PTR2[i2]) - ((double) PTR2[i1]); \
                            if (dv < 0.0)                               \
                                dv = -1.0;                              \
                            else if (dv > 0.0)                          \
                                dv = 1.0;

#define CLOSE_SIGN_LOOP                                                 \
                            *otg = dv;                                  \
                            otg = &otg[nepre];                          \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
            break;

#define SIGN_LOOP_FL(DTYPE, PTR1, PTR2)                                 \
    OPEN_SIGN_LOOP_FL(DTYPE, PTR1, PTR2)                                \
    CLOSE_SIGN_LOOP
#define SIGN_LOOP_INT(DTYPE, PTR1, PTR2)                                \
    OPEN_SIGN_LOOP_INT(DTYPE, PTR1, PTR2)                               \
    CLOSE_SIGN_LOOP

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int *dim;
    int odim[64];
    const double *dbl, *idbl;
    const float *sng, *isng;
    const signed char *si8, *isi8;
    const unsigned char *ui8, *iui8;
    const signed short *si16, *isi16;
    const unsigned short *ui16, *iui16;
    const signed int *si32, *isi32;
    const unsigned int *ui32, *iui32;
    int td, cid, nd, ne, len, nepre, inc, tinc, c1, c2, tc2, i1, i2;
    double *odbl, *otg;
    double dv;
    VARS_FOR_ISINFNAN

	if (nrhs != 2 || nlhs != 1)
		mexErrMsgTxt("Bad number of input/output arguments.");
    INIT_INF_NAN_BAD_VAL()
    if ((mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ||
        (mxGetNumberOfElements(prhs[1]) != 1) ||
        (!mxIsNumeric(*prhs)))
        mexErrMsgTxt("Invalid dim argument dims/class.");
    dv = *((const double *) mxGetPr(prhs[1]));
    IF_IS_BAD_VAL(dv) {
        mexErrMsgTxt("Bad dim argument.");
    }
    if ((dv < 1) ||
        (dv > mxGetNumberOfDimensions(*prhs)))
        mexErrMsgTxt("Bad dim argument.");
    td = ((int) dv) - 1;
    cid = mxGetClassID(*prhs);
    nd = mxGetNumberOfDimensions(*prhs);
    ne = mxGetNumberOfElements(*prhs);
    dim = mxGetDimensions(*prhs);
    len = dim[td];
    nepre = inc = tinc = 1;
    c2 = (int) (0.5 * ((double) len) * (((double) len) - 1.0));
    for (c1 = 0; c1 <= nd; ++c1) {
        odim[c1] = dim[c1];
        if (c1 < td)
            nepre *= dim[c1];
        else if (c1 == td) {
            odim[c1] = c2;
            inc = nepre * dim[c1];
            tinc = nepre * c2;
        }
    }
    len *= nepre;
    *plhs = mxCreateNumericArray(nd, odim, mxDOUBLE_CLASS, mxREAL);

    /* return early if no ties possible */
    if (c2 < 1)
        return;

    /* get output pointer */
    odbl = (double *) mxGetPr(*plhs);

    /* switch over class */
    switch (cid) {
        case mxDOUBLE_CLASS: SIGN_LOOP_FL(double, idbl, dbl)
        case mxSINGLE_CLASS: SIGN_LOOP_FL(float, isng, sng)
        case mxINT8_CLASS:   SIGN_LOOP_INT(signed char, isi8, si8)
        case mxUINT8_CLASS:  SIGN_LOOP_INT(unsigned char, iui8, ui8)
        case mxINT16_CLASS:  SIGN_LOOP_INT(signed short, isi16, si16)
        case mxUINT16_CLASS: SIGN_LOOP_INT(unsigned short, iui16, ui16)
        case mxINT32_CLASS:  SIGN_LOOP_INT(signed int, isi32, si32)
        case mxUINT32_CLASS: SIGN_LOOP_INT(unsigned int, iui32, ui32)
        default:
            mexErrMsgTxt("Invalid input class.");
    }
}

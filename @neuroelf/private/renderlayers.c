/*

rendering of transimg layers into a HxWx3 uint8 image

FORMAT:       r = renderlayers(tio)

Input fields:

      tio         1x1 struct of transimg content

Output fields:

      r           HxWx3 uint8 image

% Version:  v0.9b
% Build:    12012119
% Date:     Aug-07 2010, 6:02 PM EST
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* dimensions */
    int h = 0, w = 0, hw = 0, hw3 = 0, nl = 0, lc = 0,
        lfield = -1, tfield = -1, pfield = -1, afield = -1,
        rpfield = -1, rafield = -1, trfield = -1, ndp = 0, nda = 0,
        pc = 0, sx = 0, sy = 0, tx = 0, ty = 0, sxy = 0, sxy3 = 0;
    int od[3] = {0, 0, 3};
    const mxArray *layer = NULL, *intype = NULL, *inpp = NULL, *inpa = NULL;
    char typec = 0;

    /* pointers */
    const double *dp = NULL;
    const float *ap = NULL, *fpp = NULL, *fpp2 = NULL, *fpp3 = NULL;
    const unsigned char *pp = NULL, *pp2 = NULL, *pp3 = NULL;
    const unsigned int *lpdim = NULL;
    unsigned char *opp = NULL;
    const mxChar *typestr = NULL;

    /* background color */
    float bg1 = 0.0, bg2 = 0.0, bg3 = 0.0, al = 1.0, alm = 0.0;

    /* buffer */
    float *buffer = NULL, *bf1 = NULL, *bf2 = NULL, *bf3 = NULL,
          *cbf1 = NULL, *cbf2 = NULL, *cbf3 = NULL;

    /* variable output string */
    /* char vstr[256]; */

    /* check number, type, fields of in/out arguments */
	if ((nrhs != 1) ||
        (nlhs > 1))
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (!mxIsStruct(*prhs))
        mexErrMsgTxt("Input must be of type struct.");
    lfield = mxGetFieldNumber(*prhs, "Layer");
    if ((mxGetFieldNumber(*prhs, "Background") < 0) ||
        (mxGetFieldNumber(*prhs, "Height") < 0) ||
        (lfield < 0) ||
        (mxGetFieldNumber(*prhs, "Width") < 0))
        mexErrMsgTxt("Required field missing.");
    layer = mxGetFieldByNumber(*prhs, 0, lfield);
    if (!mxIsStruct(layer))
        mexErrMsgTxt(".Layer field must be of type struct.");
    tfield = mxGetFieldNumber(layer, "Type");
    pfield = mxGetFieldNumber(layer, "Pixel");
    afield = mxGetFieldNumber(layer, "Alpha");
    rpfield = mxGetFieldNumber(layer, "RPixel");
    rafield = mxGetFieldNumber(layer, "RAlpha");
    trfield = mxGetFieldNumber(layer, "Trans");
    if ((tfield < 0) ||
        (pfield < 0) ||
        (afield < 0) ||
        (rpfield < 0) ||
        (rafield < 0) ||
        (trfield < 0))
        mexErrMsgTxt(".Layer field must have required subfields.");
    nl = mxGetNumberOfElements(layer);
    inpa = (const mxArray *) mxGetField(*prhs, 0, "Background");
    if (mxGetNumberOfElements(inpa) != 3)
        mexErrMsgTxt("Invalid .Background field.");
    if (!mxIsUint8(inpa))
        mexErrMsgTxt("Invalid .Background field datatype.");
    pp = (const unsigned char *) mxGetData(inpa);
    bg1 = (float) *pp++;
    bg2 = (float) *pp++;
    bg3 = (float) *pp;
    inpa = (const mxArray *) mxGetField(*prhs, 0, "Height");
    if (mxGetNumberOfElements(inpa) != 1)
        mexErrMsgTxt("Invalid .Height field.");
    if (!mxIsDouble(inpa))
        mexErrMsgTxt("Invalid .Height field datatype.");
    dp = (const double *) mxGetData(inpa);
    h = (int) *dp;
    inpa = (const mxArray *) mxGetField(*prhs, 0, "Width");
    if (mxGetNumberOfElements(inpa) != 1)
        mexErrMsgTxt("Invalid .Width field.");
    if (!mxIsDouble(inpa))
        mexErrMsgTxt("Invalid .Width field datatype.");
    dp = (const double *) mxGetData(inpa);
    w = (int) *dp;
    *od = h;
    od[1] = w;
    hw = h * w;
    hw3 = 3 * hw;

    /* create output */
    *plhs = (mxArray *) mxCreateNumericArray(3, od, mxUINT8_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output mxArray.");
    opp = (unsigned char *) mxGetData(*plhs);
    if (opp == NULL)
        mexErrMsgTxt("Error getting data pointer to output mxArray.");

    /* create float representation of output for temp computation */
    buffer = (float *) mxCalloc(hw3, sizeof(float));
    if (buffer == NULL) {
        mxDestroyArray(*plhs);
        mexErrMsgTxt("Error allocating temp buffer.");
    }

    /* set background color */
    for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], pc = hw; pc > 0; --pc) {
        *bf1++ = bg1;
        *bf2++ = bg2;
        *bf3++ = bg3;
    }

    /* parse layers */
    for (lc = 0; lc < nl; ++lc) {

        /* get type information */
        intype = mxGetFieldByNumber(layer, lc, tfield);
        if (mxGetNumberOfElements(intype) != 1) {
            mxFree(buffer);
            mxDestroyArray(*plhs);
            mexErrMsgTxt("Invalid .Type field.");
        }
        typestr = (const mxChar *) mxGetData(intype);
        typec = (char) *typestr;

        /* full information */
        if ((typec == 'f') ||
            (typec == 't')) {

            /* get pixel and alpha information */
            if (typec == 'f') {
                inpp = mxGetFieldByNumber(layer, lc, pfield);
                inpa = mxGetFieldByNumber(layer, lc, afield);
            } else {
                inpp = mxGetFieldByNumber(layer, lc, rpfield);
                inpa = mxGetFieldByNumber(layer, lc, rafield);
            }
            ndp = mxGetNumberOfDimensions(inpp);
            nda = mxGetNumberOfDimensions(inpa);
            if (((!mxIsUint8(inpp)) &&
                 (!mxIsSingle(inpp))) ||
                (!mxIsSingle(inpa)) ||
                (ndp > 3) ||
                (nda > 2)) {
                mxFree(buffer);
                mxDestroyArray(*plhs);
                mexErrMsgTxt("Invalid .Pixel/.Alpha field types/dims.");
            }
            ndp = mxGetNumberOfElements(inpp);
            nda = mxGetNumberOfElements(inpa);
            if (((ndp != hw) &&
                 (ndp != hw3)) ||
                ((nda != 1) &&
                 (nda != hw))) {
                mxFree(buffer);
                mxDestroyArray(*plhs);
                mexErrMsgTxt("Invalid .Pixel/.Alpha field sizes.");
            }
            ap = (const float *) mxGetData(inpa);

            /* skip invisible layers */
            if ((nda == 1) &&
                (*ap == 0.0))
                continue;

            /* depending on type -> uint8 */
            if (mxIsUint8(inpp)) {

                /* get data (as uint8 -> unsigned char) */
                pp = (const unsigned char *) mxGetData(inpp);

                /* depending on size of alpha */
                if (nda == 1) {

                    /* get alpha value */
                    al = (float) *ap;
                    alm = 1.0 - al;

                    /* depending on size of pixels */
                    if (ndp == hw3) {

                        /* iterate over pixel */
                        for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], pp2 = &pp[hw], pp3 = &pp2[hw], pc = hw; pc > 0; --pc) {
                            *bf1 = alm * *bf1 + al * *pp++; ++bf1;
                            *bf2 = alm * *bf2 + al * *pp2++; ++bf2;
                            *bf3 = alm * *bf3 + al * *pp3++; ++bf3;
                        }
                    } else {
                        for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], pc = hw; pc > 0; --pc) {
                            *bf1 = alm * *bf1 + al * *pp; ++bf1;
                            *bf2 = alm * *bf2 + al * *pp; ++bf2;
                            *bf3 = alm * *bf3 + al * *pp++; ++bf3;
                        }
                    }

                } else {

                    if (ndp == hw3) {
                        for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], pp2 = &pp[hw], pp3 = &pp2[hw], pc = hw; pc > 0; --pc) {
                            al = (float) *ap++;
                            alm = 1.0 - al;
                            *bf1 = alm * *bf1 + al * *pp++; ++bf1;
                            *bf2 = alm * *bf2 + al * *pp2++; ++bf2;
                            *bf3 = alm * *bf3 + al * *pp3++; ++bf3;
                        }
                    } else {
                        for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], pc = hw; pc > 0; --pc) {
                            al = (float) *ap++;
                            alm = 1.0 - al;
                            *bf1 = alm * *bf1 + al * *pp; ++bf1;
                            *bf2 = alm * *bf2 + al * *pp; ++bf2;
                            *bf3 = alm * *bf3 + al * *pp++; ++bf3;
                        }
                    }
                }

            /* type -> single */
            } else {

                /* get data (as single -> float) */
                fpp = (const float *) mxGetData(inpp);

                /* depending on size of alpha */
                if (nda == 1) {

                    /* get alpha value */
                    al = (float) *ap;
                    alm = 1.0 - al;

                    /* depending on size of pixels */
                    if (ndp == hw3) {

                        /* iterate over pixel */
                        for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], fpp2 = &fpp[hw], fpp3 = &fpp2[hw], pc = hw; pc > 0; --pc) {
                            *bf1 = alm * *bf1 + al * *fpp++; ++bf1;
                            *bf2 = alm * *bf2 + al * *fpp2++; ++bf2;
                            *bf3 = alm * *bf3 + al * *fpp3++; ++bf3;
                        }
                    } else {
                        for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], pc = hw; pc > 0; --pc) {
                            *bf1 = alm * *bf1 + al * *fpp; ++bf1;
                            *bf2 = alm * *bf2 + al * *fpp; ++bf2;
                            *bf3 = alm * *bf3 + al * *fpp++; ++bf3;
                        }
                    }

                } else {

                    if (ndp == hw3) {
                        for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], fpp2 = &fpp[hw], fpp3 = &fpp2[hw], pc = hw; pc > 0; --pc) {
                            al = (float) *ap++;
                            alm = 1.0 - al;
                            *bf1 = alm * *bf1 + al * *fpp++; ++bf1;
                            *bf2 = alm * *bf2 + al * *fpp2++; ++bf2;
                            *bf3 = alm * *bf3 + al * *fpp3++; ++bf3;
                        }
                    } else {
                        for (bf1 = buffer, bf2 = &bf1[hw], bf3 = &bf2[hw], pc = hw; pc > 0; --pc) {
                            al = (float) *ap++;
                            alm = 1.0 - al;
                            *bf1 = alm * *bf1 + al * *fpp; ++bf1;
                            *bf2 = alm * *bf2 + al * *fpp; ++bf2;
                            *bf3 = alm * *bf3 + al * *fpp++; ++bf3;
                        }
                    }
                }
            }

        /* partial data (and data created from shapes/text) */
        } else {

            /* get pixel and alpha information */
            inpp = mxGetFieldByNumber(layer, lc, pfield);
            inpa = mxGetFieldByNumber(layer, lc, afield);
            ndp = mxGetNumberOfDimensions(inpp);
            nda = mxGetNumberOfDimensions(inpa);
            if (((!mxIsUint8(inpp)) &&
                 (!mxIsSingle(inpp))) ||
                (!mxIsSingle(inpa)) ||
                (ndp > 3) ||
                (nda > 2)) {
                mxFree(buffer);
                mxDestroyArray(*plhs);
                mexErrMsgTxt("Invalid .Pixel/.Alpha field types/dims.");
            }
            ap = (const float *) mxGetData(inpa);

            /* skip invisible layers */
            if ((nda == 1) &&
                (*ap == 0.0))
                continue;

            /* get offset and size */
            intype = mxGetFieldByNumber(layer, lc, trfield);
            if (!mxIsDouble(intype) ||
                (mxGetNumberOfElements(intype) != 2)) {
                mxFree(buffer);
                mxDestroyArray(*plhs);
                mexErrMsgTxt("Invalid .Trans field.");
            }
            dp = (const double *) mxGetData(intype);
            tx = (int) *dp++;
            ty = (int) *dp;
            lpdim = (const unsigned int *) mxGetDimensions(inpp);
            sx = *lpdim++;
            sy = *lpdim;
            if ((tx < 0) ||
                (ty < 0) ||
                ((tx + sx) > h) ||
                ((ty + sy) > w)) {
                mxFree(buffer);
                mxDestroyArray(*plhs);
                mexErrMsgTxt(".Pixel must fit within buffer.");
            }
            lpdim = (const unsigned int *) mxGetDimensions(inpa);
            if (((*lpdim * lpdim[1]) != 1) &&
                ((*lpdim != sx) ||
                 (lpdim[1] != sy))) {
                mxFree(buffer);
                mxDestroyArray(*plhs);
                mexErrMsgTxt(".Pixel and .Alpha fields must match in size.");
            }
            ndp = mxGetNumberOfElements(inpp);
            nda = mxGetNumberOfElements(inpa);
            sxy = sx * sy;
            sxy3 = 3 * sxy;

            /* depending on type -> uint8 */
            if (mxIsUint8(inpp)) {

                /* get data (as uint8 -> unsigned char) */
                pp = (const unsigned char *) mxGetData(inpp);

                /* depending on size of alpha */
                if (nda == 1) {

                    /* get alpha value */
                    al = (float) *ap;
                    alm = 1.0 - al;

                    /* depending on size of pixels */
                    if (ndp == sxy3) {

                        /* iterate over columns */
                        for (cbf1 = &buffer[tx + ty * h], cbf2 = &cbf1[hw], cbf3 = &cbf2[hw], pp2 = &pp[sxy], pp3 = &pp2[sxy]; sy > 0; --sy, cbf1 = &cbf1[h], cbf2 = &cbf2[h], cbf3 = &cbf3[h]) {

                            /* iterate over pixel */
                            for (bf1 = cbf1, bf2 = cbf2, bf3 = cbf3, pc = sx; pc > 0; --pc) {
                                *bf1 = alm * *bf1 + al * *pp++; ++bf1;
                                *bf2 = alm * *bf2 + al * *pp2++; ++bf2;
                                *bf3 = alm * *bf3 + al * *pp3++; ++bf3;
                            }
                        }
                    } else {
                        for (cbf1 = &buffer[tx + ty * h], cbf2 = &cbf1[hw], cbf3 = &cbf2[hw]; sy > 0; --sy, cbf1 = &cbf1[h], cbf2 = &cbf2[h], cbf3 = &cbf3[h]) {
                            for (bf1 = cbf1, bf2 = cbf2, bf3 = cbf3, pc = sx; pc > 0; --pc) {
                                *bf1 = alm * *bf1 + al * *pp; ++bf1;
                                *bf2 = alm * *bf2 + al * *pp; ++bf2;
                                *bf3 = alm * *bf3 + al * *pp++; ++bf3;
                            }
                        }
                    }

                } else {

                    if (ndp == hw3) {
                        for (cbf1 = &buffer[tx + ty * h], cbf2 = &cbf1[hw], cbf3 = &cbf2[hw], pp2 = &pp[sxy], pp3 = &pp2[sxy]; sy > 0; --sy, cbf1 = &cbf1[h], cbf2 = &cbf2[h], cbf3 = &cbf3[h]) {
                            for (bf1 = cbf1, bf2 = cbf2, bf3 = cbf3, pc = sx; pc > 0; --pc) {
                                al = (float) *ap++;
                                alm = 1.0 - al;
                                *bf1 = alm * *bf1 + al * *pp++; ++bf1;
                                *bf2 = alm * *bf2 + al * *pp2++; ++bf2;
                                *bf3 = alm * *bf3 + al * *pp3++; ++bf3;
                            }
                        }
                    } else {
                        for (cbf1 = &buffer[tx + ty * h], cbf2 = &cbf1[hw], cbf3 = &cbf2[hw]; sy > 0; --sy, cbf1 = &cbf1[h], cbf2 = &cbf2[h], cbf3 = &cbf3[h]) {
                            for (bf1 = cbf1, bf2 = cbf2, bf3 = cbf3, pc = sx; pc > 0; --pc) {
                                al = (float) *ap++;
                                alm = 1.0 - al;
                                *bf1 = alm * *bf1 + al * *pp; ++bf1;
                                *bf2 = alm * *bf2 + al * *pp; ++bf2;
                                *bf3 = alm * *bf3 + al * *pp++; ++bf3;
                            }
                        }
                    }
                }

            /* type -> single */
            } else {

                /* get data (as single -> float) */
                fpp = (const float *) mxGetData(inpp);

                /* depending on size of alpha */
                if (nda == 1) {

                    /* get alpha value */
                    al = (float) *ap;
                    alm = 1.0 - al;

                    /* depending on size of pixels */
                    if (ndp == sxy3) {

                        /* iterate over columns */
                        for (cbf1 = &buffer[tx + ty * h], cbf2 = &cbf1[hw], cbf3 = &cbf2[hw], fpp2 = &fpp[sxy], fpp3 = &fpp2[sxy]; sy > 0; --sy, cbf1 = &cbf1[h], cbf2 = &cbf2[h], cbf3 = &cbf3[h]) {

                            /* iterate over pixel */
                            for (bf1 = cbf1, bf2 = cbf2, bf3 = cbf3, pc = sx; pc > 0; --pc) {
                                *bf1 = alm * *bf1 + al * *fpp++; ++bf1;
                                *bf2 = alm * *bf2 + al * *fpp2++; ++bf2;
                                *bf3 = alm * *bf3 + al * *fpp3++; ++bf3;
                            }
                        }
                    } else {
                        for (cbf1 = &buffer[tx + ty * h], cbf2 = &cbf1[hw], cbf3 = &cbf2[hw]; sy > 0; --sy, cbf1 = &cbf1[h], cbf2 = &cbf2[h], cbf3 = &cbf3[h]) {
                            for (bf1 = cbf1, bf2 = cbf2, bf3 = cbf3, pc = sx; pc > 0; --pc) {
                                *bf1 = alm * *bf1 + al * *fpp; ++bf1;
                                *bf2 = alm * *bf2 + al * *fpp; ++bf2;
                                *bf3 = alm * *bf3 + al * *fpp++; ++bf3;
                            }
                        }
                    }

                } else {

                    if (ndp == sxy3) {
                        for (cbf1 = &buffer[tx + ty * h], cbf2 = &cbf1[hw], cbf3 = &cbf2[hw], fpp2 = &fpp[sxy], fpp3 = &fpp2[sxy]; sy > 0; --sy, cbf1 = &cbf1[h], cbf2 = &cbf2[h], cbf3 = &cbf3[h]) {
                            for (bf1 = cbf1, bf2 = cbf2, bf3 = cbf3, pc = sx; pc > 0; --pc) {
                                al = (float) *ap++;
                                alm = 1.0 - al;
                                *bf1 = alm * *bf1 + al * *fpp++; ++bf1;
                                *bf2 = alm * *bf2 + al * *fpp2++; ++bf2;
                                *bf3 = alm * *bf3 + al * *fpp3++; ++bf3;
                            }
                        }
                    } else {
                        for (cbf1 = &buffer[tx + ty * h], cbf2 = &cbf1[hw], cbf3 = &cbf2[hw]; sy > 0; --sy, cbf1 = &cbf1[h], cbf2 = &cbf2[h], cbf3 = &cbf3[h]) {
                            for (bf1 = cbf1, bf2 = cbf2, bf3 = cbf3, pc = sx; pc > 0; --pc) {
                                al = (float) *ap++;
                                alm = 1.0 - al;
                                *bf1 = alm * *bf1 + al * *fpp; ++bf1;
                                *bf2 = alm * *bf2 + al * *fpp; ++bf2;
                                *bf3 = alm * *bf3 + al * *fpp++; ++bf3;
                            }
                        }
                    }
                }
            }
        }
    }

    /* put into output */
    for (bf1 = buffer, pc = hw3; pc > 0; --pc)
        *opp++ = (unsigned char) (((float) 0.5) + *bf1++);

    /* free buffer */
    mxFree(buffer);

}

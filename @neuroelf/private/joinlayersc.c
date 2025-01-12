/*

joining of transimg layers into a HxWx3 uint8 image + HxW double alpha map

FORMAT:       [j, a] = joinlayersc(tio, l [, o])

Input fields:

      tio         1x1 struct of transimg content
      l           layer spec (numbers)
      o           override flag (default: false)

Output fields:

      j           HxWx3 single image [0.0 .. 255.0]
      a           HxW single alpha map [0.0 .. 1.0]

% Version:  v1.1b
% Build:    25011122
% Date:     Jan-11 2025, 10:18 PM EST
% Author:   Jochen Weber, NeuroElf
% URL/Info: http://neuroelf.net/

Copyright (c) 2010 - 2014, 2025, Jochen Weber
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
    int h = 0, w = 0, hw = 0, hw3 = 0, nl = 0, ln = 0, lc = 0,
        lfield = -1, pfield = -1, afield = -1,
        ndp = 0, nda = 0, pc = 0, nlist = 0;
    mwSize od[3] = {0, 0, 3};
    const mxArray *layer = NULL, *inpp = NULL, *inpa = NULL;

    /* pointers */
    const double *dp = NULL, *llist = NULL;
    const unsigned char *pp = NULL, *pp2 = NULL, *pp3 = NULL;
    const float *ap = NULL, *fpp = NULL, *fpp2 = NULL, *fpp3 = NULL;
    float *app = NULL, *appp = NULL, *opp = NULL;

    /* background color */
    float al = 1.0, alm = 0.0, oal = 1.0;

    /* override flag */
    bool override = 0;

    /* output pointers */
    float *bf1 = NULL, *bf2 = NULL, *bf3 = NULL;

    /* variable output string */
    /* char vstr[256]; */

    /* check number, type, fields of in/out arguments */
	if ((nrhs < 2) ||
        (nlhs != 2))
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (!mxIsStruct(*prhs))
        mexErrMsgTxt("First input must be of type struct.");
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("Second input must be of type double.");
    nlist = mxGetNumberOfElements(prhs[1]);
    llist = (const double *) mxGetData(prhs[1]);
    lfield = mxGetFieldNumber(*prhs, "Layer");
    if ((mxGetFieldNumber(*prhs, "Height") < 0) ||
        (lfield < 0) ||
        (mxGetFieldNumber(*prhs, "Width") < 0))
        mexErrMsgTxt("Required field missing.");
    layer = mxGetFieldByNumber(*prhs, 0, lfield);
    if (!mxIsStruct(layer))
        mexErrMsgTxt(".Layer field must be of type struct.");
    pfield = mxGetFieldNumber(layer, "Pixel");
    afield = mxGetFieldNumber(layer, "Alpha");
    if ((pfield < 0) ||
        (afield < 0))
        mexErrMsgTxt(".Layer field must have subfields .Pixel and .Alpha.");
    nl = mxGetNumberOfElements(layer);
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

    /* check override flag */
    if ((nrhs > 2) &&
        (mxIsLogical(prhs[2])) &&
        (mxGetNumberOfElements(prhs[2]) == 1)) {

        /* set to true if true */
        if ((* ((const unsigned char *) mxGetData(prhs[2]))) > 0)
            override = 1;
    }

    /* create output */
    *plhs = (mxArray *) mxCreateNumericArray(3, od, mxSINGLE_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating j output.");
    opp = (float *) mxGetData(*plhs);
    if (opp == NULL)
        mexErrMsgTxt("Error getting data pointer to j output mxArray.");
    plhs[1] = (mxArray *) mxCreateNumericArray(2, od, mxSINGLE_CLASS, mxREAL);
    if (plhs[1] == NULL) {
        mxDestroyArray(*plhs);
        mexErrMsgTxt("Error creating a output.");
    }
    app = (float *) mxGetData(plhs[1]);
    if (app == NULL) {
        mxDestroyArray(*plhs);
        mexErrMsgTxt("Error getting data pointer to a output mxArray.");
    }

    /* parse layers */
    for (lc = 0; lc < nlist; ++lc) {

        /* get pixel and alpha information */
        ln = (int) *llist++;
        if ((ln < 1) ||
            (ln > nl))
            continue;
        --ln;
        inpp = mxGetFieldByNumber(layer, ln, pfield);
        inpa = mxGetFieldByNumber(layer, ln, afield);
        ndp = mxGetNumberOfDimensions(inpp);
        nda = mxGetNumberOfDimensions(inpa);
        if (((!mxIsUint8(inpp)) &&
             (!mxIsSingle(inpp))) ||
            (!mxIsSingle(inpa)) ||
            (ndp > 3) ||
            (nda > 2)) {
            mxDestroyArray(*plhs);
            mxDestroyArray(plhs[1]);
            mexErrMsgTxt("Invalid .Pixel/.Alpha field types/dims.");
        }
        ndp = mxGetNumberOfElements(inpp);
        nda = mxGetNumberOfElements(inpa);
        if (((ndp != hw) &&
             (ndp != hw3)) ||
            ((nda != 1) &&
             (nda != hw))) {
            mxDestroyArray(*plhs);
            mxDestroyArray(plhs[1]);
            mexErrMsgTxt("Invalid .Pixel/.Alpha field sizes.");
        }
        ap = (const float *) mxGetData(inpa);

        /* skip invisible layers */
        if ((nda == 1) &&
            (*ap == 0.0))
            continue;

        /* depending on override flag -> override on (not default!) */
        if (override) {

            /* depending on type -> uint8 */
            if (mxIsUint8(inpp)) {

                /* get data (as uint8 -> unsigned char) */
                pp = (const unsigned char *) mxGetData(inpp);

                /* depending on size of alpha */
                if (nda == 1) {

                    /* get alpha value */
                    al = *ap;

                    /* depending on size of pixels */
                    if (ndp == hw3) {

                        /* iterate over pixel */
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], pp2 = &pp[hw], pp3 = &pp2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            oal = *appp * (1.0 - al);
                            *bf1 = oal * *bf1 + al * ((float) *pp++); ++bf1;
                            *bf2 = oal * *bf2 + al * ((float) *pp2++); ++bf2;
                            *bf3 = oal * *bf3 + al * ((float) *pp3++); ++bf3;
                            *appp++ = oal + al;
                        }
                    } else {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            oal = *appp * (1.0 - al);
                            alm = al * ((float) *pp++);
                            *bf1 = oal * *bf1 + alm; ++bf1;
                            *bf2 = oal * *bf2 + alm; ++bf2;
                            *bf3 = oal * *bf3 + alm; ++bf3;
                            *appp++ = oal + al;
                        }
                    }

                } else {

                    if (ndp == hw3) {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], pp2 = &pp[hw], pp3 = &pp2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            al = *ap++;
                            if (al > 0.0) {
                                oal = *appp * (1.0 - al);
                                *bf1 = oal * *bf1 + al * ((float) *pp++); ++bf1;
                                *bf2 = oal * *bf2 + al * ((float) *pp2++); ++bf2;
                                *bf3 = oal * *bf3 + al * ((float) *pp3++); ++bf3;
                                *appp++ = oal + al;
                            } else {
                                ++bf1;
                                ++bf2;
                                ++bf3;
                                ++appp;
                                ++pp;
                                ++pp2;
                                ++pp3;
                            }
                        }
                    } else {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            al = *ap++;
                            if (al > 0.0) {
                                oal = *appp * (1.0 - al);
                                alm = al * ((float) *pp++);
                                *bf1 = oal * *bf1 + alm; ++bf1;
                                *bf2 = oal * *bf2 + alm; ++bf2;
                                *bf3 = oal * *bf3 + alm; ++bf3;
                                *appp++ = oal + al;
                            } else {
                                ++bf1;
                                ++bf2;
                                ++bf3;
                                ++appp;
                                ++pp;
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
                    al = *ap;

                    /* depending on size of pixels */
                    if (ndp == hw3) {

                        /* iterate over pixel */
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], fpp2 = &fpp[hw], fpp3 = &fpp2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            oal = *appp * (1.0 - al);
                            *bf1 = oal * *bf1 + al * *fpp++; ++bf1;
                            *bf2 = oal * *bf2 + al * *fpp2++; ++bf2;
                            *bf3 = oal * *bf3 + al * *fpp3++; ++bf3;
                            *appp++ = oal + al;
                        }
                    } else {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            oal = *appp * (1.0 - al);
                            alm = al * *fpp++;
                            *bf1 = oal * *bf1 + alm; ++bf1;
                            *bf2 = oal * *bf2 + alm; ++bf2;
                            *bf3 = oal * *bf3 + alm; ++bf3;
                            *appp++ = oal + al;
                        }
                    }

                } else {

                    if (ndp == hw3) {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], fpp2 = &fpp[hw], fpp3 = &fpp2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            al = *ap++;
                            if (al > 0.0) {
                                oal = *appp * (1.0 - al);
                                *bf1 = oal * *bf1 + al * *fpp++; ++bf1;
                                *bf2 = oal * *bf2 + al * *fpp2++; ++bf2;
                                *bf3 = oal * *bf3 + al * *fpp3++; ++bf3;
                                *appp++ = oal + al;
                            } else {
                                ++bf1;
                                ++bf2;
                                ++bf3;
                                ++appp;
                                ++fpp;
                                ++fpp2;
                                ++fpp3;
                            }
                        }
                    } else {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            al = *ap++;
                            if (al > 0.0) {
                                oal = *appp * (1.0 - al);
                                alm = al *  *fpp++;
                                *bf1 = oal * *bf1 + alm; ++bf1;
                                *bf2 = oal * *bf2 + alm; ++bf2;
                                *bf3 = oal * *bf3 + alm; ++bf3;
                                *appp++ = oal + al;
                            } else {
                                ++bf1;
                                ++bf2;
                                ++bf3;
                                ++appp;
                                ++fpp;
                            }
                        }
                    }
                }
            }

        /* override off (default!) */
        } else {

            /* depending on type -> uint8 */
            if (mxIsUint8(inpp)) {

                /* get data (as uint8 -> unsigned char) */
                pp = (const unsigned char *) mxGetData(inpp);

                /* depending on size of alpha */
                if (nda == 1) {

                    /* get alpha value */
                    al = *ap;

                    /* depending on size of pixels */
                    if (ndp == hw3) {

                        /* iterate over pixel */
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], pp2 = &pp[hw], pp3 = &pp2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            *bf1++ += al * ((float) *pp++);
                            *bf2++ += al * ((float) *pp2++);
                            *bf3++ += al * ((float) *pp3++);
                            *appp++ += al;
                        }
                    } else {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            alm = al * ((float) *pp++);
                            *bf1++ += alm;
                            *bf2++ += alm;
                            *bf3++ += alm;
                            *appp++ += al;
                        }
                    }

                } else {

                    if (ndp == hw3) {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], pp2 = &pp[hw], pp3 = &pp2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            al = *ap++;
                            *bf1++ += al * ((float) *pp++);
                            *bf2++ += al * ((float) *pp2++);
                            *bf3++ += al * ((float) *pp3++);
                            *appp++ += al;
                        }
                    } else {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            al = *ap++;
                            alm = al * ((float) *pp++);
                            *bf1++ += alm;
                            *bf2++ += alm;
                            *bf3++ += alm;
                            *appp++ += al;
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
                    al = *ap;

                    /* depending on size of pixels */
                    if (ndp == hw3) {

                        /* iterate over pixel */
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], fpp2 = &fpp[hw], fpp3 = &fpp2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            *bf1++ += al * *fpp++;
                            *bf2++ += al * *fpp2++;
                            *bf3++ += al * *fpp3++;
                            *appp++ += al;
                        }
                    } else {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            alm = al * *fpp++;
                            *bf1++ += alm;
                            *bf2++ += alm;
                            *bf3++ += alm;
                            *appp++ += al;
                        }
                    }

                } else {

                    if (ndp == hw3) {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], fpp2 = &fpp[hw], fpp3 = &fpp2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            al = *ap++;
                            *bf1++ += al * *fpp++;
                            *bf2++ += al * *fpp2++;
                            *bf3++ += al * *fpp3++;
                            *appp++ += al;
                        }
                    } else {
                        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], appp = app, pc = hw; pc > 0; --pc) {
                            al = *ap++;
                            alm = al * *fpp++;
                            *bf1++ += alm;
                            *bf2++ += alm;
                            *bf3++ += alm;
                            *appp++ += al;
                        }
                    }
                }
            }
        }
    }

    /* max scale alpha to 1.0 */
    if (override)
        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], pc = hw; pc > 0; --pc, ++app, ++bf1, ++bf2, ++bf3) {
            al = *app;
            if (al > 0.0) {
                *bf1 /= al;
                *bf2 /= al;
                *bf3 /= al;
            } else {
                *bf1 = 0.0;
                *bf2 = 0.0;
                *bf3 = 0.0;
            }
        }
    else
        for (bf1 = opp, bf2 = &bf1[hw], bf3 = &bf2[hw], pc = hw; pc > 0; --pc, ++app, ++bf1, ++bf2, ++bf3) {
            al = *app;
            if (al > 0.0) {
                *bf1 /= al;
                *bf2 /= al;
                *bf3 /= al;
            }
            if (al > 1.0)
                *app = 1.0;
        }
}

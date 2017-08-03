/*

regularly divide triangles of a (closed surface) mesh 4:1

FORMAT:       [tri, coord] = mesh_tridivide(tri [, coord [, dtmcomp]])

Input fields:

      tri         Tx3 triangle vertices (one-based)
      coord       Cx3 coordinates (only needed for second output)
      dtmcomp     1x1 boolean flag, dividetrimesh.m compatibility (true)

Output fields:

      tri         4Tx3 divided triangles
      coord       UCx3 new coordinates with triangles


% Version:  v0.9d
% Build:    14071115
% Date:     Jul-11 2014, 3:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2014, Jochen Weber
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

#define MAXNRNEI 30
#define MAXNRNEI2 60

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* input pointers for triangles and coordinates */
    const double *t1, *t2, *t3, *crd1, *crd2, *crd3;
    const int *di = NULL;

    /* number of triangles and coordinates */
    int nt, nt3, nt4, nc = 0, nv = 0, nnv, sp1, sp2, sp3, tc, tc1, tc1m, tc2, tc2m, tc3, tc3m;
    double ti;
    char vc;

    /* flags */
    bool rc      = 0, /* recompute coordinates */
         split   = 0, /* edge already split */
         dtmcomp = 1; /* dividetrimesh.m compatibility */

    /* output pointer */
    double *ot1, *ot1n, *ot2, *ot2n, *ot3, *ot3n;

    /* backup pointer of triangles, lookup arrays */
    double *ot = NULL, *otb1, *otb2, *otb3;
    int *vl = NULL, *vl2;
    int *sp = NULL, *spx, *spi = NULL;

    /* variable output string */
    /* char vstr[256]; */

    /* check number of in/out arguments */
    if ((nrhs < 1) || (nlhs > 2) || ((nrhs < 2) && (nlhs > 1)))
        mexErrMsgTxt("Bad number of input/output arguments.");

    /* check argument types and sizes, and get data pointers */
    if (!mxIsDouble(*prhs) || ((nrhs > 1) && (!mxIsDouble(prhs[1]))))
        mexErrMsgTxt("Input arguments must be of type double.");
    if ((mxGetNumberOfElements(*prhs) < 3) || ((nrhs > 1) && (mxGetNumberOfElements(prhs[1]) < 3)))
        mexErrMsgTxt("Inputs must be Tx3 (and Cx3) in size.");
    di = (const int *) mxGetDimensions(*prhs);
    if (di[1] != 3)
        mexErrMsgTxt("Triangles must be Tx3 in size.");
    nt = *di;
    nt3 = 3 * nt;
    nt4 = 4 * nt;
    if (nrhs > 1) {
        di = (const int *) mxGetDimensions(prhs[1]);
        if (di[1] != 3)
            mexErrMsgTxt("Coordinates must be Cx3 in size.");
        nc = *di;
        if ((nrhs > 2) &&
            (mxGetClassID(prhs[2]) == mxLOGICAL_CLASS) &&
            (mxGetNumberOfElements(prhs[2]) == 1))
            dtmcomp = (*( (const char *) mxGetData(prhs[2])) != 0) ? 1 : 0;
    }
    t1 = (const double *) mxGetData(*prhs);
    t2 = &t1[nt];
    t3 = &t2[nt];

    /* create first output argument */
    *plhs = mxCreateDoubleMatrix(4 * nt, 3, mxREAL);
    ot1 = (double *) mxGetData(*plhs);
    if (ot1 == NULL)
        mexErrMsgTxt("Error creating first output array.");
    ot2 = &ot1[4 * nt];
    ot3 = &ot2[4 * nt];

    /* recompute coordinates */
    if (nlhs > 1)
        rc = 1;

    /* copy triangles, and inquire what is the last indexed vertex */
    ot = (double *) mxCalloc(nt3, sizeof(double));
    if (ot == NULL)
        mexErrMsgTxt("Error allocating triangles copy buffer.");
    otb1 = ot;
    otb2 = &otb1[nt];
    otb3 = &otb2[nt];
    for (tc = 0; tc < nt; ++tc) {
        ti = t1[tc];
        *otb1++ = ti;
        if (ti > nv)
            nv = ti;
        ti = t2[tc];
        *otb2++ = ti;
        if (ti > nv)
            nv = ti;
        ti = t3[tc];
        *otb3++ = ti;
        if (ti > nv)
            nv = ti;
    }
    if ((rc == 1) &&
        (nv > nc)) {
        mxFree(ot);
        mexErrMsgTxt("Triangles reference unknown coordinates.");
    }
    otb1 = ot;
    otb2 = &otb1[nt];
    otb3 = &otb2[nt];
    if (dtmcomp) {
        ot1n = &ot1[nt];
        ot2n = &ot2[nt];
        ot3n = &ot3[nt];
    }

    /* create lookup arrays -> for vertex recomputation */
    if (rc == 1) {
        vl = (int *) mxCalloc(8 * nt + 2, sizeof(int));
        if (vl == NULL) {
            mxFree(ot);
            mexErrMsgTxt("Error allocating lookup array.");
        }
        vl2 = &vl[nt4 + 1];
    }

    /* create lookup arrays for splits */
    sp = (int *) mxCalloc(nv * MAXNRNEI2, sizeof(int));
    if (sp == NULL) {
        if (rc == 1)
            mxFree(vl);
        mxFree(ot);
        mexErrMsgTxt("Error allocating splits array.");
    }
    spi = (int *) mxCalloc(2 * nv, sizeof(int));
    if (spi == NULL) {
        mxFree(sp);
        if (rc == 1)
            mxFree(vl);
        mxFree(ot);
        mexErrMsgTxt("Error allocating splits indexing array.");
    }
    
    /* initialize counter */
    nnv = nv;

    /* for each triangle */
    for (tc = 0; tc < nt; ++tc) {

        /* get the three vertices */
        tc1 = ((int) *otb1++);
        tc1m = tc1 - 1;
        tc2 = ((int) *otb2++);
        tc2m = tc2 - 1;
        tc3 = ((int) *otb3++);
        tc3m = tc3 - 1;

        /* first pair (P1, P2) not already split? */
        for (vc = 0, split = 0, spx = &sp[tc1m * MAXNRNEI2]; vc <= spi[tc1m]; ++vc) {
            if (*spx == tc2) {
                split = 1;
                break;
            }
            ++spx;
        }
        if (split)
            sp1 = spx[MAXNRNEI];
        else {
            ++nnv;
            sp1 = (double) nnv;
            --spx;
            *spx = tc2;
            spx[MAXNRNEI] = sp1;
            ++spi[tc1m];
            sp[tc2m * MAXNRNEI2 + spi[tc2m]] = tc1;
            sp[tc2m * MAXNRNEI2 + spi[tc2m] + MAXNRNEI] = sp1;
            ++spi[tc2m];
        }

        /* repeat for second and third pairs */
        for (vc = 0, split = 0, spx = &sp[tc2m * MAXNRNEI2]; vc <= spi[tc2m]; ++vc) {
            if (*spx == tc3) {
                split = 1;
                break;
            }
            ++spx;
        }
        if (split)
            sp2 = spx[MAXNRNEI];
        else {
            ++nnv;
            sp2 = (double) nnv;
            --spx;
            *spx = tc3;
            spx[MAXNRNEI] = sp2;
            ++spi[tc2m];
            sp[tc3m * MAXNRNEI2 + spi[tc3m]] = tc2;
            sp[tc3m * MAXNRNEI2 + spi[tc3m] + MAXNRNEI] = sp2;
            ++spi[tc3m];
        }
        for (vc = 0, split = 0, spx = &sp[tc3m * MAXNRNEI2]; vc <= spi[tc3m]; ++vc) {
            if (*spx == tc1) {
                split = 1;
                break;
            }
            ++spx;
        }
        if (split)
            sp3 = spx[MAXNRNEI];
        else {
            ++nnv;
            sp3 = (double) nnv;
            --spx;
            *spx = tc1;
            spx[MAXNRNEI] = sp3;
            ++spi[tc3m];
            sp[tc1m * MAXNRNEI2 + spi[tc1m]] = tc3;
            sp[tc1m * MAXNRNEI2 + spi[tc1m] + MAXNRNEI] = sp3;
            ++spi[tc1m];
        }

        /* replace original triangle with new one */
        *ot1++ = (double) tc1;
        *ot2++ = (double) sp1;
        *ot3++ = (double) sp3;

        /* and add three new triangles */
        if (dtmcomp) {
            *ot1n++ = (double) tc2;
            *ot2n++ = (double) sp2;
            *ot3n++ = (double) sp1;
            *ot1n++ = (double) tc3;
            *ot2n++ = (double) sp3;
            *ot3n++ = (double) sp2;
            *ot1n++ = (double) sp1;
            *ot2n++ = (double) sp2;
            *ot3n++ = (double) sp3;
        } else {
            *ot1++ = (double) tc2;
            *ot2++ = (double) sp2;
            *ot3++ = (double) sp1;
            *ot1++ = (double) tc3;
            *ot2++ = (double) sp3;
            *ot3++ = (double) sp2;
            *ot1++ = (double) sp1;
            *ot2++ = (double) sp2;
            *ot3++ = (double) sp3;
        }

        /* add recomputation vectors */
        if (rc) {
            vl[sp1] = tc1m;
            vl2[sp1] = tc2m;
            vl[sp2] = tc2m;
            vl2[sp2] = tc3m;
            vl[sp3] = tc3m;
            vl2[sp3] = tc1m;
        }
    }
    mxFree(sp);
    mxFree(ot);

    /* recompute vertices */
    if (rc) {
        plhs[1] = mxCreateDoubleMatrix(nnv, 3, mxREAL);
        if (plhs[1] == NULL) {
            mxFree(vl);
            mexErrMsgTxt("Error allocating output array for coordinates.");
        }
        ot1 = (double *) mxGetData(plhs[1]);
        ot2 = &ot1[nnv];
        ot3 = &ot2[nnv];
        crd1 = (const double *) mxGetData(prhs[1]);
        crd2 = &crd1[nc];
        crd3 = &crd2[nc];
        for (tc = 0; tc < nv; ++tc) {
            *ot1++ = *crd1++;
            *ot2++ = *crd2++;
            *ot3++ = *crd3++;
        }
        crd1 = (const double *) mxGetData(prhs[1]);
        crd2 = &crd1[nc];
        crd3 = &crd2[nc];
        for (++tc; tc <= nnv; ++tc) {
            *ot1++ = 0.5 * (crd1[vl[tc]] + crd1[vl2[tc]]);
            *ot2++ = 0.5 * (crd2[vl[tc]] + crd2[vl2[tc]]);
            *ot3++ = 0.5 * (crd3[vl[tc]] + crd3[vl2[tc]]);
        }
    }
    mxFree(vl);
}

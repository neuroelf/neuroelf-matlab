/*

using triangles of mesh to create neighbors list

FORMAT:       [nei, bn [, trl]] = mesh_trianglestoneighbors(nc, tri)

Input fields:

      nc          1x1 number of vertices
      tri         Cx3 triangles (one-based)

Output fields:

      nei         Nx2 neighbors list (1-based)
      bn          1xN neighbors list with warnings
      tri         Nx1 triangle list (1-based)


% Version:  v0.9d
% Build:    14062814
% Date:     Jun-28 2014, 2:34 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2010, 2014, Jochen Weber
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

#define MAXNRNEI 30
#define MAXNRNE2 60
#define MAXNRNEM 59

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* number of coordinates and counter (including double + triple) */
    double dnc, dn1, dn2, dn3;
    double *dnei, *snei;
    double dneib[MAXNRNEI], dneis[MAXNRNEI];
    int nc, cc, nt, nt2, sc;

    /* neighbors list pointer */
    int *nei, *neip;
    unsigned char *nnei, *nneic;
    unsigned char cinei, cnnei, tnnei;
    int n1, n2, n3;
    mxArray *cnei;
    mxArray *onnei;
    mxArray *onei;

    /* double array for bad neighbors */
    bool trackbn = 0;
    double *bnd = NULL;
    signed long bnc = 0;
    unsigned long *bnl = NULL;

    /* number of dimensions and dim size */
    int nd;
    const int *di;
    int cdi[2] = {1, 2};
    int ndi[2] = {1, 1};

    /* variable output string */
    char vstr[256];

    /* check number of in/out arguments */
    if ((nrhs != 2) | (nlhs > 3))
            mexErrMsgTxt("Bad number of input/output arguments.");

    /* check argument types */
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || (mxGetNumberOfElements(prhs[0]) != 1))
        mexErrMsgTxt("Input arguments must be of type double.");
    dnc = *((double *) mxGetPr(prhs[0]));

    /* good number */
    if (mxIsInf(dnc) || mxIsNaN(dnc) || dnc < 1 || dnc > 10000000 || dnc != ((double) ((int) dnc)))
        mexErrMsgTxt("Bad number of coordinates given.");
    nc = (int) dnc;

    /* check dims */
    nd = mxGetNumberOfDimensions(prhs[1]);
    di = mxGetDimensions(prhs[1]);
    if (nd != 2 || di[1] != 3)
        mexErrMsgTxt("Bad size of triangles argument");
    nt = di[0];
    nt2 = nt * 2;
    dnei = mxGetPr(prhs[1]);

    /* create internal neighbors list */
    nei = (int *) mxCalloc(MAXNRNE2 * nc, sizeof(int));
    if (nei == NULL)
        mexErrMsgTxt("Error allocating memory to store neighbors list.");
    nnei = (unsigned char *) mxCalloc(nc, sizeof(unsigned char));
    if (nnei == NULL)
        mexErrMsgTxt("Error allocating memory to store number of neighbors list.");
    nneic = (unsigned char *) mxCalloc(nc, sizeof(unsigned char));
    if (nneic == NULL)
        mexErrMsgTxt("Error allocating memory to store number of neighbors list.");

    /* create output argument */
    cdi[0] = nc;
    plhs[0] = mxCreateCellArray(2, cdi);
    if (plhs[0] == NULL)
        mexErrMsgTxt("Error allocating memory for cell array.");

    /* create second output */
    if (nlhs > 1) {
        bnl = (unsigned long *) mxCalloc(nc, sizeof(unsigned long));
        if (bnl == NULL)
            mexErrMsgTxt("Error allocating space for bad neighbors list.");
        trackbn = 1;
    }

    /* create third output */
    if (nlhs > 2) {
        cdi[1] = 1;
        plhs[2] = mxCreateCellArray(2, cdi);
        if (plhs[2] == NULL)
            mexErrMsgTxt("Error allocating memory for supplementary cell array.");
    }

    /* iterate over triangles */
    for (cc = 0; cc < nt; ++cc) {

        /* get three vertix numbers */
        dn1 = dnei[cc];
        dn2 = dnei[cc + nt];
        dn3 = dnei[cc + nt2];

        /* bail out if invalid */
        if (mxIsInf(dn1) || mxIsInf(dn2) || mxIsInf(dn3) ||
            mxIsNaN(dn1) || mxIsNaN(dn2) || mxIsNaN(dn3) ||
            dn1 < 1 || dn1 > nc || dn2 < 1 || dn2 > nc || dn3 < 1 || dn3 > nc)
            mexErrMsgTxt("Bad triangles argument.");

        /* get integer version */
        n1 = (int) (dn1 - 1);
        n2 = (int) (dn2 - 1);
        n3 = (int) (dn3 - 1);

        /* get offset to store neighbors of first vertex */
        sc = MAXNRNE2 * n1 + 2 * (nnei[n1]);

        /* check if number of neighbors is already too many */
        if (++nnei[n1] > MAXNRNEI) {
            sprintf(vstr, "Too many triangles at vertex %d.", n1 + 1);
            mexErrMsgTxt(vstr);
        }

        /* store neighbors */
        nei[sc++] = n2 + 1;
        nei[sc] = n3 + 1;

        /* repeat for second vertex */
        sc = MAXNRNE2 * n2 + 2 * (nnei[n2]);
        if (++nnei[n2] > MAXNRNEI) {
            sprintf(vstr, "Too many triangles at vertex %d.", n2 + 1);
            mexErrMsgTxt(vstr);
        }
        nei[sc++] = n3 + 1;
        nei[sc] = n1 + 1;

        /* and third vertex also */
        sc = MAXNRNE2 * n3 + 2 * (nnei[n3]);
        if (++nnei[n3] > MAXNRNEI) {
            sprintf(vstr, "Too many triangles at vertex %d.", n3 + 1);
            mexErrMsgTxt(vstr);
        }
        nei[sc++] = n1 + 1;
        nei[sc] = n2 + 1;
    }

    /* create third output first */
    if (nlhs > 2) {

        /* create necessary cell arrays */
        ndi[0] = 1;
        cnei = plhs[2];
        for (cc = 0; cc < nc; ++cc) {
            ndi[1] = nnei[cc];
            onei = mxCreateNumericArray(2, ndi, mxDOUBLE_CLASS, mxREAL);
            if (onei == NULL) {
                sprintf(vstr, "Error allocating triangle list array for vertex %d.", cc + 1);
                mexErrMsgTxt(vstr);
            }
            mxSetCell(cnei, cc, onei);
        }

        /* iterate over triangles */
        for (cc = 0; cc < nt; ++cc) {

            /* set good triangle number */
            dn2 = (double) (cc + 1);

            /* get three vertix numbers and store in according array */
            dn1 = dnei[cc];
            n1 = (int) (dn1 - 1);
            snei = (double *) mxGetPr(mxGetCell(cnei, n1));
            snei[nneic[n1]++] = dn2;

            dn1 = dnei[cc + nt];
            n1 = (int) (dn1 - 1);
            snei = (double *) mxGetPr(mxGetCell(cnei, n1));
            snei[nneic[n1]++] = dn2;

            dn1 = dnei[cc + nt2];
            n1 = (int) (dn1 - 1);
            snei = (double *) mxGetPr(mxGetCell(cnei, n1));
            snei[nneic[n1]++] = dn2;
        }
    }

    /* fill output array */
    cnei = plhs[0];
    cdi[0] = 1;
    ndi[1] = 1;
    for (cc = 0; cc < nc; ++cc) {

        /* set value */
        tnnei = nnei[cc];
        if (tnnei == 1)
            tnnei = 2;

        /* get offset into neighbors list array */
        neip = &nei[MAXNRNE2 * cc];

        /* get first neighbor */
        n1 = *neip;

        /* store into dnei buffer */
        *dneib = (double) n1;

        /* initialize counter to 1 */
        cinei = 1;

        /* get second neighbor and repeat until same as first */
        for (n2 = neip[1]; n2 != n1; ) {

            /* store next neighbor */
            dneib[cinei] = (double) n2;

            /* find next pair */
            for (cnnei = 0; cnnei < MAXNRNE2; cnnei += 2) {
                if (neip[cnnei] == n2) {
                    n2 = neip[++cnnei];
                    ++cinei;
                    break;
                }
            }

            /* if we haven't found all neighbors, search backwards with orphan */
            if ((cnnei > MAXNRNEM) || (cinei > tnnei)) {

                /* we must increase the number of neighbors first */
                sprintf(vstr, "Trying to resolve neighbors for vertex %d...", cc + 1);
                mexWarnMsgTxt(vstr);
                tnnei = nnei[cc] + 1;

                /* this time start with the last one (orphan), tracking back */
                dneis[0] = n2;
                cinei = 1;

                /* find next pair */
                for (cnnei = 2; cnnei <= MAXNRNE2; cnnei += 2) {
                    if (neip[cnnei-1] == n2) {
                        n2 = neip[cnnei-2];
                        dneis[cinei++] = n2;
                        cnnei = 0;
                        if (cinei >= MAXNRNEI)
                            break;
                    }
                }

                /* now, cinei must be the correct number of neighbors */
                if (cinei != tnnei) {
                    sprintf(vstr, "Unresolvable neighbors for vertex %d (%d, %d).", cc + 1, tnnei, cinei);
                    mexErrMsgTxt(vstr);
                }

                /* now, we must reverse the order ! */
                for (cnnei = 0; cnnei < tnnei; ++cnnei) {
                    dneib[cnnei] = dneis[tnnei-(1+cnnei)];
                }

                /* and not continue with the first for loop ! */
                break;
            }

            /* check total counter */
            if (cinei > tnnei) {
                sprintf(vstr, "Too many neighbors for vertex %d (%d, %d).", cc + 1, tnnei, cinei);
                mexErrMsgTxt(vstr);
            }

        }

        if (cinei != tnnei) {
            if (trackbn)
                bnl[bnc++] = cc + 1;
            sprintf(vstr, "Invalid neighborhood for vertex %d; cutting back.", cc + 1);
            mexWarnMsgTxt(vstr);
            tnnei = cinei;
        }

        /* create number of neighbors 1x1 array */
        onnei = mxCreateNumericArray(2, ndi, mxDOUBLE_CLASS, mxREAL);
        if (onnei == NULL)
            mexErrMsgTxt("Cannot allocate number of neighbors scalar.");
        *(mxGetPr(onnei)) = (double) tnnei;
        nnei[cc] = tnnei;

        /* create matching output array and copy */
        cdi[1] = tnnei;
        onei = mxCreateNumericArray(2, cdi, mxDOUBLE_CLASS, mxREAL);
        if (onei == NULL)
            mexErrMsgTxt("Cannot allocate list of neighbors array.");
        dnei = (double *) mxGetPr(onei);
        if (dnei == NULL)
            mexErrMsgTxt("Error getting double pointer to store neighbors.");
        for (cnnei = 0; cnnei < tnnei; ++cnnei)
            *dnei++ = dneib[cnnei];

        /* set to cell */
        mxSetCell(cnei, cc, onnei);
        mxSetCell(cnei, cc + nc, onei);
    }

    /* return bad neighbors list */
    if (trackbn) {
        if (bnc == 0)
            ndi[0] = 0;
        ndi[1] = bnc;
        plhs[1] = mxCreateNumericArray(2, ndi, mxDOUBLE_CLASS, mxREAL);
        if (plhs[1] == NULL)
            mexErrMsgTxt("Error allocating output array for bad neighbors list.");
        bnd = (double*) mxGetPr(plhs[1]);
        for (--bnc; bnc >= 0; --bnc)
            bnd[bnc] = (double) bnl[bnc];
        mxFree(bnl);
    }

}

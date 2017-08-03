/*

creating a CxN neighbors list

FORMAT:       na = mesh_neighborsarray(nei [, dp])

Input fields:

      nei         Cx2 neighbors array
      dp          default value (either {0} or 1; if 1, use own point)

Output fields:

      na          CxN neighbors array


% Version:  v0.9d
% Build:    14061514
% Date:     Jun-15 2014, 2:57 PM EST
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
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* input dimensions */
    const int *idim;
    const double *val;
    const mxArray *cell;

    /* counter and number of points */
    int c, cc, mc, nc;

    /* max number of neighbors */
    signed char maxnei = 0, sc;

    /* default point */
    bool dpself = 0;

    /* output array */
    double *out = NULL, *out2 = NULL, dc;

    /* temporary array */
    int *neigh = NULL, *na;

    /* char vstr[256]; */

    /* check number of in/out arguments */
    if ((nrhs < 1) | (nlhs > 3))
            mexErrMsgTxt("Bad number of input/output arguments.");

    /* check argument type */
    if (!mxIsCell(*prhs))
        mexErrMsgTxt("Input argument must be of type cell.");

    /* and dims */
    if (mxGetNumberOfDimensions(*prhs) > 2)
        mexErrMsgTxt("Input must be 2D.");
    idim = mxGetDimensions(*prhs);
    if (idim[1] != 2)
        mexErrMsgTxt("Input must be Cx2 in size.");
    nc = *idim;

    /* default point */
    if ((nrhs > 1) &&
        (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfElements(prhs[1]) == 1)) {

        val = (const double *) mxGetData(prhs[1]);
        if (!mxIsInf(*val) &&
            !mxIsNaN(*val) &&
            (*val == 1.0))
            dpself = 1;
    }

    /* one output argument */
    if (nlhs == 1) {

        /* empty set */
        if (nc == 0) {

            /* create output */
            *plhs = mxCreateDoubleMatrix(nc, 12, mxREAL);
            return;
        }
        mc = 12 * nc;

        /* allocate temporary memory */
        neigh = mxCalloc(12 * nc, sizeof(int));
        if (neigh == NULL)
            mexErrMsgTxt("Error allocating temporary neighbors array.");

        /* loop over neighbors (second column!) */
        for (c = 0; c < nc; ++c) {

            /* fill neighbors with default value? */
            if (dpself)
                for (cc = c; cc < mc; cc += nc)
                    neigh[cc] = c + 1;
            else
                for (cc = c; cc < mc; cc += nc)
                    neigh[cc] = 0;

            /* get cell */
            cell = mxGetCell(*prhs, c + nc);

            /* check cell */
            if (cell == NULL)
                continue;
            cc = mxGetNumberOfElements(cell);
            if ((mxGetClassID(cell) != mxDOUBLE_CLASS) ||
                (cc > 12)) {
                mxFree(neigh);
                mexErrMsgTxt("Invalid neighborhood data.");
            }

            /* get cell contents */
            val = mxGetData(cell);

            /* check data */
            if (val == NULL)
                continue;

            /* max number */
            sc = cc;
            if (maxnei < sc)
                maxnei = sc;

            /* copy data */
            for (cc = c; sc > 0; --sc, cc += nc)
                neigh[cc] = (int) *val++;
        }

        /* create final output */
        *plhs = mxCreateDoubleMatrix(nc, (int) maxnei, mxREAL);
        if (*plhs == NULL) {
            mxFree(neigh);
            mexErrMsgTxt("Error allocating output array.");
        }
        out = (double *) mxGetData(*plhs);
        if (out == NULL) {
            mxFree(neigh);
            mexErrMsgTxt("Error addressing output array.");
        }

        /* copy as many elements as needed */
        for (na = neigh, cc = ((int) maxnei) * nc; cc > 0; --cc)
            *out++ = (double) *na++;

        /* release neighbor array */
        mxFree(neigh);

    /* more than one output */
    } else {

        /* two outputs */
        if (nlhs == 2) {

            /* first count the number of total neighbors */
            for (c = nc-1, mc = 0; c >= 0; --c) {

                /* get data */
                cell = (const mxArray *) mxGetCell(*prhs, c + nc);

                /* then count */
                mc += mxGetNumberOfElements(cell);
            }

        /* three outputs */
        } else {

            /* first, create third output */
            plhs[2] = mxCreateDoubleMatrix(nc, 1, mxREAL);
            if (plhs[2] == NULL)
                mexErrMsgTxt("Error creating third output array.");
            out2 = (double *) mxGetData(plhs[2]);
            if (out2 == NULL)
                mexErrMsgTxt("Error getting third output pointer.");

            /* then loop */
            c = nc-1;
            for (mc = 0, out2 = &out2[c]; c >= 0; --c) {

                /* get data */
                cell = (const mxArray *) mxGetCell(*prhs, c + nc);

                /* then count */
                cc = mxGetNumberOfElements(cell);
                mc += cc;

                /* and store info */
                *out2-- = (double) cc;
            }
        }

        /* then allocate outputs */
        if (dpself)
            *plhs = mxCreateDoubleMatrix(mc + nc, 1, mxREAL);
        else
            *plhs = mxCreateDoubleMatrix(mc, 1, mxREAL);
        if (*plhs == NULL)
            mexErrMsgTxt("Error creating first output array.");
        out = (double *) mxGetData(*plhs);
        if (out == NULL)
            mexErrMsgTxt("Error getting first output pointer.");
        if (dpself)
            plhs[1] = mxCreateDoubleMatrix(mc + nc, 1, mxREAL);
        else
            plhs[1] = mxCreateDoubleMatrix(mc, 1, mxREAL);
        if (plhs[1] == NULL)
            mexErrMsgTxt("Error creating second output array.");
        out2 = (double *) mxGetData(plhs[1]);
        if (out2 == NULL)
            mexErrMsgTxt("Error getting second output pointer.");

        /* then work backwards (again) */
        if (dpself) {
            mc += nc - 1;
            for (c = nc, out = &out[mc], out2 = &out2[mc]; c > 0; --c) {

                /* get contents again */
                cell = (const mxArray *) mxGetCell(*prhs, c + nc - 1);

                /* get pointer */
                val = (const double *) mxGetData(cell);
                for (cc = mxGetNumberOfElements(cell) - 1, dc = (double) c; cc >= 0; --cc) {
                    *out-- = dc;
                    *out2-- = val[cc];
                }

                /* add self */
                *out-- = dc;
                *out2-- = dc;
            }
        } else {
            --mc;
            for (c = nc, out = &out[mc], out2 = &out2[mc]; c > 0; --c) {

                /* get contents again */
                cell = (const mxArray *) mxGetCell(*prhs, c + nc - 1);

                /* get pointer */
                val = (const double *) mxGetData(cell);
                for (cc = mxGetNumberOfElements(cell) - 1, dc = (double) c; cc >= 0; --cc) {
                    *out-- = (double) c;
                    *out2-- = val[cc];
                }
            }
        }
    }
}

/*

set sparse double array values

FORMAT:       sp = setsparseval(v)

FORMAT:       sp = setsparseval(v, sp)

FORMAT:       sp = setsparseval(v, m, n, i, j)

Input fields:

      v           Vx1 values, if single argument, create VxV diagonal matrix
      sp          sparse double array with exactly Vx1 non-zero elements
      m, n        1x1 arguments for size of newly to be created matrix
      i, j        row and column values for sparse array

Output fields:

      sp          MxN (or original size) sparse array

Note: in the second-argument case, the content of sp will be overwritten!


% Version:  v0.9d
% Build:    14071115
% Date:     Jul-11 2014, 3:32 PM EST
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* counters */
    mwIndex c, r, nv, oc, or, vc, tc;
    mwSize m, n;
    double iv;

    /* pointers */
    mwIndex *io, *jo;
    const double *v, *i, *j;
    double *vo;

    /* input check */
    if ((nrhs != 1) &&
        (nrhs != 2) &&
        (nrhs != 5))
        mexErrMsgTxt("Requires 1, 2, or 5 input arguments.");
    if ((mxGetClassID(*prhs) != mxDOUBLE_CLASS) ||
         mxIsSparse(*prhs))
        mexErrMsgTxt("First input argument must be of type full double.");
    nv = mxGetNumberOfElements(*prhs);
    if ((nrhs == 2) &&
        ((!mxIsSparse(prhs[1])) ||
         (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ||
         (mxGetNzmax(prhs[1]) < nv)))
        mexErrMsgTxt("For two-argument call, second input argument must be of type sparse double with at least V values.");
    else if ((nrhs == 5) &&
        ((mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ||
         (mxIsSparse(prhs[1])) ||
         (mxGetNumberOfElements(prhs[1]) != 1) ||
         (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ||
         (mxIsSparse(prhs[2])) ||
         (mxGetNumberOfElements(prhs[2]) != 1) ||
         (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) ||
         (mxIsSparse(prhs[3])) ||
         (mxGetNumberOfElements(prhs[3]) != nv) ||
         (mxGetClassID(prhs[4]) != mxDOUBLE_CLASS) ||
         (mxIsSparse(prhs[4])) ||
         (mxGetNumberOfElements(prhs[4]) != nv)))
        mexErrMsgTxt("For five-argument call, arguments must be Vx1 v, 1x1 m, 1x1 n, Vx1 i, and Vx1 j.");

    /* get data */
    v = (const double *) mxGetData(*prhs);

    /* one argument only */
    if (nrhs == 1) {

        /* create output */
        *plhs = mxCreateSparse(nv, nv, nv, mxREAL);
        if (*plhs == NULL)
            mexErrMsgTxt("Error creating sparse double matrix.");
        vo = (double *) mxGetData(*plhs);
        if (vo == NULL)
            mexErrMsgTxt("Error retrieving double data pointer.");
        io = (mwIndex *) mxGetIr(*plhs);
        if (io == NULL)
            mexErrMsgTxt("Error retrieving Ir index pointer.");
        jo = (mwIndex *) mxGetJc(*plhs);
        if (jo == NULL)
            mexErrMsgTxt("Error retrieving Jc index pointer.");

        /* loop */
        for (vc = 0; vc < nv; ++vc) {
            *vo++ = *v++;
            *io++ = vc;
            *jo++ = vc;
        }
        *jo = vc;

    /* two arguments, overwrite existing sparse */
    } else if (nrhs == 2) {

        /* copy pointer */
        *plhs = (mxArray *) ((void *) prhs[1]);

        /* get data pointer */
        vo = (double *) mxGetData(*plhs);

        /* loop */
        for (vc = nv; vc > 0; --vc)
            *vo++ = *v++;

    /* five arguments, create new sparse from complex data */
    } else {

        /* get size */
        i = (const double *) mxGetData(prhs[1]);
        m = (int) *i;
        i = (const double *) mxGetData(prhs[2]);
        n = (int) *i;
        i = (const double *) mxGetData(prhs[3]);
        j = (const double *) mxGetData(prhs[4]);

        /* create sparse matrix */
        *plhs = mxCreateSparse(m, n, nv, mxREAL);
        if (*plhs == NULL)
            mexErrMsgTxt("Error creating sparse double matrix.");
        vo = (double *) mxGetData(*plhs);
        if (vo == NULL)
            mexErrMsgTxt("Error retrieving double data pointer.");
        io = (mwIndex *) mxGetIr(*plhs);
        if (io == NULL)
            mexErrMsgTxt("Error retrieving Ir index pointer.");
        jo = (mwIndex *) mxGetJc(*plhs);
        if (jo == NULL)
            mexErrMsgTxt("Error retrieving Jc index pointer.");

        /* loop */
        for (vc = nv, or = 0, oc = 0, tc = 0; vc > 0; --vc) {

            /* get target index */
            r = (mwIndex) *i++;
            c = (mwIndex) *j++;
            iv = *v++;

            /* only allow if valid */
            if ((c < oc) ||
                ((c == oc) && (r <= or)) ||
                (r < 1) ||
                (r > m) ||
                (c < 1) ||
                (c > n))
                continue;

            /* update column ? */
            while (c > oc) {
                *jo++ = tc;
                ++oc;
            }

            /* write into arrays */
            *io++ = (r - 1);
            *vo++ = iv;

            /* keep track of written values */
            or = r;
            ++tc;
        }

        /* extend column */
        while (oc <= n) {
            *jo++ = tc;
            ++oc;
        }
    }
}

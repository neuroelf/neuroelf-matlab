/*

   psetdists.c -- return the pair-wise distances between points of two sets

% Version:  v0.9b
% Build:    11051113
% Date:     May-11 2011, 12:41 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2011, Jochen Weber
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int ndim, c1, c2, cc, nr, nc, nd;
    const int *idim = NULL;
    const double *p1 = NULL, *p2 = NULL;
    double *out = NULL, cdiff, tdiff;

    /* argument check */
    if ((nrhs < 1) ||
        (nrhs > 2) ||
        (nlhs > 1))
		mexErrMsgTxt("Bad number of input/output arguments.");

    /* check dims and types */
    ndim = mxGetNumberOfDimensions(*prhs);
    if ((ndim > 2) ||
        (mxGetClassID(*prhs) != mxDOUBLE_CLASS))
        mexErrMsgTxt("Invalid first input argument.");
    idim = mxGetDimensions(*prhs);
    nr = *idim++;
    nd = *idim;
    p1 = (const double *) mxGetData(*prhs);
    if (nrhs > 1) {
        ndim = mxGetNumberOfDimensions(prhs[1]);
        idim = mxGetDimensions(prhs[1]);
        if ((ndim > 2) ||
            (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ||
            (idim[1] != nd))
            mexErrMsgTxt("Invalid second input argument.");
        nc = *idim;
        p2 = (const double *) mxGetData(prhs[1]);
    } else {
        nc = nr;
        p2 = p1;
    }

    /* create output */
    *plhs = mxCreateDoubleMatrix(nr, nc, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output matrix.");
    out = (double *) mxGetData(*plhs);
    if (out == NULL)
        mexErrMsgTxt("Error getting output matrix data pointer.");

    /* loops */
    for (c2 = 0; c2 < nc; ++c2) {
        for (c1 = 0; c1 < nr; ++c1) {
            for (cc = 0, tdiff = 0.0; cc < nd; ++cc) {
                cdiff = p1[c1 + cc * nr] - p2[c2 + cc * nc];
                tdiff += cdiff * cdiff;
            }
            *out++ = sqrt(tdiff);
        }
    }
}

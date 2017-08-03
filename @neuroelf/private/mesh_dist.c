/*

mesh_dist  - computing distances for pairs of vertices

FORMAT:       d = mesh_dist(c, v1, v2)

Input fields:

      c           Cx3 coordinates
      v1, v2      two Vx1 (or one Vx2) vertex indices

Output fields:

      d           Vx1 vertex pair distances (lengths)


% Version:  v0.9d
% Build:    14071115
% Date:     Jul-11 2014, 3:29 PM EST
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
#include "math.h"
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* inputs */
    const int *idim;
    const double *c1, *c2, *c3, *v1, *v2;

    /* outputs */
    double *od;

    /* temp variables */
    double l1, l2, l3;

    /* Loop variable */
    mwIndex i, iv1, iv2, nc, nl = 0;

    /* Check for proper number of arguments. */
    if ((nrhs < 2) ||
        (nrhs > 3))
        mexErrMsgTxt("2 or 3 inputs are required.");
    idim = (const int *) mxGetDimensions(*prhs);
    if (nlhs > 1)
        mexErrMsgTxt("Up to 1 output supported.");
    if (mxGetClassID(*prhs) != mxDOUBLE_CLASS)
        mexErrMsgTxt("First argument must be of type double.");
    if (mxGetNumberOfElements(*prhs) < 3)
        mexErrMsgTxt("First argument must be Cx3.");
    if (idim[1] != 3)
        mexErrMsgTxt("First argument must be Cx3.");
    nc = *idim;
    if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgTxt("Second argument must be of type double.");
    if (nrhs > 2) {
        idim = (const int *) mxGetDimensions(prhs[2]);
        nl = *idim;
        if ((idim[1] > 1) ||
            (mxGetNumberOfElements(prhs[2]) != nl))
            mexErrMsgTxt("Third argument must be Vx1.");
        v2 = (const double *) mxGetData(prhs[2]);
    }
    idim = (const int*) mxGetDimensions(prhs[1]);
    v1 = (const double *) mxGetData(prhs[1]);
    if (nl == 0) {
        nl = *idim;
        if (mxGetNumberOfElements(prhs[1]) != (2 * nl))
            mexErrMsgTxt("Second argument must be Vx2.");
        v2 = &v1[nl];
    } else if (nl != *idim)
        mexErrMsgTxt("Second and third argument must match in first dim.");
    c1 = (const double *) mxGetData(*prhs);
    c2 = &c1[nc];
    c3 = &c2[nc];

    /* generate output */
    *plhs = mxCreateDoubleMatrix(nl, 1, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output array.");
    od = (double *) mxGetData(*plhs);
    if (od == NULL)
        mexErrMsgTxt("Error getting output data pointer.");

    /* loop and do the work */
    for (i = nl; i > 0; --i) {

        /* get indices */
        iv1 = ((mwIndex) *v1++) - 1;
        iv2 = ((mwIndex) *v2++) - 1;

        /* get lengths */
        l1 = c1[iv1] - c1[iv2];
        l2 = c2[iv1] - c2[iv2];
        l3 = c3[iv1] - c3[iv2];

        /* compute length */
        *od++ = sqrt(l1 * l1 + l2 * l2 + l3 * l3);
    }
}

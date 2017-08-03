/*

xffsrfwriteneighborsc.c  - create stream of neighbors for writing

FORMAT:       neistream = xffsrfwriteneighborsc(neighbors)

Input fields:

      neighbors   Nx2 neighbors cell array

Output fields:

      neistream   1xS uint32 stream with neighbor information

% Version:  v0.9d
% Build:    14072111
% Date:     Jul-21 2014, 11:32 AM EST
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int nv = 0, nv2 = 0, ne = 0, nn, vc;
    int *nnei, *nneip, *out;
    const int *ind;
    const double *din;
    const mxArray *cell;

    /* argument check */
    if (nrhs != 1 || nlhs != 1)
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (!mxIsCell(*prhs) ||
        (mxGetNumberOfDimensions(*prhs) != 2))
		mexErrMsgTxt("Input argument must be of type MxN cell.");
    ind = mxGetDimensions(*prhs);
    if (ind[1] != 2)
        mexErrMsgTxt("Cell array must be Vx2 sized.");
    nv = *ind;

    /* allocate temporary memory for number of neighbors per vertex */
    nneip = nnei = (int *) mxCalloc(nv, sizeof(int));
    if (nneip == NULL)
        mexErrMsgTxt("Error allocating temporary array.");

    /* parse over second column of cell array */
    for (nv2 = nv + nv, vc = nv; vc < nv2; ++vc) {
        cell = (const mxArray*) mxGetCell(*prhs, vc);
        if (cell == NULL) {
            mxFree(nnei);
            mexErrMsgTxt("Invalid neighbors cell array.");
        }
        if (mxGetClassID(cell) != mxDOUBLE_CLASS) {
            mxFree(nnei);
            mexErrMsgTxt("Invalid neighbors cell array.");
        }
        nn = (int) mxGetNumberOfElements(cell);
        if (nn > 2097151) {
            mxFree(nnei);
            mexErrMsgTxt("Invalid neighbors cell array.");
        }
        *nneip++ = nn;
        ne += nn + 1;
    }

    /* prepare output */
    *plhs = mxCreateNumericMatrix(1, ne, mxUINT32_CLASS, mxREAL);
    if (*plhs == NULL) {
        mxFree(nnei);
        mexErrMsgTxt("Error creating neighbors stream.");
    }
    out = (int *) mxGetPr(*plhs);
    if (out == NULL) {
        mxFree(nnei);
        mexErrMsgTxt("Error allocating neighbors stream memory.");
    }

    /* parse over second column of cell array again */
    nneip = nnei;
    for (vc = nv; vc < nv2; ++vc) {

        /* get number of neighbors for vertex */
        nn = *nneip++;

        /* put into output first */
        *out++ = nn;
        
        /* nothing else to do */
        if (nn == 0)
            continue;

        /* get pointer to cell content */
        din = mxGetPr((const mxArray*) mxGetCell(*prhs, vc));

        /* copy numbers */
        for (; nn > 0; --nn)
            *out++ = (int) (*din++ - 1);
    }

    /* free memory */
    mxFree(nnei);
}

/*

xffsrfparseneighborsc.c  - parse neighbors more quickly

FORMAT:       [neighbors, skip] = xffsrfparseneighborsc(neid, numvtx)

Input fields:

      neid        neighbor data
      numvtx      number of vertices

Output fields:

      neighbors   Nx2 cell array with neighbors
      skip        number of elements to skip in input

% Version:  v0.9d
% Build:    14072111
% Date:     Jul-21 2014, 11:40 AM EST
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

    int nv = 0, ne = 0, skip = 0, nn, vn = 0, vc;
    const int *ind;
    mxArray *outa;
    double *outd;
    const double *dind;

    /* argument check */
    if (nrhs != 2 || nlhs != 2)
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (!mxIsUint32(*prhs) ||
        (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS))
		mexErrMsgTxt("Input arguments must be of type uint32/double.");
    if (mxGetNumberOfElements(prhs[1]) != 1)
		mexErrMsgTxt("Number of vertices must be scalar.");
    dind = (const double*) mxGetPr(prhs[1]);
    if (mxIsInf(*dind) ||
        mxIsNaN(*dind) ||
        (*dind < 3.0))
		mexErrMsgTxt("Invalid number of vertices.");
    skip = nv = (int) *dind;
    ne = mxGetNumberOfElements(*prhs);
    if (ne < nv)
		mexErrMsgTxt("Neighbors stream too short.");
    ind = (const int*) mxGetPr(*prhs);

    /* prepare output */
    *plhs = mxCreateCellMatrix(nv, 2);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating Neighbors cell array.");
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    if (plhs[1] == NULL)
        mexErrMsgTxt("Error creating skip indicator.");
    if (((double*) mxGetPr(plhs[1])) == NULL)
        mexErrMsgTxt("Error allocating skip indicator memory.");

    /* loop over vertices */
    for (vc = nv; vc > 0; --vc) {

        /* get number of neighbors for vertex */
        nn = *ind++;
        skip += nn;

        /* check skip vs ne ! */
        if (skip >= ne)
            mexErrMsgTxt("Too little data in stream.");

        /* set number of neighbors value */
        outa = mxCreateDoubleMatrix(1, 1, mxREAL);
        if (outa == NULL)
            mexErrMsgTxt("Error creating neighbor list.");
        outd = (double *) mxGetPr(outa);
        if (outd == NULL)
            mexErrMsgTxt("Error allocating memory for neighbor list.");
        *outd = (double) nn;
        mxSetCell(*plhs, vn, outa);
        
        /* nothing else to do */
        if (nn == 0) {
            ++vn;
            continue;
        }

        /* create double matrix */
        outa = mxCreateDoubleMatrix(1, nn, mxREAL);
        if (outa == NULL)
            mexErrMsgTxt("Error creating neighbor list.");
        outd = (double *) mxGetPr(outa);
        if (outd == NULL)
            mexErrMsgTxt("Error allocating memory for neighbor list.");

        /* copy numbers */
        for (; nn > 0; --nn)
            *outd++ = (double) (*ind++ + 1);

        /* set into cell */
        mxSetCell(*plhs, nv + vn++, outa);
    }

    /* set skip number */
    outd = (double *) mxGetPr(plhs[1]);
    *outd = (double) skip;

}

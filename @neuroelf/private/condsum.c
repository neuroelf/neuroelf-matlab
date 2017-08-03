/*

condsum  - computing sum conditionally

FORMAT:       cs = condsum(n, c, v)

Input fields:

      n           1x1 target size (T)
      c           Cx1 condition (target index) value (double)
      v           Cx1 values

Output fields:

      cs          Tx1 conditional sum


% Version:  v0.9d
% Build:    14071115
% Date:     Jul-11 2014, 3:28 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "mex.h"
#include "isinfnan.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* inputs */
    const double *c, *v;

    /* outputs */
    double *od, ii, iv;

    /* loop variable */
    mwIndex i, nv, ti;

    /* isinfnan stuff */
    VARS_FOR_ISINFNAN;
    INIT_INF_NAN_BAD_VAL();

    /* Check for proper number of arguments. */
    if (nrhs != 3)
        mexErrMsgTxt("3 inputs are required.");
    if (nlhs > 1)
        mexErrMsgTxt("Up to 1 output supported.");
    if (mxGetClassID(*prhs) != mxDOUBLE_CLASS)
        mexErrMsgTxt("First argument must be of type double.");
    if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgTxt("Second argument must be of type double.");
    if (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS)
        mexErrMsgTxt("Third argument must be of type double.");
    if (mxGetNumberOfElements(*prhs) < 1)
        mexErrMsgTxt("First argument must not be empty.");
    if (mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]))
        mexErrMsgTxt("Second and third argument must have equal number of elements.");

    /* output size and maximal index */
    c = (const double *) mxGetData(*prhs);
    nv = (mwIndex) *c;

    /* input arguments */
    c = (const double *) mxGetData(prhs[1]);
    v = (const double *) mxGetData(prhs[2]);

    /* generate output */
    *plhs = mxCreateDoubleMatrix(nv, 1, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output array.");
    od = (double *) mxGetData(*plhs);
    if (od == NULL)
        mexErrMsgTxt("Error getting output data pointer.");

    /* loop and do the work */
    for (i = mxGetNumberOfElements(prhs[1]); i > 0; --i) {

        /* get values */
        ii = (*c++);
        iv = *v++;

        /* no invalid values */
        IF_IS_BAD_VAL(ii) {
            continue;
        }
        IF_IS_BAD_VAL(iv) {
            continue;
        }
        ti = ((mwIndex) ii) - 1;
        if (ti >= nv)
            continue;

        /* sum up */
        od[ti] += iv;
    }
}

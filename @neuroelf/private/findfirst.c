/*

find first element != 0 in input

FORMAT:         ff = findfirst(v [, spos]);

Input fields:

      v             vector
      spos          start position for search (default: 1),
                    if negative do search from end to begin

Output fields:

      ff            first occurrance of "true" in v (as of spos)

% Version:  v0.9d
% Build:    14062016
% Date:     Jun-20 2014, 4:20 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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

const static char *ff_errbadindex = "Invalid start index given.";
const static int ff_emptyarray[2] = {0, 0};
const static int ff_scalararray[2] = {1, 1};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxClassID cid;
	const int *dim;
	const unsigned char *uc;
	const unsigned short *us;
	const unsigned int *ul;
	const double *dbl;
	bool forward = 1;
	int nd, ne, c = 0, cm, si;
	int odim[2] = {1, 2};
	double *odbl;

	if (nrhs < 1 || nrhs > 2 || nlhs > 2)
		mexErrMsgTxt("Bad number of input/output arguments.");
	cid = mxGetClassID(*prhs);
	dim = mxGetDimensions(*prhs);
	nd = mxGetNumberOfDimensions(*prhs);
	ne = mxGetNumberOfElements(*prhs);
	si = 0;
	if ((nrhs > 1) &&
	    mxIsDouble(prhs[1]) &&
	    (mxGetNumberOfElements(prhs[1]) == 1)) {
		si = (int) *((double*) mxGetPr(prhs[1]));
		if (si > 0) {
			--si;
			if (si >= ne)
				mexErrMsgTxt(ff_errbadindex);
		} else if (si < 0) {
			forward = 0;
			si = ne + si;
			if (si < 0)
				mexErrMsgTxt(ff_errbadindex);
		} else
			mexErrMsgTxt(ff_errbadindex);
	}

	switch (cid) {
		case mxLOGICAL_CLASS:
		case mxINT8_CLASS:
		case mxUINT8_CLASS:
			uc = (const unsigned char*) mxGetData(*prhs);
			uc = &uc[si];
			if (forward) {
				for (c = si; c < ne; ++c)
					if (*uc++ != 0)
						break;
			} else {
				for (c = si; c >= 0; --c)
					if (*uc-- != 0)
						break;
			}
			break;

		case mxINT16_CLASS:
		case mxUINT16_CLASS:
			us = (const unsigned short*) mxGetData(*prhs);
			us = &us[si];
			if (forward) {
				for (c = si; c < ne; ++c)
					if (*us++ != 0)
						break;
			} else {
				for (c = si; c >= 0; --c)
					if (*us-- != 0)
						break;
			}
			break;

		case mxINT32_CLASS:
		case mxUINT32_CLASS:
		case mxSINGLE_CLASS:
			ul = (const unsigned int*) mxGetData(*prhs);
			ul = &ul[si];
			if (forward) {
				for (c = si; c < ne; ++c)
					if (*ul++ != 0)
						break;
			} else {
				for (c = si; c >= 0; --c)
					if (*ul-- != 0)
						break;
			}
			break;

		case mxDOUBLE_CLASS:
			dbl = (const double*) mxGetData(*prhs);
			dbl = &dbl[si];
			if (forward) {
				for (c = si; c < ne; ++c)
					if (*dbl++ != 0.0)
						break;
			} else {
				for (c = si; c >= 0; --c)
					if (*dbl-- != 0.0)
						break;
			}
			break;

		default: mexErrMsgTxt("Invalid input class.");
			 break;
	}

	if (c < 0 || c >= ne) {

		*plhs = mxCreateNumericArray(2, ff_emptyarray, mxDOUBLE_CLASS, mxREAL);
		if (nlhs > 1)
			plhs[1] = mxCreateNumericArray(2, ff_emptyarray, mxDOUBLE_CLASS, mxREAL);

	} else {

		*plhs = mxCreateNumericArray(2, ff_scalararray, mxDOUBLE_CLASS, mxREAL);
		if (*plhs == NULL)
			mexErrMsgTxt("Error allocating 1x1 output argument.");
		odbl = (double *) mxGetData(*plhs);
		if (odbl == NULL)
			mexErrMsgTxt("Error addressing 1x1 output argument.");
		*odbl = (double) (c + 1);

		if (nlhs > 1) {
			odim[1] = nd;
			plhs[1] = mxCreateNumericArray(2, odim, mxDOUBLE_CLASS, mxREAL);
			if (plhs[1] == NULL)
				mexErrMsgTxt("Error allocating 1xN output argument.");
			odbl = (double *) mxGetData(plhs[1]);
			if (odbl == NULL)
				mexErrMsgTxt("Error addressing 1xN output argument.");
			for ( ; nd > 0; --nd) {
				cm = c % *dim;
				*odbl++ = (double) (cm + 1);
				c = (int) (((double) c) / ((double) *dim++));
			}
		}
	}
}

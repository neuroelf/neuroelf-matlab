/*

gluetostringc see gluetostring.m

FORMAT:       glued = gluetostring(toglue, glue, endterm);

Input fields:

      toglue      cell array with snippets
      glue        optional glue string (default: LF)
      endterm     flag, keep end termination

Output fields:
      glued       1xN char array/string with glued content

% Version:  v0.9a
% Build:    11050511
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2010, Jochen Weber
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

#define CHAR_LF 10

/* here comes the main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int ne, ncells;
    const int *dim;

    mxChar *glued, gluechar;
    const mxChar *toglue;
    const mxChar *glue;
    const mxChar defglue[1] = {CHAR_LF};
    const mxArray *snippet, *togluec;

    int odim[2] = {1, 0};
    signed long gluelen = 1, copypos = 0, totallen = 0, sc = 0;
    bool endterm = 0;

    /* check argument count and content */
	if (nrhs < 1 || nrhs > 3 || nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");
    if (mxGetClassID(*prhs) != mxCELL_CLASS)
		mexErrMsgTxt("First input must be a cell array.");
    togluec = *prhs;
    ncells = mxGetNumberOfElements(togluec);
    if (ncells == 0) {
        *odim = 0;
        *plhs = mxCreateCharArray(2, odim);
        return;
    }
    glue = defglue;
    if ((nrhs > 1) &&
        (mxGetClassID(prhs[1]) == mxCHAR_CLASS)) {
        dim = mxGetDimensions(prhs[1]);
        if ((dim[1] > 0) &&
            (dim[1] == mxGetNumberOfElements(prhs[1]))) {
            glue = (const mxChar *) mxGetData(prhs[1]);
            gluelen = dim[1];
        }
    }
    if (nrhs > 2) {
        if ((mxGetClassID(prhs[2]) == mxLOGICAL_CLASS) &&
            (mxGetNumberOfElements(prhs[2]) == 1))
            endterm = (*((unsigned char *) mxGetData(prhs[2])) != 0);
    }

    /* count real snippets and length */
    for (; sc < ncells; ++sc) {
        snippet = mxGetCell(togluec, sc);
        if (mxGetClassID(snippet) != mxCHAR_CLASS)
            mexErrMsgTxt("All cells must be of type char.");
        totallen += mxGetNumberOfElements(snippet);
    }
    totallen += gluelen * ((endterm) ? (ncells) : (ncells - 1));

    /* create output */
    odim[1] = totallen;
    *plhs = mxCreateCharArray(2, odim);
    glued = (mxChar *) mxGetData(*plhs);

    /* depending on gluelen */
    gluechar = *glue;
    --ncells;
    switch (gluelen) {
        case 0:
            for (sc = 0; sc < ncells; ++sc) {
                snippet = mxGetCell(togluec, sc);
                toglue = (const mxChar*) mxGetData(snippet);
                ne = mxGetNumberOfElements(snippet);
                for (copypos = ne; copypos > 0; --copypos)
                    *glued++ = *toglue++;
            }
            break;
        case 1:
            for (sc = 0; sc < ncells; ++sc) {
                snippet = mxGetCell(togluec, sc);
                toglue = (const mxChar*) mxGetData(snippet);
                ne = mxGetNumberOfElements(snippet);
                for (copypos = ne; copypos > 0; --copypos)
                    *glued++ = *toglue++;
                *glued++ = gluechar;
            }
            break;
        case 2:
            for (sc = 0; sc < ncells; ++sc) {
                snippet = mxGetCell(togluec, sc);
                toglue = (const mxChar*) mxGetData(snippet);
                ne = mxGetNumberOfElements(snippet);
                for (copypos = ne; copypos > 0; --copypos)
                    *glued++ = *toglue++;
                *glued++ = gluechar;
                *glued++ = glue[1];
            }
            break;
        default:
            for (sc = 0; sc < ncells; ++sc) {
                snippet = mxGetCell(togluec, sc);
                toglue = (const mxChar*) mxGetData(snippet);
                ne = mxGetNumberOfElements(snippet);
                for (copypos = ne; copypos > 0; --copypos)
                    *glued++ = *toglue++;
                *glued++ = gluechar;
                for (copypos = 1; copypos < gluelen; ++copypos)
                    *glued++ = glue[copypos];
            }
            break;
    }
    snippet = mxGetCell(togluec, sc);
    ne = mxGetNumberOfElements(snippet);
    if (ne > 0) {
        toglue = (const mxChar*) mxGetData(snippet);
        for (copypos = ne; copypos > 0; --copypos)
            *glued++ = *toglue++;
    }
    if (endterm)
        for (copypos = 0; copypos < gluelen; ++copypos)
            *glued++ = *glue++;
}

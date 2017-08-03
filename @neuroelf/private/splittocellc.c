/*

splittocellc see splittocell.m

FORMAT:       split = splittocellc(tosplit, delim, multi, anyd, quotes);

Input fields:

      tosplit     character line with content to split
      delim       char array with delimiter(s)
      multi       boolean flag, contract multiple delimiters
      anyd        boolean flag, delim is a set of chars, each matching
      quotes      boolean flag, guess for quotes

Output fields:

      s           1xN cell array with split content
      ns          number of splits

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

#define CHAR_DQUOTE 34
#define CHAR_SQUOTE 39
#define CHAR_BSLASH 92

/* here comes the main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int *dim;
    const mxChar *tosplit, *delim;
    mxChar *copy;
    mxChar delimval = 0, cval = 0, quotec = 63;
    const mxChar defdelim[1] = {9};
    bool multi = 0, anyd = 0, quotes = 0;
    signed long tosplitlen = 0, delimlen = 1, delimlenm = 0, lastfrom = 0;

    int odim[2] = {1, 1};
    int odim0[2] = {0, 0};
    signed long *frompos, *topos;
    signed long cpos = 0, dpos = 0, spos = 0;
    unsigned long nrofsplits = 0, realnrofsplits = 0;
    bool inquote = 0, wbslash = 0;
    mxArray *split, *psplit;
    double *nsplit;



    /* check argument count and content */
	if (nrhs < 1 || nrhs > 5 || nlhs > 2)
		mexErrMsgTxt("Bad number of input/output arguments.");
    dim = mxGetDimensions(*prhs);
    if ((mxGetClassID(*prhs) != mxCHAR_CLASS) ||
        (dim[1] != mxGetNumberOfElements(*prhs)))
        mexErrMsgTxt("First argument must be a 1xN char array.");
    tosplit = (mxChar *) mxGetData(*prhs);
    tosplitlen = dim[1];
    if (tosplitlen == 0) {
        *plhs = mxCreateCellMatrix(1, 0);
        if (nlhs > 1)
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        return;
    }
    delim = defdelim;
    if ((nrhs > 1) &&
        (mxGetClassID(prhs[1]) == mxCHAR_CLASS)) {
        dim = mxGetDimensions(prhs[1]);
        if ((dim[1] > 0) &&
            (dim[1] == mxGetNumberOfElements(prhs[1]))) {
            delim = (mxChar *) mxGetData(prhs[1]);
            delimlen = dim[1];
        }
    }
    if (nrhs > 2) {
        if ((mxGetClassID(prhs[2]) == mxLOGICAL_CLASS) &&
            (mxGetNumberOfElements(prhs[2]) == 1))
            multi = (*((unsigned char *) mxGetData(prhs[2])) != 0);
    }
    if (nrhs > 3) {
        if ((mxGetClassID(prhs[3]) == mxLOGICAL_CLASS) &&
            (mxGetNumberOfElements(prhs[3]) == 1))
            anyd = (*((unsigned char *) mxGetData(prhs[3])) != 0);
    }
    if (nrhs > 4) {
        if ((mxGetClassID(prhs[4]) == mxLOGICAL_CLASS) &&
            (mxGetNumberOfElements(prhs[4]) == 1))
            quotes = (*((unsigned char *) mxGetData(prhs[4])) != 0);
    }

    /* create necessary temporary storage */
    frompos = (signed long *) mxCalloc(tosplitlen, sizeof(signed long));
    if (frompos == NULL)
        mexErrMsgTxt("Error allocating temporary array.");
    topos = (signed long *) mxCalloc(tosplitlen, sizeof(signed long));
    if (topos == NULL) {
        mxFree(frompos);
        mexErrMsgTxt("Error allocating temporary array.");
    }

    /* perform rudimentary search */
    delimval = *delim;
    delimlenm = delimlen - 1;

    /* if only one character, anyd doesn't matter */
    if (!quotes) {
        if (delimlen == 1) {
            for (; cpos < tosplitlen; ++cpos) {
                if (tosplit[cpos] == delimval) {
                    frompos[nrofsplits] = lastfrom;
                    topos[nrofsplits++] = cpos - 1;
                    lastfrom = cpos + 1;
                }
            }
            --cpos;
            if (tosplit[cpos] != delimval) {
                frompos[nrofsplits] = lastfrom;
                topos[nrofsplits++] = cpos;
            }

        /* if not anyd, compare the whole delim string */
        } else if (!anyd) {
            for (; cpos < tosplitlen; ++cpos) {
                if (tosplit[cpos] == delimval) {
                    for (spos = cpos + delimlenm, dpos = delimlenm; dpos > 0; --dpos)
                        if (tosplit[spos--] != delim[dpos])
                            break;
                    if (dpos == 0) {
                        frompos[nrofsplits] = lastfrom;
                        topos[nrofsplits++] = cpos - 1;
                        lastfrom = cpos + delimlen;
                    }
                }
            }
            --cpos;
            if (tosplit[cpos] != delim[delimlenm]) {
                frompos[nrofsplits] = lastfrom;
                topos[nrofsplits++] = cpos;
            }

        /* otherwise each matches on its own */
        } else {
            for (; cpos < tosplitlen; ++cpos) {
                cval = tosplit[cpos];
                for (dpos = delimlenm; dpos >= 0; --dpos)
                    if (cval == delim[dpos])
                        break;
                if (dpos >= 0) {
                    frompos[nrofsplits] = lastfrom;
                    topos[nrofsplits++] = cpos - 1;
                    lastfrom = cpos + 1;
                }
            }
            --cpos;
            for (dpos = delimlenm; dpos >= 0; --dpos) {
                if (tosplit[cpos] == delim[dpos])
                    break;
            }
            if (dpos < 0) {
                frompos[nrofsplits] = lastfrom;
                topos[nrofsplits++] = cpos;
            }
        }

    /* heed quotes */
    } else {
        if (delimlen == 1) {
            for (; cpos < tosplitlen; ++cpos) {
                cval = tosplit[cpos];
                switch (cval) {
                    case CHAR_BSLASH:
                        wbslash = !wbslash;
                        continue;
                    case CHAR_DQUOTE:
                    case CHAR_SQUOTE:
                        if (wbslash) {
                            wbslash = 0;
                            continue;
                        }
                        else if (!inquote)
                            quotec = cval;
                        else if (cval != quotec)
                            continue;
                        inquote = !inquote;
                        continue;
                }
                wbslash = 0;
                if (inquote)
                    continue;
                if (cval == delimval) {
                    frompos[nrofsplits] = lastfrom;
                    topos[nrofsplits++] = cpos - 1;
                    lastfrom = cpos + 1;
                }
            }
            --cpos;
            if (tosplit[cpos] != delimval) {
                frompos[nrofsplits] = lastfrom;
                topos[nrofsplits++] = cpos;
            }

        /* if not anyd, compare the whole delim string */
        } else if (!anyd) {
            for (; cpos < tosplitlen; ++cpos) {
                cval = tosplit[cpos];
                switch (cval) {
                    case CHAR_BSLASH:
                        wbslash = !wbslash;
                        continue;
                    case CHAR_DQUOTE:
                    case CHAR_SQUOTE:
                        if (wbslash) {
                            wbslash = 0;
                            continue;
                        }
                        else if (!inquote)
                            quotec = cval;
                        else if (cval != quotec)
                            continue;
                        inquote = !inquote;
                        continue;
                }
                wbslash = 0;
                if (inquote)
                    continue;
                if (cval == delimval) {
                    for (spos = cpos + delimlenm, dpos = delimlenm; dpos > 0; --dpos)
                        if (tosplit[spos--] != delim[dpos])
                            break;
                    if (dpos == 0) {
                        frompos[nrofsplits] = lastfrom;
                        topos[nrofsplits++] = cpos - 1;
                        lastfrom = cpos + delimlen;
                    }
                }
            }
            --cpos;
            if (tosplit[cpos] != delim[delimlenm]) {
                frompos[nrofsplits] = lastfrom;
                topos[nrofsplits++] = cpos;
            }

        /* otherwise each matches on its own */
        } else {
            for (; cpos < tosplitlen; ++cpos) {
                cval = tosplit[cpos];
                switch (cval) {
                    case CHAR_BSLASH:
                        wbslash = !wbslash;
                        continue;
                    case CHAR_DQUOTE:
                    case CHAR_SQUOTE:
                        if (wbslash) {
                            wbslash = 0;
                            continue;
                        }
                        else if (!inquote)
                            quotec = cval;
                        else if (cval != quotec)
                            continue;
                        inquote = !inquote;
                        continue;
                }
                wbslash = 0;
                if (inquote)
                    continue;
                for (dpos = delimlenm; dpos >= 0; --dpos)
                    if (cval == delim[dpos])
                        break;
                if (dpos >= 0) {
                    frompos[nrofsplits] = lastfrom;
                    topos[nrofsplits++] = cpos - 1;
                    lastfrom = cpos + 1;
                }
            }
            --cpos;
            for (dpos = delimlenm; dpos >= 0; --dpos) {
                if (tosplit[cpos] == delim[dpos])
                    break;
            }
            if (dpos < 0) {
                frompos[nrofsplits] = lastfrom;
                topos[nrofsplits++] = cpos;
            }
        }
    }

    /* keep track of max index */
    realnrofsplits = nrofsplits--;

    /* heed quotes: remove bracketing quotes as well */
    if (quotes)
        for (cpos = nrofsplits; cpos >= 0; --cpos) {
            if (topos[cpos] > (frompos[cpos] + 1)) {
                cval = tosplit[frompos[cpos]];
                if ((cval == tosplit[topos[cpos]]) &&
                    ((cval == CHAR_DQUOTE) ||
                     (cval == CHAR_SQUOTE))) {
                    ++frompos[cpos];
                    --topos[cpos];
                }
            }
        }

    /* any clean up on cells */
    if (multi)
        for (cpos = nrofsplits; cpos >= 0; --cpos)
            if (topos[cpos] < frompos[cpos]) {
                frompos[cpos] = -1;
                --realnrofsplits;
            }

    *plhs = mxCreateCellMatrix(1, realnrofsplits);
    split = *plhs;
    spos = 0;
    for (cpos = 0; cpos <= nrofsplits; ++cpos) {
        delimlenm = frompos[cpos];
        if (delimlenm < 0)
            continue;
        delimlen = topos[cpos] - delimlenm;
        if (++delimlen > 0) {
            delim = &tosplit[delimlenm];
            odim[1] = delimlen;
            psplit = mxCreateCharArray(2, odim);
            copy = (mxChar *) mxGetData(psplit);
            for ( ; delimlen > 0; --delimlen)
                *copy++ = *delim++;
        } else
            psplit = mxCreateCharArray(2, odim0);
        mxSetCell(split, spos++, psplit);
    }

    /* clean up temp arrays */
    mxFree(frompos);
    mxFree(topos);

    /* second output ? */
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        nsplit = (double *) mxGetData(plhs[1]);
        *nsplit = (double) realnrofsplits;
    }
}

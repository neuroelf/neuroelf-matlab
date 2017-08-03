/*

flexible data masking

FORMAT:       dout = flexmask(din [, t , mt , ot, odim, offs])

Input fields:

      din         2D/3D data (logical or numeric)
      t           threshold (masking operand; default: 0)
      mt          type of masking (-1: <=, 0: ==, 1: >=, 2: none; default 1)
      ot          output type (0: as input, 1: boolean; default 0)
      odim        1x3 size, put result in array of size; default as input
      offs        1x3 offset, e.g. [1, 1, 1] skips one index in each dim
                  default: center in new array (floor offset)

Output fields:

      dout        masked data

% Version:  v0.9a
% Build:    11043013
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

/* functions defined in header */
#include "flexmask.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* input and output dimensions */
    int ndim = 3;
    const int *idmc = NULL;
    mxClassID icls;
    const void *dp = NULL;
    unsigned long ne = 0, ns = 0;
    int idim[3], odim[3];
    unsigned long offs[3] = {0, 0, 0};
    int offm[3] = {0, 0, 0};
    mxArray *oa = NULL;
    void *ovp = NULL;
    unsigned char *o8p = NULL;
    double *odp = NULL;

    /* double access */
    const double *da = NULL;

    /* switches */
    signed char masktype = 1;
    signed char outptype = 0;
    bool copyonly = 0;
    bool guessoff = 0;
    bool rslargen = 0;
    bool rsshrink = 0;

    /* threshold value */
    double td   = 0.0;

    /* check input */
    if (nrhs < 1 ||
        (mxGetNumberOfDimensions(*prhs) > 3) ||
        (!mxIsLogical(*prhs) &&
         !mxIsNumeric(*prhs)))
        mexErrMsgTxt("Missing or bad sized/class first argument.");

    /* get dims, etc. */
    ndim = mxGetNumberOfDimensions(*prhs);
    idmc = mxGetDimensions(*prhs);
    icls = mxGetClassID(*prhs);
    dp = (const void*) mxGetData(*prhs);
    idim[0] = *idmc++;
    idim[1] = *idmc++;
    if (ndim > 2)
        idim[2] = *idmc;
    else
        idim[2] = 1;
    odim[0] = idim[0];
    odim[1] = idim[1];
    odim[2] = idim[2];
    ne = idim[0] * idim[1] * idim[2];

    /* threshold must be given as double */
    if (nrhs > 1 &&
        mxIsDouble(prhs[1]) &&
        (mxGetNumberOfElements(prhs[1]) == 1)) {
        da = (const double*) mxGetPr(prhs[1]);
        if (!mxIsInf(*da) &&
            !mxIsNaN(*da))
            td = *da;
    }

    /* type of masking */
    if (nrhs > 2 &&
        mxIsDouble(prhs[2]) &&
        (mxGetNumberOfElements(prhs[2]) == 1)) {
        da = (const double*) mxGetPr(prhs[2]);
        if (mxIsInf(*da) ||
            mxIsNaN(*da)) {
            copyonly = 1;
        } else {
            if ((*da >= -1.0) &&
                (*da <=  2.0))
                masktype = (signed char) *da;
        }
    }

    /* output type */
    if (nrhs > 3 &&
        mxIsDouble(prhs[3]) &&
        (mxGetNumberOfElements(prhs[3]) == 1)) {
        da = (const double*) mxGetPr(prhs[3]);
        if (!mxIsInf(*da) &&
            !mxIsNaN(*da)) {
            if (*da == 1.0)
                outptype = 1;
            else if (*da == 2.0)
                outptype = 2;
        }
    }

    /* output dimensions */
    if (nrhs > 4 &&
        mxIsDouble(prhs[4]) &&
        (mxGetNumberOfElements(prhs[4]) == 3)) {
        da = (const double*) mxGetPr(prhs[4]);
        if (!mxIsInf(da[0]) &&
            !mxIsNaN(da[0]) &&
            !mxIsInf(da[1]) &&
            !mxIsNaN(da[1]) &&
            !mxIsInf(da[2]) &&
            !mxIsNaN(da[2]) &&
            (da[0] >= 1.0) &&
            (da[1] >= 1.0) &&
            (da[2] >= 1.0)) {
            odim[0] = (int) *da++;
            odim[1] = (int) *da++;
            odim[2] = (int) *da;
            if ((odim[0] >= idim[0]) &&
                (odim[1] >= idim[1]) &&
                (odim[2] >= idim[2])) {
                offm[0] = odim[0] - idim[0];
                offm[1] = odim[1] - idim[1];
                offm[2] = odim[2] - idim[2];
                rslargen = 1;
            } else if ((odim[0] <= idim[0]) &&
                (odim[1] <= idim[1]) &&
                (odim[2] <= idim[2])) {
                offm[0] = idim[0] - odim[0];
                offm[1] = idim[1] - odim[1];
                offm[2] = idim[2] - odim[2];
                rsshrink = 1;
            } else {
                mexErrMsgTxt("Array cannot be shrunk/enlarged at once.");
            }
            guessoff = 1;
        }
    }

    /* output offset */
    if (nrhs > 5 &&
        mxIsDouble(prhs[5]) &&
        (mxGetNumberOfElements(prhs[5]) == 3)) {
        da = (const double*) mxGetPr(prhs[5]);
        if (!mxIsInf(da[0]) &&
            !mxIsNaN(da[0]) &&
            !mxIsInf(da[1]) &&
            !mxIsNaN(da[1]) &&
            !mxIsInf(da[2]) &&
            !mxIsNaN(da[2]) &&
            (da[0] >= 0.0) &&
            (da[0] <= (double) offm[0]) &&
            (da[1] >= 0.0) &&
            (da[1] <= (double) offm[1]) &&
            (da[2] >= 0.0) &&
            (da[2] <= (double) offm[2])) {
            guessoff = 0;
            offs[0] = (int) *da++;
            offs[1] = (int) *da++;
            offs[2] = (int) *da;
        }
    }

    /* guess offset */
    if (guessoff) {
        offs[0] = (int) (((double) offm[0]) / 2.0);
        offs[1] = (int) (((double) offm[1]) / 2.0);
        offs[2] = (int) (((double) offm[2]) / 2.0);
    }

    /* threshold depends on class */
    switch (icls) {
        case mxLOGICAL_CLASS:
            if ((td != 0.0) &&
                (td != 1.0))
                mexErrMsgTxt(flexmask_err_badthresh);
            break;
        case mxINT8_CLASS:
            if ((td < -128.0) ||
                (td >  127.5))
                mexErrMsgTxt(flexmask_err_badthresh);
            break;
        case mxINT16_CLASS:
            if ((td < -32768.0) ||
                (td >  32767.5))
                mexErrMsgTxt(flexmask_err_badthresh);
            break;
        case mxINT32_CLASS:
            if ((td < -2147483648.0) ||
                (td >  2147483647.5))
                mexErrMsgTxt(flexmask_err_badthresh);
            break;
        case mxUINT8_CLASS:
            if ((td < 0.0) ||
                (td > 255.5))
                mexErrMsgTxt(flexmask_err_badthresh);
            break;
        case mxUINT16_CLASS:
            if ((td < 0.0) ||
                (td > 65535.5))
                mexErrMsgTxt(flexmask_err_badthresh);
            break;
        case mxUINT32_CLASS:
            if ((td < 0.0) ||
                (td > 4294967295.5))
                mexErrMsgTxt(flexmask_err_badthresh);
            break;
        case mxSINGLE_CLASS:
            if (mxIsInf((double) ((float) td)) ||
                mxIsNaN((double) ((float) td)))
                mexErrMsgTxt(flexmask_err_badthresh);
            break;
        case mxDOUBLE_CLASS:
            break;
        default:
            mexErrMsgTxt("Unsupported input class.");
    }

    /* function used depends on output array type, so switch on that */
    /* first */
    switch (outptype) {

        /* as input */
        case 0:

            /* create array as input type */
            oa = mxCreateNumericArray(3, odim, icls, mxREAL);
            if (oa == NULL)
                mexErrMsgTxt(flexmask_err_memalloc);
            ovp = (void *) mxGetData(oa);
            if (ovp == NULL)
                mexErrMsgTxt(flexmask_err_memhandle);

            /* further depends on size operation */
            /* output larger */
            if (rslargen) {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_COPYOP_PART(copyle_intopart) break;
                    case  0: FLXMSK_SWCASE_COPYOP_PART(copyeq_intopart) break;
                    case  1: FLXMSK_SWCASE_COPYOP_PART(copyge_intopart) break;
                    case  2: FLXMSK_SWCASE_COPY_PART(copy_intopart) break;
                }

            /* output smaller (part) */
            } else if (rsshrink) {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_COPYOP_PART(copyle_frompart) break;
                    case  0: FLXMSK_SWCASE_COPYOP_PART(copyeq_frompart) break;
                    case  1: FLXMSK_SWCASE_COPYOP_PART(copyge_frompart) break;
                    case  2: FLXMSK_SWCASE_COPY_PART(copy_frompart) break;
                }

            /* output as big as input */
            } else {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_COPYOP(copyle) break;
                    case  0: FLXMSK_SWCASE_COPYOP(copyeq) break;
                    case  1: FLXMSK_SWCASE_COPYOP(copyge) break;
                    case  2: FLXMSK_SWCASE_COPY(copy) break;
                }
            }

            /* end of case "same as input" */
            break;

        /* as logical */
        case 1:

            /* create logical array */
            FLXMSK_CREATEARRAY(LOGICAL, odim, oa, unsigned char, o8p)

            /* further depends on size operation */
            /* output larger */
            if (rslargen) {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_SET1OP_PART(set1le_intopart) break;
                    case  0: FLXMSK_SWCASE_SET1OP_PART(set1eq_intopart) break;
                    case  1: FLXMSK_SWCASE_SET1OP_PART(set1ge_intopart) break;
                    default: mxDestroyArray(oa); mexErrMsgTxt(flexmask_err_nologcopy); break;
                }

            /* output smaller (part) */
            } else if (rsshrink) {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_SET1OP_PART(set1le_frompart) break;
                    case  0: FLXMSK_SWCASE_SET1OP_PART(set1eq_frompart) break;
                    case  1: FLXMSK_SWCASE_SET1OP_PART(set1ge_frompart) break;
                    default: mxDestroyArray(oa); mexErrMsgTxt(flexmask_err_nologcopy); break;
                }

            /* output as big as input */
            } else {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_SET1OP(set1le) break;
                    case  0: FLXMSK_SWCASE_SET1OP(set1eq) break;
                    case  1: FLXMSK_SWCASE_SET1OP(set1ge) break;
                    default: mxDestroyArray(oa); mexErrMsgTxt(flexmask_err_nologcopy); break;
                }
            }

            /* end of case "to logical" */
            break;

        /* as double */
        case 2:

            /* create double array */
            FLXMSK_CREATEARRAY(DOUBLE, odim, oa, double, odp)

            /* further depends on size operation */
            /* output larger */
            if (rslargen) {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_COPYOP_TODBL_PART(copyle_todbl_intopart) break;
                    case  0: FLXMSK_SWCASE_COPYOP_TODBL_PART(copyeq_todbl_intopart) break;
                    case  1: FLXMSK_SWCASE_COPYOP_TODBL_PART(copyge_todbl_intopart) break;
                    case  2: FLXMSK_SWCASE_COPY_TODBL_PART(copy_todbl_intopart) break;
                }

            /* output smaller (part) */
            } else if (rsshrink) {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_COPYOP_TODBL_PART(copyle_todbl_frompart) break;
                    case  0: FLXMSK_SWCASE_COPYOP_TODBL_PART(copyeq_todbl_frompart) break;
                    case  1: FLXMSK_SWCASE_COPYOP_TODBL_PART(copyge_todbl_frompart) break;
                    case  2: FLXMSK_SWCASE_COPY_TODBL_PART(copy_todbl_frompart) break;
                }

            /* output as big as input */
            } else {

                /* and depends on masking type */
                switch (masktype) {
                    case -1: FLXMSK_SWCASE_COPYOP_TODBL(copyle_todbl) break;
                    case  0: FLXMSK_SWCASE_COPYOP_TODBL(copyeq_todbl) break;
                    case  1: FLXMSK_SWCASE_COPYOP_TODBL(copyge_todbl) break;
                    case  2: FLXMSK_SWCASE_COPY_TODBL(copy_todbl) break;
                }
            }

            /* end of case "to double" */
            break;
    }
    plhs[0] = oa;

    /* second argument ? */
    if (nlhs > 1) {
        plhs[1] = mxCreateNumericArray(2, flexmask_scalar_dim, mxDOUBLE_CLASS, mxREAL);
        if (plhs[1] == NULL)
            mexErrMsgTxt(flexmask_err_memalloc);

        /* return number of (set) elements */
        *((double *) mxGetPr(plhs[1])) = (double) ns;
    }
}

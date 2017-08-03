/*

extract (linearly interpolate) data from slice

FORMAT:       [islice, aslice] = renderv3dxia(slice, aslice, x, y, perc, c)

Input fields:

      slice       XxYx1 or XxYx3 slice data
      aslice      XxYx1 optional alpha data (otherwise returned blank)
      x, y        coordinates from which to extract data (and + 1)
      perc        1x4 weighting vector
      c           slicing direction (see renderv3d code)

No output fields.

% Version:  v0.9d
% Build:    14071115
% Date:     Jul-11 2014, 3:31 PM EST
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
#include "isinfnan.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* variables */
    int a, c, d3, x1, x2, xc, xs, xy, xyo, y1, y2, yc, ys, ndim;
    int odim[3] = {0, 0, 1};
    const int *idim;
    double p1, p2, p3, p4, v1, v2, v3, v4, o1, o2, o3;
    double *outaslice, *outslice, *outsliceg, *outsliceb, *outarow, *outrow, *outrowg, *outrowb;
    const double *aslice, *slice, *sliceg, *sliceb, *inval, *inval2, *inval3, *inval4;

    /* debug */
    /* char dbgmsg[256]; */

    /* for IS_BAD_VAL */
    VARS_FOR_ISINFNAN
    INIT_INF_NAN_BAD_VAL();

    /* check number, type, fields of in/out arguments */
	if (nrhs != 6)
		mexErrMsgTxt("Requires exactly six input arguments.");
    if (nlhs != 2)
        mexErrMsgTxt("Requires to receive two outputs.");
    if (!mxIsDouble(*prhs) ||
       (mxGetNumberOfElements(*prhs) < 1))
        mexErrMsgTxt("First argument must be a non-empty double.");
    if (!mxIsDouble(prhs[2]) ||
       (mxGetNumberOfElements(prhs[2]) < 1))
        mexErrMsgTxt("Third argument must be a non-empty double.");
    if (!mxIsDouble(prhs[3]) ||
       (mxGetNumberOfElements(prhs[3]) < 1))
        mexErrMsgTxt("Forth argument must be a non-empty double.");
    if (!mxIsDouble(prhs[4]) ||
       (mxGetNumberOfElements(prhs[4]) != 4))
        mexErrMsgTxt("Fifth argument must be a 1x4 double.");
    if (!mxIsDouble(prhs[5]) ||
       (mxGetNumberOfElements(prhs[5]) != 1))
        mexErrMsgTxt("Sixth argument must be a 1x1 double.");

    /* get dims */
    ndim = mxGetNumberOfDimensions(*prhs);
    idim = (const int*) mxGetDimensions(*prhs);
    if (ndim > 3)
        mexErrMsgTxt("Slice data cannot be more than 3D.");
    xs = *idim++;
    ys = *idim++;
    xy = xs * ys;
    if ((ndim > 2) &&
        (*idim != 1) &&
        (*idim != 3))
        mexErrMsgTxt("Slice data must have 1 or 3 planes of data.");

    /* get slice buffer */
    slice = (const double*) mxGetData(*prhs);
    d3 = (((ndim > 2) && (*idim > 1)) ? 1 : 0);
    if (d3 == 1) {
        sliceg = &slice[xy];
        sliceb = &sliceg[xy];
        odim[2] = 3;
    }

    /* alpha slice */
    a = 0;
    if (mxGetNumberOfElements(prhs[1]) > 0) {
        if (!mxIsDouble(prhs[1]))
            mexErrMsgTxt("Alpha slice must be of type double.");
        idim = (const int*) mxGetDimensions(prhs[1]);
        if ((*idim != xs) ||
            (idim[1] != ys))
            mexErrMsgTxt("Alpha slice must match intensity slice in X/Y size.");
        aslice = (const double*) mxGetData(prhs[1]);
        a = 1;
    }

    /* get slicing direction */
    inval = (const double*) mxGetData(prhs[5]);
    if (mxIsInf(*inval) ||
        mxIsNaN(*inval))
        mexErrMsgTxt("Invalid slicing direction value.");
    c = (int) *inval;
    if ((c < 1) ||
        (c > 6))
        mexErrMsgTxt("Invalid slicing direction value (out of range).");

    /* positions */
    xc = mxGetNumberOfElements(prhs[2]);
    if ((c != 2) &&
        (c != 5))
        *odim = xc;
    else
        odim[1] = xc;
    inval = (const double*) mxGetData(prhs[2]);
    xyo = xc;
    --xc;
    if (mxIsInf(*inval) ||
        mxIsNaN(*inval) ||
        mxIsInf(inval[xc]) ||
        mxIsNaN(inval[xc]))
        mexErrMsgTxt("Invalid x-coordinate values.");
    x1 = (int) *inval;
    x2 = (int) inval[xc];
    if ((xc > 0) &&
        (x2 == inval[xc-1]))
        --x2;
    yc = mxGetNumberOfElements(prhs[3]);
    if ((c != 2) &&
        (c != 5))
        odim[1] = yc;
    else
        *odim = yc;
    inval = (const double*) mxGetData(prhs[3]);
    xyo *= yc;
    --yc;
    if (mxIsInf(*inval) ||
        mxIsNaN(*inval) ||
        mxIsInf(inval[yc]) ||
        mxIsNaN(inval[yc]))
        mexErrMsgTxt("Invalid y-coordinate values.");
    y1 = (int) *inval;
    y2 = (int) inval[yc];
    if ((yc > 0) &&
        (y2 == inval[yc-1]))
        --y2;

    /* get sizes */
    if ((x1 < 1) ||
        (x2 < x1) ||
        (y1 < 1) ||
        (y2 < y1) ||
        (x2 >= xs) ||
        (y2 >= ys))
        mexErrMsgTxt("Access x/y into slice out of range.");

    /* get interpolation factors */
    inval = (const double*) mxGetData(prhs[4]);
    p1 = *inval++;
    if (mxIsInf(p1) ||
        mxIsNaN(p1))
        mexErrMsgTxt("Invalid perc(1) value.");
    if ((p1 < 0.0) ||
        (p1 > 1.0))
        mexErrMsgTxt("Invalid perc(1) value (out of range).");
    p2 = *inval++;
    if (mxIsInf(p2) ||
        mxIsNaN(p2))
        mexErrMsgTxt("Invalid perc(2) value.");
    if ((p2 < 0.0) ||
        (p2 > 1.0))
        mexErrMsgTxt("Invalid perc(2) value (out of range).");
    p3 = *inval++;
    if (mxIsInf(p3) ||
        mxIsNaN(p3))
        mexErrMsgTxt("Invalid perc(3) value.");
    if ((p3 < 0.0) ||
        (p3 > 1.0))
        mexErrMsgTxt("Invalid perc(3) value (out of range).");
    p4 = *inval;
    if (mxIsInf(p4) ||
        mxIsNaN(p4))
        mexErrMsgTxt("Invalid perc(4) value.");
    if ((p4 < 0.0) ||
        (p4 > 1.0))
        mexErrMsgTxt("Invalid perc(4) value (out of range).");

    /* create output arrays */
    if (d3 == 0)
        *plhs = (mxArray*) mxCreateNumericArray(2, odim, mxDOUBLE_CLASS, mxREAL);
    else
        *plhs = (mxArray*) mxCreateNumericArray(3, odim, mxDOUBLE_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating intensity output array.");
    outslice = (double*) mxGetData(*plhs);
    if (d3 == 1) {
        outsliceg = &outslice[xyo];
        outsliceb = &outsliceg[xyo];
    }
    if (a == 1) {
            plhs[1] = (mxArray*) mxCreateNumericArray(2, odim, mxDOUBLE_CLASS, mxREAL);
        if (plhs[1] == NULL)
            mexErrMsgTxt("Error creating alpha output array.");
        outaslice = (double*) mxGetData(plhs[1]);
    } else
        plhs[1] = (mxArray*) mxCreateDoubleMatrix(0, 0, mxREAL);

    /* for slicing directions other than 2 and 5 */
    --x1;
    --y1;
    ndim = xs + 1;
    if ((c != 2) &&
        (c != 5)) {

        /* 1D (grayscale) */
        if (d3 == 0) {

            /* no alpha */
            if (a == 0) {

                /* loop over y, then x coordinates */
                for (yc = y1; yc < y2; ++yc) {
                    sliceg = &slice[yc * xs + x1];
                    sliceb = &sliceg[1];
                    inval = &sliceg[xs];
                    inval2 = &inval[1];
                    for (xc = x1; xc < x2; ++xc, ++outslice, ++sliceg, ++sliceb, ++inval, ++inval2) {
                        v1 = *sliceg;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = *inval;
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = *sliceb;
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = *inval2;
                        IF_IS_BAD_VAL(v4)
                            continue;
                        *outslice = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                    }
                }

            /* with alpha */
            } else {

                /* loop over y, then x coordinates */
                for (yc = y1; yc < y2; ++yc) {
                    xc = yc * xs + x1;
                    sliceg = &slice[xc];
                    sliceb = &sliceg[xs];
                    inval = &aslice[xc];
                    inval2 = &inval[xs];
                    for (xc = x1; xc < x2; ++xc, ++outslice, ++outaslice) {
                        v1 = *sliceg++;
                        v2 = *sliceb++;
                        IF_IS_BAD_VAL(v1) {
                            ++inval;
                            ++inval2;
                            continue;
                        }
                        IF_IS_BAD_VAL(v2) {
                            ++inval;
                            ++inval2;
                            continue;
                        }
                        v3 = *sliceg;
                        IF_IS_BAD_VAL(v3) {
                            ++inval;
                            ++inval2;
                            continue;
                        }
                        v4 = *sliceb;
                        IF_IS_BAD_VAL(v4) {
                            ++inval;
                            ++inval2;
                            continue;
                        }
                        *outslice = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval++;
                        v2 = *inval2++;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = *inval;
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = *inval2;
                        IF_IS_BAD_VAL(v4)
                            continue;
                        *outaslice = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                    }
                }
            }

        /* RGB slice */
        } else {

            /* no alpha */
            if (a == 0) {

                /* loop over y, then x coordinates */
                for (yc = y1; yc < y2; ++yc) {
                    xc = yc * xs + x1;
                    inval = &slice[xc];
                    inval2 = &sliceg[xc];
                    inval3 = &sliceb[xc];
                    for (xc = x1; xc < x2; ++xc, ++outslice, ++outsliceg, ++outsliceb, ++inval, ++inval2, ++inval3) {
                        v1 = *inval;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval[xs];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval[1];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o1 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval2;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval2[xs];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval2[1];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval2[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o2 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval3;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval3[xs];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval3[1];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval3[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        *outslice = o1;
                        *outsliceg = o2;
                        *outsliceb = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                    }
                }

            /* with alpha */
            } else {

                /* loop over y, then x coordinates */
                for (yc = y1; yc < y2; ++yc) {
                    xc = yc * xs + x1;
                    inval = &slice[xc];
                    inval2 = &sliceg[xc];
                    inval3 = &sliceb[xc];
                    inval4 = &aslice[xc];
                    for (xc = x1; xc < x2; ++xc, ++outslice, ++outsliceg, ++outsliceb, ++outaslice, ++inval, ++inval2, ++inval3, ++inval4) {
                        v1 = *inval4;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval4[xs];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval4[1];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval4[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o3 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        if (o3 <= 0.0)
                            continue;
                        v1 = *inval;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval[xs];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval[1];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o1 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval2;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval2[xs];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval2[1];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval2[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o2 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval3;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval3[xs];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval3[1];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval3[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        *outslice = o1;
                        *outsliceg = o2;
                        *outsliceb = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        *outaslice = o3;
                    }
                }
            }
        }

    /* for slicing directions 2 or 5 */
    } else {

        /* row increment */
        xyo = *odim;

        /* 1D */
        if (d3 == 0) {

            /* no alpha */
            if (a == 0) {

                /* loop over y, then x coordinates */
                for (yc = y1; yc < y2; ++yc) {
                    sliceg = &slice[yc * xs + x1];
                    sliceb = &sliceg[1];
                    inval = &sliceg[xs];
                    inval2 = &inval[1];
                    outrow = outslice++;
                    for (xc = x1; xc < x2; ++xc, outrow = &outrow[xyo], ++sliceg, ++sliceb, ++inval, ++inval2) {
                        v1 = *sliceg;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = *sliceb;
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = *inval;
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = *inval2;
                        IF_IS_BAD_VAL(v4)
                            continue;
                        *outrow = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                    }
                }

            /* with alpha */
            } else {

                /* loop over y, then x coordinates */
                for (yc = y1; yc < y2; ++yc) {
                    xc = yc * xs + x1;
                    sliceg = &slice[xc];
                    sliceb = &sliceg[xs];
                    inval = &aslice[xc];
                    inval2 = &inval[xs];
                    outrow = outslice++;
                    outarow = outaslice++;
                    for (xc = x1; xc < x2; ++xc, outrow = &outrow[xyo], outarow = &outarow[xyo]) {
                        v1 = *sliceg++;
                        v3 = *sliceb++;
                        IF_IS_BAD_VAL(v1) {
                            ++inval;
                            ++inval2;
                            continue;
                        }
                        IF_IS_BAD_VAL(v3) {
                            ++inval;
                            ++inval2;
                            continue;
                        }
                        v2 = *sliceg;
                        IF_IS_BAD_VAL(v2) {
                            ++inval;
                            ++inval2;
                            continue;
                        }
                        v4 = *sliceb;
                        IF_IS_BAD_VAL(v4) {
                            ++inval;
                            ++inval2;
                            continue;
                        }
                        *outrow = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval++;
                        v3 = *inval2++;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v2 = *inval;
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v4 = *inval2;
                        IF_IS_BAD_VAL(v4)
                            continue;
                        *outarow = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                    }
                }
            }

        /* RGB image slice */
        } else {

            /* no alpha */
            if (a == 0) {

                /* loop over y, then x coordinates */
                for (yc = y1; yc < y2; ++yc) {
                    xc = yc * xs + x1;
                    inval = &slice[xc];
                    inval2 = &sliceg[xc];
                    inval3 = &sliceb[xc];
                    outrow = outslice++;
                    outrowg = outsliceg++;
                    outrowb = outsliceb++;
                    for (xc = x1; xc < x2; ++xc, outrow = &outrow[xyo], outrowg = &outrowg[xyo], outrowb = &outrowb[xyo], ++inval, ++inval2, ++inval3) {
                        v1 = *inval;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval[1];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval[xs];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o1 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval2;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval2[1];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval2[xs];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval2[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o2 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval3;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval3[1];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval3[xs];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval3[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        *outrow = o1;
                        *outrowg = o2;
                        *outrowb = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                    }
                }

            /* with alpha */
            } else {

                /* loop over y, then x coordinates */
                for (yc = y1; yc < y2; ++yc) {
                    xc = yc * xs + x1;
                    inval = &slice[xc];
                    inval2 = &sliceg[xc];
                    inval3 = &sliceb[xc];
                    inval4 = &aslice[xc];
                    outrow = outslice++;
                    outrowg = outsliceg++;
                    outrowb = outsliceb++;
                    outarow = outaslice++;
                    for (xc = x1; xc < x2; ++xc, outrow = &outrow[xyo], outrowg = &outrowg[xyo], outrowb = &outrowb[xyo], outarow = &outarow[xyo], ++inval, ++inval2, ++inval3, ++inval4) {
                        v1 = *inval4;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval4[1];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval4[xs];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval4[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o3 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval[1];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval[xs];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o1 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval2;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval2[1];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval2[xs];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval2[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        o2 = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        v1 = *inval3;
                        IF_IS_BAD_VAL(v1)
                            continue;
                        v2 = inval3[1];
                        IF_IS_BAD_VAL(v2)
                            continue;
                        v3 = inval3[xs];
                        IF_IS_BAD_VAL(v3)
                            continue;
                        v4 = inval3[ndim];
                        IF_IS_BAD_VAL(v4)
                            continue;
                        *outrow = o1;
                        *outrowg = o2;
                        *outrowb = p1 * v1 + p2 * v2 + p3 * v3 + p4 * v4;
                        *outarow = o3;
                    }
                }
            }
        }
    }
}

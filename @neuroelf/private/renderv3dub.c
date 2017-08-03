/*

updating the internal render buffer (Ibuffer)

FORMAT:       renderv3dub(data, vc)

Input fields:

      data        struct with required fields
      vc          lookup for which volume

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
    int i3, isloc, isaloc, fieldnum, nvol, vc, ndim, inumc, numx, numy,
        numxy, bufx, bufy, bufxy, bufxy2, x1, x2, xc, y1, y2, yc;
    const int *idim;
    double imin, imaxmin, aval, avali, numa, numc, val, valg, valb;
    double *Ibuffer;
    const double *atable, *ctable, *ctableg, *ctableb, *inval, *inaloc, *inloc, *inlocg, *inlocb;
    mxArray *outmatrix;
    const mxArray *atable_cell, *inmatrix;

    /* debug
    char dbgmsg[256]; */

    /* for IS_BAD_VAL */
    VARS_FOR_ISINFNAN
    INIT_INF_NAN_BAD_VAL();

    /* check number, type, fields of in/out arguments */
	if (nrhs != 2)
		mexErrMsgTxt("Requires exactly two input arguments.");
    if (!mxIsStruct(*prhs) ||
       (mxGetNumberOfElements(*prhs) != 1))
        mexErrMsgTxt("First argument must be a 1x1 struct.");
    if (!mxIsDouble(prhs[1]) ||
       (mxGetNumberOfElements(prhs[1]) != 1))
        mexErrMsgTxt("Second argument must be a 1x1 double.");

    /* get Ibuffer field */
    fieldnum = mxGetFieldNumber(*prhs, "Ibuffer");
    if (fieldnum < 0)
        mexErrMsgTxt("Field Ibuffer not found.");
    outmatrix = (mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    ndim = mxGetNumberOfDimensions(outmatrix);
    idim = (const int*) mxGetDimensions(outmatrix);
    bufx = *idim++;
    bufy = *idim++;
    i3 = ((*idim == 3) ? 1 : 0);
    bufxy = bufx * bufy;
    bufxy2 = 2 * bufxy;
    Ibuffer = (double*) mxGetData(outmatrix);

    /* get indices into buffer */
    fieldnum = mxGetFieldNumber(*prhs, "px");
    if (fieldnum < 0)
        mexErrMsgTxt("Field px not found.");
    inmatrix = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    if (!mxIsDouble(inmatrix))
        mexErrMsgTxt("Field px must be of type double.");
    ndim = mxGetNumberOfElements(inmatrix);
    if (ndim < 1)
        mexErrMsgTxt("Field px must not be empty.");
    inval = (const double*) mxGetData(inmatrix);
    x1 = (int) *inval;
    x2 = (int) inval[ndim-1];
    fieldnum = mxGetFieldNumber(*prhs, "py");
    if (fieldnum < 0)
        mexErrMsgTxt("Field py not found.");
    inmatrix = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    if (!mxIsDouble(inmatrix))
        mexErrMsgTxt("Field py must be of type double.");
    ndim = mxGetNumberOfElements(inmatrix);
    if (ndim < 1)
        mexErrMsgTxt("Field py must not be empty.");
    inval = (const double*) mxGetData(inmatrix);
    y1 = (int) *inval;
    y2 = (int) inval[ndim-1];

    /* get atable and number of volumes from there */
    fieldnum = mxGetFieldNumber(*prhs, "atable");
    if (fieldnum < 0)
        mexErrMsgTxt("Field atable not found.");
    atable_cell = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    if (!mxIsCell(atable_cell) ||
       (mxGetNumberOfElements(atable_cell) < 1))
        mexErrMsgTxt("Field atable must be of type cell and not empty.");
    nvol = mxGetNumberOfElements(atable_cell);

    /* check volume counter */
    inval = (const double*) mxGetData(prhs[1]);
    val = *inval;
    if (mxIsInf(val) ||
        mxIsNaN(val))
        mexErrMsgTxt("Second argument must not be Inf/NaN.");
    if ((val < 1) ||
        (val > ((double) nvol)))
        mexErrMsgTxt("Second argument out of range.");
    vc = (int) (val - 1.0);

    /* get ctable */
    fieldnum = mxGetFieldNumber(*prhs, "ctable");
    if (fieldnum < 0)
        mexErrMsgTxt("Field ctable not found.");
    inmatrix = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    if (!mxIsCell(inmatrix) ||
        (mxGetNumberOfElements(inmatrix) < (vc + 1)))
        mexErrMsgTxt("Field ctable not of type cell or too few elements.");
    inmatrix = (const mxArray*) mxGetCell(inmatrix, vc);
    if (!mxIsDouble(inmatrix))
        mexErrMsgTxt("Field ctable{vc} element must be of type double.");
    ndim = mxGetNumberOfDimensions(inmatrix);
    idim = (const int*) mxGetDimensions(inmatrix);
    if ((ndim != 2) ||
        (idim[1] != 3))
        mexErrMsgTxt("Field ctable{vc} element must be Cx3 in size.");
    inumc = *idim;
    numc = (double) (inumc - 1);
    ctable = (const double*) mxGetData(inmatrix);
    ctableg = &ctable[inumc];
    ctableb = &ctableg[inumc];

    /* get intensity field and check type and size */
    isloc = 0;
    fieldnum = mxGetFieldNumber(*prhs, "intensity_rgb");
    if (fieldnum < 0)
        mexErrMsgTxt("Field intensity_rgb not found.");
    inmatrix = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    if (mxGetNumberOfElements(inmatrix) == 0) {
        isloc = 1;
        fieldnum = mxGetFieldNumber(*prhs, "intensity_loc");
        if (fieldnum < 0)
            mexErrMsgTxt("Field intensity_loc not found.");
        inmatrix = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    }
    if (!mxIsDouble(inmatrix))
        mexErrMsgTxt("Field intensity_loc or intensity_rgb must be of type double.");
    ndim = mxGetNumberOfDimensions(inmatrix);
    if ((ndim < 2) ||
        (ndim > 3))
        mexErrMsgTxt("Field intensity_loc or intensity_rgb must be 2D or 3D.");
    idim = (const int*) mxGetDimensions(inmatrix);
    numx = *idim++;
    numy = *idim++;
    if ((isloc == 0) &&
        (*idim != 3))
        mexErrMsgTxt("Field intensity_rgb must be XxYx3 in size.");
    numxy = numx * numy;
    inloc = (const double*) mxGetData(inmatrix);
    if (isloc == 0) {
        inlocg = &inloc[numxy];
        inlocb = &inlocg[numxy];
    }
    for (val = 0.0, fieldnum = 0; fieldnum < numxy; ++fieldnum)
        val += inloc[fieldnum];

    /* get alpha loc field, and check if it's empty */
    isaloc = 0;
    fieldnum = mxGetFieldNumber(*prhs, "alpha_loc");
    if (fieldnum >= 0) {
        inmatrix = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
        if (mxIsDouble(inmatrix) &&
           (mxGetNumberOfElements(inmatrix) == numxy)) {
            isaloc = 1;
            inaloc = (const double*) mxGetData(inmatrix);
        }
    }

    /* check atable content */
    if (isaloc == 0) {
        inmatrix = (const mxArray*) mxGetCell(atable_cell, vc);
        if (!mxIsDouble(inmatrix) ||
           (mxGetNumberOfElements(inmatrix) < 32))
            mexErrMsgTxt("Invalid atable content for selected volume.");
        atable = (const double*) mxGetData(inmatrix);
        numa = (double) (mxGetNumberOfElements(inmatrix) - 1);
    }

    /* get other fields, and check contents for size */
    fieldnum = mxGetFieldNumber(*prhs, "imin");
    if (fieldnum < 0)
        mexErrMsgTxt("Field imin not found.");
    inmatrix = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    if (!mxIsDouble(inmatrix) ||
       (mxGetNumberOfElements(inmatrix) != nvol))
        mexErrMsgTxt("Field imin mismatches number of volumes in size.");
    inval = (const double*) mxGetData(inmatrix);
    imin = inval[vc];
    fieldnum = mxGetFieldNumber(*prhs, "imaxmin");
    if (fieldnum < 0)
        mexErrMsgTxt("Field imaxmin not found.");
    inmatrix = (const mxArray*) mxGetFieldByNumber(*prhs, 0, fieldnum);
    if (!mxIsDouble(inmatrix) ||
       (mxGetNumberOfElements(inmatrix) != nvol))
        mexErrMsgTxt("Field imaxmin mismatches number of volumes in size.");
    inval = (const double*) mxGetData(inmatrix);
    imaxmin = 1.0 / inval[vc];


    /* with alpha volume */
    --x1;
    --y1;
    if (isaloc == 1) {

        /* with intensity data (not RGB) */
        if (isloc == 1) {

            /* calculation in two loops */
            for (yc = y1; yc < y2; ++yc) {
                for (xc = x1; xc < x2; ++xc, ++inloc, ++inaloc) {

                    /* get values alpha value */
                    aval = *inaloc;
                    IF_IS_BAD_VAL(aval)
                        continue;

                    /* continue if <= 0 */
                    if (aval <= 0.0)
                        continue;

                    /* replace with 1.0 for > 1.0 */
                    else if (aval > 1.0)
                        aval = 1.0;

                    /* convert to index */
                    val = *inloc;
                    IF_IS_BAD_VAL(val)
                        continue;

                    val = (imaxmin * (val - imin));
                    if (val < 0.0)
                        val = 0.0;
                    else if (val > 1.0)
                        val = 1.0;
                    inumc = (int) (numc * val);

                    /* get color for value */
                    val  = ctable[inumc];
                    valg = ctableg[inumc];
                    valb = ctableb[inumc];

                    /* update */
                    avali = 1.0 - aval;
                    Ibuffer[xc + bufx * yc]          = avali * Ibuffer[xc + bufx * yc]          + aval * val;
                    Ibuffer[xc + bufx * yc + bufxy]  = avali * Ibuffer[xc + bufx * yc + bufxy]  + aval * valg;
                    Ibuffer[xc + bufx * yc + bufxy2] = avali * Ibuffer[xc + bufx * yc + bufxy2] + aval * valb;
                }
            }

        /* RGB image with alpha table */
        } else {

            /* follow similar logic */
            for (yc = y1; yc < y2; ++yc) {
                for (xc = x1; xc < x2; ++xc, ++inloc, ++inlocg, ++inlocb, ++inaloc) {

                    /* get values alpha value */
                    aval = *inaloc;
                    IF_IS_BAD_VAL(aval)
                        continue;

                    /* continue if <= 0 */
                    if (aval <= 0.0)
                        continue;

                    /* replace with 1.0 for > 1.0 */
                    else if (aval > 1.0)
                        aval = 1.0;

                    /* get values directly */
                    val = *inloc;
                    IF_IS_BAD_VAL(val)
                        continue;
                    valg = *inlocg;
                    IF_IS_BAD_VAL(valg)
                        continue;
                    valb = *inlocb;
                    IF_IS_BAD_VAL(valb)
                        continue;

                    /* if values are ok */
                    val = (imaxmin * (val - imin));
                    if (val < 0.0)
                        val = 0.0;
                    else if (val > 1.0)
                        val = 1.0;
                    valg = (imaxmin * (valg - imin));
                    if (valg < 0.0)
                        valg = 0.0;
                    else if (valg > 1.0)
                        valg = 1.0;
                    valb = (imaxmin * (valb - imin));
                    if (valb < 0.0)
                        valb = 0.0;
                    else if (valb > 1.0)
                        valb = 1.0;

                    /* update */
                    avali = 1.0 - aval;
                    Ibuffer[xc + bufx * yc]          = avali * Ibuffer[xc + bufx * yc]          + aval * val;
                    Ibuffer[xc + bufx * yc + bufxy]  = avali * Ibuffer[xc + bufx * yc + bufxy]  + aval * valg;
                    Ibuffer[xc + bufx * yc + bufxy2] = avali * Ibuffer[xc + bufx * yc + bufxy2] + aval * valb;
                }
            }
        }

    /* without alpha volume */
    } else {

        /* with intensity data (not RGB) */
        if (isloc == 1) {

            /* calculation in two loops */
            for (yc = y1; yc < y2; ++yc) {
                for (xc = x1; xc < x2; ++xc, ++inloc) {

                    /* get values alpha value */
                    val = *inloc;
                    IF_IS_BAD_VAL(val)
                        continue;
                    val -= imin;

                    /* continue if <= 0 */
                    if (val < 0.0)
                        continue;

                    /* get alpha value */
                    val *= imaxmin;
                    if (val > 1.0)
                        val = 1.0;
                    aval = atable[(int) (numa * val)];

                    /* replace with 1.0 for > 1.0 */
                    if (aval <= 0.0)
                        continue;

                    /* convert to index */
                    if (val < 0.0)
                        val = 0.0;
                    else if (val > 1.0)
                        val = 1.0;
                    inumc = (int) (numc * val);

                    /* get color for value */
                    val  = ctable[inumc];
                    valg = ctableg[inumc];
                    valb = ctableb[inumc];

                    /* update */
                    avali = 1.0 - aval;
                    Ibuffer[xc + bufx * yc]          = avali * Ibuffer[xc + bufx * yc]          + aval * val;
                    Ibuffer[xc + bufx * yc + bufxy]  = avali * Ibuffer[xc + bufx * yc + bufxy]  + aval * valg;
                    Ibuffer[xc + bufx * yc + bufxy2] = avali * Ibuffer[xc + bufx * yc + bufxy2] + aval * valb;
                }
            }

        /* RGB image without alpha volume */
        } else {

            /* calculation in two loops */
            for (yc = y1; yc < y2; ++yc) {
                for (xc = x1; xc < x2; ++xc, ++inloc, ++inlocg, ++inlocb) {

                    /* get values and compute alpha value */
                    val = *inloc;
                    IF_IS_BAD_VAL(val)
                        continue;
                    valg = *inlocg;
                    IF_IS_BAD_VAL(valg)
                        continue;
                    valb = *inlocb;
                    IF_IS_BAD_VAL(valb)
                        continue;
                    val = imaxmin * (val - imin);
                    valg = imaxmin * (valg - imin);
                    valb = imaxmin * (valb - imin);

                    /* get alpha value */
                    aval = 0.309 * val + 0.461 * valg + 0.230 * valb;

                    /* fix between 0.0 and 1.0 */
                    if (val < 0.0)
                        val = 0.0;
                    else if (val > 1.0)
                        val = 1.0;
                    if (valg < 0.0)
                        valg = 0.0;
                    else if (valg > 1.0)
                        valg = 1.0;
                    if (valb < 0.0)
                        valb = 0.0;
                    else if (valb > 1.0)
                        valb = 1.0;
                    if (aval < 0.0)
                        aval = 0.0;
                    else if (aval > 1.0)
                        aval = 1.0;

                    /* get alpha value */
                    aval = atable[(int) (numa * aval)];

                    /* replace with 1.0 for > 1.0 */
                    if (aval <= 0.0)
                        continue;

                    /* update */
                    avali = 1.0 - aval;
                    Ibuffer[xc + bufx * yc]          = avali * Ibuffer[xc + bufx * yc]          + aval * val;
                    Ibuffer[xc + bufx * yc + bufxy]  = avali * Ibuffer[xc + bufx * yc + bufxy]  + aval * valg;
                    Ibuffer[xc + bufx * yc + bufxy2] = avali * Ibuffer[xc + bufx * yc + bufxy2] + aval * valb;
                }
            }
        }
    }
}

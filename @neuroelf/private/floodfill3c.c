/*

floodfill3c see floodfill3.m

FORMAT:       [v, n] = floodfill3c(v, s, m);

Input fields:

      v           input/output volume
      s           start index (1-based)
      m           method
                  1 face-connectivity (3D)
                  2 edge-connectivity (3D)
                  3 vertex-connectivity (3D)
                  4 face-connectivity (XY-slice)
                  5 face-connectivity (XY-slice, zstart = :)
                  6 edge-connectivity (XY-slice)
                  7 edge-connectivity (XY-slice, zstart = :)

% Version:  v0.9d
% Build:    14070916
% Date:     Jul-09 2014, 4:39 PM EST
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
#include "math.h"
#include <stdio.h>

#define CREATE_SCALAR(SNAME, SVAL) { \
    SNAME = mxCreateNumericArray(2, scalardim, mxDOUBLE_CLASS, mxREAL); \
    if (SNAME == NULL) \
        mexErrMsgTxt("Error allocating memory."); \
    *((double *) mxGetPr(SNAME)) = (double) SVAL; \
}

/* functions used for copying */
void ff_copy_intopart(unsigned long x1, unsigned long y1, unsigned long z1, const unsigned char *sd, const int *tnd, unsigned char *td, const unsigned long *off);
unsigned long ff_set1_frompart_eq(unsigned long x1, unsigned long y1, const unsigned char *sd, const int *tnd, unsigned char *td, const unsigned long *off, unsigned char thr);

/* up to 13 directions (and their reverse) */
static const signed long dlx[13] = { 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1 };
static const signed long dly[13] = { 0, 1, 0, 1,-1, 0, 0, 1, 1, 1, 1,-1,-1 };
static const signed long dlz[13] = { 0, 0, 1, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1 };

/* scalar dims */
static const int scalardim[2] = {1, 1};

/* here comes the main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* general purpose pointers */
    const double *rdp = NULL;

    /* sizes and position */
    const int *dim = NULL;
    int idim[3] = {1, 1, 1};
    int tdim[3] = {1, 1, 1};
    unsigned long toff[3] = {1, 1, 1};
    unsigned long dx = 0, dy = 0, dz = 0, dxy = 0, sx = 0, sy = 0, sz = 0, np = 0;

    /* data pointer */
    const unsigned char *data = NULL;
    unsigned char *copy = NULL;
    signed long *next1 = NULL;
    signed long *next2 = NULL;
    signed long *nextb = NULL;
    unsigned long nextc = 0;
    unsigned long nextcb = 0;
    signed long nextp = 0, nextpi = 0;
    bool donext = 1;

    /* direction increments */
    signed long dlxyz[26] =
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    /* counters */
    unsigned long c = 0;
    signed char dc = 0;
    unsigned long mi = 65536;

    /* max directions */
    unsigned long md = 0;

    /* max next poiters */
    unsigned long mn = 0;

    /* variable output text */
    /* char vstr[256]; */

    /* nr. of arguments check */
    if ((nrhs < 3) ||
        (nrhs > 4) ||
        (nlhs > 2))
        mexErrMsgTxt("Bad number of input/output arguments.");

    /* argument class and size check */
    if (!mxIsLogical(prhs[0]) ||
        (mxGetNumberOfElements(prhs[0]) < 8) ||
        (mxGetNumberOfDimensions(prhs[0]) > 3) ||
        !mxIsDouble(prhs[1]) ||
        (mxGetNumberOfElements(prhs[1]) != 3) ||
        !mxIsDouble(prhs[2]) ||
        (mxGetNumberOfElements(prhs[2]) != 1))
        mexErrMsgTxt("Invalid argument type or dimension.");

    /* get & check argument content */
    dim = mxGetDimensions(prhs[0]);
    idim[0] = dim[0];
    idim[1] = dim[1];
    if (mxGetNumberOfDimensions(prhs[0]) > 2)
        idim[2] = dim[2];
    dx = idim[0];
    dy = idim[1];
    dz = idim[2];
    dxy = dx * dy;

    /* check start position */
    rdp = (const double *) mxGetPr(prhs[1]);
    if (mxIsInf(*rdp) ||
        mxIsNaN(*rdp) ||
        mxIsInf(rdp[1]) ||
        mxIsNaN(rdp[1]) ||
        mxIsInf(rdp[2]) ||
        mxIsNaN(rdp[2]) ||
        (*rdp < 1) ||
        (*rdp > ((double) dx)) ||
        (rdp[1] < 1) ||
        (rdp[1] > ((double) dy)) ||
        (rdp[2] < 1) ||
        (rdp[2] > ((double) dz)))
        mexErrMsgTxt("Bad start position given.");
    sx = (unsigned long) *rdp;
    sy = (unsigned long) rdp[1];
    sz = (unsigned long) rdp[2];

    /* max iteration override */
    if ((nrhs > 3) &&
        (mxIsDouble(prhs[3])) &&
        (mxGetNumberOfElements(prhs[3]) == 1)) {
        rdp = (const double *) mxGetPr(prhs[3]);
        if (!mxIsInf(*rdp) &&
            !mxIsNaN(*rdp) &&
            (*rdp >= 1.0))
            mi = (unsigned long) *rdp;
    }

    /* check connectivity flag */
    rdp = (const double *) mxGetPr(prhs[2]);
    if (mxIsInf(*rdp) ||
        mxIsNaN(*rdp) ||
        (*rdp < 1.0) ||
        (*rdp > 7.0))
        mexErrMsgTxt("Invalid connectivity flag.");
    if (*rdp == 1.0)
        md = 3;
    else if (*rdp == 2.0)
        md = 9;
    else if (*rdp == 3.0)
        md = 13;
    else
        md = 2;

    /* get largest possible plane size */
    /* dy MUST be in product */
    if (dx < dy) {

        /* product of dy and dz if dx is smallest */
        if (dx < dz)
            mn = dy * dz;

        /* otherwise ... */
        else
            mn = dy * dx;

    /* dx MUST be in product */
    } else {

        /* product of dx and dz if dy is smallest */
        if  (dy < dz)
            mn = dx * dz;

        /* otherwise ... */
        else
            mn = dx * dy;
    }
    mn *= 9;

    /* create output argument */
    plhs[0] = mxCreateNumericArray(3, idim, mxLOGICAL_CLASS, mxREAL);
    if (plhs[0] == NULL)
        mexErrMsgTxt("Error allocation memory for output argument.");

    /* check whether sx, sy, sz is not set in volume */
    data = (const unsigned char *) mxGetData(prhs[0]);
    if (data[sx + sy * dx + sz * dxy - (1 + dx + dxy)] == 0) {

        /* possible create second arg, too */
        if (nlhs > 1)
            CREATE_SCALAR(plhs[1], 0.0)

         /* early return */
        return;
    }

    /* reserve space next iteration pointers */
    next1 = (signed long *) mxCalloc(mn, sizeof(signed long));
    next2 = (signed long *) mxCalloc(mn, sizeof(signed long));
    if ((next1 == NULL) ||
        (next2 == NULL))
        mexErrMsgTxt("Error allocating temporary memory.");

    /* increase numbers accordingly */
    if (idim[0] > 1)
        dx += 2;
    else {
        --sx;
        toff[0] = 0;
    }
    if (idim[1] > 1)
        dy += 2;
    else {
        --sy;
        toff[1] = 0;
    }
    if (idim[2] > 1)
        dz += 2;
    else {
        --sz;
        toff[2] = 0;
    }
    dxy = dx * dy;
    *next1 = sx + sy * dx + sz * dxy;
    tdim[0] = dx;
    tdim[1] = dy;
    tdim[2] = dz;

    /* fill direction increment array */
    nextb = dlxyz;
    np = 0;
    for (dc = 0; dc < md; ++dc) {
        if (((idim[0] > 1) ||
             (dlx[dc] == 0)) &&
            ((idim[1] > 1) ||
             (dly[dc] == 0)) &&
            ((idim[2] > 1) ||
             (dlz[dc] == 0))) {
            nextb[np] = dlx[dc] + dly[dc] * dx + dlz[dc] * dxy;
            nextb[np+1] = -nextb[np];
            np += 2;
        }
    }
    if ((*rdp > 5.0) &&
        (idim[0] > 1) &&
        (idim[1] > 1)) {
        for (dc = 3; dc < 5; ++dc) {
            nextb[np] = dlx[dc] + dly[dc] * dx + dlz[dc] * dxy;
            nextb[np+1] = -nextb[np];
            np += 2;
        }
    }
    md = np - 1;

    /* copy into temporary array */
    copy = (unsigned char *) mxCalloc(dxy * dz, sizeof(unsigned char));
    if (copy == NULL)
        mexErrMsgTxt("Error allocating temporary memory");
    ff_copy_intopart(idim[0], idim[1], idim[2], data, tdim, copy, toff);
    if ((*rdp != 5.0) &&
        (*rdp != 7.0)) {
        copy[sx + sy * dx + sz * dxy] = 2;
        nextc = 1;
    } else {
        for (np = 0, nextp = sx + sy * dx + (tdim[2] - 1) * dxy; nextp > 0; nextp -= dxy) {
            if (copy[nextp] == 1) {
                copy[nextp] = 2;
                next1[np++] = nextp;
            }
        }
        nextc = np;
        if (nextc == 0)
            donext = 0;
    }

    /* iterate while new points found */
    nextb = next1;
    for (; donext; ) {

        /* iterate over indices to check */
        for (nextcb = nextc, c = nextcb, nextc = 0; c > 0; --c) {

            /* get voxel position and iterate over directions */
            nextp = *nextb++;
            for (dc = md;  dc >= 0; --dc) {

                /* check voxel with increment */
                nextpi = nextp + dlxyz[dc];
                if (copy[nextpi] == 1) {

                        /* set data to 2 and add index to list */
                        copy[nextpi] = 2;
                        next2[nextc++] = nextpi;
                }
            }
        }

        /* swap pointers */
        nextb = next2;
        next2 = next1;
        next1 = nextb;

        /* if no more found, leave loop */
        --mi;
        if ((nextc == 0) || (mi == 0))
            donext = 0;
    }

    /* free lists */
    mxFree(next1);
    mxFree(next2);

    /* copy to output copy */
    np = ff_set1_frompart_eq(dx, dy, copy, idim, (unsigned char*) mxGetData(plhs[0]), toff, 2);
    mxFree(copy);

    /* if second output needed create a scalar */
    if (nlhs > 1)
        CREATE_SCALAR(plhs[1], np);
}

void ff_copy_intopart(unsigned long x1, unsigned long y1, unsigned long z1, const unsigned char *sd, const int *tnd, unsigned char *td, const unsigned long *off)
{
    unsigned long xc, yc, zc, x2, y2, z2, xy, to, yo, zo;
    x2 = *tnd++;
    y2 = *tnd++;
    z2 = *tnd;
    xy = x2 * y2;
    to = off[0] + off[1] * x2 + off[2] * xy - 1;
    for (zc = 0; zc < z1; ++zc) {
        zo = to + zc * xy;
        for (yc = 0; yc < y1; ++yc) {
            yo = zo + yc * x2;
            for (xc = 0; xc < x1; ++xc)
                td[++yo] = *sd++;
        }
    }
}

unsigned long ff_set1_frompart_eq(unsigned long x1, unsigned long y1, const unsigned char *sd, const int *tnd, unsigned char *td, const unsigned long *off, unsigned char thr)
{
    unsigned long rv, xc, yc, zc, xy, x2, y2, z2, to, yo, zo;
    rv = 0;
    xy = x1 * y1;
    x2 = *tnd++;
    y2 = *tnd++;
    z2 = *tnd;
	to = off[0] + off[1] * x1 + off[2] * xy - 1;
	for (zc = 0; zc < z2; ++zc) {
        zo = to + zc * xy;
        for (yc = 0; yc < y2; ++yc) {
            yo = zo + yc * x1;
            for (xc = 0; xc < x2; ++xc) {
                if (sd[++yo] == thr) {
                    *td = 1;
                    ++rv;
                }
                ++td;
            }
        }
    }
    return rv;
}

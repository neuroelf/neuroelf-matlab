/*

mesh_reconstruct -  getting a list of points and triangles from voxels

FORMAT:       [p, t, n] = mesh_reconstruct(vox [, tps])

Input fields:

      vox         3D uint8/logical data, 0:= background, >0:= volume
      tps         triangles per face ({2} or 4)

Output fields:

      p           Cx3 coordinates
      t           Tx3 triangles
      n           optional Cx2 neighbors array (for BV formats)


% Version:  v0.9b
% Build:    11050511
% Date:     Jul-28 2010, 5:30 PM EST
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
#include <stdio.h>

/* voxel surfaces:
 * 1: -X
 * 2: +X
 * 3: -Y
 * 4: +Y
 * 5: -Z
 * 6: +Z
 *
 * voxel points:
 * [0]: [0, 0, 0]
 * [1]: [1, 0, 0]
 * [2]: [0, 1, 0]
 * [3]: [1, 1, 0]
 * [4]: [0, 0, 1]
 * [5]: [1, 0, 1]
 * [6]: [0, 1, 1]
 * [7]: [1, 1, 1]
 *
 * voxel face edges:
 * [00]: s1, p[0] - p[4] : [0, 0, 0] - [0, 0, 1], vp[00], vp[01]
 * [01]: s1, p[2] - p[6] : [0, 1, 0] - [0, 1, 1], vp[02], vp[03]
 * [02]: s1, p[0] - p[2] : [0, 0, 0] - [0, 1, 0], vp[00], vp[02]
 * [03]: s1, p[4] - p[6] : [0, 0, 1] - [0, 1, 1], vp[01], vp[03]
 * [04]: s2, p[1] - p[5] : [1, 0, 0] - [1, 0, 1], vp[04], vp[05]
 * [05]: s2, p[3] - p[7] : [1, 1, 0] - [1, 1, 1], vp[06], vp[07]
 * [06]: s2, p[1] - p[3] : [1, 0, 0] - [1, 1, 0], vp[04], vp[06]
 * [07]: s2, p[5] - p[7] : [1, 0, 1] - [1, 1, 1], vp[05], vp[07]
 * [08]: s3, p[0] - p[4] : [0, 0, 0] - [0, 0, 1], vp[08], vp[09]
 * [09]: s3, p[1] - p[5] : [1, 0, 0] - [1, 0, 1], vp[10], vp[11]
 * [10]: s3, p[0] - p[1] : [0, 0, 0] - [1, 0, 0], vp[08], vp[10]
 * [11]: s3, p[4] - p[5] : [0, 0, 1] - [1, 0, 1], vp[09], vp[11]
 * [12]: s4, p[2] - p[6] : [0, 1, 0] - [0, 1, 1], vp[12], vp[13]
 * [13]: s4, p[3] - p[7] : [1, 1, 0] - [1, 1, 1], vp[14], vp[15]
 * [14]: s4, p[2] - p[3] : [0, 1, 0] - [1, 1, 0], vp[12], vp[14]
 * [15]: s4, p[6] - p[7] : [0, 1, 1] - [1, 1, 1], vp[13], vp[15]
 * [16]: s5, p[0] - p[2] : [0, 0, 0] - [0, 1, 0], vp[16], vp[17]
 * [17]: s5, p[1] - p[3] : [1, 0, 0] - [1, 1, 0], vp[18], vp[19]
 * [18]: s5, p[0] - p[1] : [0, 0, 0] - [1, 0, 0], vp[16], vp[18]
 * [19]: s5, p[2] - p[3] : [0, 1, 0] - [1, 1, 0], vp[17], vp[19]
 * [20]: s6, p[4] - p[6] : [0, 0, 1] - [0, 1, 1], vp[20], vp[21]
 * [21]: s6, p[5] - p[7] : [1, 0, 1] - [1, 1, 1], vp[22], vp[23]
 * [22]: s6, p[4] - p[5] : [0, 0, 1] - [1, 0, 1], vp[20], vp[22]
 * [23]: s6, p[6] - p[7] : [0, 1, 1] - [1, 1, 1], vp[21], vp[23]
 */
static const signed char mesh_pts_crd[72] = {
    0, 0, 0,
    0, 0, 1,
    0, 1, 0,
    0, 1, 1,
    1, 0, 0,
    1, 0, 1,
    1, 1, 0,
    1, 1, 1,
    0, 0, 0,
    0, 0, 1,
    1, 0, 0,
    1, 0, 1,
    0, 1, 0,
    0, 1, 1,
    1, 1, 0,
    1, 1, 1,
    0, 0, 0,
    0, 1, 0,
    1, 0, 0,
    1, 1, 0,
    0, 0, 1,
    0, 1, 1,
    1, 0, 1,
    1, 1, 1};

/* neighbors or points, 3 kind triplets x 24 points array:
 * 1. index in triplet: same voxel (adjacent face vertex)
 * 2. index in triplet: planar neighboring voxel (same orientation face)
 * 3. index in triplet: diagonal neighboring voxel face vertex
 * 4. neighboring point (on same voxel) for neighborhood information
 * ALL INDICES ARE GIVEN IN ANTI-CLOCKWISE DIRECTION! */
static const signed char mesh_neigh_pts[96] = {
     8,  2, 14,  1,
    20,  0, 18,  3,
    17,  3, 23,  0,
    13,  1, 11,  2,
    18,  5, 20,  6,
    11,  7, 13,  4,
    14,  4,  8,  7,
    23,  6, 17,  5,
    16,  9, 21, 10,
     1, 11,  7,  8,
     4,  8,  2, 11,
    22, 10, 19,  9,
     2, 14,  4, 13,
    21, 12, 16, 15,
    19, 15, 22, 12,
     7, 13,  1, 14,
     0, 18,  5, 17,
    12, 16,  9, 19,
    10, 19, 15, 16,
     6, 17,  3, 18,
     9, 21, 12, 22,
     3, 23,  6, 20,
     5, 20,  0, 23,
    15, 22, 10, 21};

/* search dim definition: 5-by-24 series,
   main dim+dir, neighbor search dim+dir, face point 3rd dim dir*/
static const signed char mesh_neigh_sdim[120] = {
     1, -1,  2, -1, -1,
     1, -1,  3,  1, -1,
     1, -1,  3, -1,  1,
     1, -1,  2,  1,  1,
     1,  1,  3, -1, -1,
     1,  1,  2, -1,  1,
     1,  1,  2,  1, -1,
     1,  1,  3,  1,  1,
     2, -1,  3, -1, -1,
     2, -1,  1, -1,  1,
     2, -1,  1,  1, -1,
     2, -1,  3,  1,  1,
     2,  1,  1, -1, -1,
     2,  1,  3,  1, -1,
     2,  1,  3, -1,  1,
     2,  1,  1,  1,  1,
     3, -1,  1, -1, -1,
     3, -1,  2,  1, -1,
     3, -1,  2, -1,  1,
     3, -1,  1,  1,  1,
     3,  1,  2, -1, -1,
     3,  1,  1, -1,  1,
     3,  1,  1,  1, -1,
     3,  1,  2,  1,  1};

/* function to check whether/set a point of an edge/is already set */
void mesh_checkaddpoint(int *voxfpts, const int *voxidx, const unsigned char *vox,
        int d1, int d2, int d3, int d12,
        signed short c1, signed short c2, signed short c3,
        int vi, signed char pi,
        int *pnp, double **pp1, double **pp2, double **pp3,
        int *neigh) {

    /* voxel+point chain */
    int vin;
    signed char pin;
    signed char maxnei = 12;
    bool keepsearching = 1;

    /* search dimensions and limits */
    signed char sd1, sd2, sd3;
    signed short dc1, dc2, dc3;
    signed short ld1, ld2, ld3, ad1, ad2;
    signed short df1 = 0, df2 = 0, df3 = 0;
    int fd1, fd2, fd3, vc;

    char vstr[256]; /* debugging */

    /* return if point already set */
    if (voxfpts[24 * vi + pi] >= 0)
        return;

    /* add point to array */
    **pp1 = ((double) c1 + mesh_pts_crd[3 * pi]) + 0.5;
    **pp2 = ((double) c2 + mesh_pts_crd[3 * pi + 1]) + 0.5;
    **pp3 = ((double) c3 + mesh_pts_crd[3 * pi + 2]) + 0.5;
    ++(*pp1);
    ++(*pp2);
    ++(*pp3);

    /* set next element to first element */
    vin = vi;
    pin = pi;

    /* arrange neighbor pointer */
    if (neigh != NULL)
        neigh = &neigh[12 * *pnp];

    /* repeat until all done */
    while (keepsearching) {

        /* all tries up? */
        if (maxnei < 0) {
            /* sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1 - 1) + (c2 - 1) * d2 + (c3 - 1) * d12],
                voxidx[(c1 - 1) + (c2    ) * d2 + (c3 - 1) * d12],
                voxidx[(c1 - 1) + (c2 + 1) * d2 + (c3 - 1) * d12]);
            mexWarnMsgTxt(vstr);
            sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1    ) + (c2 - 1) * d2 + (c3 - 1) * d12],
                voxidx[(c1    ) + (c2    ) * d2 + (c3 - 1) * d12],
                voxidx[(c1    ) + (c2 + 1) * d2 + (c3 - 1) * d12]);
            mexWarnMsgTxt(vstr);
            sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1 + 1) + (c2 - 1) * d2 + (c3 - 1) * d12],
                voxidx[(c1 + 1) + (c2    ) * d2 + (c3 - 1) * d12],
                voxidx[(c1 + 1) + (c2 + 1) * d2 + (c3 - 1) * d12]);
            mexWarnMsgTxt(vstr);
            mexWarnMsgTxt(" ");
            sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1 - 1) + (c2 - 1) * d2 + (c3    ) * d12],
                voxidx[(c1 - 1) + (c2    ) * d2 + (c3    ) * d12],
                voxidx[(c1 - 1) + (c2 + 1) * d2 + (c3    ) * d12]);
            mexWarnMsgTxt(vstr);
            sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1    ) + (c2 - 1) * d2 + (c3    ) * d12],
                voxidx[(c1    ) + (c2    ) * d2 + (c3    ) * d12],
                voxidx[(c1    ) + (c2 + 1) * d2 + (c3    ) * d12]);
            mexWarnMsgTxt(vstr);
            sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1 + 1) + (c2 - 1) * d2 + (c3    ) * d12],
                voxidx[(c1 + 1) + (c2    ) * d2 + (c3    ) * d12],
                voxidx[(c1 + 1) + (c2 + 1) * d2 + (c3    ) * d12]);
            mexWarnMsgTxt(vstr);
            mexWarnMsgTxt(" ");
            sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1 - 1) + (c2 - 1) * d2 + (c3 + 1) * d12],
                voxidx[(c1 - 1) + (c2    ) * d2 + (c3 + 1) * d12],
                voxidx[(c1 - 1) + (c2 + 1) * d2 + (c3 + 1) * d12]);
            mexWarnMsgTxt(vstr);
            sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1    ) + (c2 - 1) * d2 + (c3 + 1) * d12],
                voxidx[(c1    ) + (c2    ) * d2 + (c3 + 1) * d12],
                voxidx[(c1    ) + (c2 + 1) * d2 + (c3 + 1) * d12]);
            mexWarnMsgTxt(vstr);
            sprintf(vstr, " %6d %6d %6d",
                voxidx[(c1 + 1) + (c2 - 1) * d2 + (c3 + 1) * d12],
                voxidx[(c1 + 1) + (c2    ) * d2 + (c3 + 1) * d12],
                voxidx[(c1 + 1) + (c2 + 1) * d2 + (c3 + 1) * d12]);
            mexWarnMsgTxt(vstr);
            mexWarnMsgTxt(" "); */
            sprintf(vstr,
                "Error looking up neighborhood information for voxel %d (%d, %d, %d), point %d.",
                vi, c1 + 1, c2 + 1, c3 + 1, pi);
            mexErrMsgTxt(vstr);
        }
        --maxnei;

        /* take note of neighborhood information */
        if (neigh != NULL)
            *neigh++ = 24 * vin + ((int) mesh_neigh_pts[4 * pin + 3]);

        /* set point */
        voxfpts[24 * vin + pin] = *pnp;
        /* sprintf(vstr, "setting v(%d).p(%d) to %d", vin, pin, *pnp);
        mexWarnMsgTxt(vstr); */

        /* which dimensions and direction to use? */
        sd1 = mesh_neigh_sdim[5 * pin];
        sd2 = mesh_neigh_sdim[5 * pin + 2];
        sd3 = 6 - (sd1 + sd2);
        if (sd1 == 1) {

            /* d1 = 1, d2 = 2 */
            if (sd2 == 2) {
                ld1 = d1 - 1;
                fd1 = 1;
                dc1 = c1 + df1;
                ld2 = d2 - 1;
                fd2 = d1;
                dc2 = c2 + df2;
                ld3 = d3 - 1;
                fd3 = d12;
                dc3 = c3 + df3;

            /* d1 = 1, d2 = 3 */
            } else {
                ld1 = d1 - 1;
                fd1 = 1;
                dc1 = c1 + df1;
                ld2 = d3 - 1;
                fd2 = d12;
                dc2 = c3 + df3;
                ld3 = d2 - 1;
                fd3 = d1;
                dc3 = c2 + df2;
            }

        } else if (sd1 == 2) {

            /* d1 = 2, d2 = 1 */
            if (sd2 == 1) {
                ld1 = d2 - 1;
                fd1 = d1;
                dc1 = c2 + df2;
                ld2 = d1 - 1;
                fd2 = 1;
                dc2 = c1 + df1;
                ld3 = d3 - 1;
                fd3 = d12;
                dc3 = c3 + df3;

            /* d1 = 2, d2 = 3 */
            } else {
                ld1 = d2 - 1;
                fd1 = d1;
                dc1 = c2 + df2;
                ld2 = d3 - 1;
                fd2 = d12;
                dc2 = c3 + df3;
                ld3 = d1 - 1;
                fd3 = 1;
                dc3 = c1 + df1;
            }

        } else {

            /* d1 = 3, d2 = 1 */
            if (sd2 == 1) {
                ld1 = d3 - 1;
                fd1 = d12;
                dc1 = c3 + df3;
                ld2 = d1 - 1;
                fd2 = 1;
                dc2 = c1 + df1;
                ld3 = d2 - 1;
                fd3 = d1;
                dc3 = c2 + df2;

            /* d1 = 3, d2 = 2 */
            } else {
                ld1 = d3 - 1;
                fd1 = d12;
                dc1 = c3 + df3;
                ld2 = d2 - 1;
                fd2 = d1;
                dc2 = c2 + df2;
                ld3 = d1 - 1;
                fd3 = 1;
                dc3 = c1 + df1;
            }
        }
        /* compute base coordinate */
        vc = ((int) dc1) * ((int) fd1) + ((int) dc2) * ((int) fd2) + ((int) dc3) * ((int) fd3);

        /* the limits depend on their direction */
        if (mesh_neigh_sdim[5 * pin + 1] < 0) {
            ld1 = 0;
            fd1 = -fd1;
        }
        if (mesh_neigh_sdim[5 * pin + 3] < 0) {
            ld2 = 0;
            fd2 = -fd2;
        }
        if (mesh_neigh_sdim[5 * pin + 4] < 0) {
            ld3 = 0;
            fd3 = -fd3;
        }

        /* sprintf(vstr, "checking point %d on voxel %d (%d, %d, %d) in dim (%d,%d,%d) dir(%d,%d,%d => %d,%d,%d)",
            pin, vin, dc1, dc2, dc3, sd1, sd2, sd3, (int) (ld1 > 0), (int) (ld2 > 0), (int) (ld3 > 0), fd1, fd2, fd3);
        mexWarnMsgTxt(vstr); */

        /* neighbor in 1st+2nd direction (potential diagonal neighbor) */
        if ((dc1 != ld1) &&
            (dc2 != ld2) &&
            (voxidx[vc + fd1 + fd2] >= 0) &&

            /* and connected by the neighbor in 2nd dim directly */
            ((vox[vc + fd2] > 0))) {

            /* or via one of two side paths along the 3rd dimension
             ((dc3 != ld3) &&
              (vox[vc + fd3] > 0) &&
              (vox[vc + fd1 + fd2 + fd3] > 0) */

            /* strict checking:
              &&
              ((vox[vc + fd1 + fd3] > 0) ||
               (vox[vc + fd2 + fd3] > 0)) */


            /* next voxel is in both directions */
            vin = voxidx[vc + fd1 + fd2];
            ad1 = mesh_neigh_sdim[5 * pin + 1];
            ad2 = mesh_neigh_sdim[5 * pin + 3];

            /* next point to check is 3rd of the triplet */
            pin = mesh_neigh_pts[4 * pin + 2];

            /* sprintf(vstr, "diag: %d, %d", vin, pin);
            mexWarnMsgTxt(vstr); */

        /* no neighbor in 2nd direction */
        } else if ((dc2 == ld2) || (voxidx[vc + fd2] < 0)) {

            /* it's the same voxel */
            ad1 = 0;
            ad2 = 0;

            /* and the 1st of the triplet */
            pin = mesh_neigh_pts[4 * pin];

            /* sprintf(vstr, "same: %d, %d", vin, pin);
            mexWarnMsgTxt(vstr); */

        /* otherwise => no diagonal neighbor but direct neighbor in 2nd dim */
        } else {

            /* next voxel is in 2nd direction */
            vin = voxidx[vc + fd2];
            ad1 = 0;
            ad2 = mesh_neigh_sdim[5 * pin + 3];

            /* and the second of the triplet */
            pin = mesh_neigh_pts[4 * pin + 1];

            /* sprintf(vstr, "plan: %d, %d", vin, pin);
            mexWarnMsgTxt(vstr); */
        }

        /* add to correct coordinate */
        if (ad1 != 0) {
            if (sd1 == 1)
                df1 += ad1;
            else if(sd1 == 2)
                df2 += ad1;
            else
                df3 += ad1;
        }
        if (ad2 != 0) {
            if (sd2 == 1)
                df1 += ad2;
            else if (sd2 == 2)
                df2 += ad2;
            else
                df3 += ad2;
        }

        /* check values */
        if ((df1 < -1) ||
            (df1 > 1) ||
            (df2 < -1) ||
            (df2 > 1) ||
            (df3 < -1) ||
            (df3 > 1))
            mexErrMsgTxt("Invalid neighborhood chain.");

        /* loop complete? */
        if ((vin == vi) &&
            (pin == pi))
            keepsearching = 0;
    }
    /* mexWarnMsgTxt("found"); */

    /* increase number of points */
    ++(*pnp);
}

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* input data pointer */
    const unsigned char *vox = NULL;

    /* size arguments */
    const int *idim = NULL;
    signed short dim1, dim2, dim3, lim1, lim2, lim3;
    int dim12;

    /* triangles per vertex option data */
    const double *dtps = NULL;
    int tps = 2;

    /* 3D voxel->index copy (same size as vox) */
    int *voxidx = NULL;

    /* general counter */
    int c;

    /* coordinate counters */
    signed short c1, c2, c3;

    /* voxel coordinate (sum of c1 + c2 * dim1 + c3 * dim12) */
    int vc;

    /* number of found voxels, surfaces (-> triangles), estimated points */
    int nv, nf, nt, np, pc;

    /* 2D 24-by-number of found voxel face-points buffer */
    int *voxfpts = NULL, tv;

    /* storage array pointers */
    double *tri1 = NULL, *tri2, *tri3,
           *pts = NULL, *pts1, *pts2, *pts3;

    /* central coordinate index */
    double cpi = 0.0;

    /* neighbors cell array arguments */
    int *neigh = NULL, *pneigh;
    mxArray *cneigh = NULL;
    int cdi[2] = {0, 2};
    int neighvi, neighvn, lneigh[12];
    signed char nneigh, maxnei;

    /* char vstr[256]; */

    /* check number of in/out arguments */
    if ((nrhs < 1) ||
        (nlhs < 2) ||
        (nlhs > 3))
        mexErrMsgTxt("Bad number of input/output arguments.");

    /* check type, dims, and elements of input volume */
    if (((mxGetClassID(*prhs) != mxUINT8_CLASS) &&
         (mxGetClassID(*prhs) != mxLOGICAL_CLASS)) ||
        (mxGetNumberOfDimensions(*prhs) != 3) ||
        (mxGetNumberOfElements(*prhs) == 0))
        mexErrMsgTxt("Bad voxel data input argument.");

    /* tps argument given? */
    if ((nrhs > 1) &&
        (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfElements(prhs[1]) == 1)) {

        /* check for value := 4 */
        dtps = (const double*) mxGetData(prhs[1]);
        if (!mxIsInf(*dtps) &&
            !mxIsNaN(*dtps) &&
            (*dtps == 4.0))

            /* overwrite default tps value */
            tps = 4;
    }

    /* get data pointer and dimensions */
    vox = (const unsigned char*) mxGetData(*prhs);
    idim = mxGetDimensions(*prhs);

    /* check dimensions */
    if (*idim > 32767)
        mexErrMsgTxt("Each dimension is restricted to 32767 voxels.");
    dim1 = (signed short) *idim++;
    if (*idim > 32767)
        mexErrMsgTxt("Each dimension is restricted to 32767 voxels.");
    dim2 = (signed short) *idim++;
    if (*idim > 32767)
        mexErrMsgTxt("Each dimension is restricted to 32767 voxels.");
    dim3 = (signed short) *idim;

    /* multiply to shorthand dim12 */
    dim12 = ((int) dim1) * ((int) dim2);

    /* also assign limit values for voxel search */
    lim1 = dim1 - 1;
    lim2 = dim2 - 1;
    lim3 = dim3 - 1;

    /* allocate temporary buffer to voxel index volume */
    voxidx = (int *) mxCalloc(dim12 * ((int) dim3), sizeof(int));
    if (voxidx == NULL)
        mexErrMsgTxt("Error allocating voxel pointer volume.");

    /* fill with -1's */
    for (c = (dim12 * ((int) dim3)) - 1; c >= 0; --c)
        voxidx[c] = -1;

    /* count voxels and faces */
    nf = nv = 0;
    for (c3 = 0; c3 < dim3; ++c3)
        for (c2 = 0; c2 < dim2; ++c2)
            for (c1 = 0; c1 < dim1; ++c1) {

                /* get voxel coordinate */
                vc = ((int) c3) * dim12 + ((int) c2) * ((int) dim1) + ((int) c1);

                /* test if voxel has any face to draw */
                if ((vox[vc] > 0) &&
                    (((c1 == 0) || (vox[vc - 1] == 0)) ||
                     ((c1 == lim1) || (vox[vc + 1] == 0)) ||
                     ((c2 == 0) || (vox[vc - dim1] == 0)) ||
                     ((c2 == lim2) || (vox[vc + dim1] == 0)) ||
                     ((c3 == 0) || (vox[vc - dim12] == 0)) ||
                     ((c3 == lim3) || (vox[vc + dim12] == 0)))) {

                    /* mark voxel as border voxel */
                    voxidx[vc] = nv;
                    ++nv;

                    /* test each surface separately to count faces */
                    if ((c1 == 0) || (vox[vc - 1] == 0))
                        ++nf;
                    if ((c1 == lim1) || (vox[vc + 1] == 0))
                        ++nf;
                    if ((c2 == 0) || (vox[vc - dim1] == 0))
                        ++nf;
                    if ((c2 == lim2) || (vox[vc + dim1] == 0))
                        ++nf;
                    if ((c3 == 0) || (vox[vc - dim12] == 0))
                        ++nf;
                    if ((c3 == lim3) || (vox[vc + dim12] == 0))
                        ++nf;
                }
            }

    /* empty set */
    if (nv == 0) {

        /* create empty outputs */
        *plhs = mxCreateDoubleMatrix(0, 3, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(0, 3, mxREAL);
        return;
    }

    /* compute number of triangles */
    nt = tps * nf;
    plhs[1] = mxCreateDoubleMatrix(nt, 3, mxREAL);
    if (plhs[1] == NULL) {
        mxFree(voxidx);
        mexErrMsgTxt("Error allocating triangles array.");
    }

    /* get pointers */
    tri1 = (double *) mxGetData(plhs[1]);
    if (tri1 == NULL) {
        mxFree(voxidx);
        mexErrMsgTxt("Error retrieving triangles array data pointer.");
    }
    tri2 = &tri1[nt];
    tri3 = &tri2[nt];

    /* allocate voxel face-points buffer (used to collapse true neighbors) */
    voxfpts = (int *) mxCalloc(24 * nv, sizeof(int));
    if (voxfpts == NULL) {
        mxFree(voxidx);
        mexErrMsgTxt("Error allocating temporary collapsing array.");
    }

    /* compute max. estimated number of points */
    if (tps == 2)
        np = 4 * nf;
    else
        np = 5 * nf;
    /* sprintf(vstr, "%d voxels with %d surfaces -> %d triangles with max. %d points",
        nv, nf, nt, np);
    mexWarnMsgTxt(vstr); */

    /* allocate temporary points buffer */
    pts = (double *) mxCalloc(3 * np, sizeof(double));
    if (pts == NULL) {
        mxFree(voxfpts);
        mxFree(voxidx);
        mexErrMsgTxt("Error allocating temporary points array.");
    }

    /* and if required also an array for neighbors */
    if (nlhs > 2) {
        neigh = (int *) mxCalloc(12 * np, sizeof(int));
        if (neigh == NULL) {
            mxFree(pts);
            mxFree(voxfpts);
            mxFree(voxidx);
            mexErrMsgTxt("Error allocating memory for neighborhood information.");
        }

        /* pre-set points/neighbor information */
        for (c = 12 * np - 1; c >= 0; --c)
            neigh[c] = -1;
    }

    /* get shorthand pointers */
    pts1 = pts;
    pts2 = &pts1[np];
    pts3 = &pts2[np];

    /* pre-set points/neighbor information */
    for (c = 24 * nv - 1; c >= 0; --c)
        voxfpts[c] = -1;

    /* initialize counters */
    pc = 0;

    /* iterate over volume again */
    for (c3 = 0; c3 < dim3; ++c3)
        for (c2 = 0; c2 < dim2; ++c2)
            for (c1 = 0; c1 < dim1; ++c1) {

                /* get voxel coordinate */
                vc = c3 * dim12 + c2 * dim1 + c1;

                /* voxel needs attention */
                if (voxidx[vc] >= 0) {

                    /* get voxel number */
                    nv = voxidx[vc];

                    /* we need a surface in -X direction (s1) */
                    if ((c1 == 0) || (vox[vc - 1] == 0)) {

                        /* index to points */
                        tv = 24 * nv;

                        /* if we need a central voxel */
                        if (tps == 4) {

                            /* add it */
                            *pts1++ = ((double) c1) + 0.5;
                            *pts2++ = (double) (c2 + 1);
                            *pts3++ = (double) (c3 + 1);

                            /* and add to neighbors for that point */
                            if (neigh != NULL) {
                                pneigh = &neigh[12 * pc];
                                *pneigh++ = tv;
                                *pneigh++ = tv + 2;
                                *pneigh++ = tv + 3;
                                *pneigh++ = tv + 1;
                            }

                            /* and note the number of the point! */
                            cpi = (double) ++pc;
                        }

                        /* ensure neighbor points for face points 0 - 3 */
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 0, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 1, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 2, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 3, &pc, &pts1, &pts2, &pts3, neigh);

                        /* two triangles */
                        if (tps == 2) {
                            *tri1++ = (double) (voxfpts[tv] + 1);
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv + 3] + 1);
                            *tri1++ = (double) (voxfpts[tv + 3] + 1);
                            *tri2++ = (double) (voxfpts[tv + 1] + 1);
                            *tri3++ = (double) (voxfpts[tv] + 1);

                        /* four triangles */
                        } else {
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv + 3] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 1] + 1);
                            *tri3++ = (double) (voxfpts[tv] + 1);
                        }
                    }

                    /* we need a surface in +X direction (s2) */
                    if ((c1 == lim1) || (vox[vc + 1] == 0)) {

                        /* index to points */
                        tv = 24 * nv + 4;

                        /* if we need a central voxel */
                        if (tps == 4) {

                            /* add it */
                            *pts1++ = ((double) c1) + 1.5;
                            *pts2++ = (double) (c2 + 1);
                            *pts3++ = (double) (c3 + 1);

                            /* and add to neighbors for that point */
                            if (neigh != NULL) {
                                pneigh = &neigh[12 * pc];
                                *pneigh++ = tv;
                                *pneigh++ = tv + 1;
                                *pneigh++ = tv + 3;
                                *pneigh++ = tv + 2;
                            }

                            /* and note the number of the point! */
                            cpi = (double) ++pc;
                        }

                        /* ensure neighbor points for face points 4 - 7 */
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 4, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 5, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 6, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 7, &pc, &pts1, &pts2, &pts3, neigh);

                        /* two triangles */
                        if (tps == 2) {
                            *tri1++ = (double) (voxfpts[tv] + 1);
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);
                            *tri1++ = (double) (voxfpts[tv + 3] + 1);
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);

                        /* four triangles */
                        } else {
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 1] + 1);
                            *tri3++ = (double) (voxfpts[tv + 3] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv] + 1);
                        }
                    }

                    /* we need a surface in -Y direction (s3) */
                    if ((c2 == 0) || (vox[vc - dim1] == 0)) {

                        /* index to points */
                        tv = 24 * nv + 8;

                        /* if we need a central voxel */
                        if (tps == 4) {

                            /* add it */
                            *pts1++ = (double) (c1 + 1);
                            *pts2++ = ((double) c2) + 0.5;
                            *pts3++ = (double) (c3 + 1);

                            /* and add to neighbors for that point */
                            if (neigh != NULL) {
                                pneigh = &neigh[12 * pc];
                                *pneigh++ = tv;
                                *pneigh++ = tv + 1;
                                *pneigh++ = tv + 3;
                                *pneigh++ = tv + 2;
                            }

                            /* and note the number of the point! */
                            cpi = (double) ++pc;
                        }

                        /* ensure neighbor points for face points 8 - 11 */
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 8, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 9, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 10, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 11, &pc, &pts1, &pts2, &pts3, neigh);

                        /* two triangles */
                        if (tps == 2) {
                            *tri1++ = (double) (voxfpts[tv] + 1);
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);
                            *tri1++ = (double) (voxfpts[tv + 3] + 1);
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);

                        /* four triangles */
                        } else {
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 1] + 1);
                            *tri3++ = (double) (voxfpts[tv + 3] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv] + 1);
                        }
                    }

                    /* we need a surface in +Y direction (s4) */
                    if ((c2 == lim2) || (vox[vc + dim1] == 0)) {

                        /* index to points */
                        tv = 24 * nv + 12;

                        /* if we need a central voxel */
                        if (tps == 4) {

                            /* add it */
                            *pts1++ = (double) (c1 + 1);
                            *pts2++ = ((double) c2) + 1.5;
                            *pts3++ = (double) (c3 + 1);

                            /* and add to neighbors for that point */
                            if (neigh != NULL) {
                                pneigh = &neigh[12 * pc];
                                *pneigh++ = tv;
                                *pneigh++ = tv + 2;
                                *pneigh++ = tv + 3;
                                *pneigh++ = tv + 1;
                            }

                            /* and note the number of the point! */
                            cpi = (double) ++pc;
                        }

                        /* ensure neighbor points for face points 12 - 15 */
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 12, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 13, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 14, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 15, &pc, &pts1, &pts2, &pts3, neigh);

                        /* two triangles */
                        if (tps == 2) {
                            *tri1++ = (double) (voxfpts[tv] + 1);
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);
                            *tri1++ = (double) (voxfpts[tv + 2] + 1);
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);

                        /* four triangles */
                        } else {
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv + 3] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 1] + 1);
                            *tri3++ = (double) (voxfpts[tv] + 1);
                        }
                    }

                    /* we need a surface in -Z direction (s5) */
                    if ((c3 == 0) || (vox[vc - dim12] == 0)) {

                        /* index to points */
                        tv = 24 * nv + 16;

                        /* if we need a central voxel */
                        if (tps == 4) {

                            /* add it */
                            *pts1++ = (double) (c1 + 1);
                            *pts2++ = (double) (c2 + 1);
                            *pts3++ = ((double) c3) + 0.5;

                            /* and add to neighbors for that point */
                            if (neigh != NULL) {
                                pneigh = &neigh[12 * pc];
                                *pneigh++ = tv;
                                *pneigh++ = tv + 2;
                                *pneigh++ = tv + 3;
                                *pneigh++ = tv + 1;
                            }

                            /* and note the number of the point! */
                            cpi = (double) ++pc;
                        }

                        /* ensure neighbor points for face points 16 - 19 */
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 16, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 17, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 18, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 19, &pc, &pts1, &pts2, &pts3, neigh);

                        /* two triangles */
                        if (tps == 2) {
                            *tri1++ = (double) (voxfpts[tv] + 1);
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);
                            *tri1++ = (double) (voxfpts[tv + 3] + 1);
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);

                        /* four triangles */
                        } else {
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv + 3] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 1] + 1);
                            *tri3++ = (double) (voxfpts[tv] + 1);
                        }
                    }

                    /* we need a surface in +Z direction (s6) */
                    if ((c3 == lim3) || (vox[vc + dim12] == 0)) {

                        /* index to points */
                        tv = 24 * nv + 20;

                        /* if we need a central voxel */
                        if (tps == 4) {

                            /* add it */
                            *pts1++ = (double) (c1 + 1);
                            *pts2++ = (double) (c2 + 1);
                            *pts3++ = ((double) c3) + 1.5;

                            /* and add to neighbors for that point */
                            if (neigh != NULL) {
                                pneigh = &neigh[12 * pc];
                                *pneigh++ = tv;
                                *pneigh++ = tv + 1;
                                *pneigh++ = tv + 3;
                                *pneigh++ = tv + 2;
                            }

                            /* and note the number of the point! */
                            cpi = (double) ++pc;
                        }

                        /* ensure neighbor points for face points 20 - 23 */
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 20, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 21, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 22, &pc, &pts1, &pts2, &pts3, neigh);
                        mesh_checkaddpoint(voxfpts, voxidx, vox, dim1, dim2, dim3, dim12, c1, c2, c3, nv, 23, &pc, &pts1, &pts2, &pts3, neigh);

                        /* two triangles */
                        if (tps == 2) {
                            *tri1++ = (double) (voxfpts[tv + 2] + 1);
                            *tri2++ = (double) (voxfpts[tv + 1] + 1);
                            *tri3++ = (double) (voxfpts[tv + 3] + 1);
                            *tri1++ = (double) (voxfpts[tv + 1] + 1);
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv] + 1);

                        /* four triangles */
                        } else {
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv] + 1);
                            *tri3++ = (double) (voxfpts[tv + 1] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 1] + 1);
                            *tri3++ = (double) (voxfpts[tv + 3] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 3] + 1);
                            *tri3++ = (double) (voxfpts[tv + 2] + 1);
                            *tri1++ = cpi;
                            *tri2++ = (double) (voxfpts[tv + 2] + 1);
                            *tri3++ = (double) (voxfpts[tv] + 1);
                        }
                    }
                }
            }

    /* clear first temporary array (no longer needed) */
    mxFree(voxidx);

    /* neighbors also? */
    if (nlhs > 2) {

        /* allocate cell array */
        *cdi = pc;
        plhs[2] = mxCreateCellArray(2, cdi);
        *cdi = 1;

        /* iterate over the points' neighbors lists! */
        for (c = 0; c < pc; ++c) {

            /* get shorthand */
            pneigh = &neigh[12 * c];

            /* then copy until -1, but only 12 max! */
            for (maxnei = 12, neighvi = -1, nneigh = 0; *pneigh >= 0; ) {
                neighvn = voxfpts[*pneigh++];
                if (neighvn != neighvi) {
                    neighvi = neighvn;
                    lneigh[nneigh++] = neighvi;
                }
                --maxnei;
                if (maxnei == 0)
                    break;
            }

            /* ensure no double point is given */
            if ((nneigh > 1) &&
                (*lneigh == lneigh[nneigh - 1]))
                --nneigh;

            /* create required cell arrays */
            cneigh = mxCreateDoubleMatrix(1, 1, mxREAL);
            if (cneigh == NULL) {
                mxFree(neigh);
                mxFree(voxfpts);
                mxFree(pts);
                mexErrMsgTxt("Error allocating 1x1 number of neighbors scalar.");
            }

            /* simply set (without too much ado) */
            *((double *) mxGetData(cneigh)) = (double) nneigh;
            mxSetCell(plhs[2], c, cneigh);

            /* then create the actual list of neighbors array */
            cneigh = mxCreateDoubleMatrix(1, (int) nneigh, mxREAL);
            if (cneigh == NULL) {
                mxFree(neigh);
                mxFree(voxfpts);
                mxFree(pts);
                mexErrMsgTxt("Error allocating 1xN neighbors list.");
            }

            /* get pointer and copy */
            tri1 = mxGetData(cneigh);
            for (--nneigh ; nneigh >= 0; --nneigh)
                tri1[nneigh] = (double) (lneigh[nneigh] + 1);

            /* then set cell */
            mxSetCell(plhs[2], pc + c, cneigh);
        }

        /* free neighbor memory */
        mxFree(neigh);
    }

    /* free one more temporary array */
    mxFree(voxfpts);

    /* create final points output */
    *plhs = mxCreateDoubleMatrix(pc, 3, mxREAL);
    if (*plhs == NULL) {
        mxFree(pts);
        mexErrMsgTxt("Error allocating points output.");
    }
    /* sprintf(vstr, "%d estimated points -> %d real points", np, pc);
    mexWarnMsgTxt(vstr); */

    tri1 = (double *) mxGetData(*plhs);
    if (tri1 == NULL) {
        mxFree(pts);
        mexErrMsgTxt("Error getting points output data pointer.");
    }
    tri2 = &tri1[pc];
    tri3 = &tri2[pc];

    /* copy points */
    pts1 = pts;
    pts2 = &pts1[np];
    pts3 = &pts2[np];
    for (c = 0; c < pc; ++c) {
        *tri1++ = *pts1++;
        *tri2++ = *pts2++;
        *tri3++ = *pts3++;
    }

    /* free temporary points */
    mxFree(pts);
}

/*

applyspmsnc.c -- internal matrix multiplication for SPM norm

FORMAT:       nc = applyspmsnc(c, tp, tdim [, itrf [, otrf]])

Input fields:

      c           Cx3 coordinates or 4x3 coordinate range
      tp          transformation parameters (sn mat-struct field Tr)
      tdim        spatial dimension of template (sn mat-struct field
                  VG(1).dim
      itrf        input transformation (sn mat-struct field
                  inv(VG(1).mat)
      otrf        output transformation (sn mat-struct field Affine,
                  potentially pre-multiplied with
                  inv(dataset.mat) * normimage.mat)

% Version:  v0.9d
% Build:    14062015
% Date:     Jun-20 2014, 3:52 PM EST
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
#include "math.h"
#include "quattrans.h"
#define TR_PI 3.14159265358979323846

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    const int *idim, *tdim;
    int c1, c2, c3, cc, d1, d2, d3, d12, d123, nc, r1 = 0, r2 = 0, r3 = 0, s1, s2, s3;
    double cvx, cvy, cvz, tv;
    const double *cx, *cy, *cz, *tb, *tx, *ty, *tz;
    double *ocx, *ocy, *ocz,
           *bx = NULL, *by = NULL, *bz = NULL, *bxb, *byb, *bzb, *bzbb,
           *xx = NULL, *xy = NULL, *xz = NULL, *xxb, *xyb, *xzb,
            ds1, ds2, ds3, rsh1, rsh2, rsh3,
            rf1 = 0.0, rf2 = 0.0, rf3 = 0.0,
            rp1 = 0.0, rp2 = 0.0, rp3 = 0.0,
            rs1 = 0.0, rs2 = 0.0, rs3 = 0.0;
    int mnout[2] = {0, 3};
    bool q1 = 0, q2 = 0, cr = 0;
    /* char vstr[256]; */

    QUATVARS
    QUATVARS2

    /* argument check */
    if (nrhs < 3 || nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");
    if ((mxGetNumberOfDimensions(*prhs) != 2) ||
        (mxGetClassID(*prhs) != mxDOUBLE_CLASS) ||
        (mxGetNumberOfDimensions(prhs[1]) != 4) ||
        (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ||
        (mxGetNumberOfElements(prhs[1]) == 0) ||
        (mxGetNumberOfDimensions(prhs[2]) != 2) ||
        (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ||
        (mxGetNumberOfElements(prhs[2]) != 3))
        mexErrMsgTxt("Invalid input argument.");
    idim = mxGetDimensions(*prhs);
    if (idim[1] != 3)
        mexErrMsgTxt("c argument must be of Cx3 size.");
    tdim = mxGetDimensions(prhs[1]);
    if (tdim[3] != 3)
        mexErrMsgTxt("tr argument must be of XxYxZx3 size.");

    /* get size of TR volume */
    cx = (const double*) mxGetPr(prhs[2]);
    if (cx == NULL)
        mexErrMsgTxt("Error getting tr volume size data pointer.");
    if (mxIsInf(*cx) ||
        mxIsNaN(*cx) ||
        (*cx < 1.0) ||
        (*cx > 1024.0))
        mexErrMsgTxt("Invalid 1st dim for tr volume size.");
    s1 = (int) (*cx++);
    if (mxIsInf(*cx) ||
        mxIsNaN(*cx) ||
        (*cx < 1.0) ||
        (*cx > 1024.0))
        mexErrMsgTxt("Invalid 2nd dim for tr volume size.");
    s2 = (int) (*cx++);
    if (mxIsInf(*cx) ||
        mxIsNaN(*cx) ||
        (*cx < 1.0) ||
        (*cx > 1024.0))
        mexErrMsgTxt("Invalid 3rd dim for tr volume size.");
    s3 = (int) *cx;

    /* get number of coordinates and tr dims */
    nc = *idim;
    d1 = *tdim++;
    d2 = *tdim++;
    d3 = *tdim;

    /* compute helper dims */
    d12 = d1 * d2;
    d123 = d12 * d3;

    /* transformation matrices? */
    if ((nrhs > 3) &&
        (mxGetClassID(prhs[3]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfDimensions(prhs[3]) == 2) &&
        (mxGetM(prhs[3]) == 4) &&
        (mxGetN(prhs[3]) == 4)) {
        cx = (const double*) mxGetPr(prhs[3]);
        if (cx != NULL) {
            for (c1 = 0; c1 < 16; ++c1)
                if (mxIsInf(cx[c1]) ||
                    mxIsNaN(cx[c1]))
                    break;
            if ((cx[3] != 0.0) ||
                (cx[7] != 0.0) ||
                (cx[11] != 0.0) ||
                (cx[15] != 1.0))
                c1 = 0;
            if (c1 > 15) {
                GET_QUATVARS(cx)
                q1 = 1;
            }
        }
    }
    if ((nrhs > 4) &&
        (mxGetClassID(prhs[4]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfDimensions(prhs[4]) == 2) &&
        (mxGetM(prhs[4]) == 4) &&
        (mxGetN(prhs[4]) == 4)) {
        cx = (const double*) mxGetPr(prhs[4]);
        if (cx != NULL) {
            for (c1 = 0; c1 < 16; ++c1)
                if (mxIsInf(cx[c1]) ||
                    mxIsNaN(cx[c1]))
                    break;
            if ((cx[3] != 0.0) ||
                (cx[7] != 0.0) ||
                (cx[11] != 0.0) ||
                (cx[15] != 1.0))
                c1 = 0;
            if (c1 > 15) {
                GET_QUATVARS2(cx)
                q2 = 1;
            }
        }
    }

    /* get data pointers */
    cx = (const double*) mxGetPr(*prhs);
    if (cx == NULL)
        mexErrMsgTxt("Error getting coordinate data pointer.");
    tx = (const double*) mxGetPr(prhs[1]);
    if (tx == NULL)
        mexErrMsgTxt("Error getting tr parameters data pointer.");
    cy = &cx[nc];
    cz = &cy[nc];
    ty = &tx[d123];
    tz = &ty[d123];

    /* check for coordinate range */
    if ((nc == 4) &&
         mxIsInf(*cx) &&
         mxIsInf(*cy) &&
         mxIsInf(*cz)) {

         /* get range spec */
         rf1 = cx[1];
         rf2 = cy[1];
         rf3 = cz[1];
         rs1 = cx[2];
         rs2 = cy[2];
         rs3 = cz[2];
         rp1 = cx[3];
         rp2 = cy[3];
         rp3 = cz[3];

         /* handle bad range */
         if (mxIsInf(rf1) ||
             mxIsNaN(rf1) ||
             mxIsInf(rf2) ||
             mxIsNaN(rf2) ||
             mxIsInf(rf3) ||
             mxIsNaN(rf3) ||
             mxIsInf(rs1) ||
             mxIsNaN(rs1) ||
             mxIsInf(rs2) ||
             mxIsNaN(rs2) ||
             mxIsInf(rs3) ||
             mxIsNaN(rs3) ||
             mxIsInf(rp1) ||
             mxIsNaN(rp1) ||
             mxIsInf(rp2) ||
             mxIsNaN(rp2) ||
             mxIsInf(rp3) ||
             mxIsNaN(rp3) ||
             (rs1 == 0.0) ||
             (rs2 == 0.0) ||
             (rs3 == 0.0)) {
                 *plhs = mxCreateNumericArray(2, mnout, mxDOUBLE_CLASS, mxREAL);
                 return;
         }

         /* get range correct */
         if (rs1 > 0)
             for (r1 = 0; rf1 <= rp1; ++r1) rf1 += rs1;
         else
             for (r1 = 0; rf1 >= rp1; ++r1) rf1 += rs1;
         if (rs2 > 0)
             for (r2 = 0; rf2 <= rp2; ++r2) rf2 += rs2;
         else
             for (r2 = 0; rf2 >= rp2; ++r2) rf2 += rs2;
         if (rs2 > 0)
             for (r3 = 0; rf3 <= rp3; ++r3) rf3 += rs3;
         else
             for (r3 = 0; rf3 >= rp3; ++r3) rf3 += rs3;

         /* handle empty range */
         if ((r1 == 0) ||
             (r2 == 0) ||
             (r3 == 0)) {
                 *plhs = mxCreateNumericArray(2, mnout, mxDOUBLE_CLASS, mxREAL);
                 return;
         }

         /* re-get start positions */
         rf1 = cx[1];
         rf2 = cy[1];
         rf3 = cz[1];

         /* enable range */
         cr = 1;

         /* number of coordinates */
         nc = r1 * r2 * r3;
    }

    /* create output matrices */
    *mnout = nc;
    *plhs = mxCreateNumericArray(2, mnout, mxDOUBLE_CLASS, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output argument. Out of memory.");
    ocx = (double *) mxGetData(*plhs);
    if (ocx == NULL)
        mexErrMsgTxt("Error getting data pointer.");
    ocy = &ocx[nc];
    ocz = &ocy[nc];

    /* compute helping parameters */
    ds1 = (double) (2 * s1);
    ds2 = (double) (2 * s2);
    ds3 = (double) (2 * s3);
    rsh1 = sqrt(2.0 / ((double) s1));
    rsh2 = sqrt(2.0 / ((double) s2));
    rsh3 = sqrt(2.0 / ((double) s3));

    /* general coordinates (not a range) */
    if (!cr) {

        /* create temporary vectors */
        bx = (double *) mxCalloc(d1, sizeof(double));
        by = (double *) mxCalloc(d2, sizeof(double));
        bz = (double *) mxCalloc(d3, sizeof(double));
        if ((bx == NULL) ||
            (by == NULL) ||
            (bz == NULL)) {
            if (bx != NULL)
                mxFree(bx);
            if (by != NULL)
                mxFree(by);
            if (bz != NULL)
                mxFree(bz);
            mexErrMsgTxt("Error allocating temporary vectors.");
        }

        /* set first weight fix! */
        *bx = 1.0 / sqrt((double) s1);
        *by = 1.0 / sqrt((double) s2);
        *bz = 1.0 / sqrt((double) s3);

        /* loop over coordinates */
        for (cc = nc; cc > 0; --cc) {

            /* get coordinate value */
            cvx = *cx++;
            cvy = *cy++;
            cvz = *cz++;

            /* initial transform? */
            if (q1) {
                QUAT_MULT(cvx, cvy, cvz);
            }

            /* compute DCTs for X/Y/Z */
            for (c1 = 1; c1 < d1; ++c1)
                bx[c1] = rsh1 * cos(TR_PI * (2.0 * cvx - 1.0) * ((double) c1) / ds1);
            for (c2 = 1; c2 < d2; ++c2)
                by[c2] = rsh2 * cos(TR_PI * (2.0 * cvy - 1.0) * ((double) c2) / ds2);
            for (c3 = 1; c3 < d3; ++c3)
                bz[c3] = rsh3 * cos(TR_PI * (2.0 * cvz - 1.0) * ((double) c3) / ds3);

            /* for X coordinate -> set transformation value to 0 and 3D multiply */
            for (tv = 0.0, tb = tx, c3 = 0, bzb = bz; c3 < d3; ++c3, ++bzb)
                for (c2 = 0, byb = by; c2 < d2; ++c2, ++byb)
                    for (c1 = 0, bxb = bx; c1 < d1; ++c1)
                        tv += *tb++ * *bxb++ * *byb * *bzb;
            cvx += tv;

            /* repeat for Y/Z */
            for (tv = 0.0, tb = ty, c3 = 0, bzb = bz; c3 < d3; ++c3, ++bzb)
                for (c2 = 0, byb = by; c2 < d2; ++c2, ++byb)
                    for (c1 = 0, bxb = bx; c1 < d1; ++c1)
                        tv += *tb++ * *bxb++ * *byb * *bzb;
            cvy += tv;
            for (tv = 0.0, tb = tz, c3 = 0, bzb = bz; c3 < d3; ++c3, ++bzb)
                for (c2 = 0, byb = by; c2 < d2; ++c2, ++byb)
                    for (c1 = 0, bxb = bx; c1 < d1; ++c1)
                        tv += *tb++ * *bxb++ * *byb * *bzb;
            cvz += tv;

            /* final transform? */
            if (q2) {
                QUAT_MULT2(cvx, cvy, cvz);
            }

            /* set in output */
            *ocx++ = cvx;
            *ocy++ = cvy;
            *ocz++ = cvz;
        }

    /* coordinate range */
    } else {

        /* DCTs cacheable? */
        if ((!q1) ||
            ((qt21 > -0.000000000001) &&
             (qt21 <  0.000000000001) &&
             (qt31 > -0.000000000001) &&
             (qt31 <  0.000000000001) &&
             (qt12 > -0.000000000001) &&
             (qt12 <  0.000000000001) &&
             (qt32 > -0.000000000001) &&
             (qt32 <  0.000000000001) &&
             (qt13 > -0.000000000001) &&
             (qt13 <  0.000000000001) &&
             (qt23 > -0.000000000001) &&
             (qt23 <  0.000000000001))) {

            /* apply transformation on range! */
            QUAT_MULT(rf1, rf2, rf3)
            rs1 *= qt11;
            rs2 *= qt22;
            rs3 *= qt33;

            /* create temporary vectors */
            bx = (double *) mxCalloc(d1 * r1, sizeof(double));
            by = (double *) mxCalloc(d2 * r2, sizeof(double));
            bz = (double *) mxCalloc(d3 * r3, sizeof(double));
            xx = (double *) mxCalloc(d1 * d2, sizeof(double));
            xy = (double *) mxCalloc(d1 * d2, sizeof(double));
            xz = (double *) mxCalloc(d1 * d2, sizeof(double));
            if ((bx == NULL) ||
                (by == NULL) ||
                (bz == NULL) ||
                (xx == NULL) ||
                (xy == NULL) ||
                (xz == NULL)) {
                if (bx != NULL)
                    mxFree(bx);
                if (by != NULL)
                    mxFree(by);
                if (bz != NULL)
                    mxFree(bz);
                if (xx != NULL)
                    mxFree(xx);
                if (xy != NULL)
                    mxFree(xy);
                if (xz != NULL)
                    mxFree(xz);
                mexErrMsgTxt("Error allocating temporary vectors.");
            }

            /* cache DCTs */
            for (bxb = bx, rp1 = rf1, c1 = 0; c1 < r1; ++c1, rp1 += rs1) {
                *bxb++ = 1.0 / sqrt((double) s1);
                for (c2 = 1; c2 < d1; ++c2)
                    *bxb++ = rsh1 * cos(TR_PI * (2.0 * rp1 - 1.0) * ((double) c2) / ds1);
            }
            for (byb = by, rp2 = rf2, c1 = 0; c1 < r2; ++c1, rp2 += rs2) {
                *byb++ = 1.0 / sqrt((double) s2);
                for (c2 = 1; c2 < d2; ++c2)
                    *byb++ = rsh2 * cos(TR_PI * (2.0 * rp2 - 1.0) * ((double) c2) / ds2);
            }
            for (bzb = bz, rp3 = rf3, c1 = 0; c1 < r3; ++c1, rp3 += rs3) {
                *bzb++ = 1.0 / sqrt((double) s3);
                for (c2 = 1; c2 < d3; ++c2)
                    *bzb++ = rsh3 * cos(TR_PI * (2.0 * rp3 - 1.0) * ((double) c2) / ds3);
            }

            /* loop over three dims */
            for (rp3 = rf3, s3 = 0; s3 < r3; ++s3, rp3 += rs3) {

                /* set bzb position */
                bzb = &bz[s3 * d3];

                /* cache temporary arrays for faster multiplication */
                for (tb = tx, xxb = xx, c2 = 0; c2 < d2; ++c2)
                    for (c1 = 0; c1 < d1; ++c1) {
                        for (tv = 0.0, bzbb = bzb, c3 = 0; c3 < d3; ++c3)
                            tv += tb[c1+c2*d1+c3*d12] * *bzbb++;
                        *xxb++ = tv;
                    }
                for (tb = ty, xyb = xy, c2 = 0; c2 < d2; ++c2)
                    for (c1 = 0; c1 < d1; ++c1) {
                        for (tv = 0.0, bzbb = bzb, c3 = 0; c3 < d3; ++c3)
                            tv += tb[c1+c2*d1+c3*d12] * *bzbb++;
                        *xyb++ = tv;
                    }
                for (tb = tz, xzb = xz, c2 = 0; c2 < d2; ++c2)
                    for (c1 = 0; c1 < d1; ++c1) {
                        for (tv = 0.0, bzbb = bzb, c3 = 0; c3 < d3; ++c3)
                            tv += tb[c1+c2*d1+c3*d12] * *bzbb++;
                        *xzb++ = tv;
                    }

                /* loop over plane dims */
                for (rp2 = rf2, s2 = 0; s2 < r2; ++s2, rp2 += rs2)
                    for (rp1 = rf1, s1 = 0; s1 < r1; ++s1, rp1 += rs1) {
                        cvx = rp1;
                        cvy = rp2;
                        cvz = rp3;

                        /* for each coordinate, set transformation value to 0 and 3D multiply */
                        for (tv = 0.0, xxb = xx, c2 = 0, byb = &by[s2 * d2]; c2 < d2; ++c2, ++byb)
                            for (c1 = 0, bxb = &bx[s1 * d1]; c1 < d1; ++c1)
                                tv += *xxb++ * *bxb++ * *byb;
                        cvx = rp1 + tv;
                        for (tv = 0.0, xyb = xy, c2 = 0, byb = &by[s2 * d2]; c2 < d2; ++c2, ++byb)
                            for (c1 = 0, bxb = &bx[s1 * d1]; c1 < d1; ++c1)
                                tv += *xyb++ * *bxb++ * *byb;
                        cvy = rp2 + tv;
                        for (tv = 0.0, xzb = xz, c2 = 0, byb = &by[s2 * d2]; c2 < d2; ++c2, ++byb)
                            for (c1 = 0, bxb = &bx[s1 * d1]; c1 < d1; ++c1)
                                tv += *xzb++ * *bxb++ * *byb;
                        cvz = rp3 + tv;

                        /* final transform? */
                        if (q2) {
                            QUAT_MULT2(cvx, cvy, cvz);
                        }

                        /* set in output */
                        *ocx++ = cvx;
                        *ocy++ = cvy;
                        *ocz++ = cvz;
                    }
            }

            /* free additional pointers */
            mxFree(xx);
            mxFree(xy);
            mxFree(xz);

        /* DCTs not cacheable */
        } else {

            /* create temporary vectors */
            bx = (double *) mxCalloc(d1, sizeof(double));
            by = (double *) mxCalloc(d2, sizeof(double));
            bz = (double *) mxCalloc(d3, sizeof(double));
            if ((bx == NULL) ||
                (by == NULL) ||
                (bz == NULL)) {
                if (bx != NULL)
                    mxFree(bx);
                if (by != NULL)
                    mxFree(by);
                if (bz != NULL)
                    mxFree(bz);
                mexErrMsgTxt("Error allocating temporary vectors.");
            }

            /* set first weight fix! */
            *bx = 1.0 / sqrt((double) s1);
            *by = 1.0 / sqrt((double) s2);
            *bz = 1.0 / sqrt((double) s3);

            /* loop over three dims */
            for (rp3 = rf3, s3 = 0; s3 < r3; ++s3, rp3 += rs3)
                for (rp2 = rf2, s2 = 0; s2 < r2; ++s2, rp2 += rs2)
                    for (rp1 = rf1, s1 = 0; s1 < r1; ++s1, rp1 += rs1) {
                        cvx = rp1;
                        cvy = rp2;
                        cvz = rp3;

                        /* initial transform? */
                        if (q1) {
                            QUAT_MULT(cvx, cvy, cvz);
                        }

                        /* compute DCTs for X/Y/Z */
                        for (c1 = 1; c1 < d1; ++c1)
                            bx[c1] = rsh1 * cos(TR_PI * (2.0 * cvx - 1.0) * ((double) c1) / ds1);
                        for (c2 = 1; c2 < d2; ++c2)
                            by[c2] = rsh2 * cos(TR_PI * (2.0 * cvy - 1.0) * ((double) c2) / ds2);
                        for (c3 = 1; c3 < d3; ++c3)
                            bz[c3] = rsh3 * cos(TR_PI * (2.0 * cvz - 1.0) * ((double) c3) / ds3);

                        /* for X coordinate -> set transformation value to 0 and 3D multiply */
                        for (tv = 0.0, tb = tx, c3 = 0, bzb = bz; c3 < d3; ++c3, ++bzb)
                            for (c2 = 0, byb = by; c2 < d2; ++c2, ++byb)
                                for (c1 = 0, bxb = bx; c1 < d1; ++c1)
                                    tv += *tb++ * *bxb++ * *byb * *bzb;
                        cvx += tv;

                        /* repeat for Y/Z */
                        for (tv = 0.0, tb = ty, c3 = 0, bzb = bz; c3 < d3; ++c3, ++bzb)
                            for (c2 = 0, byb = by; c2 < d2; ++c2, ++byb)
                                for (c1 = 0, bxb = bx; c1 < d1; ++c1)
                                    tv += *tb++ * *bxb++ * *byb * *bzb;
                        cvy += tv;
                        for (tv = 0.0, tb = tz, c3 = 0, bzb = bz; c3 < d3; ++c3, ++bzb)
                            for (c2 = 0, byb = by; c2 < d2; ++c2, ++byb)
                                for (c1 = 0, bxb = bx; c1 < d1; ++c1)
                                    tv += *tb++ * *bxb++ * *byb * *bzb;
                        cvz += tv;

                        /* final transform? */
                        if (q2) {
                            QUAT_MULT2(cvx, cvy, cvz);
                        }

                        /* set in output */
                        *ocx++ = cvx;
                        *ocy++ = cvy;
                        *ocz++ = cvz;
                    }
        }
    }

    /* free temporary space */
    mxFree(bx);
    mxFree(by);
    mxFree(bz);
}

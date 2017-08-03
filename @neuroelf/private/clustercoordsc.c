/*

clustercoordsc see clustercoods.m

FORMAT:       [cs, cv, l, c] = clustercoordsc(v [, m [, t]]);

Input fields:

      v           binary input volume
      m           method, default: 2
                  1 face-connectivity (3D)
                  2 edge-connectivity (3D)
                  3 vertex-connectivity (3D)
                  4 face-connectivity (XY-slice)
                  5 edge-connectivity (XY-slice)
      t           cluster-size threshold, default: 1

Output fields:
      cs          list of cluster sizes
      cv          clustered volume (uint32 with 0's for false voxels
                  and 1...n for clustered voxels)
      l           Vx4 list of cluster voxels [x(:), y(:), z(:), n(:)]
      c           Cx1 cell array with lists of coordinates

% Version:  v0.9d
% Build:    14062015
% Date:     Jun-20 2014, 3:53 PM EST
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
#include <stdio.h>

#define CREATE_MxN(SNAME, MDIM, NDIM) \
    SNAME = mxCreateNumericMatrix(MDIM, NDIM, mxDOUBLE_CLASS, mxREAL);
#define CREATE_SCALAR(SNAME, SVAL) { \
    SNAME = mxCreateNumericArray(2, scalardim, mxDOUBLE_CLASS, mxREAL); \
    if (SNAME == NULL) \
        mexErrMsgTxt("Error allocating memory."); \
    *((double *) mxGetPr(SNAME)) = (double) SVAL; \
}

/* functions used for copying */
void ff_copy_intopart(unsigned int x1, unsigned int y1, unsigned int z1, const unsigned char *sd, const int *tnd, unsigned char *td, const unsigned int *off);
void ff_copy_frompart(unsigned int x1, unsigned int y1, const unsigned int *sd, const int *tnd, unsigned int *td, const unsigned int *off);

/* up to 13 directions (and their reverse) */
static const signed int dlx[13] = { 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1 };
static const signed int dly[13] = { 0, 1, 0, 1,-1, 0, 0, 1, 1, 1, 1,-1,-1 };
static const signed int dlz[13] = { 0, 0, 1, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1 };

/* scalar dims */
static const int scalardim[2] = {1, 1};

/* here comes the main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* default settings */
    /* connectivity pattern: 2 */
    const double DEF_clconn = 2.0;

    /* threshold 1 (consider also single voxel clusters) */
    const double DEF_kthr   = 1.0;

    /* setting pointers (for argument read-out) */
    const double *clconn = NULL;
    const double *kthr   = NULL;

    /* settings (internal) */
    /* max directions */
    unsigned int md = 2;

    /* cluster threshold */
    unsigned int kthrl = 1;



    /* sizes and position */
    const int *dim = NULL;
    int idim[3] = {1, 1, 1};
    int tdim[3] = {1, 1, 1};

    /* offset (must be variable, as copy might not require additional slices) */
    unsigned int toff[3] = {1, 1, 1};

    /* dimX, dimY, dimZ, dimXY, next position */
    unsigned int dx = 0, dy = 0, dz = 0, dxy = 0, dxyz = 0, np = 0;

    /* data pointers */
    /* (enlarged) copy */
    unsigned char *copy = NULL;

    /* list of cluster voxels */
    unsigned int clvoxcnt;
    unsigned int clnxtvox;
    signed int *clvox = NULL;

    /* list of cluster sizes and total number of clustered voxels */
    unsigned int mxclnum = 0;
    unsigned int clnum = 0;
    unsigned int *clsiz = NULL;
    unsigned int cltvox = 0, cltvox2 = 0, cltvox3 = 0;

    /* clustered volume (internal) */
    unsigned int *clvol = NULL;

    /* direction increments (as index differences) */
    signed int dlxyz[26] =
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    /* next positions (signed!) */
    signed int nextp = 0, nextpi = 0;

    /* counters */
    signed char dc = 0;
    unsigned int cx = 0, cy = 0, cz = 0;

    /* output pointer(s) */
    mxArray *ma = NULL;
    double *op = NULL, *sop = NULL;

    /* variable output text */
    /* char vstr[256]; */



    /* code starts here */



    /* nr. of arguments check */
    if ((nrhs < 1) ||
        (nlhs > 4))
        mexErrMsgTxt("Bad number of input/output arguments.");

    /* first argument class and size check */
    if (!mxIsLogical(prhs[0]) ||
        (mxGetNumberOfElements(prhs[0]) < 8) ||
        (mxGetNumberOfDimensions(prhs[0]) > 3))
        mexErrMsgTxt("Invalid argument type or dimension.");

    /* check connectivity flag */
    if ((nrhs < 2) ||
        !mxIsDouble(prhs[1]) ||
        (mxGetNumberOfElements(prhs[1]) != 1))
	clconn = &DEF_clconn;
    else {
	clconn = (const double*) mxGetPr(prhs[1]);
	if (mxIsNaN(*clconn) ||
	    mxIsInf(*clconn) ||
	   ((*clconn) < 1.0) ||
	   ((*clconn) > 5.0))
	    clconn = &DEF_clconn;
    }

    /* check threshold */
    if ((nrhs < 3) ||
        !mxIsDouble(prhs[2]) ||
        (mxGetNumberOfElements(prhs[2]) != 1))
        kthr = &DEF_kthr;
    else {
        kthr = (const double*) mxGetPr(prhs[2]);
        if (mxIsNaN(*kthr) ||
            mxIsInf(*kthr) ||
           ((*kthr) < 1.0) ||
           ((*kthr) > (0.5 * mxGetNumberOfElements(prhs[0]))))
            kthr = &DEF_kthr;
    }
    kthrl = (unsigned int) *kthr;

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
    dxyz = dxy * dz;
    mxclnum = (unsigned int) (((double) (dxyz + 3)) / 4.0);

    /* set connectivity flag -> max directions (opposite direction later) */
    /* face connectivity -> three main directions */
    if ((*clconn >= 1.0) && (*clconn < 1.5))
        md = 3;

    /* edge connectivity -> add 6 edge diagonals */
    else if ((*clconn >= 1.5) && (*clconn < 2.5))
        md = 9;

    /* vertex connectivity -> add 4 3D diagonals */
    else if ((*clconn >= 2.5) && (*clconn < 3.5))
        md = 13;
    /* otherwise, leave at 2, we stay "in-plane" */

    /* reserve space for next cluster list and cluster sizes */
    clvox = (signed int *) mxCalloc(dxyz, sizeof(signed int));
    if (clvox == NULL)
	mexErrMsgTxt("Error allocating next cluster voxel list.");
    clsiz = (unsigned int *) mxCalloc(mxclnum, sizeof(signed int));
    if (clsiz == NULL) {
        mxFree(clvox);
	mexErrMsgTxt("Error allocating next cluster voxel list.");
    }

    /* increase numbers accordingly */
    if (idim[0] > 1)
        dx += 2;
    else
        toff[0] = 0;
    if (idim[1] > 1)
        dy += 2;
    else
        toff[1] = 0;
    if (idim[2] > 1)
        dz += 2;
    else
        toff[2] = 0;

    /* compute new plane/vol numel */
    dxy = dx * dy;
    dxyz = dxy * dz;

    /* reserve space for clustered volume (if needed) */
    if (nlhs > 1) {
        clvol = (unsigned int *) mxCalloc(dxyz, sizeof(unsigned int));
        if (clvol == NULL) {
            mxFree(clsiz);
            mxFree(clvox);
            mexErrMsgTxt("Error allocating clustered volume.");
        }
    }

    /* set temp dim array */
    tdim[0] = dx;
    tdim[1] = dy;
    tdim[2] = dz;

    /* fill direction increment array */
    np = 0;
    for (dc = 0; dc < md; ++dc) {
        if (((idim[0] > 1) ||
             (dlx[dc] == 0)) &&
            ((idim[1] > 1) ||
             (dly[dc] == 0)) &&
            (((idim[2] > 1) &&
	      (*clconn < 3.5)) ||
             (dlz[dc] == 0))) {
            dlxyz[np] = dlx[dc] + dly[dc] * dx + dlz[dc] * dxy;
            dlxyz[np+1] = -dlxyz[np];
            np += 2;
        }
    }
    if ((*clconn >= 4.5) &&
        (idim[0] > 1) &&
        (idim[1] > 1)) {
        for (dc = 3; dc < 5; ++dc) {
            dlxyz[np] = dlx[dc] + dly[dc] * dx + dlz[dc] * dxy;
            dlxyz[np+1] = -dlxyz[np];
            np += 2;
        }
    }

    /* re-set */
    md = np - 1;

    /* create temp. copy array */
    copy = (unsigned char *) mxCalloc(dxyz, sizeof(unsigned char));
    if (copy == NULL) {
        mxFree(clsiz);
        mxFree(clvox);
        mexErrMsgTxt("Error allocating temporary memory");
    }

    /* copy into temporary array */
    ff_copy_intopart(idim[0], idim[1], idim[2], (unsigned char *) mxGetPr(prhs[0]), tdim, copy, toff);



    /* preparations complete, now we can "search"... */



    /* iterate while new points found */
    for (np = 0; np < dxyz; ++np) {

	/* check voxel */
        if (copy[np] == 1) {

            /* start new cluster */
            clvoxcnt = 1;
            clvox[0] = np;

            /* remove from search volume */
            copy[np] = 0;

            /* continue until no more neighbors found */
            for (clnxtvox = 0; clnxtvox < clvoxcnt; ++clnxtvox) {

                /* check neighbors */
                nextp = clvox[clnxtvox];
                for (dc = md;  dc >= 0; --dc) {

                    /* check voxel with increment */
                    nextpi = nextp + dlxyz[dc];
                    if (copy[nextpi] == 1) {

                            /* equally remove from search volume */
                            copy[nextpi] = 0;

                            /* and add to list */
                            clvox[clvoxcnt++] = nextpi;
                    }
                }
            }

            /* rest of code only if surpasses threshold */
            if (clvoxcnt >= kthrl) {

                /* write to list of sizes and increase cluster number count */
                clsiz[clnum++] = clvoxcnt;
                cltvox += clvoxcnt;

                /* write to output volume (if needed) */
                if (nlhs > 1) {
                    for (nextp = clvoxcnt - 1; nextp >= 0; --nextp)
                        clvol[clvox[nextp]] = clnum;
                }
            }
        }
    }

    /* we're done with this */
    mxFree(copy);

    /* create output argument #1 */
    CREATE_MxN(plhs[0], 1, clnum);
    if (plhs[0] == NULL) {
        if (clvol != NULL)
            mxFree(clvol);
        mxFree(clsiz);
        mxFree(clvox);
        mexErrMsgTxt("Error creating output argument cs.");
    }

    /* if requested create other outputs first */
    if (nlhs > 1) {
        plhs[1] = mxCreateNumericArray(3, idim, mxUINT32_CLASS, mxREAL);
        if (plhs[1] == NULL) {
            mxDestroyArray(plhs[0]);
            mxFree(clvol);
            mxFree(clsiz);
            mxFree(clvox);
            mexErrMsgTxt("Error creating output argument cv.");
        }
    }

    /* if requested create other outputs first */
    if (nlhs > 2) {
        CREATE_MxN(plhs[2], cltvox, 4);
        if (plhs[2] == NULL) {
            mxDestroyArray(plhs[1]);
            mxDestroyArray(plhs[0]);
            mxFree(clvol);
            mxFree(clsiz);
            mxFree(clvox);
            mexErrMsgTxt("Error creating output argument l.");
        }
    }

    /* if requested create other outputs first */
    if (nlhs > 3) {
        plhs[3] = mxCreateCellMatrix(1, clnum);
        if (plhs[3] == NULL) {
            mxDestroyArray(plhs[2]);
            mxDestroyArray(plhs[1]);
            mxDestroyArray(plhs[0]);
            mxFree(clvol);
            mxFree(clsiz);
            mxFree(clvox);
            mexErrMsgTxt("Error creating output argument c.");
        }
        /* put numeric array into each cell */
        for (nextp = 0; nextp < clnum; ++nextp) {
            CREATE_MxN(ma, clsiz[nextp], 3);
            if (ma == NULL) {
                mxDestroyArray(plhs[3]);
                mxDestroyArray(plhs[2]);
                mxDestroyArray(plhs[1]);
                mxDestroyArray(plhs[0]);
                mxFree(clvol);
                mxFree(clsiz);
                mxFree(clvox);
                mexErrMsgTxt("Error creating output cell contents in argument c.");
            }
            mxSetCell(plhs[3], nextp, ma);
        }
    }

    /* copy cluster sizes into output */
    op = (double *) mxGetPr(plhs[0]);
    for (nextp = 0; nextp < clnum; ++nextp)
        *op++ = (double) clsiz[nextp];

    /* copy volume to output */
    if (nlhs > 1) {
        ff_copy_frompart(dx, dy, clvol, idim, (unsigned int*) mxGetData(plhs[1]), toff);
        mxFree(clvol);
        clvol = (unsigned int *) mxGetPr(plhs[1]);
    }

    /* create list of voxels */
    if (nlhs > 2) {

        /* get output pointer */
        op = (double *) mxGetPr(plhs[2]);

        /* 2/3 * total count for fast indexing */
        cltvox2 = 2 * cltvox;
        cltvox3 = 3 * cltvox;

        /* if required re-init counter */
        if (nlhs > 3) {
            for (nextp = clnum - 1; nextp >= 0; --nextp)
                clvox[nextp] = 0;
        }

        dx = idim[0];
        dy = idim[1];
        dz = idim[2];
        dxy = dx * dy;
        for (cz = 0; cz < dz; ++cz) {
            dxyz = dxy * cz;
            for (cy = 0; cy < dy; ++cy) {
                np = dxyz + cy * dx;
                for (cx = 0; cx < dx; ++cx) {
                    clnum = *clvol++;
                    if (clnum > 0) {
                        op[cltvox3] = (double) clnum;
                        op[cltvox2] = (double) (cz + 1);
                        op[cltvox]  = (double) (cy + 1);
                        *op++ = (double) (cx + 1);
                        if (nlhs > 3) {
                            --clnum;
                            sop = (double *) mxGetPr(mxGetCell(plhs[3], clnum));
                            sop[2 * clsiz[clnum] + clvox[clnum]] = (double) (cz + 1);
                            sop[clsiz[clnum] + clvox[clnum]] = (double) (cy + 1);
                            sop[clvox[clnum]++] = (double) (cx + 1);
                        }
                    }
                }
            }
        }
    }

    /* free other temp arrays */
    mxFree(clsiz);
    mxFree(clvox);

}

void ff_copy_intopart(unsigned int x1, unsigned int y1, unsigned int z1, const unsigned char *sd, const int *tnd, unsigned char *td, const unsigned int *off)
{
    unsigned int xc, yc, zc, x2, y2, z2, xy, to, yo, zo;
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

void ff_copy_frompart(unsigned int x1, unsigned int y1, const unsigned int *sd, const int *tnd, unsigned int *td, const unsigned int *off)
{
    unsigned int rv, xc, yc, zc, xy, x2, y2, z2, to, yo, zo;
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
            for (xc = 0; xc < x2; ++xc)
                *td++ = sd[++yo];
        }
    }
}

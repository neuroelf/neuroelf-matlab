/*

morphing of mesh coordinates with neighbor information

FORMAT:       [c, d] = mesh_morph(c, n, tri, opts)

Input fields:

      c           Cx3 coordinate list
      n           Cx2 cell array of SRF xff object (1-based !)
      tri         Tx3 triangle list (1-based !)
      opts        mandatory struct with settings
       .force     1x1 double smoothing force
       .niter     1x1 double number of iterations
                - optionally provided settings
       .areac     if 1x1 double := 1, keep area constant
                  (from initial state, requires .tri to be set!)
       .areaw     if given and [1], perform smoothing with area-based
                  weighting (default: false)
       .distc     if given and between 0 .. 1, perform distortion corr
       .distw     if given and [1], perform smoothing with distance-based
                  weighting (default: false)
       .equalw    apply equal weighting (assuming 6 neighbors, false)
       .norm      force along normal vector (default: 0)
       .normramp  if 1x1 double := 1, ramp along-normal force: 0 ... .norm 
       .sphere    to-sphere force
       .type      1xN char type, currently only 'smooth' supported

Output fields:

      c           morphed coordinates
      d           map values giving the triangle size/vertex after
                  morphing (optional, requires areac/tri to be set!) with
                  formula log(sqrt(TriangleSizeAtVertex / MeanTriangleSize));

Note: for smoothing the neighbor information is used, for distortion
      correction the triangle data is used!

% Version:  v0.9d
% Build:    14070113
% Date:     Jul-01 2014, 1:37 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2010, 2011, 2014, Jochen Weber
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

#define MESH_PI 3.14159265358979323846
#define O_AREAMEAN             0
#define O_AREATOTAL            1
#define O_CENTEROFGRAVITY      2
#define O_DENSITYSMPDATA       3
#define O_NORMALVECTORS        4
#define O_ORIGAREAMEAN         5
#define O_ORIGAREATOTAL        6
#define O_ORIGDENSITYSMPDATA   7
#define O_RADIUSMEAN           8
#define O_XNUMBEROFFIELDS      9

double subFunc_area(int numc, int numt, const double *c, const unsigned long *t);
double subFunc_area_a(int numc, int numt, const double *c, const unsigned long *t, double *ca);
double subFunc_area_c(int numc, int numt, const double *c, const unsigned long *t, double *ca, const unsigned char *cad);
void subFunc_center(int numc, const double *c, double *cx, double *cy, double *cz);
void subFunc_normals(int numc, int numt, const double *c, const unsigned long *t, double *nrm);
double subFunc_radius(int numc, const double *c, double *cr);
void subFunc_vertexnrareas(int numc, int numt, const unsigned long *t, unsigned char *cad);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* number of coordinates (including double + triple), counter */
    int nc, nc2, nc3;
    int cc = 0;
    double dnc;

    /* number of triangles (including triple) */
    int tc, tc3;

    /* coordinate pointer (and backup pointer) */
    double *c, *c2, *c3, *bc, *bc2, *bc3;

    /* coordinate backup pointer */
    double *cb, *cb2, *cb3, *bcb, *bcb2, *bcb3;

    /* neighbors/triangle list pointer */
    unsigned long *n = NULL;
    unsigned long *tn = NULL;

    /* number of neighbors list pointer */
    unsigned char *nn = NULL;

    /* triangle list */
    const double *tria = NULL;
    unsigned long *tri = NULL;

    /* vertex area list, original vertex area list and number of areas per vertex list */
    double *ca = NULL;
    double *oca = NULL;
    unsigned char *cad = NULL;

    /* smoothing force (with default value) */
    double force = 0.07, normforce = 0.0, rampforce = 0.0, rampfrom = 0.0, selfforce;

    /* distortion correction values */
    double distc = 0.5;
    bool distcs = 0;

    /* distance weighting */
    bool areaw = 0, distw = 0, distwl = 0, distwsq = 0, equalw = 0, norm = 0, normramp = 0;

    /* neighbor area, mean area, smoothing values */
    double na[36];
    double *bna;
    double ma, dsv, cav;

    /* to-sphere force and radius values */
    double tosph = 0.0, e1x, e1y, e1z, e2x, e2y, e2z, nsx, nsy, nsz, nx, ny, nz, vlength;
    bool sphere = 0;

    /* mean radius and radius correction, radius list */
    double mr;
    double *cr, *bcr;

    /* number of iterations */
    double niterd = 1.0;
    int niter = 1;
    int iterc;

    /* morphtype (from options) */
    char mtypec[11];
    int mtype = 0;

    /* number of dimensions and dim size */
    int nd;
    const int *di;

    /* field number, field array and array pointer */
    int fn;
    mxArray *faout;
    const mxArray *fa;
    double *fapout;
    const double *fap;

    /* coordinate and smoothing values, number of neighbors (and counter)<s */
    double cvx, cvy, cvz, svx, svy, svz, dnn;
    unsigned char cnn, cnc;

    /* number of neighbor list elements, element counter, neighbor list offset */
    int nel, ec, no;

    /* area constantness */
    bool areac = 0;
    double areav, nareav, tarea, ocx = 0.0, ocy = 0.0, ocz = 0.0, ncx, ncy, ncz, scale;

    /* double triangle value and sides */
    double tv;

    /* compute output flag */
    bool compdm = 0;
    int compds[2] = {1, 1};
    int compdmnumfields = O_XNUMBEROFFIELDS;
    const char *compdmfields[] = {
        "AreaMean",
        "AreaTotal",
        "CenterOfGravity",
        "DensitySMPData",
        "NormalVectors",
        "OrigAreaMean",
        "OrigAreaTotal",
        "OrigDensitySMPData",
        "RadiusMean"
    };

    /* variable output string */
    char vstr[256]; /* debugging */

    /* check number of in/out arguments */
	if ((nrhs != 4) || (nlhs > 2))
		mexErrMsgTxt("Bad number of input/output arguments.");

    /* density map requested -> requires triangles */
    if (nlhs > 1)
        compdm = 1;

    /* check input argument types */
    if (!mxIsDouble(prhs[0]) || !mxIsCell(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsStruct(prhs[3]))
        mexErrMsgTxt("Bad argument type. Must be (double, cell, struct).");

    /* check coordinates argument ndims */
	nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd != 2)
        mexErrMsgTxt("Coordinates must be 2-D array.");

    /* check coordinates argument size (and get number of coordinates from this) */
	di = mxGetDimensions(prhs[0]);
    nc = di[0];
    nc2 = nc * 2;
    nc3 = nc * 3;
    dnc = (double) nc;
    if (di[1] != 3)
        mexErrMsgTxt("Coordinates must be Cx3 matrix.");

    /* check neighbors argument ndims */
    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd != 2)
        mexErrMsgTxt("Neighbors must be 2-D array");

    /* check neighbors argument size */
    di = mxGetDimensions(prhs[1]);
    if (di[0] != nc || di[1] != 2)
        mexErrMsgTxt("Neighbors must be a Cx2 matrix.");

    /* check triangles argument ndims */
    nd = mxGetNumberOfDimensions(prhs[2]);
    if (nd != 2)
        mexErrMsgTxt("Triangles must be 2-D array");

    /* check triangles argument size */
    di = mxGetDimensions(prhs[2]);
    if (di[1] != 3)
        mexErrMsgTxt("Triangles must be a Tx3 matrix.");
    tc = di[0];
    tc3 = tc * 3;
    tria = (double *) mxGetPr(prhs[2]);

    /* check options argument ndims */
    nd = mxGetNumberOfDimensions(prhs[3]);
    if (nd != 2)
        mexErrMsgTxt("Options must be a 2-D struct array.");

    /* check options argument size */
    if (mxGetNumberOfElements(prhs[3]) != 1)
        mexErrMsgTxt("Options must be a 1x1 struct array.");

    /* check options argument force field */
    fn = mxGetFieldNumber(prhs[3], "force");
    if (fn < 0)
        mexErrMsgTxt("Options does not contain force information.");

    /* check force field type */
    fa = mxGetFieldByNumber(prhs[3], 0, fn);
    if (!mxIsDouble(fa) || (mxGetNumberOfElements(fa) != 1))
        mexErrMsgTxt("Options force field must be a 1x1 double.");

    /* get and check force value */
    force = *((double *) mxGetPr(fa));
    if (mxIsInf(force) || mxIsNaN(force) || force < -16.0 || force > 16.0 || force == 0)
        mexErrMsgTxt("Options force field contains invalid force value.");
    selfforce = 1.0 - force;

    /* check options argument niter field */
    fn = mxGetFieldNumber(prhs[3], "niter");
    if (fn < 0)
        mexErrMsgTxt("Options does not contain niter information.");

    /* check niter field type */
    fa = mxGetFieldByNumber(prhs[3], 0, fn);
    if (!mxIsDouble(fa) || (mxGetNumberOfElements(fa) != 1))
        mexErrMsgTxt("Options niter field must be a 1x1 double.");

    /* get and check niter value */
    niterd = *((double *) mxGetPr(fa));
    if (mxIsInf(niterd) || mxIsNaN(niterd) || niterd < 0 || niterd > 50000)
        mexErrMsgTxt("Options niter field contains invalid niter value.");
    niter = (int) niterd;

    /* check options argument type field */
    fn = mxGetFieldNumber(prhs[3], "type");
    if (fn >= 0) {
        mtype = -1;
        /* check type field type */
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        nd = mxGetNumberOfDimensions(fa);
        di = mxGetDimensions(fa);
        if (!mxIsChar(fa) || (nd != 2) || (di[0] != 1) || (di[1] < 1) || (di[1] > 10))
            mexErrMsgTxt("Options type field must be a 1xC char array.");

        /* get and check type */
        mxGetString(fa, mtypec, 10);
        switch (mtypec[0]) {
            case 's':
            case 'S':
                mtype = 0;
                break;
        }
    }
    if (mtype == -1)
        mexErrMsgTxt("Options type field contains invalid value.");

    /* check for distance weighting */
    fn = mxGetFieldNumber(prhs[3], "areac");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            fap = (double *) mxGetPr(fa);
            if (fap[0] == 1)
                areac = 1;
        }
    }

    /* check for area-based weighting */
    fn = mxGetFieldNumber(prhs[3], "areaw");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            distc = *((double *) mxGetPr(fa));
            if (!mxIsInf(distc) && !mxIsNaN(distc) && distc == 1.0)
                areaw = 1;
        }
    }

    /* check for normal-force ramping */
    fn = mxGetFieldNumber(prhs[3], "normramp");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            distc = *((double *) mxGetPr(fa));
            if (!mxIsInf(distc) && !mxIsNaN(distc) && distc == 1.0)
                normramp = 1;
        }
    }

    /* check for distance-based weighting */
    fn = mxGetFieldNumber(prhs[3], "distw");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            distc = *((double *) mxGetPr(fa));
            if (!mxIsInf(distc) && !mxIsNaN(distc) && distc == 1.0)
                distw = 1;
        }
    }
    fn = mxGetFieldNumber(prhs[3], "distwl");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            distc = *((double *) mxGetPr(fa));
            if (!mxIsInf(distc) && !mxIsNaN(distc) && distc == 1.0) {
                distwl = 1;
                if (distwl)
                    distw = 1;
            }
        }
    }

    fn = mxGetFieldNumber(prhs[3], "distwsq");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            distc = *((double *) mxGetPr(fa));
            if (!mxIsInf(distc) && !mxIsNaN(distc) && distc == 1.0) {
                distwsq = 1;
                if (distwsq)
                    distw = 1;
            }
        }
    }


    /* check for equal-based weighting */
    fn = mxGetFieldNumber(prhs[3], "equalw");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            distc = *((double *) mxGetPr(fa));
            if (!mxIsInf(distc) && !mxIsNaN(distc) && distc == 1.0)
                equalw = 1;
        }
    }

    /* check for distortion correction */
    fn = mxGetFieldNumber(prhs[3], "distc");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            distc = *((double *) mxGetPr(fa));
            if (!mxIsInf(distc) && !mxIsNaN(distc) && distc > 0.0 && distc <= 6.0)
                distcs = 1;
        }
    }

    /* apply along-normal force */
    fn = mxGetFieldNumber(prhs[3], "rampfrom");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            normforce = *((double *) mxGetPr(fa));
            if (!mxIsInf(normforce) && !mxIsNaN(normforce) && normforce >= -1.0 && normforce <= 1.0 && normforce != 0.0)
                rampfrom = normforce;
            normforce = 0.0;
        }
    }
    fn = mxGetFieldNumber(prhs[3], "norm");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            normforce = *((double *) mxGetPr(fa));
            if (!mxIsInf(normforce) && !mxIsNaN(normforce) && normforce >= -1.0 && normforce <= 1.0 && normforce != 0.0)
                norm = 1;
        }
    }

    /* apply to-sphere force */
    fn = mxGetFieldNumber(prhs[3], "sphere");
    if (fn >= 0) {
        fa = mxGetFieldByNumber(prhs[3], 0, fn);
        if (mxIsDouble(fa) && (mxGetNumberOfElements(fa) == 1)) {
            fap = (double *) mxGetPr(fa);
            if (!mxIsInf(*fap) && !mxIsNaN(*fap) && *fap > 0.0 && *fap <= 1.0) {
                tosph = *fap;
                sphere = 1;
            }
        }
    }

    /* create coordinates backup array */
    cb = (double *) mxCalloc(3 * nc, sizeof(double));
    if (cb == NULL)
        mexErrMsgTxt("Error allocating coordinate backup buffer.");

    /* prepare neighbors list */
    n = (unsigned long*) mxCalloc(36 * nc, sizeof(unsigned long));
    if (n == NULL)
        mexErrMsgTxt("Error allocating neighbors list in memory.");
    nn = (unsigned char*) mxCalloc(nc, sizeof(unsigned char));
    if (nn == NULL)
        mexErrMsgTxt("Error allocating number of neighbors list in memory.");

    /* create triangle lookup array */
    tri = (unsigned long *) mxCalloc(tc3, sizeof(unsigned long));
    if (tri == NULL)
        mexErrMsgTxt("Unable to allocate triangle list memory.");

    /* parse triangle data */
    if (tria == NULL)
        mexErrMsgTxt("Triangle data required but not given.");
    for (cc = 0; cc < tc3; ++cc) {
        tv = tria[cc];
        if (mxIsInf(tv) || mxIsNaN(tv) || tv < 1.0 || tv > dnc)
            mexErrMsgTxt("Invalid tri option (each triangle vertex must be between 1 .. C).");
        tri[cc] = ((unsigned long) tv) - 1;
    }

    /* parse neighbors into list */
    for (cc = 0; cc < nc; ++cc) {
        fa = mxGetCell(prhs[1], nc + cc);
        if ((fa == NULL) || !mxIsDouble(fa))
            mexErrMsgTxt("Invalid neighbors argument (empty neighbor).");
        nel = mxGetNumberOfElements(fa);
        if (nel > 36)
            mexErrMsgTxt("Invalid neighbors argument (vertex must not have more than 36 neighbors).");
        nn[cc] = (unsigned char) nel;
        fap = (double *) mxGetPr(fa);
        no = 36 * cc;
        for (ec = 0; ec < nel; ++ec) {
            tv = fap[ec];
            if (mxIsInf(tv) || mxIsInf(tv) || tv < 1.0 || tv > dnc)
                mexErrMsgTxt("Invalid neighbors argument (each neighbor must be between 1 .. C).");
            n[no + ec] = ((unsigned long) tv) - 1;
        }
    }

    /* create output */
    di = mxGetDimensions(prhs[0]);
    plhs[0] = mxCreateNumericArray(2, di, mxDOUBLE_CLASS, mxREAL);
    c = (double *) mxGetPr(plhs[0]);

    /* create second argument as well */
    if (compdm) {
        plhs[1] = mxCreateStructArray(2, compds, compdmnumfields, compdmfields);
        if (plhs[1] == NULL)
            mexErrMsgTxt("Error creating second output argument.");
    }
    /* copy coordinates from source */
    fap = (double *) mxGetPr(prhs[0]);
    bc = c;
    for (cc = 0; cc < nc3; ++cc)
        *c++ = *fap++;
    c = bc;

    /* get original area per vertex and total area */
    ca = (double *) mxCalloc(nc, sizeof(double));
    if (ca == NULL)
        mexErrMsgTxt("Error allocating memory for coordinate area array.");
    cad = (unsigned char *) mxCalloc(nc, sizeof(unsigned char));
    if (cad == NULL)
        mexErrMsgTxt("Error allocating memory for number of areas per coordinate array.");
    cr = (double *) mxCalloc(nc, sizeof(double));
    if (cr == NULL)
        mexErrMsgTxt("Error allocating memory for coordinate area array.");
    oca = (double *) mxCalloc(nc, sizeof(double));
    if (oca == NULL)
        mexErrMsgTxt("Error allocating memory for coordinate area array.");
    subFunc_vertexnrareas(nc, tc, tri, cad);
    areav = subFunc_area(nc, tc, c, tri);
    ma = areav / dnc;
    if (compdm) {
        compds[0] = 1;
        compds[1] = 1;
        faout = mxCreateNumericArray(2, compds, mxDOUBLE_CLASS, mxREAL);
        if (faout == NULL)
            mexErrMsgTxt("Error creating AreaMean field in memory.");
        fapout = (double *) mxGetPr(faout);
        if (fapout == NULL)
            mexErrMsgTxt("Error getting double pointer to AreaMean field.");
        *fapout = ma;
        mxSetFieldByNumber(plhs[1], 0, O_AREAMEAN, faout);
        faout = mxCreateNumericArray(2, compds, mxDOUBLE_CLASS, mxREAL);
        if (faout == NULL)
            mexErrMsgTxt("Error creating AreaTotal field in memory.");
        fapout = (double *) mxGetPr(faout);
        if (fapout == NULL)
            mexErrMsgTxt("Error getting double pointer to AreaTotal field.");
        *fapout = areav;
        mxSetFieldByNumber(plhs[1], 0, O_AREATOTAL, faout);
    }

    /* make backup of pointers */
    bc = c;
    bc2 = &bc[nc];
    bc3 = &bc[nc2];
    bcb = cb;
    bcb2 = &bcb[nc];
    bcb3 = &bcb[nc2];
    bcr = cr;

    /* get mean */
    if (areac || sphere || norm) {

        /* get and subtract center */
        subFunc_center(nc, bc, &ocx, &ocy, &ocz);
        if (compdm) {
            compds[0] = 1;
            compds[1] = 3;
            faout = mxCreateNumericArray(2, compds, mxDOUBLE_CLASS, mxREAL);
            if (faout == NULL)
                mexErrMsgTxt("Error creating CenterOfGravity field in memory.");
            fapout = (double *) mxGetPr(faout);
            if (fapout == NULL)
                mexErrMsgTxt("Error getting double pointer to CenterOfGravity field.");
            *fapout++ = ocx;
            *fapout++ = ocy;
            *fapout++ = ocz;
            mxSetFieldByNumber(plhs[1], 0, O_CENTEROFGRAVITY, faout);
        }
        c = bc;
        c2 = bc2;
        c3 = bc3;
        for (cc = 0; cc < nc; ++cc) {
            *c++ -= ocx;
            *c2++ -= ocy;
            *c3++ -= ocz;
        }
    }

    /* no rampimg of normal force */
    if (!normramp)
        rampforce = normforce;

    /* get radius */
    if (sphere) {
        mr = subFunc_radius(nc, bc, bcr);
        if (compdm) {
            compds[0] = 1;
            compds[1] = 1;
            faout = mxCreateNumericArray(2, compds, mxDOUBLE_CLASS, mxREAL);
            if (faout == NULL)
                mexErrMsgTxt("Error creating RadiusMean field in memory.");
            fapout = (double *) mxGetPr(faout);
            if (fapout == NULL)
                mexErrMsgTxt("Error getting double pointer to RadiusMean field.");
            *fapout = mr;
            mxSetFieldByNumber(plhs[1], 0, O_RADIUSMEAN, faout);
        }

        /* make "unit sphere like" vectors */
        c = bc;
        cr = bcr;
        for (cc = 0; cc < nc3; ++cc)
            *c++ /= mr;
    }

    /* further processing depends on morph type */
    switch (mtype) {

        /* smoothing */
        case 0:

            /* for equal weighting */
            if ((equalw) && (!areaw) && (!distw))
                force /= 6.0;

            /* iterated morph */
            for (iterc = 0; iterc < niter; ++iterc) {

                /* check for NaNs */
                if (mxIsNaN(*bc)) {
                    sprintf(vstr, "Coordinates turned to NaN at %d. iteration.", iterc + 1);
                    mexErrMsgTxt(vstr);
                }

                /* ramping of normal force */
                if (norm && normramp)
                    rampforce = rampfrom + (normforce - rampfrom) * ((double) iterc / ((double) (niter - 1)));

                /* to-sphere force */
                if (sphere) {

                    /* subtract any center */
                    subFunc_center(nc, bc, &ncx, &ncy, &ncz);
                    c  = bc;
                    c2 = bc2;
                    c3 = bc3;
                    for (cc = 0; cc < nc; ++cc) {
                        *c++ -= ncx;
                        *c2++ -= ncy;
                        *c3++ -= ncz;
                    }

                    /* pull vertices towards radius 1 */
                    c  = bc;
                    c2 = bc2;
                    c3 = bc3;
                    for (cc = 0; cc < nc; ++cc) {
                        cvx = *c;
                        cvy = *c2;
                        cvz = *c3;
                        mr = sqrt(cvx * cvx + cvy * cvy + cvz * cvz);
                        scale = 1.0 - tosph * (mr - 1) / mr;
                        *c++ *= scale;
                        *c2++ *= scale;
                        *c3++ *= scale;
                    }
                }

                /* copy coordinates */
                c  = bc;
                cb = bcb;
                for (cc = 0; cc < nc3; ++cc)
                    *cb++ = *c++;

                /* apply smoothing */
                c  = bc;
                c2 = bc2;
                c3 = bc3;
                cb = bcb;
                cb2 = bcb2;
                cb3 = bcb3;

                /* with area-based weighting */
                if (areaw) {

                    /* get average triangle-area and area associated with each vertex */
                    tarea = subFunc_area_a(nc, tc, bc, tri, ca) / dnc;

                    /* iterate over vertices */
                    for (cc = 0; cc < nc; ++cc) {

                        /* get number of neighbors and offset */
                        cnn = nn[cc];
                        dnn = (double) cnn;
                        tn = &n[36 * cc];

                        /* sum neighbor coordinates */
                        for (svx = svy = svz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                            svx += bcb[*tn];
                            svy += bcb2[*tn];
                            svz += bcb3[*tn++];
                        }

                        /* subtract coordinate from mean and multiply by force */
                        cvx = *cb++;
                        cvy = *cb2++;
                        cvz = *cb3++;
                        scale = (ca[cc] / tarea);
                        if (scale > 2.0)
                            scale = 2.0;
                        scale *= force;
                        *c++  = cvx + (svx / dnn - cvx) * scale;
                        *c2++ = cvy + (svy / dnn - cvy) * scale;
                        *c3++ = cvz + (svz / dnn - cvz) * scale;
                    }

                /* with distance-based weighting */
                } else if (distw) {

                    /* square of distance */
                    if (distwsq) {

                        /* with normal force */
                        if (norm) {

                            /* iterate over vertices */
                            for (cc = 0; cc < nc; ++cc) {

                                /* get number of neighbors and offset */
                                cnn = nn[cc];
                                tn = &n[36 * cc];

                                /* sum neighbor coordinates */
                                cvx = *cb;
                                cvy = *cb2;
                                cvz = *cb3;
                                e2x = bcb[tn[cnn - 1]] - cvx;
                                e2y = bcb2[tn[cnn - 1]] - cvy;
                                e2z = bcb3[tn[cnn - 1]] - cvz;
                                vlength = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 2.3e-16;
                                e2x /= vlength;
                                e2y /= vlength;
                                e2z /= vlength;
                                for (scale = svx = svy = svz = nsx = nsy = nsz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                    e1x = e2x;
                                    e1y = e2y;
                                    e1z = e2z;
                                    ncx = bcb[*tn];
                                    ncy = bcb2[*tn];
                                    ncz = bcb3[*tn++];
                                    e2x = ncx - cvx;
                                    e2y = ncy - cvy;
                                    e2z = ncz - cvz;
                                    dnn = e2x * e2x + e2y * e2y + e2z * e2z;
                                    vlength = sqrt(dnn) + 2.3e-16;
                                    scale += dnn;
                                    svx += dnn * ncx;
                                    svy += dnn * ncy;
                                    svz += dnn * ncz;
                                    e2x /= vlength;
                                    e2y /= vlength;
                                    e2z /= vlength;
                                    nx = e2y * e1z - e2z * e1y;
                                    ny = e2z * e1x - e2x * e1z;
                                    nz = e2x * e1y - e2y * e1x;
                                    vlength = sqrt(nx * nx + ny * ny + nz * nz) + 2.3e-16;
                                    nsx += nx / vlength;
                                    nsy += ny / vlength;
                                    nsz += nz / vlength;
                                }
                                dnn = (double) cnn;
                                nsx /= dnn;
                                nsy /= dnn;
                                nsz /= dnn;
                                vlength = rampforce / (sqrt(nsx * nsx + nsy * nsy + nsz * nsz) + 2.3e-16);

                                /* subtract coordinate from mean and multiply by force */
                                dnn = force / scale;
                                *c++  = selfforce * (*cb++)  + dnn * svx + nsx * vlength;
                                *c2++ = selfforce * (*cb2++) + dnn * svy + nsy * vlength;
                                *c3++ = selfforce * (*cb3++) + dnn * svz + nsz * vlength;
                            }

                        /* without normal force */
                        } else {

                            /* iterate over vertices */
                            for (cc = 0; cc < nc; ++cc) {

                                /* get number of neighbors and offset */
                                cnn = nn[cc];
                                tn = &n[36 * cc];

                                /* sum neighbor coordinates */
                                for (scale = svx = svy = svz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                    cvx = bcb[*tn] - *cb;
                                    cvy = bcb2[*tn] - *cb2;
                                    cvz = bcb3[*tn] - *cb3;
                                    dnn = cvx * cvx + cvy * cvy + cvz * cvz;
                                    scale += dnn;
                                    svx += dnn * bcb[*tn];
                                    svy += dnn * bcb2[*tn];
                                    svz += dnn * bcb3[*tn++];
                                }

                                /* subtract coordinate from mean and multiply by force */
                                dnn = force / scale;
                                *c++  = selfforce * (*cb++)  + dnn * svx;
                                *c2++ = selfforce * (*cb2++) + dnn * svy;
                                *c3++ = selfforce * (*cb3++) + dnn * svz;
                            }
                        }

                    /* log of (1 + dist) */
                    } else if (distwl) {

                        /* with normal force */
                        if (norm) {

                            /* iterate over vertices */
                            for (cc = 0; cc < nc; ++cc) {

                                /* get number of neighbors and offset */
                                cnn = nn[cc];
                                tn = &n[36 * cc];

                                /* sum neighbor coordinates */
                                cvx = *cb;
                                cvy = *cb2;
                                cvz = *cb3;
                                e2x = bcb[tn[cnn - 1]] - cvx;
                                e2y = bcb2[tn[cnn - 1]] - cvy;
                                e2z = bcb3[tn[cnn - 1]] - cvz;
                                vlength = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 2.3e-16;
                                e2x /= vlength;
                                e2y /= vlength;
                                e2z /= vlength;
                                for (scale = svx = svy = svz = nsx = nsy = nsz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                    e1x = e2x;
                                    e1y = e2y;
                                    e1z = e2z;
                                    ncx = bcb[*tn];
                                    ncy = bcb2[*tn];
                                    ncz = bcb3[*tn++];
                                    e2x = ncx - cvx;
                                    e2y = ncy - cvy;
                                    e2z = ncz - cvz;
                                    dnn = sqrt(e2x * e2x + e2y * e2y + e2z * e2z);
                                    vlength = dnn + 2.3e-16;
                                    dnn = log(1.0 + dnn) + 2.3e-16;
                                    scale += dnn;
                                    svx += dnn * ncx;
                                    svy += dnn * ncy;
                                    svz += dnn * ncz;
                                    e2x /= vlength;
                                    e2y /= vlength;
                                    e2z /= vlength;
                                    nx = e2y * e1z - e2z * e1y;
                                    ny = e2z * e1x - e2x * e1z;
                                    nz = e2x * e1y - e2y * e1x;
                                    vlength = sqrt(nx * nx + ny * ny + nz * nz) + 2.3e-16;
                                    nsx += nx / vlength;
                                    nsy += ny / vlength;
                                    nsz += nz / vlength;
                                }
                                dnn = (double) cnn;
                                nsx /= dnn;
                                nsy /= dnn;
                                nsz /= dnn;
                                vlength = rampforce / (sqrt(nsx * nsx + nsy * nsy + nsz * nsz) + 2.3e-16);

                                /* subtract coordinate from mean and multiply by force */
                                dnn = force / scale;
                                *c++  = selfforce * (*cb++)  + dnn * svx + nsx * vlength;
                                *c2++ = selfforce * (*cb2++) + dnn * svy + nsy * vlength;
                                *c3++ = selfforce * (*cb3++) + dnn * svz + nsz * vlength;
                            }

                        /* without normal force */
                        } else {

                            /* iterate over vertices */
                            for (cc = 0; cc < nc; ++cc) {

                                /* get number of neighbors and offset */
                                cnn = nn[cc];
                                tn = &n[36 * cc];

                                /* sum neighbor coordinates */
                                for (scale = svx = svy = svz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                    cvx = bcb[*tn] - *cb;
                                    cvy = bcb2[*tn] - *cb2;
                                    cvz = bcb3[*tn] - *cb3;
                                    dnn = log(1.0 + sqrt(cvx * cvx + cvy * cvy + cvz * cvz)) + 2.3e-16;
                                    scale += dnn;
                                    svx += dnn * bcb[*tn];
                                    svy += dnn * bcb2[*tn];
                                    svz += dnn * bcb3[*tn++];
                                }

                                /* subtract coordinate from mean and multiply by force */
                                dnn = force / scale;
                                *c++  = selfforce * (*cb++)  + dnn * svx;
                                *c2++ = selfforce * (*cb2++) + dnn * svy;
                                *c3++ = selfforce * (*cb3++) + dnn * svz;
                            }
                        }

                    /* regular distance weighting */
                    } else {

                        /* with normal force */
                        if (norm) {

                            /* iterate over vertices */
                            for (cc = 0; cc < nc; ++cc) {

                                /* get number of neighbors and offset */
                                cnn = nn[cc];
                                tn = &n[36 * cc];

                                /* sum neighbor coordinates */
                                cvx = *cb;
                                cvy = *cb2;
                                cvz = *cb3;
                                e2x = bcb[tn[cnn - 1]] - cvx;
                                e2y = bcb2[tn[cnn - 1]] - cvy;
                                e2z = bcb3[tn[cnn - 1]] - cvz;
                                vlength = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 2.3e-16;
                                e2x /= vlength;
                                e2y /= vlength;
                                e2z /= vlength;
                                for (scale = svx = svy = svz = nsx = nsy = nsz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                    e1x = e2x;
                                    e1y = e2y;
                                    e1z = e2z;
                                    ncx = bcb[*tn];
                                    ncy = bcb2[*tn];
                                    ncz = bcb3[*tn++];
                                    e2x = ncx - cvx;
                                    e2y = ncy - cvy;
                                    e2z = ncz - cvz;
                                    dnn = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 2.3e-16;
                                    scale += dnn;
                                    svx += dnn * ncx;
                                    svy += dnn * ncy;
                                    svz += dnn * ncz;
                                    e2x /= dnn;
                                    e2y /= dnn;
                                    e2z /= dnn;
                                    nx = e2y * e1z - e2z * e1y;
                                    ny = e2z * e1x - e2x * e1z;
                                    nz = e2x * e1y - e2y * e1x;
                                    vlength = sqrt(nx * nx + ny * ny + nz * nz) + 2.3e-16;
                                    nsx += nx / vlength;
                                    nsy += ny / vlength;
                                    nsz += nz / vlength;
                                }
                                dnn = (double) cnn;
                                nsx /= dnn;
                                nsy /= dnn;
                                nsz /= dnn;
                                vlength = rampforce / (sqrt(nsx * nsx + nsy * nsy + nsz * nsz) + 2.3e-16);

                                /* subtract coordinate from mean and multiply by force */
                                dnn = force / scale;
                                *c++  = selfforce * (*cb++)  + dnn * svx + nsx * vlength;
                                *c2++ = selfforce * (*cb2++) + dnn * svy + nsy * vlength;
                                *c3++ = selfforce * (*cb3++) + dnn * svz + nsz * vlength;
                            }

                        /* without normal force */
                        } else {

                            /* iterate over vertices */
                            for (cc = 0; cc < nc; ++cc) {

                                /* get number of neighbors and offset */
                                cnn = nn[cc];
                                tn = &n[36 * cc];

                                /* sum neighbor coordinates */
                                for (scale = svx = svy = svz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                    cvx = bcb[*tn] - *cb;
                                    cvy = bcb2[*tn] - *cb2;
                                    cvz = bcb3[*tn] - *cb3;
                                    dnn = sqrt(cvx * cvx + cvy * cvy + cvz * cvz + 2.3e-16);
                                    scale += dnn;
                                    svx += dnn * bcb[*tn];
                                    svy += dnn * bcb2[*tn];
                                    svz += dnn * bcb3[*tn++];
                                }

                                /* subtract coordinate from mean and multiply by force */
                                dnn = force / scale;
                                *c++  = selfforce * (*cb++)  + dnn * svx;
                                *c2++ = selfforce * (*cb2++) + dnn * svy;
                                *c3++ = selfforce * (*cb3++) + dnn * svz;
                            }
                        }
                    }

                /* without equal weighting */
                } else if (equalw) {

                    /* with normal force */
                    if (norm) {

                        /* iterate over vertices */
                        for (cc = 0; cc < nc; ++cc) {

                            /* get number of neighbors and offset */
                            cnn = nn[cc];
                            tn = &n[36 * cc];

                            /* sum neighbor coordinates */
                            cvx = *cb;
                            cvy = *cb2;
                            cvz = *cb3;
                            e2x = bcb[tn[cnn - 1]] - cvx;
                            e2y = bcb2[tn[cnn - 1]] - cvy;
                            e2z = bcb3[tn[cnn - 1]] - cvz;
                            vlength = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 2.3e-16;
                            e2x /= vlength;
                            e2y /= vlength;
                            e2z /= vlength;
                            for (svx = svy = svz = nsx = nsy = nsz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                e1x = e2x;
                                e1y = e2y;
                                e1z = e2z;
                                ncx = bcb[*tn];
                                ncy = bcb2[*tn];
                                ncz = bcb3[*tn++];
                                e2x = ncx - cvx;
                                e2y = ncy - cvy;
                                e2z = ncz - cvz;
                                svx += e2x;
                                svy += e2y;
                                svz += e2z;
                                vlength = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 2.3e-16;
                                e2x /= vlength;
                                e2y /= vlength;
                                e2z /= vlength;
                                nx = e2y * e1z - e2z * e1y;
                                ny = e2z * e1x - e2x * e1z;
                                nz = e2x * e1y - e2y * e1x;
                                vlength = sqrt(nx * nx + ny * ny + nz * nz) + 2.3e-16;
                                nsx += nx / vlength;
                                nsy += ny / vlength;
                                nsz += nz / vlength;
                            }
                            dnn = (double) cnn;
                            nsx /= dnn;
                            nsy /= dnn;
                            nsz /= dnn;
                            vlength = rampforce / (sqrt(nsx * nsx + nsy * nsy + nsz * nsz) + 2.3e-16);

                            /* subtract coordinate from mean and multiply by force */
                            *c++  = cvx + svx * force + nsx * vlength;
                            *c2++ = cvy + svy * force + nsy * vlength;
                            *c3++ = cvz + svz * force + nsz * vlength;
                        }

                    /* without normal force */
                    } else {

                        /* iterate over vertices */
                        for (cc = 0; cc < nc; ++cc) {

                            /* get number of neighbors and offset */
                            cnn = nn[cc];
                            tn = &n[36 * cc];

                            /* sum neighbor coordinates */
                            cvx = *cb++;
                            cvy = *cb2++;
                            cvz = *cb3++;
                            for (svx = svy = svz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                svx += bcb[*tn] - cvx;
                                svy += bcb2[*tn] - cvy;
                                svz += bcb3[*tn++] - cvz;
                            }

                            /* subtract coordinate from mean and multiply by force */
                            *c++  = cvx + svx * force;
                            *c2++ = cvy + svy * force;
                            *c3++ = cvz + svz * force;
                        }
                    }

                /* without weighting */
                } else {

                    /* with normal force */
                    if (norm) {

                        /* iterater over vertices */
                        for (cc = 0; cc < nc; ++cc) {

                            /* get number of neighbors and offset */
                            cnn = nn[cc];
                            dnn = (double) cnn;
                            tn = &n[36 * cc];

                            /* sum neighbor coordinates */
                            cvx = *cb++;
                            cvy = *cb2++;
                            cvz = *cb3++;
                            e2x = bcb[tn[cnn - 1]] - cvx;
                            e2y = bcb2[tn[cnn - 1]] - cvy;
                            e2z = bcb3[tn[cnn - 1]] - cvz;
                            vlength = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 2.3e-16;
                            e2x /= vlength;
                            e2y /= vlength;
                            e2z /= vlength;
                            for (svx = svy = svz = nsx = nsy = nsz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                e1x = e2x;
                                e1y = e2y;
                                e1z = e2z;
                                ncx = bcb[*tn];
                                ncy = bcb2[*tn];
                                ncz = bcb3[*tn++];
                                svx += ncx;
                                svy += ncy;
                                svz += ncz;
                                e2x = ncx - cvx;
                                e2y = ncy - cvy;
                                e2z = ncz - cvz;
                                vlength = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 2.3e-16;
                                e2x /= vlength;
                                e2y /= vlength;
                                e2z /= vlength;
                                nx = e2y * e1z - e2z * e1y;
                                ny = e2z * e1x - e2x * e1z;
                                nz = e2x * e1y - e2y * e1x;
                                vlength = sqrt(nx * nx + ny * ny + nz * nz) + 2.3e-16;
                                nsx += nx / vlength;
                                nsy += ny / vlength;
                                nsz += nz / vlength;
                            }
                            nsx /= dnn;
                            nsy /= dnn;
                            nsz /= dnn;
                            vlength = rampforce / (sqrt(nsx * nsx + nsy * nsy + nsz * nsz) + 2.3e-16);

                            /* subtract coordinate from mean and multiply by force */
                            *c++  = cvx + (svx / dnn - cvx) * force + nsx * vlength;
                            *c2++ = cvy + (svy / dnn - cvy) * force + nsy * vlength;
                            *c3++ = cvz + (svz / dnn - cvz) * force + nsz * vlength;
                        }

                    /* without normal force */
                    } else {

                        /* iterate over vertices */
                        for (cc = 0; cc < nc; ++cc) {

                            /* get number of neighbors and offset */
                            cnn = nn[cc];
                            dnn = (double) cnn;
                            tn = &n[36 * cc];

                            /* sum neighbor coordinates */
                            for (svx = svy = svz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                                svx += bcb[*tn];
                                svy += bcb2[*tn];
                                svz += bcb3[*tn++];
                            }

                            /* subtract coordinate from mean and multiply by force */
                            cvx = *cb++;
                            cvy = *cb2++;
                            cvz = *cb3++;
                            *c++  = cvx + (svx / dnn - cvx) * force;
                            *c2++ = cvy + (svy / dnn - cvy) * force;
                            *c3++ = cvz + (svz / dnn - cvz) * force;
                        }
                    }
                }

                /* distortion correction */
                if (distcs) {

                    /* get area */
                    subFunc_area_c(nc, tc, bc, tri, ca, cad);

                    /* iterate over vertices again */
                    c  = bc;
                    c2 = bc2;
                    c3 = bc3;
                    cb = bcb;
                    cb2 = bcb2;
                    cb3 = bcb3;
                    for (cc = 0; cc < nc; ++cc) {

                        /* get number of neighbors and offset */
                        tn = &n[36 * cc];
                        cnn = nn[cc];
                        dnn = (double) cnn;
                        cav = ca[cc];

                        /* get neighbor area and fill forces area */
                        dsv = 0.0;
                        bna = na;
                        for (cnc = 0; cnc < cnn; ++cnc) {
                            scale = -(ca[*tn++] / cav);
                            if (scale > 0.0)
                                *bna = scale * scale;
                            else
                                *bna = -scale * scale;
                            dsv += *bna++;
                        }
                        dsv *= distc;

                        /* get coordinate value */
                        cvx = *cb++;
                        cvy = *cb2++;
                        cvz = *cb3++;

                        /* add neighbor differences to correction vector */
                        tn = &n[36 * cc];
                        bna = na;
                        for (svx = svy = svz = 0.0, cnc = 0; cnc < cnn; ++cnc) {
                            scale = *bna++;
                            svx += scale * (bcb[*tn] - cvx);
                            svy += scale * (bcb2[*tn] - cvy);
                            svz += scale * (bcb3[*tn++] - cvz);
                        }

                        /* add weighted correction vector */
                        *c++  = cvx + svx / dsv;
                        *c2++ = cvy + svy / dsv;
                        *c3++ = cvz + svz / dsv;
                    }
                }
            }

            /* area constantness */
            if (areac) {

                /* get new area and center */
                nareav = subFunc_area(nc, tc, bc, tri);
                subFunc_center(nc, bc, &ncx, &ncy, &ncz);

                /* compute scale */
                scale = sqrt(areav / nareav);

                /* iterate over coordinates */
                c = bc;
                c2 = bc2;
                c3 = bc3;
                for (cc = 0; cc < nc; ++cc) {
                    *c++ = (bc[cc] - ncx) * scale + ocx;
                    *c2++ = (bc2[cc] - ncy) * scale + ocy;
                    *c3++ = (bc3[cc] - ncz) * scale + ocz;
                }

            /* normal-vector force was given */
            } else if (norm) {
                c = bc;
                c2 = bc2;
                c3 = bc3;
                for (cc = 0; cc < nc; ++cc) {
                    *c++ += ocx;
                    *c2++ += ocy;
                    *c3++ += ocz;
                }
            }
            break;
    }

    /* build density map output */
    if (compdm) {
        compds[0] = nc;
        compds[1] = 1;
        faout = mxCreateNumericArray(2, compds, mxDOUBLE_CLASS, mxREAL);
        fapout = (double *) mxGetPr(faout);
        nareav = subFunc_area_c(nc, tc, bc, tri, fapout, cad);
        mxSetFieldByNumber(plhs[1], 0, O_DENSITYSMPDATA, faout);
    }

    /* free memory */
    mxFree(bcb);
    mxFree(n);
    mxFree(nn);
    mxFree(tri);
    mxFree(ca);
    mxFree(cad);
    mxFree(cr);
    mxFree(oca);
}

double subFunc_area(int numc, int numt, const double *c, const unsigned long *t)
{
    /* return value */
    double rv = 0.0;

    /* helper pointers */
    const double *c2, *c3;
    const unsigned long *t2, *t3;
    unsigned long it1, it2, it3;

    /* coordinates / sides */
    double c1x, c1y, c1z, c2x, c2y, c2z, s1, cosval;

    /* triangle counter */
    int tc;

    /* get helper pointers */
    c2 = &c[numc];
    c3 = &c[numc * 2];
    t2 = &t[numt];
    t3 = &t[numt * 2];

    /* iterate over triangles */
    for (tc = 0; tc < numt; ++tc) {

        /* get coordinate indices */
        it1 = t[tc];
        it2 = t2[tc];
        it3 = t3[tc];

        /* compute orthogonal sides */
        c1x  = - c[it1];
        c1y  = -c2[it1];
        c1z  = -c3[it1];
        c2x  =  c[it2] + c1x;
        c2y  = c2[it2] + c1y;
        c2z  = c3[it2] + c1z;
        c1x +=  c[it3];
        c1y += c2[it3];
        c1z += c3[it3];
        s1 = (c1x * c1x + c1y * c1y + c1z * c1z);
        cosval = (c1x * c2x + c1y * c2y + c1z * c2z) / s1;
        c2x -= cosval * c1x;
        c2y -= cosval * c1y;
        c2z -= cosval * c1z;
        cosval = 0.5 * sqrt(s1 * (c2x * c2x + c2y * c2y + c2z * c2z));

        /* add to area */
        rv += cosval;
    }
    return rv;
}

double subFunc_area_a(int numc, int numt, const double *c, const unsigned long *t, double *ca)
{
    /* return value */
    double rv = 0.0;

    /* variable output string */
    /* char vstr[256]; */

    /* helper pointers */
    const double *c2, *c3;
    const unsigned long *t2, *t3;
    unsigned long it1, it2, it3;

    /* coordinates and sides and value */
    double c1x, c1y, c1z, c2x, c2y, c2z, s1, cosval;

    /* triangle counter */
    int tc;

    /* double pointer backup */
    double *dpb = ca;

    /* get helper pointers */
    c2 = &c[numc];
    c3 = &c[numc * 2];
    t2 = &t[numt];
    t3 = &t[numt * 2];

    /* delete area values in ca */
    for (tc = 0; tc < numc; ++tc) *dpb++ = 0;

    /* iterate over triangles */
    for (tc = 0; tc < numt; ++tc) {

        /* get coordinate indices */
        it1 = t[tc];
        it2 = t2[tc];
        it3 = t3[tc];

        /* compute orthogonal sides */
        c1x  = - c[it1];
        c1y  = -c2[it1];
        c1z  = -c3[it1];
        c2x  =  c[it2] + c1x;
        c2y  = c2[it2] + c1y;
        c2z  = c3[it2] + c1z;
        c1x +=  c[it3];
        c1y += c2[it3];
        c1z += c3[it3];
        s1 = (c1x * c1x + c1y * c1y + c1z * c1z);
        cosval = (c1x * c2x + c1y * c2y + c1z * c2z) / s1;
        c2x -= cosval * c1x;
        c2y -= cosval * c1y;
        c2z -= cosval * c1z;
        cosval = 0.5 * sqrt(s1 * (c2x * c2x + c2y * c2y + c2z * c2z));

        /* add to area */
        rv += cosval;

        /* add to area of vertices */
        ca[it1] += cosval;
        ca[it2] += cosval;
        ca[it3] += cosval;
    }

    /* devide by number of triangles */
    for (tc = 0; tc < numc; ++tc)
        *ca++ /= 3.0;

    return rv;
}

double subFunc_area_c(int numc, int numt, const double *c, const unsigned long *t, double *ca, const unsigned char *cad)
{
    /* return value */
    double rv = 0.0;

    /* variable output string */
    /* char vstr[256]; */

    /* helper pointers */
    const double *c2, *c3;
    const unsigned long *t2, *t3;
    unsigned long it1, it2, it3;

    /* coordinates and sides and value */
    double c1x, c1y, c1z, c2x, c2y, c2z, s1, cosval;

    /* triangle counter */
    int tc;

    /* double pointer backup */
    double *dpb = ca;

    /* get helper pointers */
    c2 = &c[numc];
    c3 = &c[numc * 2];
    t2 = &t[numt];
    t3 = &t[numt * 2];

    /* delete area values in ca */
    for (tc = 0; tc < numc; ++tc) *dpb++ = 0;

    /* iterate over triangles */
    for (tc = 0; tc < numt; ++tc) {

        /* get coordinate indices */
        it1 = t[tc];
        it2 = t2[tc];
        it3 = t3[tc];

        /* compute orthogonal sides */
        c1x  = - c[it1];
        c1y  = -c2[it1];
        c1z  = -c3[it1];
        c2x  =  c[it2] + c1x;
        c2y  = c2[it2] + c1y;
        c2z  = c3[it2] + c1z;
        c1x +=  c[it3];
        c1y += c2[it3];
        c1z += c3[it3];
        s1 = (c1x * c1x + c1y * c1y + c1z * c1z);
        cosval = (c1x * c2x + c1y * c2y + c1z * c2z) / s1;
        c2x -= cosval * c1x;
        c2y -= cosval * c1y;
        c2z -= cosval * c1z;
        cosval = 0.5 * sqrt(s1 * (c2x * c2x + c2y * c2y + c2z * c2z));

        /* add to area */
        rv += cosval;

        /* add to area of vertices */
        ca[it1] += cosval;
        ca[it2] += cosval;
        ca[it3] += cosval;
    }

    /* devide by number of triangles */
    for (tc = 0; tc < numc; ++tc)
        *ca++ /= (double) *cad++;

    return rv;
}

void subFunc_center(int numc, const double *c, double *cx, double *cy, double *cz)
{
    /* double number of c */
    const double *c2, *c3;
    double numcd = (double) numc;

    /* vertex counter */
    int vc;

    /* set to 0 */
    double ccx = 0.0;
    double ccy = 0.0;
    double ccz = 0.0;

    /* get helper pointers */
    c2 = &c[numc];
    c3 = &c[numc * 2];

    /* loop */
    for (vc = 0; vc < numc; ++vc) {
        ccx += *c++;
        ccy += *c2++;
        ccz += *c3++;
    }
    *cx = ccx / numcd;
    *cy = ccy / numcd;
    *cz = ccz / numcd;
}

void subFunc_normals(int numc, int numt, const double *c, const unsigned long *t, double *nrm)
{
}

double subFunc_radius(int numc, const double *c, double *cr)
{
    /* double number of c and mean radius */
    int numc2 = numc * 2;
    double cx, cy, cz;
    double mv = 0.0;
    double rv;

    /* vertex counter */
    int vc;

    /* loop */
    for (vc = 0; vc < numc; ++vc) {
        cy = c[numc];
        cz = c[numc2];
        cx = *c++;
        rv = sqrt(cx * cx + cy * cy + cz * cz);
        *cr++ = rv;
        mv += rv;
    }

    /* mean radius */
    mv /= (double) numc;
    return mv;
}

void subFunc_vertexnrareas(int numc, int numt, const unsigned long *t, unsigned char *cad)
{
    /* double number of t */
    int numt2 = numt * 2;

    /* triangle counter */
    int tc;

    /* delete area values in ca */
    unsigned char *dpb = cad;
    for (tc = 0; tc < numc; ++tc) *dpb++ = 0;

    /* iterate over triangles */
    for (tc = 0; tc < numt; ++tc) {

        /* increase counters */
        ++cad[t[tc]];
        ++cad[t[tc + numt]];
        ++cad[t[tc + numt2]];
    }
}

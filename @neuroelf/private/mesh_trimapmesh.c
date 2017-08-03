/*

mapping of two meshes, barycoord-based

FORMAT:       [t, el] = mesh_trimapmesh(t, s, st, str)

Input fields:

      t           Tx3 target coordinate list
      s           Sx3 source coordinate list
      st          Tx3 source triangle list
      str         Sx1 cell array with triangle back references

Output fields:

      t           source vertex (one output) or triangle (two outputs)
      w           weights for vertices 2, 3 (w1 = 1 - (w2+w3))

Note: coordinates must be centered around [0, 0, 0] and rescaled to r=1

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

#define MESH_PIQ 0.78539816339744830962
#define MESH_PIH 1.57079632679489661923
#define MESH_PIT 2.35619449019234492885
#define MESH_PI  3.14159265358979323846
#define MESH_PID 6.28318530717958647692
#define MAXNEI 36

/* mexCallMATLAB(1, plhs, 4, flexin, "FUNCTION"); */

void spherecoords(double cx, double cy, double cz, double *sr, double *sp, double *st)
{
    *sr = sqrt(cx * cx + cy * cy + cz * cz);
    *sp = acos(cz / *sr);
    *st = atan2(cy, cx);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* counters */
    int cc = 0, cc2 = 0;

    /* number of target coordinates (incl. double) */
    int tnc = 0, tnc2 = 0;

    /* number of source coordinates (incl. double) and triangles (incl. triple) */
    int snc = 0, snc2 = 0, stc = 0, stc3 = 0;

    /* pointer to target and source coordinates and source triangles */
    const double *tc = NULL, *sc = NULL, *stri = NULL;

    /* target and source coordinate value, computational value, distance and min dist */
    double tcx = 0.0, tcy = 0.0, tcz = 0.0;
    double scx = 0.0, scy = 0.0, scz = 0.0;
    double dval = 0.0, vdist = 0.0, mdist = 0.0;
    double mvl1 = 0.0, mvl2 = 0.0, minvl1 = 0.0, minvl2 = 0.0, minvls = 0.0;

    /* triangle vectors, matrix variables, etc. */
    double tvx1 = 0.0, tvy1 = 0.0, tvz1 = 0.0, tvx2 = 0.0, tvy2 = 0.0, tvz2 = 0.0;
    double pif = 0.0, cvm11 = 0.0, cvm12 = 0.0, cvm22 = 0.0, cim11 = 0.0, cim12 = 0.0, cim22 = 0.0;

    /* const any-purpose pointer, movable target and source coordinate pointers */
    const double *capdp = NULL, *cdp1 = NULL, *cdp2 = NULL, *cdp3 = NULL, *csp2 = NULL, *csp3 = NULL;

    /* integer triangle indices list, backreference and triangle list pointers */
    int *striu = NULL, *strib = NULL, *t1v = NULL, *t2v = NULL, *t3v = NULL;

    /* any-purpose, movable int pointer */
    int *aplp = NULL;

    /* match counters and matching index, triangle vertices */
    int m1 = 0, m2 = 0, mdisti = 0, tv1 = 0, tv2 = 0, tv3 = 0;
    bool matchok = 0;

    /* quad splits, number, squared, number - 1, and squared - qs */
    int qs = 0, qss = 0, qsm = 0, qssm = 0;

    /* quad identifier for target(1) and source(1) as well as counters */
    int qc1 = 0, qc2 = 0, qc1f = 0, qc1t = 0, qc2f = 0, qc2t = 0, qc12f = 0, qc12t = 0;

    /* quad pointers */
    int *q2 = NULL, *qp = NULL, *qsrt = NULL;

    /* interger values and any-purpose pointers */
    int ival1 = 0, ival2 = 0, *apip = NULL, *apip2 = NULL;

    /* dim argument for sort call */
    int sortinsz[2] = {0, 1};
    double dqs, dqsh;
    mxArray *sortin;
    mxArray *sortres[2];
    const int *qsrti;

    /* RO cell element pointer */
    const mxArray *cellcont;

    /* output array size and pointers */
    const int *di;
    int ndi[2] = {0, 0};
    double *otl = NULL, *oll = NULL, *ol2 = NULL;

    /* sorting result */

    /* variable output string and NaN counter */
    char vstr[256];

    /* check number of in/out arguments */
	if ((nrhs != 4) || (nlhs > 2))
		mexErrMsgTxt("Bad number of input/output arguments.");

    /* check input argument types */
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsCell(prhs[3]))
        mexErrMsgTxt("Bad argument type. Must be (double, double, double, cell).");
    if ((mxGetNumberOfDimensions(prhs[0]) != 2) ||
        (mxGetNumberOfDimensions(prhs[1]) != 2) ||
        (mxGetNumberOfDimensions(prhs[2]) != 2))
        mexErrMsgTxt("Bad number of dims on coordinate/triangle lists.");

    /* check coordinates argument size (and get number of coordinates from this) */
    di = mxGetDimensions(prhs[0]);
    if (di[1] != 3)
        mexErrMsgTxt("Target coordinates must be Cx3 matrix.");
    tnc = di[0];
    tnc2 = tnc * 2;
    tc = (double *) mxGetPr(prhs[0]);

    di = mxGetDimensions(prhs[1]);
    if (di[1] != 3)
        mexErrMsgTxt("Source coordinates must be Cx3 matrix.");
    snc = di[0];
    snc2 = snc * 2;
    sc = (double *) mxGetPr(prhs[1]);

    di = mxGetDimensions(prhs[2]);
    if (di[1] != 3)
        mexErrMsgTxt("Triangles must be a Tx3 matrix.");
    stc = di[0];
    stc3 = stc * 3;
    stri = (double *) mxGetPr(prhs[2]);

    /* check options argument size */
    if (mxGetNumberOfElements(prhs[3]) != snc)
        mexErrMsgTxt("Number of triangle reference lists must match with coordinates.");

    /* find good measure for quad areas */
    qs = 2 * ((int) floor(sqrt(((double) snc) / 2.0) / 2.0));
    dqs = (double) qs;
    dqsh = dqs / 2.0;
    qss = qs * qs;
    qsm = qs - 1;
    qssm = qss - qs;

    /* */
    /* argument consistency check done, allocate memory */
    /* */
    striu = (int *) mxCalloc(3 * stc, sizeof(int));
    if (striu == NULL) {
        mexErrMsgTxt("Error allocating memory for uint32 triangles list.");
    }
    strib = (int *) mxCalloc(MAXNEI * snc, sizeof(int));
    if (strib == NULL) {
        mxFree(striu);
        mexErrMsgTxt("Error allocating memory for uint32 triangle back reference lists.");
    }
    *sortinsz = stc;
    sortin = mxCreateNumericArray(2, sortinsz, mxINT32_CLASS, mxREAL);
    if (sortin == NULL) {
        mxFree(striu);
        mxFree(strib);
        mexErrMsgTxt("Error allocating memory for quad list of source coordinates.");
    }
    q2 = (int *) mxGetData(sortin);
    if (q2 == NULL) {
        mxFree(striu);
        mxFree(strib);
        mxDestroyArray(sortin);
        mexErrMsgTxt("Error getting pointer to memory for quad list of source coordinates.");
    }
    qsrt = (int *) mxCalloc(stc, sizeof(int));
    if (qsrt == NULL) {
        mxFree(striu);
        mxFree(strib);
        mxDestroyArray(sortin);
        mexErrMsgTxt("Error allocating memory for sorted quad list of source coordinates.");
    }
    qp = (int *) mxCalloc(qss + 1, sizeof(int));
    if (qp == NULL) {
        mxFree(striu);
        mxFree(strib);
        mxDestroyArray(sortin);
        mxFree(qsrt);
        mexErrMsgTxt("Error allocating memory for quad list sort pointers.");
    }

    /* create output arrays */
    ndi[0] = tnc;
    ndi[1] = 1;
    plhs[0] = mxCreateNumericArray(2, ndi, mxDOUBLE_CLASS, mxREAL);
    if ((plhs[0] == NULL) || ((otl = (double *) mxGetPr(plhs[0])) == NULL)) {
        mxFree(striu);
        mxFree(strib);
        mxDestroyArray(sortin);
        mxFree(qsrt);
        mxFree(qp);
        mexErrMsgTxt("Error allocating output: list of triangle indices.");
    }
    ndi[1] = 2;
    plhs[1] = mxCreateNumericArray(2, ndi, mxDOUBLE_CLASS, mxREAL);
    if ((plhs[1] == NULL) || ((oll = (double *) mxGetPr(plhs[1])) == NULL)) {
        mxFree(striu);
        mxFree(strib);
        mxDestroyArray(sortin);
        mxFree(qsrt);
        mxFree(qp);
        mxDestroyArray(plhs[0]);
        mexErrMsgTxt("Error allocating output: list of triangle indices.");
    }
    ol2 = &oll[tnc];

    /* copy triangle list */
    t1v = striu; t2v = &t1v[stc]; t3v = &t2v[stc];
    for (cc = 0; cc < stc3; ++cc)
        *t1v++ = ((int) (*stri++)) - 1;
    t1v = striu;

    /* fill back reference lists */
    for (cc = 0; cc < snc; ++cc) {
        cellcont = mxGetCell(prhs[3], cc);
        if (!mxIsDouble(cellcont) || (mxGetNumberOfDimensions(cellcont) != 2)) {
            mxFree(striu);
            mxFree(strib);
            mxDestroyArray(sortin);
            mxFree(qsrt);
            mxFree(qp);
            mxDestroyArray(plhs[0]);
            mxDestroyArray(plhs[1]);
            mexErrMsgTxt("Bad cell contents (reference list).");
        }
        di = mxGetDimensions(cellcont);
        if ((di[0] != 1) || (di[1] < 3) || (di[1] >= MAXNEI)) {
            mxFree(striu);
            mxFree(strib);
            mxDestroyArray(sortin);
            mxFree(qsrt);
            mxFree(qp);
            mxDestroyArray(plhs[0]);
            mxDestroyArray(plhs[1]);
            mexErrMsgTxt("Bad cell contents (reference list entry).");
        }
        capdp = (double *) mxGetPr(cellcont);
        aplp = &strib[cc*MAXNEI];
        ival1 = di[1];
        *aplp++ = ival1;
        for (cc2 = 0; cc2 < ival1; ++cc2)
            *aplp++ = ((int) (*capdp++)) - 1;
    }

    /* computing spherical coordinates of target and source, and fill quad arrays */
    cdp1 = sc; cdp2 = &sc[snc]; cdp3 = &sc[snc2];
    apip = q2;
    for (cc = 0; cc < stc; ++cc) {
        tv1 = *t1v++;
        tv2 = *t2v++;
        tv3 = *t3v++;
        tcx = cdp1[tv1] + cdp1[tv2] + cdp1[tv3];
        tcy = cdp2[tv1] + cdp2[tv2] + cdp2[tv3];
        tcz = cdp3[tv1] + cdp3[tv2] + cdp3[tv3];
        scx = sqrt(tcx * tcx + tcy * tcy + tcz * tcz);
        ival1 = qs * ((int) floor(dqsh * (1.0 + tcz / scx)));
        ival2 = ((int) floor(dqsh * (1.0 + (atan2(tcy, tcx) / MESH_PI))));
        *apip++ = ((ival1 >= qss) ? qssm : ival1) + ((ival2 >= qs) ? qsm : ival2);
    }
    t1v = striu;
    t2v = &t1v[stc];
    t3v = &t2v[stc];

    /* sort q2 (use Matlab's sort for now) */
    mexCallMATLAB(2, sortres, 1, &sortin, "sort");

    /* copy result and denote outcome */
    qsrti = (int *) mxGetData(sortres[0]);
    capdp = (double *) mxGetPr(sortres[1]);
    apip = qsrt;
    apip2 = qp;
    *apip2++ = 0;
    ival1 = 0;
    for (cc = 0; cc < stc; ++cc) {
        *apip++ = ((int) (*capdp++)) - 1;
        ival2 = *qsrti++;
        for (; ival1 < ival2; ++ival1)
            *apip2++ = cc;
    }
    *apip2 = cc;
    mxDestroyArray(sortres[0]);
    mxDestroyArray(sortres[1]);

    /* */
    /* here begins the main loop */
    /* */
    cdp1 = tc; cdp2 = &tc[tnc]; cdp3 = &tc[tnc2]; csp2 = &sc[snc], csp3 = &sc[snc2];
    for (m1 = 0; m1 < tnc; ++m1) {

        /* target coordinate */
        tcx = *cdp1++;
        tcy = *cdp2++;
        tcz = *cdp3++;

        /* compute spherical coordinate and matching quad */
        scx = sqrt(tcx * tcx + tcy * tcy + tcz * tcz);
        scy = acos(tcz / scx);
        ival1 = ((int) floor(dqsh * (1.0 + tcz / scx)));
        ival2 = ((int) floor(dqsh * (1.0 + (atan2(tcy, tcx) / MESH_PI))));
        ival1 = (ival1 >= qs) ? qsm : ival1;
        ival2 = (ival2 >= qs) ? qsm : ival2;

        /* check which quads to look at */
        dval = (scy >= MESH_PIH) ? (MESH_PIH - scy) : scy;
        if ((ival1 == 0) || (ival1 == qsm)) {
            qc12f = 1;
            qc12t = (int) dqsh;
        } else if (dval <= 0.3) {
            qc12f = 1;
            qc12t = 2;
        } else if (dval <= 1.2) {
            qc12f = 1;
            qc12t = 1;
        } else {
            qc12f = 2;
            qc12t = 1;
        }
        qc1f = ((ival1 < qc12f) ? 0 : (ival1 - qc12f));
        qc1t = (((ival1 + qc1f) >= qs) ? (qs - 1) : (ival1 + qc12f));
        qc2f = qs + ival2 - qc12t;
        qc2t = qs + ival2 + qc12t;

        /* set min distance to 2 and index to -1 */
        mdist = 2.0;
        mdisti = -1;
        minvls = 1.05;
        matchok = 0;

        /* now check all entries in the selected quads */
        for (qc1 = qc1f; qc1 <= qc1t; ++qc1) {
            for (qc2 = qc2f; qc2 <= qc2t; ++qc2) {
                qc12f = qs * qc1 + (qc2 % qs);
                qc12t = qp[qc12f+1];
                qc12f = qp[qc12f];
                apip = &qsrt[qc12f];

                /* check  */
                for (; qc12f < qc12t; ++qc12f) {
                    m2 = *apip++;
                    tv1 = t1v[m2];
                    tv2 = t2v[m2];
                    tv3 = t3v[m2];
                    tvx1 =   sc[tv1];
                    tvy1 = csp2[tv1];
                    tvz1 = csp3[tv1];
                    scx  = tcx - tvx1;
                    scy  = tcy - tvy1;
                    scz  = tcz - tvz1;
                    vdist = sqrt(scx * scx + scy * scy + scz * scz);
                    if (vdist > (3 * mdist))
                        continue;
                    tvx2 =   sc[tv3] - tvx1;
                    tvy2 = csp2[tv3] - tvy1;
                    tvz2 = csp3[tv3] - tvz1;
                    tvx1 =   sc[tv2] - tvx1;
                    tvy1 = csp2[tv2] - tvy1;
                    tvz1 = csp3[tv2] - tvz1;

                    /* compute covariance matrix elements */
                    cvm11 = tvx1 * tvx1 + tvy1 * tvy1 + tvz1 * tvz1;
                    cvm12 = tvx1 * tvx2 + tvy1 * tvy2 + tvz1 * tvz2;
                    cvm22 = tvx2 * tvx2 + tvy2 * tvy2 + tvz2 * tvz2;

                    /* inverse matrix */
                    pif   = 1.0 / (cvm11 * cvm22 - cvm12 * cvm12);
                    cim11 = pif * cvm22;
                    cim12 = -pif * cvm12;
                    cim22 = pif * cvm11;

                    /* multiply matrix evaluation */
                    mvl1 = scx * (cim11 * tvx1 + cim12 * tvx2) + scy * (cim11 * tvy1 + cim12 * tvy2) + scz * (cim11 * tvz1 + cim12 * tvz2);
                    mvl2 = scx * (cim12 * tvx1 + cim22 * tvx2) + scy * (cim12 * tvy1 + cim22 * tvy2) + scz * (cim12 * tvz1 + cim22 * tvz2);
                    dval = ((mvl1 < 0.0) ? -mvl1 : mvl1) + ((mvl2 < 0.0) ? -mvl2 : mvl2);

                    /* if both positive, best match for now */
                    if ((mvl1 >= -0.002) && (mvl2 >= -0.002) && (dval < minvls) && (vdist < (12 * mdist))) {
                        mdisti = m2;
                        mdist = vdist;
                        minvl1 = mvl1;
                        minvl2 = mvl2;
                        minvls = dval;
                        if ((minvls < 0.999) && (minvl1 >= 0.0) && (minvl2 >= 0.0)) {
                            matchok = 1;
                            break;
                        }
                    }
                }
                if (matchok)
                    break;
            }
            if (matchok)
                break;
        }

        /* if nothing, we must extend the search radius (~sigh~) */
        if (mdisti < 0) {

            /* check which quads to look at */
            scy = acos(tcz / scx);
            dval = (scy >= MESH_PIH) ? (MESH_PIH - scy) : scy;
            if ((ival1 == 0) || (ival1 == qsm)) {
                qc12f = 3;
                qc12t = (int) dqsh;
            } else if (dval <= 0.3) {
                qc12f = 4;
                qc12t = 6;
            } else if (dval <= 1.2) {
                qc12f = 5;
                qc12t = 5;
            } else {
                qc12f = 6;
                qc12t = 4;
            }
            qc1f = ((ival1 < qc12f) ? 0 : (ival1 - qc12f));
            qc1t = (((ival1 + qc1f) >= qs) ? (qs - 1) : (ival1 + qc12f));
            qc2f = qs + ival2 - qc12t;
            qc2t = qs + ival2 + qc12t;

            /* set min distance to 2 and index to -1 */
            mdist = 2.0;
            mdisti = -1;
            minvls = 3.00 ;

            /* now check all entries in the selected quads */
            for (qc1 = qc1f; qc1 <= qc1t; ++qc1) {
                for (qc2 = qc2f; qc2 <= qc2t; ++qc2) {
                    qc12f = qs * qc1 + (qc2 % qs);
                    qc12t = qp[qc12f+1];
                    qc12f = qp[qc12f];
                    apip = &qsrt[qc12f];

                    /* check  */
                    for (; qc12f < qc12t; ++qc12f) {
                        m2 = *apip++;
                        tv1 = t1v[m2];
                        tv2 = t2v[m2];
                        tv3 = t3v[m2];
                        tvx1 =   sc[tv1];
                        tvy1 = csp2[tv1];
                        tvz1 = csp3[tv1];
                        scx  = tcx - tvx1;
                        scy  = tcy - tvy1;
                        scz  = tcz - tvz1;
                        vdist = sqrt(scx * scx + scy * scy + scz * scz);
                        if (vdist > (3 * mdist))
                            continue;
                        tvx2 =   sc[tv3] - tvx1;
                        tvy2 = csp2[tv3] - tvy1;
                        tvz2 = csp3[tv3] - tvz1;
                        tvx1 =   sc[tv2] - tvx1;
                        tvy1 = csp2[tv2] - tvy1;
                        tvz1 = csp3[tv2] - tvz1;

                        /* compute covariance matrix elements */
                        cvm11 = tvx1 * tvx1 + tvy1 * tvy1 + tvz1 * tvz1;
                        cvm12 = tvx1 * tvx2 + tvy1 * tvy2 + tvz1 * tvz2;
                        cvm22 = tvx2 * tvx2 + tvy2 * tvy2 + tvz2 * tvz2;

                        /* inverse matrix */
                        pif   = 1.0 / (cvm11 * cvm22 - cvm12 * cvm12);
                        cim11 = pif * cvm22;
                        cim12 = -pif * cvm12;
                        cim22 = pif * cvm11;

                        /* multiply matrix evaluation */
                        mvl1 = scx * (cim11 * tvx1 + cim12 * tvx2) + scy * (cim11 * tvy1 + cim12 * tvy2) + scz * (cim11 * tvz1 + cim12 * tvz2);
                        mvl2 = scx * (cim12 * tvx1 + cim22 * tvx2) + scy * (cim12 * tvy1 + cim22 * tvy2) + scz * (cim12 * tvz1 + cim22 * tvz2);
                        dval = ((mvl1 < 0.0) ? -mvl1 : mvl1) + ((mvl2 < 0.0) ? -mvl2 : mvl2);

                        /* if both positive, best match for now */
                        if ((mvl1 >= -0.001) && (mvl2 >= -0.001) && (dval < 1.002) && (vdist < (12 * mdist))) {
                            mdisti = m2;
                            mdist = vdist;
                            minvl1 = mvl1;
                            minvl2 = mvl2;
                            minvls = dval;
                            matchok = 1;
                            break;
                        } else if ((dval < minvls) && (vdist < (12 * mdist))) {
                            mdisti = m2;
                            mdist = vdist;
                            minvl1 = mvl1;
                            minvl2 = mvl2;
                            minvls = dval;
                            if ((mvl1 >= -0.001) && (mvl2 >= -0.001) && (minvls <= 0.999)) {
                                matchok = 1;
                                break;
                            }
                        }
                    }
                    if (matchok)
                        break;
                }
                if (matchok)
                    break;
            }
        }

        /* still nothing? then giving up (for now...) */
        if (mdisti < 0) {
            sprintf(vstr, "No match found for vertex %d.", m1 + 1);
            mexWarnMsgTxt(vstr);
        }

        /* set in list */
        *otl++ = (double) (mdisti + 1);
        *oll++ = minvl1;
        *ol2++ = minvl2;
    }

    /* */
    /* free memory */
    /* */
    mxFree(strib);
    mxFree(striu);
    mxDestroyArray(sortin);
    mxFree(qp);
    mxFree(qsrt);
}

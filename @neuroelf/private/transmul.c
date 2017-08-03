/*

   transmul.c -- return the multiplication of M' * M or M * M2

% Version:  v0.9b
% Build:    11050512
% Date:     Aug-26 2010, 1:53 PM EST
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    unsigned char transflag = 0;
    const int *mn, *mn2;
    int mnout[3], c, c1, c2, d3, id = 0, m = 0, mc1, n = 0, nc1, nd, nn;
    const double *a, *ac1, *ac1c, *ac2, *ac3, *b;
    double *p, s, v1, v2, v3, v4, v5,
            s11, s12, s13, s14, s15,
            s22, s23, s24, s25,
            s33, s34, s35, s44, s45, s55;

    /* argument check */
    if ((nrhs < 1) ||
        (nrhs > 3) ||
        (nlhs > 1))
		mexErrMsgTxt("Bad number of input/output arguments.");

    /* first input argument */
    nd = mxGetNumberOfDimensions(*prhs);
    if ((nd > 3) ||
        (mxGetClassID(*prhs) != mxDOUBLE_CLASS))
        mexErrMsgTxt("Invalid input argument.");

    /* get dimensions */
    mn = mxGetDimensions(*prhs);

    /* only one input argument */
    if (nrhs == 1) {

        /* the output is n * n * size(X, 3) */
        m = *mn++;
        *mnout = n = *mn++;
        mnout[1] = n;
        if (nd > 2)
            d3 = *mn;
        else
            d3 = 1;
        mnout[2] = d3;

    /* two input arguments */
    } else {

        /* check type */
        if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
            mexErrMsgTxt("Invalid second input argument.");

        /* check if second argument matches in number of dims */
        if (mxGetNumberOfDimensions(prhs[1]) != nd)
            mexErrMsgTxt("Second input must match in size.");

        /* get second dimension argument */
        mn2 = mxGetDimensions(prhs[1]);

        /* check 3rd dim of both also (if necessary) */
        if ((nd > 2) &&
            (mn[2] != mn2[2]))
            mexErrMsgTxt("3D arguments must match in 3rd dim.");

        /* check if a transpose flag is given */
        if ((nrhs > 2) &&
            (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS) &&
            (mxGetNumberOfElements(prhs[2]) == 1)) {

            /* get flag data pointer */
            b = (const double *) mxGetData(prhs[2]);
            if (!mxIsInf(*b) &&
                !mxIsNaN(*b) &&
                (*b >= 0.0) &&
                (*b <= 3.0))
                transflag = (unsigned char) *b;
        }

        /* depending on transflag */
        switch (transflag) {

            /* no transpose of either a or b */
            case 0:
                m = *mn;
                id = *mn2;
                n = mn2[1];
                if (mn[1] != id)
                    transflag = 4;
                break;

            /* transpose of a, but not b */
            case 1:
                m = mn[1];
                id = *mn;
                n = mn2[1];
                if (id != *mn2)
                    transflag = 4;
                break;

            /* transpose of b, but not a */
            case 2:
                m = *mn;
                id = mn[1];
                n = *mn2;
                if (id != mn2[1])
                    transflag = 4;
                break;

            /* transpose of a and b */
            case 3:
                m = mn[1];
                id = *mn;
                n = *mn2;
                if (id != mn2[1])
                    transflag = 4;
                break;
        }

        /* check inner product size */
        if (transflag > 3)
            mexErrMsgTxt("Inner size must match.");

        /* fill output dim */
        *mnout = m;
        mnout[1] = n;
        if (nd > 2)
            d3 = mn[2];
        else
            d3 = 1;
        mnout[2] = d3;
    }

    /* get data of first input */
    a = (const double*) mxGetData(*prhs);

    /* create output matrices */
    *plhs = mxCreateNumericArray(nd, mnout, mxDOUBLE_CLASS, mxREAL);
    p = (double *) mxGetData(*plhs);

    /* with just one input */
    if (nrhs == 1) {

        /* allow for special cases with small n */
        switch (n) {
            case 0: break;
            case 1:
                for ( ; d3 > 0; --d3) {
                    for (c = m, s = 0.0; c > 0; --c) {
                        s += *a * *a;
                        ++a;
                    }
                    *p++ = s;
                    a = &a[m];
                }
                break;
            case 2:
                for ( ; d3 > 0; --d3) {
                    ac1 = &a[m];
                    for (c = m, s11 = s12 = s22 = 0.0; c > 0; --c) {
                        v1 = *a++;
                        v2 = *ac1++;
                        s11 += v1 * v1;
                        s12 += v1 * v2;
                        s22 += v2 * v2;
                    }
                    *p++ = s11;
                    *p++ = s12;
                    *p++ = s12;
                    *p++ = s22;
                    a = ac1;
                }
                break;
            case 3:
                for ( ; d3 > 0; --d3) {
                    ac1 = &a[m];
                    ac1c = &ac1[m];
                    for (c = m, s11 = s12 = s13 = s22 = s23 = s33 = 0.0; c > 0; --c) {
                        v1 = *a++;
                        v2 = *ac1++;
                        v3 = *ac1c++;
                        s11 += v1 * v1;
                        s12 += v1 * v2;
                        s13 += v1 * v3;
                        s22 += v2 * v2;
                        s23 += v2 * v3;
                        s33 += v3 * v3;
                    }
                    *p++ = s11;
                    *p++ = s12;
                    *p++ = s13;
                    *p++ = s12;
                    *p++ = s22;
                    *p++ = s23;
                    *p++ = s13;
                    *p++ = s23;
                    *p++ = s33;
                    a = ac1c;
                }
                break;
            case 4:
                for ( ; d3 > 0; --d3) {
                    ac1 = &a[m];
                    ac1c = &ac1[m];
                    ac2 = &ac1c[m];
                    for (c = m, s11 = s12 = s13 = s14 = s22 = s23 = s24 = s33 = s34 = s44 = 0.0; c > 0; --c) {
                        v1 = *a++;
                        v2 = *ac1++;
                        v3 = *ac1c++;
                        v4 = *ac2++;
                        s11 += v1 * v1;
                        s12 += v1 * v2;
                        s13 += v1 * v3;
                        s14 += v1 * v4;
                        s22 += v2 * v2;
                        s23 += v2 * v3;
                        s24 += v2 * v4;
                        s33 += v3 * v3;
                        s34 += v3 * v4;
                        s44 += v4 * v4;
                    }
                    *p++ = s11;
                    *p++ = s12;
                    *p++ = s13;
                    *p++ = s14;
                    *p++ = s12;
                    *p++ = s22;
                    *p++ = s23;
                    *p++ = s24;
                    *p++ = s13;
                    *p++ = s23;
                    *p++ = s33;
                    *p++ = s34;
                    *p++ = s14;
                    *p++ = s24;
                    *p++ = s34;
                    *p++ = s44;
                    a = ac2;
                }
                break;
            case 5:
                for ( ; d3 > 0; --d3) {
                    ac1 = &a[m];
                    ac1c = &ac1[m];
                    ac2 = &ac1c[m];
                    ac3 = &ac2[m];
                    for (c = m, s11 = s12 = s13 = s14 = s15 = s22 = s23 = s24 = s25 = s33 = s34 = s35 = s44 = s45 = s55 = 0.0; c > 0; --c) {
                        v1 = *a++;
                        v2 = *ac1++;
                        v3 = *ac1c++;
                        v4 = *ac2++;
                        v5 = *ac3++;
                        s11 += v1 * v1;
                        s12 += v1 * v2;
                        s13 += v1 * v3;
                        s14 += v1 * v4;
                        s15 += v1 * v5;
                        s22 += v2 * v2;
                        s23 += v2 * v3;
                        s24 += v2 * v4;
                        s25 += v2 * v5;
                        s33 += v3 * v3;
                        s34 += v3 * v4;
                        s35 += v3 * v5;
                        s44 += v4 * v4;
                        s45 += v4 * v5;
                        s55 += v5 * v5;
                    }
                    *p++ = s11;
                    *p++ = s12;
                    *p++ = s13;
                    *p++ = s14;
                    *p++ = s15;
                    *p++ = s12;
                    *p++ = s22;
                    *p++ = s23;
                    *p++ = s24;
                    *p++ = s25;
                    *p++ = s13;
                    *p++ = s23;
                    *p++ = s33;
                    *p++ = s34;
                    *p++ = s35;
                    *p++ = s14;
                    *p++ = s24;
                    *p++ = s34;
                    *p++ = s44;
                    *p++ = s45;
                    *p++ = s15;
                    *p++ = s25;
                    *p++ = s35;
                    *p++ = s45;
                    *p++ = s55;
                    a = ac3;
                }
                break;
            default:
                nd = m * n;
                nn = n * n;
                /* loop over third dim */
                for ( ; d3 > 0; --d3) {
                    /* loop over columns (twice) */
                    for (c1 = 0; c1 < n; ++c1) {
                        c2 = c1;
                        mc1 = c1 * m;
                        nc1 = c1 * n;
                        ac1 = &a[mc1];
                        for (c = m, s = 0.0, ac2 = ac1; c > 0; --c) {
                            s += *ac2 * *ac2;
                            ++ac2;
                        }
                        p[nc1 + c1] = s;
                        for (++c2; c2 < n; ++c2) {
                            for (c = m, s = 0.0, ac1c = ac1; c > 0; --c)
                                s += *ac1c++ * *ac2++;
                            p[nc1 + c2] = p[c2 * n + c1] = s;
                        }
                    }
                    a = &a[nd];
                    p = &p[nn];
                }
                break;
        }

    /* two inputs */
    } else {

        /* get data of second input */
        b = (const double *) mxGetData(prhs[1]);

        /* compute advancement counters */
        nd = m * id;
        nn = id * n;

        /* depending on transflag */
        switch (transflag) {

            /* no transpose of either a or b */
            case 0:

                /* iterate over third dim */
                for ( ; d3 > 0; --d3) {

                    /* iterate over n of second argument */
                    for (c2 = 0; c2 < n; ++c2) {

                        /* iterate over m of first argument */
                        for (c1 = 0, ac1c = &b[c2 * id]; c1 < m; ++c1) {

                            /* iterate over inner dim */
                            for (c = id, s = 0.0, ac1 = ac1c, ac2 = &a[c1]; c > 0; --c) {

                                /* add to sum */
                                s += *ac1++ * *ac2;

                                /* increase counter */
                                ac2 = &ac2[m];
                            }

                            /* store in output */
                            *p++ = s;
                        }
                    }

                    /* advance pointers */
                    a = &a[nd];
                    b = &b[nn];
                }
                break;

            /* transpose of a, but not b */
            case 1:

                /* iterate over third dim */
                for ( ; d3 > 0; --d3) {

                    /* iterate over n of second argument */
                    for (c2 = 0; c2 < n; ++c2) {

                        /* iterate over m of first argument */
                        for (c1 = 0, ac1c = &b[c2 * id]; c1 < m; ++c1) {

                            /* iterate over inner dim */
                            for (c = id, s = 0.0, ac1 = ac1c, ac2 = &a[c1 * id]; c > 0; --c)

                                /* add to sum */
                                s += *ac1++ * *ac2++;

                            /* store in output */
                            *p++ = s;
                        }
                    }

                    /* advance pointers */
                    a = &a[nd];
                    b = &b[nn];
                }
                break;

            /* transpose of b, but not a */
            case 2:

                /* iterate over third dim */
                for ( ; d3 > 0; --d3) {

                    /* iterate over n of second argument */
                    for (c2 = 0; c2 < n; ++c2) {

                        /* iterate over m of first argument */
                        for (c1 = 0, ac1c = &b[c2]; c1 < m; ++c1) {

                            /* iterate over inner dim */
                            for (c = id, s = 0.0, ac1 = ac1c, ac2 = &a[c1]; c > 0; --c) {

                                /* add to sum */
                                s += *ac1 * *ac2;

                                /* increase counters */
                                ac1 = &ac1[n];
                                ac2 = &ac2[m];
                            }

                            /* store in output */
                            *p++ = s;
                        }
                    }

                    /* advance pointers */
                    a = &a[nd];
                    b = &b[nn];
                }
                break;

            /* transpose of a and b */
            case 3:

                /* iterate over third dim */
                for ( ; d3 > 0; --d3) {

                    /* iterate over n of second argument */
                    for (c2 = 0; c2 < n; ++c2) {

                        /* iterate over m of first argument */
                        for (c1 = 0, ac1c = &b[c2]; c1 < m; ++c1) {

                            /* iterate over inner dim */
                            for (c = id, s = 0.0, ac1 = ac1c, ac2 = &a[c1 * id]; c > 0; --c) {

                                /* add to sum */
                                s += *ac1 * *ac2++;

                                /* increase counters */
                                ac1 = &ac1[n];
                            }

                            /* store in output */
                            *p++ = s;
                        }
                    }

                    /* advance pointers */
                    a = &a[nd];
                    b = &b[nn];
                }
                break;
        }
    }
}

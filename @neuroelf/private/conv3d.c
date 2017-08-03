/*

convolve in 3D, see conv3d.m

FORMAT:       v = conv3d(v, k [, o, t])

Input fields:

      v           3D volume (logical or numeric)
      k           3D convolution kernel (must be double)
                  -  special 1x1 kernel values:
                     0 - only apply == threshold (values ~= thresh := 0)
                     1 - binary erode (full erosion!)
                     2 - simple binary dilate (only main directions, faces)
                     3 - complex binary dilate (edge neighbors)
                     4 - full binary dilate (vertex neighbors)
                         if input array is not logical, threshold must be
                         given, by default uses == threshold, add 64 to use
                         values >= thresh and -64 to use values <= thresh
                    32 - edge enhancer
                    48 - gradient magnitude
                    49 - mean gradient cosine/colinearity in 3x3x3 cube
                  - special 1x2 kernel values:
                    [fwhm, func] - build kernel with FWHM and func
                    func 0: gaussian
                    func 1: linear
                  - special 1x3 kernel values:
                    [ g1 ,  g2 ,  g3 ] - free spatial gradient
                    [+/-1,  0  ,  0  ] - first dim gradient
                    [ 0  , +/-1,  0  ] - second dim gradient
                    [ 0  ,  0  , +/-1] - third dim gradient
                  - special 1x4 kernel values:
                    [k1, k2, k3, func] - build kernel with FWHM
      o           output type (1x1 double, if not given or 0, as input)
                     1 - logical volume
      t           optional threshold value (see above)

Output fields:

      v           convolved and/or thresholded/computed volume

Note: for 1x1 kernel arguments 0 to 4 , a logical array with the
      requested threshold is produced

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

#include "conv3d.h"

/* create copy-from functions */
CONV3D_COPY_FPFD_TYPE_S( int8 , char  )
CONV3D_COPY_FPFD_TYPE_S( int16, short )
CONV3D_COPY_FPFD_TYPE_S( int32, int   )
CONV3D_COPY_FPFD_TYPE_U(uint8 , char  )
CONV3D_COPY_FPFD_TYPE_U(uint16, short )
CONV3D_COPY_FPFD_TYPE_U(uint32, int   )
CONV3D_COPY_FPFD_TYPE_F(single, float )
CONV3D_COPY_FPFD_TYPE_F(double, double)


/* here comes the main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* output type, type of call, type of threshold, and threshold */
    double otyped = 0.0;
    unsigned char otype = 0;
    unsigned char tot = 0;
    double othresh = 1.0;

    /* image arg mxArray class, and content pointers and dims */
    const mxArray *ia = NULL;
    const int *idim = NULL;
    mxClassID iclass;

    /* kernel arg mxArray class, and content pointers and dims and value */
    const mxArray *ka = NULL;
    unsigned long knd = 0;
    unsigned long kne = 0;
    const int *kdim = NULL;
    int kdims[3] = {1, 1, 1};
    const double *kad = NULL;
    double kdv = 0.0;
    signed char ksv = 0;

    /* self-created kernel array */
    int kdimx = 0, kdimy = 0, kdimz = 0;
    const unsigned char *ksfau8 = NULL;

    /* summation array(s) */
    mxArray *suma = NULL;
    unsigned long sumoffs[3] = {1, 1, 1};
    unsigned long soffx = 0, soffy = 0, soffz = 0;
    int sdims[3], sdimx = 0, sdimy = 0, sdimz = 0;
    unsigned char *sumau8 = NULL, *sumpu8 = NULL;
    double *sumad = NULL, *sumpd = NULL;

    /* temp array(s) */
    mxArray *tmpa = NULL;
    unsigned char *tmpau8 = NULL, *tmppu8 = NULL;
    double *tmpad = NULL, *tmppd = NULL;
    unsigned long netmp = 0;

    /* number of elements for fast convolution copy/add */
    unsigned long necpy = 0, sumpos = 0;

    /* convolution offsets (into sum/temp array) */
    signed long offxy = 0, offy = 0, offz = 0, offtmp = 0;

    /* kernel position */
    signed long krnx = 0, krny = 0, krnz = 0;

    /* output array */
    mxArray *outa = NULL;
    int odims[3] = {1, 1, 1};
    void *outav = NULL;
    unsigned char *outau8 = NULL;
    unsigned long neout = 0, necnt = 0;

    /* erosion/test threshold */
    unsigned char ertu8 = 0;
    float erts = 0.0;

    /* flexmask pointers */
    const mxArray *flexin[6];
    mxArray *flexout[2];

    /* in case we need to create arrays from scratch */
    int cadim[3] = {1, 1, 1};

    /* isinfnan.h variables */
    VARS_FOR_ISINFNAN

    /* variable output string */
    /* char vstr[256]; */


    /* the code starts here */


    /* check number of in/out arguments */
	if (nlhs > 1 ||
        nrhs  < 2 ||
        nrhs  > 5)
		mexErrMsgTxt(conv3d_err_numberofargs);


    /* check first argument */
    ia = prhs[0];
    if ((!mxIsLogical(ia) &&
         !mxIsNumeric(ia)) ||
        ( mxGetNumberOfElements(ia) == 0) ||
        ( mxGetNumberOfDimensions(ia) > 3))
        mexErrMsgTxt(conv3d_err_badfirstarg);

    /* get dims, etc. */
    idim = mxGetDimensions(ia);
    odims[0] = idim[0];
    odims[1] = idim[1];
    if (mxGetNumberOfDimensions(ia) > 2)
        odims[2] = idim[2];
    neout = odims[0] * odims[1] * odims[2];
    iclass = mxGetClassID(ia);


    /* check second argument */
    ka = prhs[1];
    if (!mxIsDouble(ka) ||
        (mxGetNumberOfElements(ka) == 0) ||
        (mxGetNumberOfDimensions(ka) > 3))
        mexErrMsgTxt(conv3d_err_badsecondarg);

    /* get dims, etc. */
    kdim = mxGetDimensions(ka);
    knd = mxGetNumberOfDimensions(ka);
    kne = mxGetNumberOfElements(ka);
    kad = (const double*) mxGetPr(ka);
    kdims[0] = kdim[0];
    kdims[1] = kdim[1];
    if (knd > 2)
        kdims[2] = kdim[2];


    /* initialize IS_* macros */
    INIT_INF_NAN_BAD_VAL()


    /* check optional argument output type (type and size) */
    if (nrhs > 2) {
        if (!mxIsDouble(prhs[2]) ||
            (mxGetNumberOfElements(prhs[2]) != 1))
            mexErrMsgTxt(conv3d_err_invalidarg_otype);
        otyped = *(const double*) mxGetPr(prhs[2]);
        if (mxIsInf(otyped) ||
            mxIsNaN(otyped) ||
           ((otyped < 0.0) ||
            (otyped > 3.0)))
            mexErrMsgTxt(conv3d_err_invalidarg_otype);

        /* which output requested */
        otype = (unsigned char) otyped;
    }


    /* check optional argument threshold */
    if (nrhs > 3) {
        if (!mxIsDouble(prhs[3]) ||
            (mxGetNumberOfElements(prhs[3]) != 1))
            mexErrMsgTxt(conv3d_err_invalidarg_thresh);

        /* get threshold value */
        othresh = *(const double*) mxGetPr(prhs[3]);
        if (mxIsInf(othresh) || mxIsNaN(othresh))
            mexErrMsgTxt(conv3d_err_invalidarg_thresh);
    }


    /* check class and threshold value */
    switch (iclass) {
        case mxINT8_CLASS:
            if (othresh < -128.0 || othresh >= 127.5)
                mexErrMsgTxt(conv3d_err_badthreshvalue);
            break;
        case mxINT16_CLASS:
            if (othresh < -32768.0 || othresh >= 32767.5)
                mexErrMsgTxt(conv3d_err_badthreshvalue);
            break;
        case mxINT32_CLASS:
            if (othresh < -2147483648.0 || othresh >= 2147483647.5)
                mexErrMsgTxt(conv3d_err_badthreshvalue);
            break;
        case mxLOGICAL_CLASS:
            if (othresh < 0.0 || othresh > 1.0)
                mexErrMsgTxt(conv3d_err_badthreshvalue);
            break;
        case mxUINT8_CLASS:
            if (othresh < 0.0 || othresh >= 255.5)
                mexErrMsgTxt(conv3d_err_badthreshvalue);
            break;
        case mxUINT16_CLASS:
            if (othresh < 0.0 || othresh >= 65535.5)
                mexErrMsgTxt(conv3d_err_badthreshvalue);
            break;
        case mxUINT32_CLASS:
            if (othresh < 0.0 || othresh > 4294967295.5)
                mexErrMsgTxt(conv3d_err_badthreshvalue);
            break;
        case mxSINGLE_CLASS:
            erts = (float) othresh;
            if (mxIsInf((double) erts) ||
                mxIsNaN((double) erts))
                mexErrMsgTxt(conv3d_err_badthreshvalue);
            break;
        case mxDOUBLE_CLASS:
            break;
        default:
            mexErrMsgTxt(conv3d_err_badinputclass);
    }


    /* checks completed, let's rumble... */


    /* type of kernel argument decides what to do */
    /* multiply with size(1) and (ndims-1) to find out special value */
    switch (kne) {


        /* 1x1 special value (erosion, dilation, gradient, ...) */
        case 1:

            /* get scalar kernel value (from double) */
            kdv = *kad;
            if (mxIsInf(kdv) ||
                mxIsNaN(kdv) ||
               (kdv < -64.0) ||
               (kdv > 127.0))
                mexErrMsgTxt(conv3d_err_badkernelvalue);
            ksv = (signed char) kdv;

            /* default thresholding (if at all) is == threshold */
            tot = CONV3D_TOT_EQ;

            /* less/equal thresholding */
            if (ksv < 0) {
                tot = CONV3D_TOT_LE;
                ksv += 64;

            /* greater/equal thresholding */
            } else if (ksv > 63) {
                tot = CONV3D_TOT_GE;
                ksv -= 64;
            }

            /* force boolean for erosion/dilation */
            if (((ksv > 0) &&
                 (ksv < 5)) ||
                (iclass == mxLOGICAL_CLASS))
                otype = 1;

            /* copy 1. arg and prepare 3 new args for call to flexin */
            flexin[0] = prhs[0];
            flexin[1] = mxCreateNumericArray(2, cadim, mxDOUBLE_CLASS, mxREAL);
            if (flexin[1] == NULL)
                mexErrMsgTxt(conv3d_err_createarrayfailed);
            flexin[2] = mxCreateNumericArray(2, cadim, mxDOUBLE_CLASS, mxREAL);
            if (flexin[2] == NULL)
                mexErrMsgTxt(conv3d_err_createarrayfailed);
            flexin[3] = mxCreateNumericArray(2, cadim, mxDOUBLE_CLASS, mxREAL);
            if (flexin[3] == NULL)
                mexErrMsgTxt(conv3d_err_createarrayfailed);

            /* second arg: thresh */
            *((double *) mxGetPr(flexin[1])) = othresh;

            /* third arg masking type */
            switch (tot) {
                case CONV3D_TOT_LE:
                    *((double *) mxGetPr(flexin[2])) = -1.0;
                    break;
                case CONV3D_TOT_GE:
                    *((double *) mxGetPr(flexin[2])) =  1.0;
                    break;
            }

            /* special case, do not apply == mask on special task */
            if ((ksv > 4) &&
                (tot == CONV3D_TOT_EQ))
                    *((double *) mxGetPr(flexin[2])) =  2.0;

            /* for all simple tasks set to boolean */
            if (ksv < 5)
                *((double *) mxGetPr(flexin[3])) = 1.0;

            /* for all other tasks to double */
            else
                *((double *) mxGetPr(flexin[3])) = 2.0;

            /* only masking/copying is requested */
            if (ksv == 0) {

                /* call flexmask */
                mexCallMATLAB(1, plhs, 4, (mxArray**) ((void**) flexin), "flexmask");

                /* destroy arrays */
                mxDestroyArray((mxArray*) ((void*) flexin[1]));
                mxDestroyArray((mxArray*) ((void*) flexin[2]));
                mxDestroyArray((mxArray*) ((void*) flexin[3]));

                /* return */
                return;
            }

            /* erosion/dilation requested */
            if (ksv < 5) {

                /* we must increase the dimensions by two in each direction */
                sdims[0] = odims[0] + 2;
                sdims[1] = odims[1] + 2;
                sdims[2] = odims[2] + 2;
                sdimx = *sdims;
                sdimy = sdims[1];
                sdimz = sdims[2];

                /* we need 5. arg, so further preparations */
                cadim[1] = 3;
                flexin[4] = mxCreateNumericArray(2, cadim, mxDOUBLE_CLASS, mxREAL);
                if (flexin[4] == NULL)
                    mexErrMsgTxt(conv3d_err_createarrayfailed);

                /* add two to the dimension of the data */
                tmpad = (double *) mxGetPr(flexin[4]);
                tmpad[0] = (double) sdimx;
                tmpad[1] = (double) sdimy;
                tmpad[2] = (double) sdimz;

                /* call flexmask */
                mexCallMATLAB(2, flexout, 5, (mxArray**) ((void**) flexin), "flexmask");

                /* destroy temporary arrays */
                mxDestroyArray((mxArray*) ((void*) flexin[1]));
                mxDestroyArray((mxArray*) ((void*) flexin[2]));
                mxDestroyArray((mxArray*) ((void*) flexin[3]));
                mxDestroyArray((mxArray*) ((void*) flexin[4]));

                /* get pointer */
                tmpau8 = (unsigned char*) mxGetData(flexout[0]);

                /* default threshold is 1 */
                ertu8 = 1;

                /* if number not set, get number */
                necpy = (unsigned long) *((double *) mxGetPr(flexout[1]));
                if (necpy == 0) {
                    tmppu8 = tmpau8;
                    for (necnt = neout; necnt > 0; --necnt)
                        if (*tmppu8++ > 0)
                            ++necpy;
                }

                /* create kernel and summation array */
                CONV3D_CREATEARRAY(UINT8, sdims, suma, unsigned char, sumau8)
                CONV3D_CREATEARRAY(LOGICAL, odims, outa, unsigned char, outau8)
                plhs[0] = outa;

                /* fill kernel array */
                switch (ksv) {

                    /* erode */
                    case 1:
                        ksfau8 = conv3d_krn_fullcb;
                        ertu8 = 27;
                        break;

                    /* dilate (along faces) */
                    case 2:
                        ksfau8 = conv3d_krn_dilate1;
                        break;

                    /* dilate (along faces and edges) */
                    case 3:
                        ksfau8 = conv3d_krn_dilate2;
                        break;

                    /* dilate (full 3x3x3 cube) */
                    case 4:
                        ksfau8 = conv3d_krn_fullcb;
                        break;
                }

                /* get offsets into array */
                offxy = sdimx * sdimy;
                offtmp = offxy + sdimx + 1;
                necnt = offxy * sdimz - 2 * offtmp;

                /* set pointers further */
                tmppu8 = &tmpau8[offtmp];
                sumpu8 = &sumau8[offtmp];

                /* cheap convolution (swapping the tmp/kernel logic) */
                if ((( (double) neout / (double) necpy ) > 6.0) ||
                    (ksv == 2)) {

                    /* for faces dilation, even cheaper */
                    if (ksv == 2) {

                        /* iterate over data and set at specified positions */
                        for ( ; necnt > 0; --necnt) {
                            if (*tmppu8++ > 0) {
                                sumpu8[-offxy] = 1;
                                sumpu8[-sdimx] = 1;
                                sumpu8[-1] = 1;
                               *sumpu8 = 1;
                                sumpu8[1] = 1;
                                sumpu8[sdimx] = 1;
                                sumpu8[offxy] = 1;
                            }
                            ++sumpu8;
                        }

                    /* another cheap solution for full kernel */
                    } else if (ksv != 3) {

                        /* iterate over data and set 3x3x3 instead */
                        for ( ; necnt > 0; --necnt) {
                            if (*tmppu8++ > 0) {
                                for (krnz = -offxy; krnz <= offxy; krnz += offxy) {
                                    for (krny = -sdimx; krny <= sdimx; krny += sdimx) {
                                        offy = krnz + krny - 2;
                                        ++sumpu8[++offy];
                                        ++sumpu8[++offy];
                                        ++sumpu8[++offy];
                                    }
                                }
                            }
                            ++sumpu8;
                        }

                    /* special convolution (must be case 3 !) */
                    } else {

                        /* iterate over data and apply kernel instead */
                        for ( ; necnt > 0; --necnt) {
                            if (*tmppu8++ > 0) {
                                ksfau8 = conv3d_krn_dilate2;
                                for (krnz = -offxy; krnz <= offxy; krnz += offxy) {
                                    for (krny = -sdimx; krny <= sdimx; krny += sdimx) {
                                        offy = krnz + krny;
                                        for (krnx = -1; krnx <= 1; ++krnx) {
                                            if (*ksfau8++ > 0)
                                                ++sumpu8[offy + krnx];
                                        }
                                    }
                                }
                            }
                            ++sumpu8;
                        }
                    }

                /* full convolution */
                } else {

                    /* iterate over kernel dimensions */
                    for (krnz = 0; krnz < 3; ++krnz) {
                        offz = offxy * krnz;
                        for (krny = 0; krny < 3; ++krny) {
                            offy = offz + sdimx * krny;
                            for (krnx = 0; krnx < 3; ++krnx) {
                                if (*ksfau8++ > 0) {
                                    tmppu8 = &tmpau8[offtmp];
                                    sumpu8 = &sumau8[offy + krnx];
                                    for (sumpos = necnt; sumpos > 0; --sumpos)
                                        *sumpu8++ += *tmppu8++;
                                }
                            }
                        }
                    }
                }

                /* set 1 in output from sum where > ert */
                conv3d_set1ge_frompart_uint8(sdims, sumau8, odims, outau8, sumoffs, ertu8);

                /* destroy temporary arrays and return early */
                mxDestroyArray(suma);
                mxDestroyArray(tmpa);
                mxDestroyArray(flexout[1]);
                mxDestroyArray(flexout[0]);
                return;
            }
            mexErrMsgTxt(conv3d_err_notyetimplemented);

            /* other special meanings here ... */
            switch (ksv) {

                /* edge enhancement */
                case 32:
                    break;

                /* gradient magnitude */
                case 48:
                    break;

                /* mean colinearity */
                case 49:
                    break;

                /* reject other cases */
                default:
                    mexErrMsgTxt(conv3d_err_invalidarg_kernel);

            }
            break;

        /* uniform smoothing */
        case 2:
            mexErrMsgTxt(conv3d_err_notyetimplemented);
            return;
            break;

        /* variable smoothing, non uniform edge enhancing, */
        /* more gradient stuff etc. */
        case 4:
            mexErrMsgTxt(conv3d_err_notyetimplemented);
            return;
            break;
    }

    /* the kernel argument *must* be odd in each dim! */
    if (((kdims[0] % 2) != 1) ||
        ((kdims[1] % 2) != 1) ||
        ((kdims[2] % 2) != 1))
        mexErrMsgTxt(conv3d_err_invalidarg_kernel);

    /* put input array in double array (use flexmask) */
    /* copy 1. arg and prepare 4 new args for call */
    flexin[0] = prhs[0];
    flexin[1] = mxCreateNumericArray(2, cadim, mxDOUBLE_CLASS, mxREAL);
    if (flexin[1] == NULL)
        mexErrMsgTxt(conv3d_err_createarrayfailed);
    flexin[2] = mxCreateNumericArray(2, cadim, mxDOUBLE_CLASS, mxREAL);
    if (flexin[2] == NULL)
        mexErrMsgTxt(conv3d_err_createarrayfailed);
    flexin[3] = mxCreateNumericArray(2, cadim, mxDOUBLE_CLASS, mxREAL);
    if (flexin[3] == NULL)
        mexErrMsgTxt(conv3d_err_createarrayfailed);
    cadim[1] = 3;
    flexin[4] = mxCreateNumericArray(2, cadim, mxDOUBLE_CLASS, mxREAL);
    if (flexin[4] == NULL)
        mexErrMsgTxt(conv3d_err_createarrayfailed);

    /* masking is off */
    *((double *) mxGetPr(flexin[2])) = 2.0;

    /* force double */
    *((double *) mxGetPr(flexin[3])) = 2.0;

    /* we must increase the dimensions by two in each direction */
    kdimx = kdims[0];
    kdimy = kdims[1];
    kdimz = kdims[2];
    soffx = kdimx - 1;
    soffy = kdimy - 1;
    soffz = kdimz - 1;
    sumoffs[0] = soffx;
    sumoffs[1] = soffy;
    sumoffs[2] = soffz;
    sdims[0] = odims[0] + 2 * kdimx - 2;
    sdims[1] = odims[1] + 2 * kdimy - 2;
    sdims[2] = odims[2] + 2 * kdimz - 2;
    sdimx = *sdims;
    sdimy = sdims[1];
    sdimz = sdims[2];

    /* add two to the dimension of the data */
    tmpad = (double *) mxGetPr(flexin[4]);
    tmpad[0] = (double) sdimx;
    tmpad[1] = (double) sdims[1];
    tmpad[2] = (double) sdims[2];

    /* call flexmask and get pointer to data */
    mexCallMATLAB(1, flexout, 5, (mxArray**) ((void**) flexin), "flexmask");
    mxDestroyArray((mxArray*) ((void*) flexin[1]));
    mxDestroyArray((mxArray*) ((void*) flexin[2]));
    mxDestroyArray((mxArray*) ((void*) flexin[3]));
    mxDestroyArray((mxArray*) ((void*) flexin[4]));
    tmpad = (double *) mxGetPr(flexout[0]);

    /* create summation array */
    CONV3D_CREATEARRAY(DOUBLE, sdims, suma, double, sumad)

    /* get offsets into array */
    offxy = sdimx * sdimy;
    offtmp = soffx * offxy + soffy * sdimx + soffz;
    netmp = offxy * sdimz;
    necnt = netmp - 2 * offtmp;

    /* expand 6 edge slices */
    for (krny = netmp - sdimx; krny >= 0; krny -= sdimx) {
        tmppd = &tmpad[krny + soffx - 1];
        sumpd = &tmpad[krny + sdimx - soffx];
        for (krnx = soffx; krnx > 0; --krnx) {
            *tmppd = tmppd[1];
            --tmppd;
            *sumpd = sumpd[-1];
            ++sumpd;
        }
    }
    for (krnz = netmp - offxy; krnz >= 0; krnz -= offxy) {
        for (krnx = sdimx-1; krnx >= 0; --krnx) {
            tmppd = &tmpad[krnz + krnx + soffy * sdimx - sdimx];
            sumpd = &tmpad[krnz + krnx + offxy - soffy * sdimx];
            for (krny = soffy; krny > 0; --krny) {
                *tmppd = tmppd[sdimx];
                tmppd = &tmppd[-sdimx];
                *sumpd = sumpd[-sdimx];
                sumpd = &sumpd[sdimx];
            }
        }
    }
    for (krnz = (soffz - 1) * offxy; krnz >= 0; krnz -= offxy) {
        tmppd = &tmpad[krnz];
        sumpd = &tmpad[krnz+offxy];
        for (krnx = offxy - 1; krnx >= 0; --krnx)
            *tmppd++ = *sumpd++;
    }
    for (krnz = (sdimz - soffz) * offxy; krnz < netmp; krnz += offxy) {
        tmppd = &tmpad[krnz];
        sumpd = &tmpad[krnz-offxy];
        for (krnx = offxy - 1; krnx >= 0; --krnx)
            *tmppd++ = *sumpd++;
    }


    /* get offsets into array */
    soffx = (unsigned long) (((double) soffx) / 2.0);
    soffy = (unsigned long) (((double) soffy) / 2.0);
    soffz = (unsigned long) (((double) soffz) / 2.0);
    offtmp = soffx * offxy + soffy * sdimx + soffz;
    necnt = netmp - 2 * offtmp;

    /* loop over kernel elements */
    for (krnz = 0; krnz < kdimz; ++krnz) {
        offz = offxy * krnz;
        for (krny = 0; krny < kdimy; ++krny) {
            offy = offz + sdimx * krny;
            for (krnx = 0; krnx < kdimx; ++krnx) {
                kdv = *kad++;
                if (kdv != 0.0) {
                    tmppd = &tmpad[offtmp];
                    sumpd = &sumad[offy + krnx];
                    for (sumpos = necnt; sumpos > 0; --sumpos) {
                        otyped = *tmppd++;
                        IF_IS_GOOD_VAL(otyped) {
                            *sumpd++ += otyped * kdv;
                        } else {
                            ++sumpd;
                        }
                    }
                }
            }
        }
    }

    /* delete temporary array */
    mxDestroyArray(flexout[0]);

    /* what kind of output requested */
    switch (otype) {

        /* as input */
        case 0:
            plhs[0] = mxCreateNumericArray(3, odims, iclass, mxREAL);
            break;

        /* thresholded */
        case 1:
            plhs[0] = mxCreateNumericArray(3, odims, mxLOGICAL_CLASS, mxREAL);
            break;

        /* force to double */
        case 2:
            plhs[0] = mxCreateNumericArray(3, odims, mxDOUBLE_CLASS, mxREAL);
            break;

        /* force to single */
        case 3:
            plhs[0] = mxCreateNumericArray(3, odims, mxSINGLE_CLASS, mxREAL);
            break;
    }
    /* check output array */
    if (plhs[0] == NULL)
        mexErrMsgTxt(conv3d_err_createarrayfailed);
    outav = (void *) mxGetData(plhs[0]);

    /* depends on class */
    switch (mxGetClassID(plhs[0])) {
        case mxLOGICAL_CLASS : conv3d_set1ge_frompart_double(      sdims, sumad, odims, (unsigned char   *) outav, sumoffs, othresh); break;
        case mxINT8_CLASS    : conv3d_copy_frompart_fromdbl_int8(  sdims, sumad, odims, (  signed char   *) outav, sumoffs);          break;
        case mxINT16_CLASS   : conv3d_copy_frompart_fromdbl_int16( sdims, sumad, odims, (  signed short  *) outav, sumoffs);          break;
        case mxINT32_CLASS   : conv3d_copy_frompart_fromdbl_int32( sdims, sumad, odims, (  signed int    *) outav, sumoffs);          break;
        case mxUINT8_CLASS   : conv3d_copy_frompart_fromdbl_uint8( sdims, sumad, odims, (unsigned char   *) outav, sumoffs);          break;
        case mxUINT16_CLASS  : conv3d_copy_frompart_fromdbl_uint16(sdims, sumad, odims, (unsigned short  *) outav, sumoffs);          break;
        case mxUINT32_CLASS  : conv3d_copy_frompart_fromdbl_uint32(sdims, sumad, odims, (unsigned int    *) outav, sumoffs);          break;
        case mxSINGLE_CLASS  : conv3d_copy_frompart_fromdbl_single(sdims, sumad, odims, (         float  *) outav, sumoffs);          break;
        case mxDOUBLE_CLASS  : conv3d_copy_frompart_fromdbl_double(sdims, sumad, odims, (         double *) outav, sumoffs);          break;
        default: /* already handled by argument check at the beginning */ break;
    }

    /* delete temp array */
    mxDestroyArray(suma);
}

unsigned long conv3d_set1ge_frompart_uint8(const int *snd, const unsigned char *sd, const int *tnd, unsigned char *td, const unsigned long *off, unsigned char t)
{
    unsigned long xc, yc, zc, xy, x1, y1, z1, x2, y2, z2, to, yo, zo;
    x1 = *snd++;
    y1 = *snd;
    xy = x1 * y1;
    z1 = 0;
    x2 = *tnd++;
    y2 = *tnd++;
    z2 = *tnd;
    to = off[0] + off[1] * x1 + off[2] * xy;
    for (zc = 0; zc < z2; ++zc) {
        zo = to + zc * xy;
        for (yc = 0; yc < y2; ++yc) {
            yo = zo + yc * x1;
            for (xc = 0; xc < x2; ++xc) {
                if (sd[yo++] >= t) {
                    *td = 1;
                    ++z1;
                }
                ++td;
            }
        }
    }
    return z1;
}
void conv3d_set1ge_frompart_double(const int *snd, const double *sd, const int *tnd, unsigned char *td, const unsigned long *off, double t)
{
    unsigned long xc, yc, zc, xy, x1, y1, x2, y2, z2, to, yo, zo;
    x1 = *snd++;
    y1 = *snd;
    xy = x1 * y1;
    x2 = *tnd++;
    y2 = *tnd++;
    z2 = *tnd;
    to = off[0] + off[1] * x1 + off[2] * xy;
    for (zc = 0; zc < z2; ++zc) {
        zo = to + zc * xy;
        for (yc = 0; yc < y2; ++yc) {
            yo = zo + yc * x1 - 1;
            for (xc = 0; xc < x2; ++xc) {
                if (sd[++yo] >= t)
                    *td = 1;
                ++td;
            }
        }
    }
}

/*
  conv3d.h

% Version:  v0.9b
% Build:    11043013
% Date:     Sep-01 2010, 11:28 AM EST
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
#include "isinfnan.h"

#define CONV3D_TOT_EQ 1
#define CONV3D_TOT_GE 2
#define CONV3D_TOT_LE 3

/* macros to define a copy frompart function */
#define CONV3D_COPY_FPFD_TYPE_F(MLtype, Ctype) \
void conv3d_copy_frompart_fromdbl_##MLtype(const int *snd, const double *sd, const int *tnd, Ctype *td, const unsigned long *off) \
{ unsigned long xc, yc, zc, xy, x1, y1, x2, y2, z2, to, yo, zo; x1 = *snd++; y1 = *snd; xy = x1 * y1; x2 = *tnd++; y2 = *tnd++; z2 = *tnd; \
  to = off[0] + off[1] * x1 + off[2] * xy - 1; \
  for (zc = 0; zc < z2; ++zc) { zo = to + zc * xy; for (yc = 0; yc < y2; ++yc) { yo = zo + yc * x1; for (xc = 0; xc < x2; ++xc) { \
    *td++ = (Ctype) sd[++yo]; \
  } } } \
}
#define CONV3D_COPY_FPFD_TYPE_S(MLtype, Ctype) \
void conv3d_copy_frompart_fromdbl_##MLtype(const int *snd, const double *sd, const int *tnd, signed Ctype *td, const unsigned long *off) \
{ unsigned long xc, yc, zc, xy, x1, y1, x2, y2, z2, to, yo, zo; double sdb = 0.0; x1 = *snd++; y1 = *snd; xy = x1 * y1; x2 = *tnd++; y2 = *tnd++; z2 = *tnd; \
  to = off[0] + off[1] * x1 + off[2] * xy - 1; \
  for (zc = 0; zc < z2; ++zc) { zo = to + zc * xy; for (yc = 0; yc < y2; ++yc) { yo = zo + yc * x1; for (xc = 0; xc < x2; ++xc) { \
    sdb = sd[++yo]; *td++ = (signed Ctype) ((sdb > 0.0) ? (sdb + 0.5) : (sdb - 0.5)); \
  } } } \
}
#define CONV3D_COPY_FPFD_TYPE_U(MLtype, Ctype) \
void conv3d_copy_frompart_fromdbl_##MLtype(const int *snd, const double *sd, const int *tnd, unsigned Ctype *td, const unsigned long *off) \
{ unsigned long xc, yc, zc, xy, x1, y1, x2, y2, z2, to, yo, zo; x1 = *snd++; y1 = *snd; xy = x1 * y1; x2 = *tnd++; y2 = *tnd++; z2 = *tnd; \
  to = off[0] + off[1] * x1 + off[2] * xy - 1; \
  for (zc = 0; zc < z2; ++zc) { zo = to + zc * xy; for (yc = 0; yc < y2; ++yc) { yo = zo + yc * x1; for (xc = 0; xc < x2; ++xc) { \
    *td++ = (unsigned Ctype) (sd[++yo] + 0.5); \
  } } } \
}

/* macro to create an mxArray and assign/check array pointer(s) */
#define CONV3D_CREATEARRAY(ATYPE, ADIM, AVAR, PTYPE, PVAR) \
    AVAR = mxCreateNumericArray(3, ADIM, mx##ATYPE##_CLASS, mxREAL); \
    if ( AVAR == NULL ) \
        mexErrMsgTxt("Error creating array."); \
	PVAR = ( PTYPE *) mxGetData( AVAR ); \
    if ( PVAR == NULL ) \
        mexErrMsgTxt("Error getting pointer to array.");

/* error messages */
static const char *conv3d_err_badfirstarg       = "Bad first input argument.";
static const char *conv3d_err_badinputclass     = "Bad input class.";
static const char *conv3d_err_badkernelvalue    = "Bad kernel value.";
static const char *conv3d_err_badsecondarg      = "Bad second argument.";
static const char *conv3d_err_badthreshvalue    = "Invalid threshold value.";
static const char *conv3d_err_createarrayfailed = "Error creating temporary array.";
static const char *conv3d_err_invalidarg_kernel = "Bad special kernel argument.";
static const char *conv3d_err_invalidarg_otype  = "Invalid argument: output type.";
static const char *conv3d_err_invalidarg_thresh = "Invalid argument: threshold.";
static const char *conv3d_err_numberofargs      = "Bad number of input/output arguments.";
static const char *conv3d_err_notyetimplemented = "Functionality not yet implemented.";

/* erosion, dilation kernels */
static const unsigned char conv3d_krn_dilate1[27] =
    {0, 0, 0,
     0, 1, 0,
     0, 0, 0,

     0, 1, 0,
     1, 1, 1,
     0, 1, 0,

     0, 0, 0,
     0, 1, 0,
     0, 0, 0};
static const unsigned char conv3d_krn_dilate2[27] =
    {0, 1, 0,
     1, 1, 1,
     0, 1, 0,

     1, 1, 1,
     1, 1, 1,
     1, 1, 1,

     0, 1, 0,
     1, 1, 1,
     0, 1, 0};
static const unsigned char conv3d_krn_fullcb[27] =
    {1, 1, 1,
     1, 1, 1,
     1, 1, 1,

     1, 1, 1,
     1, 1, 1,
     1, 1, 1,

     1, 1, 1,
     1, 1, 1,
     1, 1, 1};

/* function for binary erosion/dilation (compare to flexmask.h) */
unsigned long conv3d_set1ge_frompart_uint8(const int *snd, const unsigned char *sd, const int *tnd, unsigned char *td, const unsigned long *off, unsigned char t);
void conv3d_set1ge_frompart_double(const int *snd, const double *sd, const int *tnd, unsigned char *td, const unsigned long *off, double t);

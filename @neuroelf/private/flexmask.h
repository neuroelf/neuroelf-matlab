/*

flexmask.h (used by flexmask.c)

% Version:  v0.9a
% Build:    11102110
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

/* include Matlab stuff */
#include "mex.h"


/* function calls are defined as macros for simplification */


/* copy loop, set 1 */
#define FLXMSK_COPYLOOP1_OPEN \
    for (zc = 0; zc < z2; ++zc) { \
        zo = to + zc * xy; \
        for (yc = 0; yc < y2; ++yc) { \
            yo = zo + yc * x1; \
            for (xc = 0; xc < x2; ++xc) {

/* copy loop, set 1 */
#define FLXMSK_COPYLOOP2_OPEN \
    for (zc = 0; zc < z1; ++zc) { \
        zo = to + zc * xy; \
        for (yc = 0; yc < y1; ++yc) { \
            yo = zo + yc * x2; \
            for (xc = 0; xc < x1; ++xc) {

/* copy loop close */
#define FLXMSK_COPYLOOP_CLOSE }}}
#define FLXMSK_COPYLOOP1_CLOSE }}} return z1;
#define FLXMSK_COPYLOOP2_CLOSE }}} return z2;

/* common variables for COPY/SET */
#define FLXMSK_COPYVARS(SAD, TAD) \
    unsigned long xc, yc, zc, xy, x1, y1, z1, x2, y2, z2, to, yo, zo; \
    x1 = SAD[0]; y1 = SAD[1]; z1 = SAD[2]; \
    x2 = TAD[0]; y2 = TAD[1]; z2 = TAD[2];

/* variables for COPY/SET macros, set 1 */
#define FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
    FLXMSK_COPYVARS(SAD, TAD) \
    z1 = 0; \
    xy = x1 * y1; \
    to = OZ * xy + OY * x1 + OX;

/* variables for COPY/SET macros, set 2 */
#define FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
    FLXMSK_COPYVARS(SAD, TAD) \
    z2 = 0; \
    xy = x2 * y2; \
    to = OZ * xy + OY * x2 + OX; \

/* copy data if equal threshold value */
#define FLXMSK_COPY(CVAR, CMAX, IVAR, OVAR) \
    for (CVAR = CMAX; CVAR > 0; --CVAR) *OVAR++ = *IVAR++;
#define FLXMSK_COPY_CMD \
    unsigned long cvar; \
	FLXMSK_COPY(cvar, ne, from, to)

/* copy data from one (non type-specific) pointer to another with */
/* an offset, for use when the target is SMALLER than the source */
#define FLXMSK_COPY_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
        *TAP++ = SAP[yo++]; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_COPY_FROMPART_A1A2CMD \
    FLXMSK_COPY_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2])

/* copy data from one (non type-specific) pointer to another with */
/* an offset, for use when the target is LARGER than the source */
#define FLXMSK_COPY_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
        TAP[yo++] = *SAP++; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_COPY_INTOPART_A1A2CMD \
    FLXMSK_COPY_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2])

/* copy data if equal threshold value */
#define FLXMSK_COPYEQ(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ == TVAR) \
            OVAR[CVAR] = IVAR[-1];
#define FLXMSK_COPYEQ_CMD \
    unsigned long cvar; \
	FLXMSK_COPYEQ(cvar, ne, from, to, t)

/* copy data from one (non type-specific) pointer to another with */
/* an offset, for use when the target is SMALLER than the source, */
/* v == threshold */
#define FLXMSK_COPYEQ_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
        if (SAP[yo] == THR) \
            *TAP = SAP[yo]; \
        ++TAP; \
        ++yo; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_COPYEQ_FROMPART_A1A2CMD \
    FLXMSK_COPYEQ_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* copy data from one (non type-specific) pointer to another with */
/* an offset for use when the target is LARGER than the source, */
/* v == threshold */
#define FLXMSK_COPYEQ_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
        if (*SAP++ == THR) \
            TAP[yo] = SAP[-1]; \
        ++yo; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_COPYEQ_INTOPART_A1A2CMD \
    FLXMSK_COPYEQ_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* copy data if greater or equal than a threshold value */
#define FLXMSK_COPYGE(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ >= TVAR) \
            OVAR[CVAR] = IVAR[-1];
#define FLXMSK_COPYGE_CMD \
    unsigned long cvar; \
	FLXMSK_COPYGE(cvar, ne, from, to, t)

/* copy data if equal threshold value */
#define FLXMSK_COPYLE(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ <= TVAR) \
            OVAR[CVAR] = IVAR[-1];
#define FLXMSK_COPYLE_CMD \
    unsigned long cvar; \
	FLXMSK_COPYEQ(cvar, ne, from, to, t)

/* copy data from one (non type-specific) pointer to another with */
/* an offset for use when the target is LARGER than the source, */
/* v >= threshold */
#define FLXMSK_COPYLE_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
        if (*SAP++ <= THR) \
            TAP[yo] = SAP[-1]; \
        ++yo; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_COPYLE_INTOPART_A1A2CMD \
    FLXMSK_COPYLE_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* copy data from one (non type-specific) pointer to another with */
/* an offset and apply threshold, for use when the target is SMALLER */
/* than the source */
#define FLXMSK_COPYGE_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
        if (SAP[yo] >= THR) \
            *TAP = SAP[yo]; \
        ++TAP; \
        ++yo; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_COPYGE_FROMPART_A1A2CMD \
    FLXMSK_COPYGE_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* copy data from one (non type-specific) pointer to another with */
/* an offset and apply threshold, for use when the target is LARGER */
/* than the source */
#define FLXMSK_COPYGE_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
        if (*SAP >= THR) \
            TAP[yo] = *SAP; \
        ++SAP; \
        ++yo; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_COPYGE_INTOPART_A1A2CMD \
    FLXMSK_COPYGE_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* copy data if less or equal than a threshold value */
#define FLXMSK_COPYLE(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ <= TVAR) \
            OVAR[CVAR] = IVAR[-1];

/* copy data from one (non type-specific) pointer to another with */
/* an offset and apply threshold, for use when the target is SMALLER */
/* than the source */
#define FLXMSK_COPYLE_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
        if (SAP[yo] <= THR) \
            *TAP = SAP[yo]; \
        ++TAP; \
        ++yo; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_COPYLE_FROMPART_A1A2CMD \
    FLXMSK_COPYLE_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* create a mxArray and assign/check mxArray and data array pointers */
#define FLXMSK_CREATEARRAY(ATYPE, ADIM, AVAR, PTYPE, PVAR) \
    AVAR = mxCreateNumericArray(3, ADIM, mx##ATYPE##_CLASS, mxREAL); \
    if ( AVAR == NULL ) \
        mexErrMsgTxt("Error creating array."); \
	PVAR = ( PTYPE *) mxGetData( AVAR ); \
    if ( PVAR == NULL ) \
        mexErrMsgTxt("Error getting pointer to array.");

/* sets boolean/uint8 array to 0 where data is not equal threshold */
#define FLXMSK_SET0NE(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ != TVAR) \
            OVAR[CVAR] = 0;
#define FLXMSK_SET0NE_CMD \
    unsigned long cvar; \
	FLXMSK_SET0NE(cvar, ne, from, to, t)

/* sets boolean/uint8 to 0 from (non type-specific) pointer with */
/* an offset where value is not equal threshold */
/* for use when the target is SMALLER than the source */
#define FLXMSK_SET0NE_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
                if (SAP[yo++] != THR) \
                    *TAP = 0; \
                ++TAP; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_SET0NE_FROMPART_A1A2CMD \
    FLXMSK_SET0NE_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 to 0 from (non type-specific) pointer with */
/* an offset where value is not equal threshold */
/* for use when the target is LARGER than the source */
#define FLXMSK_SET0NE_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
                if (*SAP++ != THR) \
                    TAP[yo] = 0; \
                ++yo; \
	FLXMSK_COPYLOOP_CLOSE
#define FLXMSK_SET0NE_INTOPART_A1A2CMD \
    FLXMSK_SET0NE_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 array to 1 where data is exactly = threshold */
#define FLXMSK_SET1EQ(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ == TVAR) { \
            OVAR[CVAR] = 1; \
            ++z1; \
        } \
    return z1;
#define FLXMSK_SET1EQ_CMD \
    unsigned long cvar, z1 = 0; \
	FLXMSK_SET1EQ(cvar, ne, from, to, t)

/* sets boolean/uint8 to 1 from (non type-specific) pointer with */
/* an offset where value is equal threshold */
/* for use when the target is SMALLER than the source */
#define FLXMSK_SET1EQ_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
                if (SAP[yo++] == THR) { \
                    *TAP = 1; \
                    ++z1; \
                } \
                ++TAP; \
	FLXMSK_COPYLOOP1_CLOSE
#define FLXMSK_SET1EQ_FROMPART_A1A2CMD \
    FLXMSK_SET1EQ_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 to 1 from (non type-specific) pointer with */
/* an offset where value is equal threshold */
/* for use when the target is LARGER than the source */
#define FLXMSK_SET1EQ_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
                if (*SAP++ == THR) { \
                    TAP[yo] = 1; \
                    ++z2; \
                } \
                ++yo; \
	FLXMSK_COPYLOOP2_CLOSE
#define FLXMSK_SET1EQ_INTOPART_A1A2CMD \
    FLXMSK_SET1EQ_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 array to 1 where data is greater or equal threshold */
#define FLXMSK_SET1GE(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ >= TVAR) { \
            OVAR[CVAR] = 1; \
            ++z1; \
        } \
    return z1;
#define FLXMSK_SET1GE_CMD \
    unsigned long cvar, z1 = 0; \
	FLXMSK_SET1GE(cvar, ne, from, to, t)

/* sets boolean/uint8 to 1 from (non type-specific) pointer with */
/* an offset where value is greater or equal threshold */
/* for use when the target is SMALLER than the source */
#define FLXMSK_SET1GE_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
                if (SAP[yo++] >= THR) { \
                    *TAP = 1; \
                    ++z1; \
                } \
                ++TAP; \
	FLXMSK_COPYLOOP1_CLOSE
#define FLXMSK_SET1GE_FROMPART_A1A2CMD \
    FLXMSK_SET1GE_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 to 1 from (non type-specific) pointer with */
/* an offset where value is greater or equal threshold */
/* for use when the target is LARGER than the source */
#define FLXMSK_SET1GE_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
                if (*SAP++ >= THR) { \
                    TAP[yo] = 1; \
                    ++z2; \
                } \
                ++yo; \
	FLXMSK_COPYLOOP2_CLOSE
#define FLXMSK_SET1GE_INTOPART_A1A2CMD \
    FLXMSK_SET1GE_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 array to 1 where data is greater or equal threshold */
#define FLXMSK_SET1LE(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ <= TVAR) { \
            OVAR[CVAR] = 1; \
            ++z1; \
        } \
    return z1;
#define FLXMSK_SET1LE_CMD \
    unsigned long cvar, z1 = 0; \
	FLXMSK_SET1LE(cvar, ne, from, to, t)

/* sets boolean/uint8 to 1 from (non type-specific) pointer with */
/* an offset where value is less or equal threshold */
/* for use when the target is SMALLER than the source */
#define FLXMSK_SET1LE_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
                if (SAP[yo++] <= THR) { \
                    *TAP = 1; \
                    ++z1; \
                } \
                ++TAP; \
	FLXMSK_COPYLOOP1_CLOSE
#define FLXMSK_SET1LE_FROMPART_A1A2CMD \
    FLXMSK_SET1LE_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 to 1 from (non type-specific) pointer with */
/* an offset where value is greater or equal threshold */
/* for use when the target is LARGER than the source */
#define FLXMSK_SET1LE_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
                if (*SAP++ <= THR) { \
                    TAP[yo] = 1; \
                    ++z2; \
                } \
                ++yo; \
	FLXMSK_COPYLOOP2_CLOSE
#define FLXMSK_SET1LE_INTOPART_A1A2CMD \
    FLXMSK_SET1GE_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 array to 1 where data is not equal threshold */
#define FLXMSK_SET1NE(CVAR, CMAX, IVAR, OVAR, TVAR) \
    for (CVAR = 0; CVAR < CMAX; ++CVAR) \
        if (*IVAR++ != TVAR) { \
            OVAR[CVAR] = 1; \
            ++z1; \
        } \
    return z1;
#define FLXMSK_SET1NE_CMD \
    unsigned long cvar, z1 = 0; \
	FLXMSK_SET1NE(cvar, ne, from, to, t)

/* sets boolean/uint8 to 1 from (non type-specific) pointer with */
/* an offset where value is not equal threshold */
/* for use when the target is SMALLER than the source */
#define FLXMSK_SET1NE_FROMPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS1(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP1_OPEN \
                if (SAP[yo++] != THR) { \
                    *TAP = 1; \
                    ++z1; \
                } \
                ++TAP; \
	FLXMSK_COPYLOOP1_CLOSE
#define FLXMSK_SET1NE_FROMPART_A1A2CMD \
    FLXMSK_SET1NE_FROMPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* sets boolean/uint8 to 1 from (non type-specific) pointer with */
/* an offset where value is not equal threshold */
/* for use when the target is LARGER than the source */
#define FLXMSK_SET1NE_INTOPART_A1A2(SAD, SAP, TAD, TAP, OX, OY, OZ, THR) \
    FLXMSK_COPYVARS2(SAD, TAD, OX, OY, OZ) \
	FLXMSK_COPYLOOP2_OPEN \
                if (*SAP++ != THR) { \
                    TAP[yo] = 1; \
                    ++z2; \
                } \
                ++yo; \
	FLXMSK_COPYLOOP2_CLOSE
#define FLXMSK_SET1NE_INTOPART_A1A2CMD \
    FLXMSK_SET1NE_INTOPART_A1A2(snd, sd, tnd, td, off[0], off[1], off[2], t)

/* copy functions using macros */
#define FLXMSK_FUNC_COPY(MLtype, Ctype) \
void FLXMSK_copy_##MLtype(unsigned long ne, const Ctype *from, Ctype *to) { FLXMSK_COPY_CMD }
#define FLXMSK_FUNC_COPY_FROMPART(MLtype, Ctype) \
void FLXMSK_copy_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, Ctype *td, const unsigned long *off) { FLXMSK_COPY_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_COPY_INTOPART(MLtype, Ctype) \
void FLXMSK_copy_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, Ctype *td, const unsigned long *off) { FLXMSK_COPY_INTOPART_A1A2CMD }
#define FLXMSK_FUNC_COPY_TODBL(MLtype, Ctype) \
void FLXMSK_copy_todbl_##MLtype(unsigned long ne, const Ctype *from, double *to) { FLXMSK_COPY_CMD }
#define FLXMSK_FUNC_COPY_TODBL_FROMPART(MLtype, Ctype) \
void FLXMSK_copy_todbl_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, double *td, const unsigned long *off) { FLXMSK_COPY_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_COPY_TODBL_INTOPART(MLtype, Ctype) \
void FLXMSK_copy_todbl_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, double *td, const unsigned long *off) { FLXMSK_COPY_INTOPART_A1A2CMD }

/* copy where equal */
#define FLXMSK_FUNC_COPYEQ(MLtype, Ctype) \
void FLXMSK_copyeq_##MLtype(unsigned long ne, const Ctype *from, Ctype *to, Ctype t) { FLXMSK_COPYEQ_CMD }
#define FLXMSK_FUNC_COPYEQ_FROMPART(MLtype, Ctype) \
void FLXMSK_copyeq_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, Ctype *td, const unsigned long *off, Ctype t) { FLXMSK_COPYEQ_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_COPYEQ_INTOPART(MLtype, Ctype) \
void FLXMSK_copyeq_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, Ctype *td, const unsigned long *off, Ctype t) { FLXMSK_COPYEQ_INTOPART_A1A2CMD }
#define FLXMSK_FUNC_COPYEQ_TODBL(MLtype, Ctype) \
void FLXMSK_copyeq_todbl_##MLtype(unsigned long ne, const Ctype *from, double *to, Ctype t) { FLXMSK_COPYEQ_CMD }
#define FLXMSK_FUNC_COPYEQ_TODBL_FROMPART(MLtype, Ctype) \
void FLXMSK_copyeq_todbl_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, double *td, const unsigned long *off, Ctype t) { FLXMSK_COPYEQ_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_COPYEQ_TODBL_INTOPART(MLtype, Ctype) \
void FLXMSK_copyeq_todbl_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, double *td, const unsigned long *off, Ctype t) { FLXMSK_COPYEQ_INTOPART_A1A2CMD }

/* copy where greate or equal */
#define FLXMSK_FUNC_COPYGE(MLtype, Ctype) \
void FLXMSK_copyge_##MLtype(unsigned long ne, const Ctype *from, Ctype *to, Ctype t) { FLXMSK_COPYGE_CMD }
#define FLXMSK_FUNC_COPYGE_FROMPART(MLtype, Ctype) \
void FLXMSK_copyge_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, Ctype *td, const unsigned long *off, Ctype t) { FLXMSK_COPYGE_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_COPYGE_INTOPART(MLtype, Ctype) \
void FLXMSK_copyge_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, Ctype *td, const unsigned long *off, Ctype t) { FLXMSK_COPYGE_INTOPART_A1A2CMD }
#define FLXMSK_FUNC_COPYGE_TODBL(MLtype, Ctype) \
void FLXMSK_copyge_todbl_##MLtype(unsigned long ne, const Ctype *from, double *to, Ctype t) { FLXMSK_COPYGE_CMD }
#define FLXMSK_FUNC_COPYGE_TODBL_FROMPART(MLtype, Ctype) \
void FLXMSK_copyge_todbl_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, double *td, const unsigned long *off, Ctype t) { FLXMSK_COPYGE_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_COPYGE_TODBL_INTOPART(MLtype, Ctype) \
void FLXMSK_copyge_todbl_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, double *td, const unsigned long *off, Ctype t) { FLXMSK_COPYGE_INTOPART_A1A2CMD }

/* copy where less or equal */
#define FLXMSK_FUNC_COPYLE(MLtype, Ctype) \
void FLXMSK_copyle_##MLtype(unsigned long ne, const Ctype *from, Ctype *to, Ctype t) { FLXMSK_COPYLE_CMD }
#define FLXMSK_FUNC_COPYLE_FROMPART(MLtype, Ctype) \
void FLXMSK_copyle_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, Ctype *td, const unsigned long *off, Ctype t) { FLXMSK_COPYLE_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_COPYLE_INTOPART(MLtype, Ctype) \
void FLXMSK_copyle_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, Ctype *td, const unsigned long *off, Ctype t) { FLXMSK_COPYLE_INTOPART_A1A2CMD }
#define FLXMSK_FUNC_COPYLE_TODBL(MLtype, Ctype) \
void FLXMSK_copyle_todbl_##MLtype(unsigned long ne, const Ctype *from, double *to, Ctype t) { FLXMSK_COPYLE_CMD }
#define FLXMSK_FUNC_COPYLE_TODBL_FROMPART(MLtype, Ctype) \
void FLXMSK_copyle_todbl_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, double *td, const unsigned long *off, Ctype t) { FLXMSK_COPYLE_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_COPYLE_TODBL_INTOPART(MLtype, Ctype) \
void FLXMSK_copyle_todbl_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, double *td, const unsigned long *off, Ctype t) { FLXMSK_COPYLE_INTOPART_A1A2CMD }

/* set 1 (true) where equal */
#define FLXMSK_FUNC_SET1EQ(MLtype, Ctype) \
unsigned long FLXMSK_set1eq_##MLtype(unsigned long ne, const Ctype *from, unsigned char *to, Ctype t) { FLXMSK_SET1EQ_CMD }
#define FLXMSK_FUNC_SET1EQ_FROMPART(MLtype, Ctype) \
unsigned long FLXMSK_set1eq_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, unsigned char *td, const unsigned long *off, Ctype t) { FLXMSK_SET1EQ_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_SET1EQ_INTOPART(MLtype, Ctype) \
unsigned long FLXMSK_set1eq_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, unsigned char *td, const unsigned long *off, Ctype t) { FLXMSK_SET1EQ_INTOPART_A1A2CMD }

/* set 1 (true) where greater or equal */
#define FLXMSK_FUNC_SET1GE(MLtype, Ctype) \
unsigned long FLXMSK_set1ge_##MLtype(unsigned long ne, const Ctype *from, unsigned char *to, Ctype t) { FLXMSK_SET1GE_CMD }
#define FLXMSK_FUNC_SET1GE_FROMPART(MLtype, Ctype) \
unsigned long FLXMSK_set1ge_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, unsigned char *td, const unsigned long *off, Ctype t) { FLXMSK_SET1GE_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_SET1GE_INTOPART(MLtype, Ctype) \
unsigned long FLXMSK_set1ge_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, unsigned char *td, const unsigned long *off, Ctype t) { FLXMSK_SET1GE_INTOPART_A1A2CMD }

/* set 1 (true) where less or equal */
#define FLXMSK_FUNC_SET1LE(MLtype, Ctype) \
unsigned long FLXMSK_set1le_##MLtype(unsigned long ne, const Ctype *from, unsigned char *to, Ctype t) { FLXMSK_SET1LE_CMD }
#define FLXMSK_FUNC_SET1LE_FROMPART(MLtype, Ctype) \
unsigned long FLXMSK_set1le_frompart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, unsigned char *td, const unsigned long *off, Ctype t) { FLXMSK_SET1LE_FROMPART_A1A2CMD }
#define FLXMSK_FUNC_SET1LE_INTOPART(MLtype, Ctype) \
unsigned long FLXMSK_set1le_intopart_##MLtype(const int *snd, const Ctype *sd, const int *tnd, unsigned char *td, const unsigned long *off, Ctype t) { FLXMSK_SET1LE_INTOPART_A1A2CMD }

/* macro to define all type functions */
#define FLXMSK_TYPEFUNCS(FNNAME) \
FLXMSK_FUNC_##FNNAME( int8 , signed char) \
FLXMSK_FUNC_##FNNAME( int16, signed short) \
FLXMSK_FUNC_##FNNAME( int32, signed int) \
FLXMSK_FUNC_##FNNAME(uint8 , unsigned char) \
FLXMSK_FUNC_##FNNAME(uint16, unsigned short) \
FLXMSK_FUNC_##FNNAME(uint32, unsigned int) \
FLXMSK_FUNC_##FNNAME(single, float) \
FLXMSK_FUNC_##FNNAME(double, double) \

/* macros to make the cases for types */
#define FLXMSK_SWCASE_COPY(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: FLXMSK_##FNNAME##_double( ne , (         double *) dp,       (         double *) ovp      ); break; \
    case mxSINGLE_CLASS: FLXMSK_##FNNAME##_single( ne , (         float  *) dp,       (         float  *) ovp      ); break; \
    case mxLOGICAL_CLASS: \
    case mxINT8_CLASS:     FLXMSK_##FNNAME##_int8( ne , (  signed char   *) dp,       (  signed char   *) ovp      ); break; \
    case mxINT16_CLASS:   FLXMSK_##FNNAME##_int16( ne , (  signed short  *) dp,       (  signed short  *) ovp      ); break; \
    case mxINT32_CLASS:   FLXMSK_##FNNAME##_int32( ne , (  signed int    *) dp,       (  signed int    *) ovp      ); break; \
    case mxUINT8_CLASS:   FLXMSK_##FNNAME##_uint8( ne , (unsigned char   *) dp,       (unsigned char   *) ovp      ); break; \
    case mxUINT16_CLASS: FLXMSK_##FNNAME##_uint16( ne , (unsigned short  *) dp,       (unsigned short  *) ovp      ); break; \
    case mxUINT32_CLASS: FLXMSK_##FNNAME##_uint32( ne , (unsigned int    *) dp,       (unsigned int    *) ovp      ); break; \
    default: mexErrMsgTxt("Unsupported class."); break;                                                                      \
}
#define FLXMSK_SWCASE_COPY_PART(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: FLXMSK_##FNNAME##_double(idim, (         double *) dp, odim, (         double *) ovp, offs); break; \
    case mxSINGLE_CLASS: FLXMSK_##FNNAME##_single(idim, (         float  *) dp, odim, (         float  *) ovp, offs); break; \
    case mxLOGICAL_CLASS: \
    case mxINT8_CLASS:     FLXMSK_##FNNAME##_int8(idim, (  signed char   *) dp, odim, (  signed char   *) ovp, offs); break; \
    case mxINT16_CLASS:   FLXMSK_##FNNAME##_int16(idim, (  signed short  *) dp, odim, (  signed short  *) ovp, offs); break; \
    case mxINT32_CLASS:   FLXMSK_##FNNAME##_int32(idim, (  signed int    *) dp, odim, (  signed int    *) ovp, offs); break; \
    case mxUINT8_CLASS:   FLXMSK_##FNNAME##_uint8(idim, (unsigned char   *) dp, odim, (unsigned char   *) ovp, offs); break; \
    case mxUINT16_CLASS: FLXMSK_##FNNAME##_uint16(idim, (unsigned short  *) dp, odim, (unsigned short  *) ovp, offs); break; \
    case mxUINT32_CLASS: FLXMSK_##FNNAME##_uint32(idim, (unsigned int    *) dp, odim, (unsigned int    *) ovp, offs); break; \
    default: mexErrMsgTxt("Unsupported class."); break;                                                                      \
}
#define FLXMSK_SWCASE_COPY_TODBL(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: FLXMSK_##FNNAME##_double( ne , (         double *) dp,       odp      ); break; \
    case mxSINGLE_CLASS: FLXMSK_##FNNAME##_single( ne , (         float  *) dp,       odp      ); break; \
    case mxLOGICAL_CLASS: \
    case mxINT8_CLASS:     FLXMSK_##FNNAME##_int8( ne , (  signed char   *) dp,       odp      ); break; \
    case mxINT16_CLASS:   FLXMSK_##FNNAME##_int16( ne , (  signed short  *) dp,       odp      ); break; \
    case mxINT32_CLASS:   FLXMSK_##FNNAME##_int32( ne , (  signed int    *) dp,       odp      ); break; \
    case mxUINT8_CLASS:   FLXMSK_##FNNAME##_uint8( ne , (unsigned char   *) dp,       odp      ); break; \
    case mxUINT16_CLASS: FLXMSK_##FNNAME##_uint16( ne , (unsigned short  *) dp,       odp      ); break; \
    case mxUINT32_CLASS: FLXMSK_##FNNAME##_uint32( ne , (unsigned int    *) dp,       odp      ); break; \
    default: mexErrMsgTxt("Unsupported class."); break;                                                                      \
}
#define FLXMSK_SWCASE_COPY_TODBL_PART(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: FLXMSK_##FNNAME##_double(idim, (         double *) dp, odim, odp, offs); break; \
    case mxSINGLE_CLASS: FLXMSK_##FNNAME##_single(idim, (         float  *) dp, odim, odp, offs); break; \
    case mxLOGICAL_CLASS: \
    case mxINT8_CLASS:     FLXMSK_##FNNAME##_int8(idim, (  signed char   *) dp, odim, odp, offs); break; \
    case mxINT16_CLASS:   FLXMSK_##FNNAME##_int16(idim, (  signed short  *) dp, odim, odp, offs); break; \
    case mxINT32_CLASS:   FLXMSK_##FNNAME##_int32(idim, (  signed int    *) dp, odim, odp, offs); break; \
    case mxUINT8_CLASS:   FLXMSK_##FNNAME##_uint8(idim, (unsigned char   *) dp, odim, odp, offs); break; \
    case mxUINT16_CLASS: FLXMSK_##FNNAME##_uint16(idim, (unsigned short  *) dp, odim, odp, offs); break; \
    case mxUINT32_CLASS: FLXMSK_##FNNAME##_uint32(idim, (unsigned int    *) dp, odim, odp, offs); break; \
    default: mexErrMsgTxt("Unsupported class."); break;                                                                      \
}
#define FLXMSK_SWCASE_COPYOP(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: FLXMSK_##FNNAME##_double( ne , (         double *) dp,       (         double *) ovp,                         td); break; \
    case mxSINGLE_CLASS: FLXMSK_##FNNAME##_single( ne , (         float  *) dp,       (         float  *) ovp,       (         float ) td); break; \
    case mxLOGICAL_CLASS: \
    case mxINT8_CLASS:     FLXMSK_##FNNAME##_int8( ne , (  signed char   *) dp,       (  signed char   *) ovp,       (  signed char  ) td); break; \
    case mxINT16_CLASS:   FLXMSK_##FNNAME##_int16( ne , (  signed short  *) dp,       (  signed short  *) ovp,       (  signed short ) td); break; \
    case mxINT32_CLASS:   FLXMSK_##FNNAME##_int32( ne , (  signed int    *) dp,       (  signed int    *) ovp,       (  signed int   ) td); break; \
    case mxUINT8_CLASS:   FLXMSK_##FNNAME##_uint8( ne , (unsigned char   *) dp,       (unsigned char   *) ovp,       (unsigned char  ) td); break; \
    case mxUINT16_CLASS: FLXMSK_##FNNAME##_uint16( ne , (unsigned short  *) dp,       (unsigned short  *) ovp,       (unsigned short ) td); break; \
    case mxUINT32_CLASS: FLXMSK_##FNNAME##_uint32( ne , (unsigned int    *) dp,       (unsigned int    *) ovp,       (unsigned int   ) td); break; \
    default: mexErrMsgTxt("Unsupported class."); break;                                                                      \
}
#define FLXMSK_SWCASE_COPYOP_PART(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: FLXMSK_##FNNAME##_double(idim, (         double *) dp, odim, (         double *) ovp, offs,                   td); break; \
    case mxSINGLE_CLASS: FLXMSK_##FNNAME##_single(idim, (         float  *) dp, odim, (         float  *) ovp, offs, (         float ) td); break; \
    case mxINT8_CLASS:     FLXMSK_##FNNAME##_int8(idim, (  signed char   *) dp, odim, (  signed char   *) ovp, offs, (  signed char  ) td); break; \
    case mxINT16_CLASS:   FLXMSK_##FNNAME##_int16(idim, (  signed short  *) dp, odim, (  signed short  *) ovp, offs, (  signed short ) td); break; \
    case mxINT32_CLASS:   FLXMSK_##FNNAME##_int32(idim, (  signed int    *) dp, odim, (  signed int    *) ovp, offs, (  signed int   ) td); break; \
    case mxUINT8_CLASS:   FLXMSK_##FNNAME##_uint8(idim, (unsigned char   *) dp, odim, (unsigned char   *) ovp, offs, (unsigned char  ) td); break; \
    case mxUINT16_CLASS: FLXMSK_##FNNAME##_uint16(idim, (unsigned short  *) dp, odim, (unsigned short  *) ovp, offs, (unsigned short ) td); break; \
    case mxUINT32_CLASS: FLXMSK_##FNNAME##_uint32(idim, (unsigned int    *) dp, odim, (unsigned int    *) ovp, offs, (unsigned int   ) td); break; \
    default: mxDestroyArray(oa); mexErrMsgTxt(flexmask_err_nologmask); break; \
}
#define FLXMSK_SWCASE_COPYOP_TODBL(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: FLXMSK_##FNNAME##_double( ne , (         double *) dp,      odp,                         td); break; \
    case mxSINGLE_CLASS: FLXMSK_##FNNAME##_single( ne , (         float  *) dp,      odp,       (         float ) td); break; \
    case mxLOGICAL_CLASS: \
    case mxINT8_CLASS:     FLXMSK_##FNNAME##_int8( ne , (  signed char   *) dp,      odp,       (  signed char  ) td); break; \
    case mxINT16_CLASS:   FLXMSK_##FNNAME##_int16( ne , (  signed short  *) dp,      odp,       (  signed short ) td); break; \
    case mxINT32_CLASS:   FLXMSK_##FNNAME##_int32( ne , (  signed int    *) dp,      odp,       (  signed int   ) td); break; \
    case mxUINT8_CLASS:   FLXMSK_##FNNAME##_uint8( ne , (unsigned char   *) dp,      odp,       (unsigned char  ) td); break; \
    case mxUINT16_CLASS: FLXMSK_##FNNAME##_uint16( ne , (unsigned short  *) dp,      odp,       (unsigned short ) td); break; \
    case mxUINT32_CLASS: FLXMSK_##FNNAME##_uint32( ne , (unsigned int    *) dp,      odp,       (unsigned int   ) td); break; \
    default: mexErrMsgTxt("Unsupported class."); break;                                                                      \
}
#define FLXMSK_SWCASE_COPYOP_TODBL_PART(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: FLXMSK_##FNNAME##_double(idim, (         double *) dp, odim, odp, offs,                   td); break; \
    case mxSINGLE_CLASS: FLXMSK_##FNNAME##_single(idim, (         float  *) dp, odim, odp, offs, (         float ) td); break; \
    case mxINT8_CLASS:     FLXMSK_##FNNAME##_int8(idim, (  signed char   *) dp, odim, odp, offs, (  signed char  ) td); break; \
    case mxINT16_CLASS:   FLXMSK_##FNNAME##_int16(idim, (  signed short  *) dp, odim, odp, offs, (  signed short ) td); break; \
    case mxINT32_CLASS:   FLXMSK_##FNNAME##_int32(idim, (  signed int    *) dp, odim, odp, offs, (  signed int   ) td); break; \
    case mxUINT8_CLASS:   FLXMSK_##FNNAME##_uint8(idim, (unsigned char   *) dp, odim, odp, offs, (unsigned char  ) td); break; \
    case mxUINT16_CLASS: FLXMSK_##FNNAME##_uint16(idim, (unsigned short  *) dp, odim, odp, offs, (unsigned short ) td); break; \
    case mxUINT32_CLASS: FLXMSK_##FNNAME##_uint32(idim, (unsigned int    *) dp, odim, odp, offs, (unsigned int   ) td); break; \
    default: mxDestroyArray(oa); mexErrMsgTxt(flexmask_err_nologmask); break; \
}
#define FLXMSK_SWCASE_SET1OP(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: ns = FLXMSK_##FNNAME##_double( ne , (         double *) dp,       o8p      ,                   td); break; \
    case mxSINGLE_CLASS: ns = FLXMSK_##FNNAME##_single( ne , (         float  *) dp,       o8p      , (         float ) td); break; \
    case mxLOGICAL_CLASS: \
    case mxINT8_CLASS:   ns =   FLXMSK_##FNNAME##_int8( ne , (  signed char   *) dp,       o8p      , (  signed char  ) td); break; \
    case mxINT16_CLASS:  ns =  FLXMSK_##FNNAME##_int16( ne , (  signed short  *) dp,       o8p      , (  signed short ) td); break; \
    case mxINT32_CLASS:  ns =  FLXMSK_##FNNAME##_int32( ne , (  signed int    *) dp,       o8p      , (  signed int   ) td); break; \
    case mxUINT8_CLASS:  ns =  FLXMSK_##FNNAME##_uint8( ne , (unsigned char   *) dp,       o8p      , (unsigned char  ) td); break; \
    case mxUINT16_CLASS: ns = FLXMSK_##FNNAME##_uint16( ne , (unsigned short  *) dp,       o8p      , (unsigned short ) td); break; \
    case mxUINT32_CLASS: ns = FLXMSK_##FNNAME##_uint32( ne , (unsigned int    *) dp,       o8p      , (unsigned int   ) td); break; \
    default: mexErrMsgTxt("Unsupported class."); break;                                                                      \
}
#define FLXMSK_SWCASE_SET1OP_PART(FNNAME) \
switch (icls) { \
    case mxDOUBLE_CLASS: ns = FLXMSK_##FNNAME##_double(idim, (         double *) dp, odim, o8p, offs,               td); break; \
    case mxSINGLE_CLASS: ns = FLXMSK_##FNNAME##_single(idim, (         float  *) dp, odim, o8p, offs, (         float ) td); break; \
    case mxLOGICAL_CLASS: \
    case mxINT8_CLASS:   ns =   FLXMSK_##FNNAME##_int8(idim, (  signed char   *) dp, odim, o8p, offs, (  signed char  ) td); break; \
    case mxINT16_CLASS:  ns =  FLXMSK_##FNNAME##_int16(idim, (  signed short  *) dp, odim, o8p, offs, (  signed short ) td); break; \
    case mxINT32_CLASS:  ns =  FLXMSK_##FNNAME##_int32(idim, (  signed int    *) dp, odim, o8p, offs, (  signed int   ) td); break; \
    case mxUINT8_CLASS:  ns =  FLXMSK_##FNNAME##_uint8(idim, (unsigned char   *) dp, odim, o8p, offs, (unsigned char  ) td); break; \
    case mxUINT16_CLASS: ns = FLXMSK_##FNNAME##_uint16(idim, (unsigned short  *) dp, odim, o8p, offs, (unsigned short ) td); break; \
    case mxUINT32_CLASS: ns = FLXMSK_##FNNAME##_uint32(idim, (unsigned int    *) dp, odim, o8p, offs, (unsigned int   ) td); break; \
    default: mexErrMsgTxt("Unsupported class."); break;                                                                      \
}

/* finally create copy functions */
FLXMSK_TYPEFUNCS(COPY)
FLXMSK_TYPEFUNCS(COPY_FROMPART)
FLXMSK_TYPEFUNCS(COPY_INTOPART)
FLXMSK_TYPEFUNCS(COPY_TODBL)
FLXMSK_TYPEFUNCS(COPY_TODBL_FROMPART)
FLXMSK_TYPEFUNCS(COPY_TODBL_INTOPART)
FLXMSK_TYPEFUNCS(COPYEQ)
FLXMSK_TYPEFUNCS(COPYEQ_FROMPART)
FLXMSK_TYPEFUNCS(COPYEQ_INTOPART)
FLXMSK_TYPEFUNCS(COPYEQ_TODBL)
FLXMSK_TYPEFUNCS(COPYEQ_TODBL_FROMPART)
FLXMSK_TYPEFUNCS(COPYEQ_TODBL_INTOPART)
FLXMSK_TYPEFUNCS(COPYGE)
FLXMSK_TYPEFUNCS(COPYGE_FROMPART)
FLXMSK_TYPEFUNCS(COPYGE_INTOPART)
FLXMSK_TYPEFUNCS(COPYGE_TODBL)
FLXMSK_TYPEFUNCS(COPYGE_TODBL_FROMPART)
FLXMSK_TYPEFUNCS(COPYGE_TODBL_INTOPART)
FLXMSK_TYPEFUNCS(COPYLE)
FLXMSK_TYPEFUNCS(COPYLE_FROMPART)
FLXMSK_TYPEFUNCS(COPYLE_INTOPART)
FLXMSK_TYPEFUNCS(COPYLE_TODBL)
FLXMSK_TYPEFUNCS(COPYLE_TODBL_FROMPART)
FLXMSK_TYPEFUNCS(COPYLE_TODBL_INTOPART)

/* and also create set1 functions */
FLXMSK_TYPEFUNCS(SET1EQ)
FLXMSK_TYPEFUNCS(SET1EQ_FROMPART)
FLXMSK_TYPEFUNCS(SET1EQ_INTOPART)
FLXMSK_TYPEFUNCS(SET1GE)
FLXMSK_TYPEFUNCS(SET1GE_FROMPART)
FLXMSK_TYPEFUNCS(SET1GE_INTOPART)
FLXMSK_TYPEFUNCS(SET1LE)
FLXMSK_TYPEFUNCS(SET1LE_FROMPART)
FLXMSK_TYPEFUNCS(SET1LE_INTOPART)

/* error text(s) */
static const char *flexmask_err_badthresh = "Bad threshold value.";
static const char *flexmask_err_memalloc  = "Error allocating memory.";
static const char *flexmask_err_memhandle = "Error getting handle to allocated memory.";
static const char *flexmask_err_nologcopy = "Copying to logical arrays not supported (and useful).";
static const char *flexmask_err_nologmask = "Masking of logical arrays not supported (and useful).";

/* scalar value dims */
static const int flexmask_scalar_dim[2] = {1, 1};

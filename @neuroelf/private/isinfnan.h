/*
  isnaninf.h

% Version:  v0.9b
% Build:    11050515
% Date:     Aug-29 2010, 1:47 AM EST
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


/* this is a "cheaper" test for isinf/isnan */


/* IsInf / IsNaN combined */
#define IS_BAD_VAL()                                                    \
    (*isbadptr == posinfval ||                                          \
     *isbadptr == neginfval ||                                          \
     *isbadptr == nanval1   ||                                          \
     *isbadptr == nanval2)

#define IS_GOOD_VAL()                                                   \
    (*isbadptr != posinfval &&                                          \
     *isbadptr != neginfval &&                                          \
     *isbadptr != nanval1   &&                                          \
     *isbadptr != nanval2)

/* IsInf */
#define IS_INF()                                                        \
    (*isbadptr == posinfval  ||                                         \
     *isbadptr == neginfval)

/* IsInf (NEG only) */
#define IS_NEG_INF()                                                    \
    (*isbadptr == neginfval)

/* IsInf (POS only) */
#define IS_POS_INF()                                                    \
    (*isbadptr == posinfval)

/* IsNaN */
#define IS_NAN()                                                        \
    (*isbadptr == nanval1  ||                                           \
     *isbadptr == nanval2)

/* replacement for mxIsInf / mxIsNaN */
#define IF_IS_BAD_VAL(VAL)                                              \
    isbadval.d = (double) VAL;                                          \
    if (IS_BAD_VAL())

/* replacement for !mxIsInf && !mxIsNaN */
#define IF_IS_GOOD_VAL(VAL)                                             \
    isbadval.d = (double) VAL;                                          \
    if (IS_GOOD_VAL())

/* replacement for mxIsInf */
#define IF_IS_INF(VAL)                                                  \
    isbadval.d = (double) VAL;                                          \
    if (IS_INF())

/* replacement for mxIsInf */
#define IF_IS_NEG_INF(VAL)                                              \
    isbadval.d = (double) VAL;                                          \
    if (IS_NEG_INF())

/* replacement for mxIsInf */
#define IF_IS_POS_INF(VAL)                                              \
    isbadval.d = (double) VAL;                                          \
    if (IS_POS_INF())

/* replacement for mxIsNaN */
#define IF_IS_NAN(VAL)                                                  \
    isbadval.d = (double) VAL;                                          \
    if (IS_NAN())

/* initialization of memory for isinf/isnan test */
#define INIT_INF_NAN_BAD_VAL()                                          \
    isbadval.d = mxGetInf();                                            \
    if (isbadval.ui[0] == 0)                                            \
        isbadptr = (const unsigned int*) &isbadval.ui[1];               \
    else                                                                \
        isbadptr = (const unsigned int*) &isbadval.ui[0];               \
    posinfval = *isbadptr;                                              \
    neginfval = posinfval | 0x80000000;                                 \
    isbadval.d = mxGetNaN();                                            \
    nanval1 = *isbadptr & 0x7FFFFFFF;                                   \
    nanval2 = nanval1 | 0x80000000;

/* required variables */
#define VARS_FOR_ISINFNAN                                               \
    union u_iin_dbl_int2 {                                              \
        double d;                                                       \
        unsigned int ui[2];                                             \
        float f[2];                                                     \
    };                                                                  \
    union u_iin_dbl_int2 isbadval;                                      \
    const unsigned int *isbadptr = NULL;                                \
    unsigned int                                                        \
        posinfval = 0,                                                  \
        neginfval = 0,                                                  \
        nanval1 = 0,                                                    \
        nanval2 = 0;

/*

additional methods to index an array

FORMAT:       x = indexarray(a, i1, ...)

Input fields:

      a           N-D array (up to 4D)
      i1, ...     indexing expressions

Output fields:

      x           returned value

Note: supported indexing expressions (either double or uint32) are:
      > single argument
        - NxD number-of-values -by- number-of-input-dimensions
      > multiple arguments (as many as dimensions)
        - Nx1 (or 1xN) specific indexing for this dimension (ordered list)
        - N-D argument where as input has lower dimensionality
          e.g. indexing into a 10x10 array a with a(s, :) where s is a
          10x1x10 array with numbers
          note: for now, only first dimension can be expanded this way!!

% Version:  v0.9b
% Build:    11050511
% Date:     Aug-26 2010, 2:00 PM EST
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

/* macros */

/* create output */
#define IA_CREATE_OUTPUT(CLASS, DT, DP)                                     \
    *plhs = mxCreateNumericArray(od, odim, CLASS, mxREAL);                  \
    if (*plhs == NULL)                                                      \
        mexErrMsgTxt("Error creating output array.");                       \
    DP = (DT *) mxGetData(*plhs);                                           \
    if (DP == NULL)                                                         \
        mexErrMsgTxt("Error getting output data pointer.");

/* get input pointer */
#define IA_GET_INPUT(DT, DP)                                                \
    DP = (const DT *) mxGetData(*prhs);                                     \
    if (DP == NULL)                                                         \
        mexErrMsgTxt("Error getting input data pointer.");

/* expand 1D input */
#define IA_EXP1D(CLASS, DT, DP, ODP, IP)                                    \
    IA_GET_INPUT(DT, DP);                                                   \
    IA_CREATE_OUTPUT(CLASS, DT, ODP);                                       \
    for ( ; c1 > 0; --c1) {                                                 \
        i1 = -1 + ((int) *IP++);                                            \
        if ((i1 < 0) ||                                                     \
            (i1 >= ne))                                                     \
            mexErrMsgTxt("Indexing out of bounds.");                        \
        *ODP++ = DP[i1];                                                    \
    }

/* expand 2D input with 1st to 3rd (a.s.o.) */
#define IA_EXP2D13(CLASS, DT, DP, DPC, ODP, IPA, IPB, IPC, IPCC)            \
    IA_GET_INPUT(DT, DP);                                                   \
    IA_CREATE_OUTPUT(CLASS, DT, ODP);                                       \
    for (c3 = 0; c3 < s3; ++c3) {                                           \
        IPC = &IPA[c3 * s1];                                                \
        for (c2 = 0; c2 < s2; ++c2) {                                       \
            i2 = -1 + ((int) IPB[c2]);                                      \
            if ((i2 < 0) ||                                                 \
                (i2 >= n2))                                                 \
                mexErrMsgTxt("Indexing out of bounds.");                    \
            DPC = &DP[i2 * d1];                                             \
            IPCC = IPC;                                                     \
            for (c1 = s1; c1 > 0; --c1) {                                   \
                i1 = -1 + ((int) *IPCC++);                                  \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DPC[i1];                                           \
            }                                                               \
        }                                                                   \
    }

/* index from a single ND range argument */
#define IA_RANGE_ND(CLASS, DT, DP, ODP)                                     \
    IA_GET_INPUT(DT, DP);                                                   \
    IA_CREATE_OUTPUT(CLASS, DT, ODP);                                       \
    switch (nd) {                                                           \
        case 1:                                                             \
            if (st1 > 0.0)                                                  \
                for (fr1 = *di1; fr1 <= to1; fr1 += st1) {                  \
                    i1 = (int) (fr1 - 0.5);                                 \
                    if ((i1 < 0) ||                                         \
                        (i1 >= d1))                                         \
                        ++ODP;                                              \
                    else                                                    \
                        *ODP++ = DP[i1];                                    \
                }                                                           \
            else                                                            \
                for (fr1 = *di1; fr1 >= to1; fr1 += st1) {                  \
                    i1 = (int) (fr1 - 0.5);                                 \
                    if ((i1 < 0) ||                                         \
                        (i1 >= d1))                                         \
                        ++ODP;                                              \
                    else                                                    \
                        *ODP++ = DP[i1];                                    \
                }                                                           \
            break;                                                          \
        case 2:                                                             \
            if (st2 > 0.0) {                                                \
                if (st1 > 0.0) {                                            \
                    for (fr2 = *di2; fr2 <= to2; fr2 += st2) {              \
                        i2 = (int) (fr2 - 0.5);                             \
                        if ((i2 < 0) ||                                     \
                            (i2 >= d2))                                     \
                            ODP = &ODP[s2];                                 \
                        else {                                              \
                            i2 *= n2;                                       \
                            for (fr1 = *di1; fr1 <= to1; fr1 += st1) {      \
                                i1 = (int) (fr1 - 0.5);                     \
                                if ((i1 < 0) ||                             \
                                    (i1 >= d1))                             \
                                    ++ODP;                                  \
                                else                                        \
                                    *ODP++ = DP[i1 + i2];                   \
                            }                                               \
                        }                                                   \
                    }                                                       \
                } else {                                                    \
                    for (fr2 = *di2; fr2 <= to2; fr2 += st2) {              \
                        i2 = (int) (fr2 - 0.5);                             \
                        if ((i2 < 0) ||                                     \
                            (i2 >= d2))                                     \
                            ODP = &ODP[s2];                                 \
                        else {                                              \
                            i2 *= n2;                                       \
                            for (fr1 = *di1; fr1 >= to1; fr1 += st1) {      \
                                i1 = (int) (fr1 - 0.5);                     \
                                if ((i1 < 0) ||                             \
                                    (i1 >= d1))                             \
                                    ++ODP;                                  \
                                else                                        \
                                    *ODP++ = DP[i1 + i2];                   \
                            }                                               \
                        }                                                   \
                    }                                                       \
                }                                                           \
            } else {                                                        \
                if (st1 > 0.0) {                                            \
                    for (fr2 = *di2; fr2 >= to2; fr2 += st2) {              \
                        i2 = (int) (fr2 - 0.5);                             \
                        if ((i2 < 0) ||                                     \
                            (i2 >= d2))                                     \
                            ODP = &ODP[s2];                                 \
                        else {                                              \
                            i2 *= n2;                                       \
                            for (fr1 = *di1; fr1 <= to1; fr1 += st1) {      \
                                i1 = (int) (fr1 - 0.5);                     \
                                if ((i1 < 0) ||                             \
                                    (i1 >= d1))                             \
                                    ++ODP;                                  \
                                else                                        \
                                    *ODP++ = DP[i1 + i2];                   \
                            }                                               \
                        }                                                   \
                    }                                                       \
                } else {                                                    \
                    for (fr2 = *di2; fr2 >= to2; fr2 += st2) {              \
                        i2 = (int) (fr2 - 0.5);                             \
                        if ((i2 < 0) ||                                     \
                            (i2 >= d2))                                     \
                            ODP = &ODP[s2];                                 \
                        else {                                              \
                            i2 *= n2;                                       \
                            for (fr1 = *di1; fr1 >= to1; fr1 += st1) {      \
                                i1 = (int) (fr1 - 0.5);                     \
                                if ((i1 < 0) ||                             \
                                    (i1 >= d1))                             \
                                    ++ODP;                                  \
                                else                                        \
                                    *ODP++ = DP[i1 + i2];                   \
                            }                                               \
                        }                                                   \
                    }                                                       \
                }                                                           \
            }                                                               \
            break;                                                          \
        case 3:                                                             \
            if (st3 > 0.0) {                                                \
                if (st2 > 0.0) {                                            \
                    if (st1 > 0.0) {                                        \
                        for (fr3 = *di3; fr3 <= to3; fr3 += st3) {          \
                            i3 = (int) (fr3 - 0.5);                         \
                            if ((i3 < 0) ||                                 \
                                (i3 >= d3))                                 \
                                ODP = &ODP[s3];                             \
                            else {                                          \
                                i3 *= n3;                                   \
                                for (fr2 = *di2; fr2 <= to2; fr2 += st2) {  \
                                    i2 = (int) (fr2 - 0.5);                 \
                                    if ((i2 < 0) ||                         \
                                        (i2 >= d2))                         \
                                        ODP = &ODP[s2];                     \
                                    else {                                  \
                                        i2 *= n2;                           \
                                        i2 += i3;                           \
                                        for (fr1 = *di1; fr1 <= to1; fr1 += st1) { \
                                            i1 = (int) (fr1 - 0.5);         \
                                            if ((i1 < 0) ||                 \
                                                (i1 >= d1))                 \
                                                ++ODP;                      \
                                            else                            \
                                                *ODP++ = DP[i1 + i2];       \
                                        }                                   \
                                    }                                       \
                                }                                           \
                            }                                               \
                        }                                                   \
                    } else {                                                \
                        for (fr3 = *di3; fr3 <= to3; fr3 += st3) {          \
                            i3 = (int) (fr3 - 0.5);                         \
                            if ((i3 < 0) ||                                 \
                                (i3 >= d3))                                 \
                                ODP = &ODP[s3];                             \
                            else {                                          \
                                i3 *= n3;                                   \
                                for (fr2 = *di2; fr2 <= to2; fr2 += st2) {  \
                                    i2 = (int) (fr2 - 0.5);                 \
                                    if ((i2 < 0) ||                         \
                                        (i2 >= d2))                         \
                                        ODP = &ODP[s2];                     \
                                    else {                                  \
                                        i2 *= n2;                           \
                                        i2 += i3;                           \
                                        for (fr1 = *di1; fr1 >= to1; fr1 += st1) { \
                                            i1 = (int) (fr1 - 0.5);         \
                                            if ((i1 < 0) ||                 \
                                                (i1 >= d1))                 \
                                                ++ODP;                      \
                                            else                            \
                                                *ODP++ = DP[i1 + i2];       \
                                        }                                   \
                                    }                                       \
                                }                                           \
                            }                                               \
                        }                                                   \
                    }                                                       \
                } else {                                                    \
                    if (st1 > 0.0) {                                        \
                        for (fr3 = *di3; fr3 <= to3; fr3 += st3) {          \
                            i3 = (int) (fr3 - 0.5);                         \
                            if ((i3 < 0) ||                                 \
                                (i3 >= d3))                                 \
                                ODP = &ODP[s3];                             \
                            else {                                          \
                                i3 *= n3;                                   \
                                for (fr2 = *di2; fr2 >= to2; fr2 += st2) {  \
                                    i2 = (int) (fr2 - 0.5);                 \
                                    if ((i2 < 0) ||                         \
                                        (i2 >= d2))                         \
                                        ODP = &ODP[s2];                     \
                                    else {                                  \
                                        i2 *= n2;                           \
                                        i2 += i3;                           \
                                        for (fr1 = *di1; fr1 <= to1; fr1 += st1) { \
                                            i1 = (int) (fr1 - 0.5);         \
                                            if ((i1 < 0) ||                 \
                                                (i1 >= d1))                 \
                                                ++ODP;                      \
                                            else                            \
                                                *ODP++ = DP[i1 + i2];       \
                                        }                                   \
                                    }                                       \
                                }                                           \
                            }                                               \
                        }                                                   \
                    } else {                                                \
                        for (fr3 = *di3; fr3 <= to3; fr3 += st3) {          \
                            i3 = (int) (fr3 - 0.5);                         \
                            if ((i3 < 0) ||                                 \
                                (i3 >= d3))                                 \
                                ODP = &ODP[s3];                             \
                            else {                                          \
                                i3 *= n3;                                   \
                                for (fr2 = *di2; fr2 >= to2; fr2 += st2) {  \
                                    i2 = (int) (fr2 - 0.5);                 \
                                    if ((i2 < 0) ||                         \
                                        (i2 >= d2))                         \
                                        ODP = &ODP[s2];                     \
                                    else {                                  \
                                        i2 *= n2;                           \
                                        i2 += i3;                           \
                                        for (fr1 = *di1; fr1 >= to1; fr1 += st1) { \
                                            i1 = (int) (fr1 - 0.5);         \
                                            if ((i1 < 0) ||                 \
                                                (i1 >= d1))                 \
                                                ++ODP;                      \
                                            else                            \
                                                *ODP++ = DP[i1 + i2];       \
                                        }                                   \
                                    }                                       \
                                }                                           \
                            }                                               \
                        }                                                   \
                    }                                                       \
                }                                                           \
            } else {                                                        \
                if (st2 > 0.0) {                                            \
                    if (st1 > 0.0) {                                        \
                        for (fr3 = *di3; fr3 >= to3; fr3 += st3) {          \
                            i3 = (int) (fr3 - 0.5);                         \
                            if ((i3 < 0) ||                                 \
                                (i3 >= d3))                                 \
                                ODP = &ODP[s3];                             \
                            else {                                          \
                                i3 *= n3;                                   \
                                for (fr2 = *di2; fr2 <= to2; fr2 += st2) {  \
                                    i2 = (int) (fr2 - 0.5);                 \
                                    if ((i2 < 0) ||                         \
                                        (i2 >= d2))                         \
                                        ODP = &ODP[s2];                     \
                                    else {                                  \
                                        i2 *= n2;                           \
                                        i2 += i3;                           \
                                        for (fr1 = *di1; fr1 <= to1; fr1 += st1) { \
                                            i1 = (int) (fr1 - 0.5);         \
                                            if ((i1 < 0) ||                 \
                                                (i1 >= d1))                 \
                                                ++ODP;                      \
                                            else                            \
                                                *ODP++ = DP[i1 + i2];       \
                                        }                                   \
                                    }                                       \
                                }                                           \
                            }                                               \
                        }                                                   \
                    } else {                                                \
                        for (fr3 = *di3; fr3 >= to3; fr3 += st3) {          \
                            i3 = (int) (fr3 - 0.5);                         \
                            if ((i3 < 0) ||                                 \
                                (i3 >= d3))                                 \
                                ODP = &ODP[s3];                             \
                            else {                                          \
                                i3 *= n3;                                   \
                                for (fr2 = *di2; fr2 <= to2; fr2 += st2) {  \
                                    i2 = (int) (fr2 - 0.5);                 \
                                    if ((i2 < 0) ||                         \
                                        (i2 >= d2))                         \
                                        ODP = &ODP[s2];                     \
                                    else {                                  \
                                        i2 *= n2;                           \
                                        i2 += i3;                           \
                                        for (fr1 = *di1; fr1 >= to1; fr1 += st1) { \
                                            i1 = (int) (fr1 - 0.5);         \
                                            if ((i1 < 0) ||                 \
                                                (i1 >= d1))                 \
                                                ++ODP;                      \
                                            else                            \
                                                *ODP++ = DP[i1 + i2];       \
                                        }                                   \
                                    }                                       \
                                }                                           \
                            }                                               \
                        }                                                   \
                    }                                                       \
                } else {                                                    \
                    if (st1 > 0.0) {                                        \
                        for (fr3 = *di3; fr3 >= to3; fr3 += st3) {          \
                            i3 = (int) (fr3 - 0.5);                         \
                            if ((i3 < 0) ||                                 \
                                (i3 >= d3))                                 \
                                ODP = &ODP[s3];                             \
                            else {                                          \
                                i3 *= n3;                                   \
                                for (fr2 = *di2; fr2 >= to2; fr2 += st2) {  \
                                    i2 = (int) (fr2 - 0.5);                 \
                                    if ((i2 < 0) ||                         \
                                        (i2 >= d2))                         \
                                        ODP = &ODP[s2];                     \
                                    else {                                  \
                                        i2 *= n2;                           \
                                        i2 += i3;                           \
                                        for (fr1 = *di1; fr1 <= to1; fr1 += st1) { \
                                            i1 = (int) (fr1 - 0.5);         \
                                            if ((i1 < 0) ||                 \
                                                (i1 >= d1))                 \
                                                ++ODP;                      \
                                            else                            \
                                                *ODP++ = DP[i1 + i2];       \
                                        }                                   \
                                    }                                       \
                                }                                           \
                            }                                               \
                        }                                                   \
                    } else {                                                \
                        for (fr3 = *di3; fr3 >= to3; fr3 += st3) {          \
                            i3 = (int) (fr3 - 0.5);                         \
                            if ((i3 < 0) ||                                 \
                                (i3 >= d3))                                 \
                                ODP = &ODP[s3];                             \
                            else {                                          \
                                i3 *= n3;                                   \
                                for (fr2 = *di2; fr2 >= to2; fr2 += st2) {  \
                                    i2 = (int) (fr2 - 0.5);                 \
                                    if ((i2 < 0) ||                         \
                                        (i2 >= d2))                         \
                                        ODP = &ODP[s2];                     \
                                    else {                                  \
                                        i2 *= n2;                           \
                                        i2 += i3;                           \
                                        for (fr1 = *di1; fr1 >= to1; fr1 += st1) { \
                                            i1 = (int) (fr1 - 0.5);         \
                                            if ((i1 < 0) ||                 \
                                                (i1 >= d1))                 \
                                                ++ODP;                      \
                                            else                            \
                                                *ODP++ = DP[i1 + i2];       \
                                        }                                   \
                                    }                                       \
                                }                                           \
                            }                                               \
                        }                                                   \
                    }                                                       \
                }                                                           \
            }                                                               \
            break;                                                          \
    }

/* index from a single ND (float type) argument */
#define IA_SINGLE_NDF(CLASS, DT, DP, ODP, IPA, IPB, IPC, IPD)               \
    IA_GET_INPUT(DT, DP);                                                   \
    IA_CREATE_OUTPUT(CLASS, DT, ODP);                                       \
    switch (nd) {                                                           \
        case 1:                                                             \
            for (c1 = ne; c1 > 0; --c1) {                                   \
                i1 = (int) (*IPA++ - 0.5);                                  \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DP[i1];                                            \
            }                                                               \
            break;                                                          \
        case 2:                                                             \
            for (c1 = ne; c1 > 0; --c1) {                                   \
                i1 = (int) (*IPA++ - 0.5);                                  \
                i2 = (int) (*IPB++ - 0.5);                                  \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1) ||                                           \
                    (i2 < 0) ||                                             \
                    (i2 >= d2))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DP[i1 + i2 * n2];                                  \
            }                                                               \
            break;                                                          \
        case 3:                                                             \
            for (c1 = ne; c1 > 0; --c1) {                                   \
                i1 = (int) (*IPA++ - 0.5);                                  \
                i2 = (int) (*IPB++ - 0.5);                                  \
                i3 = (int) (*IPC++ - 0.5);                                  \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1) ||                                           \
                    (i2 < 0) ||                                             \
                    (i2 >= d2) ||                                           \
                    (i3 < 0) ||                                             \
                    (i3 >= d3))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DP[i1 + i2 * n2 + i3 * n3];                        \
            }                                                               \
            break;                                                          \
        case 4:                                                             \
            for (c1 = ne; c1 > 0; --c1) {                                   \
                i1 = (int) (*IPA++ - 0.5);                                  \
                i2 = (int) (*IPB++ - 0.5);                                  \
                i3 = (int) (*IPC++ - 0.5);                                  \
                i4 = (int) (*IPD++ - 0.5);                                  \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1) ||                                           \
                    (i2 < 0) ||                                             \
                    (i2 >= d2) ||                                           \
                    (i3 < 0) ||                                             \
                    (i3 >= d3) ||                                           \
                    (i4 < 0) ||                                             \
                    (i4 >= d4))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DP[i1 + i2 * n2 + i3 * n3 + i4 * n4];              \
            }                                                               \
            break;                                                          \
	}

/* index from a single ND (integer type) argument */
#define IA_SINGLE_NDI(CLASS, DT, DP, ODP, IPA, IPB, IPC, IPD)               \
    IA_GET_INPUT(DT, DP);                                                   \
    IA_CREATE_OUTPUT(CLASS, DT, ODP);                                       \
    switch (nd) {                                                           \
        case 1:                                                             \
            for (c1 = ne; c1 > 0; --c1) {                                   \
                i1 = -1 + ((int) *IPA++);                                   \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DP[i1];                                            \
            }                                                               \
            break;                                                          \
        case 2:                                                             \
            for (c1 = ne; c1 > 0; --c1) {                                   \
                i1 = -1 + ((int) *IPA++);                                   \
                i2 = -1 + ((int) *IPB++);                                   \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1) ||                                           \
                    (i2 < 0) ||                                             \
                    (i2 >= d2))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DP[i1 + i2 * n2];                                  \
            }                                                               \
            break;                                                          \
        case 3:                                                             \
            for (c1 = ne; c1 > 0; --c1) {                                   \
                i1 = -1 + ((int) *IPA++);                                   \
                i2 = -1 + ((int) *IPB++);                                   \
                i3 = -1 + ((int) *IPC++);                                   \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1) ||                                           \
                    (i2 < 0) ||                                             \
                    (i2 >= d2) ||                                           \
                    (i3 < 0) ||                                             \
                    (i3 >= d3))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DP[i1 + i2 * n2 + i3 * n3];                        \
            }                                                               \
            break;                                                          \
        case 4:                                                             \
            for (c1 = ne; c1 > 0; --c1) {                                   \
                i1 = -1 + ((int) *IPA++);                                   \
                i2 = -1 + ((int) *IPB++);                                   \
                i3 = -1 + ((int) *IPC++);                                   \
                i4 = -1 + ((int) *IPD++);                                   \
                if ((i1 < 0) ||                                             \
                    (i1 >= d1) ||                                           \
                    (i2 < 0) ||                                             \
                    (i2 >= d2) ||                                           \
                    (i3 < 0) ||                                             \
                    (i3 >= d3) ||                                           \
                    (i4 < 0) ||                                             \
                    (i4 >= d4))                                             \
                    mexErrMsgTxt("Indexing out of bounds.");                \
                *ODP++ = DP[i1 + i2 * n2 + i3 * n3 + i4 * n4];              \
            }                                                               \
            break;                                                          \
	}

/* MEX function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int d1 = 1, d2 = 1, d3 = 1, d4 = 1, nd, ne, c1, c2, c3, i1, i2, i3, i4, n2, n3, n4, od, s1 = 1, s2 = 1, s3 = 1;
    const int *dim, *idim;
    int odim[4];
    mxClassID cid;
    double fr1 = 0.0, fr2 = 0.0, fr3 = 0.0, st1 = 1.0, st2 = 1.0, st3 = 1.0, to1 = 1.0, to2 = 1.0, to3 = 1.0;
    const          double  *dbl, *dblc, *di1 = NULL, *di2 = NULL, *di3 = NULL, *di4 = NULL;
    const          float   *sng, *sngc;
    const   signed int     *sin, *sinc;
    const   signed short   *ssh, *sshc;
    const   signed char    *sch, *schc;
    const unsigned int     *uin, *uinc, *ui1 = NULL, *ui2 = NULL, *ui3 = NULL, *ui4 = NULL;
    const unsigned short   *ush, *ushc;
    const unsigned char    *uch, *uchc;
                   double  *odbl;
                   float   *osng;
            signed int     *osin;
            signed short   *ossh;
            signed char    *osch;
          unsigned int     *ouin;
          unsigned short   *oush;
          unsigned char    *ouch;
    /* char vstr[256]; */

    /* check number of input/output arguments */
    if (nrhs < 2 ||
        nrhs > 5 ||
        nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");

    /* check input array */
    if (mxIsComplex(*prhs) ||
        mxIsSparse(*prhs) ||
        !mxIsNumeric(*prhs))
        mexErrMsgTxt("Input array must be real, full and numeric.");
	nd = mxGetNumberOfDimensions(*prhs);
    if (nd > 4)
        mexErrMsgTxt("Input array cannot be more than 4-D.");

    /* get number of elements and dims */
    ne = mxGetNumberOfElements(*prhs);
	dim = mxGetDimensions(*prhs);

    /* if ne either = dim[0] or dim[1] -> assume 1D */
    if ((ne == *dim) ||
        (ne == dim[1]))
        nd = 1;

    /* get class */
    cid = mxGetClassID(*prhs);

    /* input dimensions */
    d1 = *dim;
    d2 = dim[1];
    if (nd == 1) {
        d1 = ne;
        d2 = 1;
    }
    if (nd > 2) {
        d3 = dim[2];
        if (nd > 3)
            d4 = dim[3];
    }
    n2 = d1;
    n3 = n2 * d2;
    n4 = n3 * d3;

    /* single index expression with NxD syntax */
    if ((nrhs == 2) &&
        ((mxIsDouble(prhs[1])) ||
         (mxIsSingle(prhs[1])) ||
         (mxIsUint32(prhs[1])) ||
         (mxIsInt32(prhs[1]))) &&
        (mxGetNumberOfDimensions(prhs[1]) < 3) &&
        (mxGetN(prhs[1]) == nd)) {

        /* special case: range specification (up to 3D) */
        if ((nd < 4) &&
            mxIsDouble(prhs[1]) &&
            (mxGetM(prhs[1]) == 4) &&
            mxIsInf(*((const double *) mxGetData(prhs[1])))) {

            /* get range pointer */
            di1 = (const double *) mxGetData(prhs[1]);

            /* check whether range is valid */
            if (((nd > 1) &&
                 !mxIsInf(di1[4])) ||
                ((nd > 2) &&
                 !mxIsInf(di1[8])))
                mexErrMsgTxt("Invalid range specification.");

            /* create pointers accordingly */
            ++di1;
            st1 = di1[1];
            if ((st1 == 0.0) ||
                mxIsInf(st1) ||
                mxIsNaN(st1))
                mexErrMsgTxt("Invalid step size for dim1.");
            to1 = di1[2];
            if (nd > 1) {
                di2 = &di1[4];
                st2 = di2[1];
                if ((st2 == 0.0) ||
                    mxIsInf(st2) ||
                    mxIsNaN(st2))
                    mexErrMsgTxt("Invalid step size for dim2.");
                to2 = di2[2];
                if (nd > 2) {
                    di3 = &di2[4];
                    st3 = di3[1];
                    if ((st3 == 0.0) ||
                        mxIsInf(st3) ||
                        mxIsNaN(st3))
                        mexErrMsgTxt("Invalid step size for dim3.");
                    to3 = di3[2];
                }
            }

            /* for each dimension compution output size */
            c3 = c2 = c1 = 0;
            if (st1 > 0.0)
                for (fr1 = *di1; fr1 <= to1; ++c1)
                    fr1 += st1;
            else
                for (fr1 = *di1; fr1 >= to1; ++c1)
                    fr1 += st1;
            if (nd > 1) {
                if (st2 > 0.0)
                    for (fr2 = *di2; fr2 <= to2; ++c2)
                        fr2 += st2;
                else
                    for (fr2 = *di2; fr2 >= to2; ++c2)
                        fr2 += st2;
                s2 = c1;
            }
            if (nd > 2) {
                if (st3 > 0.0)
                    for (fr3 = *di3; fr3 <= to3; ++c3)
                        fr3 += st3;
                else
                    for (fr3 = *di3; fr3 >= to3; ++c3)
                        fr3 += st3;
                s3 = c2 * s2;
            }

            /* force second dim to 1 if nd == 1 */
            if (nd == 1)
                c2 = 1;

            /* set output size */
            *odim = c1;
            odim[1] = c2;
            odim[2] = c3;

            /* output dimensions */
            od = (nd > 1) ? nd : 2;

            /* switch over input class */
            switch (cid) {
                case mxDOUBLE_CLASS:
                    IA_RANGE_ND( mxDOUBLE_CLASS,          double, dbl, odbl); break;
                case mxSINGLE_CLASS:
                    IA_RANGE_ND( mxSINGLE_CLASS,          float , sng, osng); break;
                case  mxUINT32_CLASS:
                    IA_RANGE_ND( mxUINT32_CLASS, unsigned int   , uin, ouin); break;
                case  mxUINT16_CLASS:
                    IA_RANGE_ND( mxUINT16_CLASS, unsigned short , ush, oush); break;
                case   mxUINT8_CLASS:
                    IA_RANGE_ND(  mxUINT8_CLASS, unsigned char  , uch, ouch); break;
                case   mxINT32_CLASS:
                    IA_RANGE_ND(  mxINT32_CLASS,   signed int   , sin, osin); break;
                case   mxINT16_CLASS:
                    IA_RANGE_ND(  mxINT16_CLASS,   signed short , ssh, ossh); break;
                case    mxINT8_CLASS:
                    IA_RANGE_ND(   mxINT8_CLASS,   signed char  , sch, osch); break;
                case mxLOGICAL_CLASS:
                    IA_RANGE_ND(mxLOGICAL_CLASS, unsigned char  , uch, ouch); break;
                default:
                    mexErrMsgTxt("Unsupported class.");
                    break;
            }

            /* return early */
            return;
        }

        /* output dimensions */
        od = 2;
        ne = mxGetM(prhs[1]);
        *odim = ne;
        odim[1] = 1;

        /* depending on index class ID */
        switch (mxGetClassID(prhs[1])) {

            /* DOUBLE index class */
            case mxDOUBLE_CLASS:

                /* get indexing array pointers */
                di1 = (const double *) mxGetData(prhs[1]);
                if (di1 == NULL)
                    mexErrMsgTxt("Error getting indexing array pointer.");
                if (nd > 1) {
                    di2 = &di1[ne];
                    if (nd > 2) {
                        di3 = &di2[ne];
                        if (nd > 3)
                            di4 = &di3[ne];
                    }
                }

                /* depending on data class ID */
                switch (cid) {
                    case  mxDOUBLE_CLASS:
                        IA_SINGLE_NDF( mxDOUBLE_CLASS,          double, dbl, odbl, di1, di2, di3, di4); break;
                    case  mxSINGLE_CLASS:
                        IA_SINGLE_NDF( mxSINGLE_CLASS,          float , sng, osng, di1, di2, di3, di4); break;
                    case  mxUINT32_CLASS:
                        IA_SINGLE_NDF( mxUINT32_CLASS, unsigned int   , uin, ouin, di1, di2, di3, di4); break;
                    case  mxUINT16_CLASS:
                        IA_SINGLE_NDF( mxUINT16_CLASS, unsigned short , ush, oush, di1, di2, di3, di4); break;
                    case   mxUINT8_CLASS:
                        IA_SINGLE_NDF(  mxUINT8_CLASS, unsigned char  , uch, ouch, di1, di2, di3, di4); break;
                    case   mxINT32_CLASS:
                        IA_SINGLE_NDF(  mxINT32_CLASS,   signed int   , sin, osin, di1, di2, di3, di4); break;
                    case   mxINT16_CLASS:
                        IA_SINGLE_NDF(  mxINT16_CLASS,   signed short , ssh, ossh, di1, di2, di3, di4); break;
                    case    mxINT8_CLASS:
                        IA_SINGLE_NDF(   mxINT8_CLASS,   signed char  , sch, osch, di1, di2, di3, di4); break;
                    case mxLOGICAL_CLASS:
                        IA_SINGLE_NDF(mxLOGICAL_CLASS, unsigned char  , uch, ouch, di1, di2, di3, di4); break;
                    default:
                        mexErrMsgTxt("Unsupported class.");
                        break;
                }
                break;

            /* UINT32 index */
            case mxUINT32_CLASS:

                /* get indexing array pointers */
                ui1 = (const unsigned int *) mxGetData(prhs[1]);
                if (ui1 == NULL)
                    mexErrMsgTxt("Error getting indexing array pointer.");
                if (nd > 1) {
                    ui2 = &ui1[ne];
                    if (nd > 2) {
                        ui3 = &ui2[ne];
                        if (nd > 3)
                            ui4 = &ui3[ne];
                    }
                }

                /* depending on data class ID */
                switch (cid) {
                    case  mxDOUBLE_CLASS:
                        IA_SINGLE_NDI( mxDOUBLE_CLASS,          double, dbl, odbl, ui1, ui2, ui3, ui4); break;
                    case  mxSINGLE_CLASS:
                        IA_SINGLE_NDI( mxSINGLE_CLASS,          float , sng, osng, ui1, ui2, ui3, ui4); break;
                    case  mxUINT32_CLASS:
                        IA_SINGLE_NDI( mxUINT32_CLASS, unsigned int   , uin, ouin, ui1, ui2, ui3, ui4); break;
                    case  mxUINT16_CLASS:
                        IA_SINGLE_NDI( mxUINT16_CLASS, unsigned short , ush, oush, ui1, ui2, ui3, ui4); break;
                    case   mxUINT8_CLASS:
                        IA_SINGLE_NDI(  mxUINT8_CLASS, unsigned char  , uch, ouch, ui1, ui2, ui3, ui4); break;
                    case   mxINT32_CLASS:
                        IA_SINGLE_NDI(  mxINT32_CLASS,   signed int   , sin, osin, ui1, ui2, ui3, ui4); break;
                    case   mxINT16_CLASS:
                        IA_SINGLE_NDI(  mxINT16_CLASS,   signed short , ssh, ossh, ui1, ui2, ui3, ui4); break;
                    case    mxINT8_CLASS:
                        IA_SINGLE_NDI(   mxINT8_CLASS,   signed char  , sch, osch, ui1, ui2, ui3, ui4); break;
                    case mxLOGICAL_CLASS:
                        IA_SINGLE_NDI(mxLOGICAL_CLASS, unsigned char  , uch, ouch, ui1, ui2, ui3, ui4); break;
                    default:
                        mexErrMsgTxt("Unsupported class.");
                        break;
                }
                break;
            default:
                mexErrMsgTxt("Unsupported class.");
                break;
        }

        /* return early */
        return;

    /* interpret indexes as multiple dims */
    }

    /* the number of indices must match the number of dimensions */
    if ((nd > 1) &&
        ((nrhs - 1) != nd))
        mexErrMsgTxt("Invalid number of indexing dimensions.");

    /* depending on the number of input dimensions */
    switch (nd) {

        /* one-dimensional input */
        case 1:

            /* check that number of input dimensions isn't too much */
            od = mxGetNumberOfDimensions(prhs[1]);
            if (od > 4)
                mexErrMsgTxt("Cannot expand beyond 4D.");

            /* copy dimensions */
            idim = mxGetDimensions(prhs[1]);
            for (c1 = 0; c1 < od; ++c1)
                odim[c1] = *idim++;

            /* get number of elements */
            c1 = mxGetNumberOfElements(prhs[1]);

            /* simply index into 1D array, depending on index class ID */
            switch (mxGetClassID(prhs[1])) {

                /* DOUBLE index class */
                case mxDOUBLE_CLASS:

                    /* get indexing array pointers */
                    di1 = (const double *) mxGetData(prhs[1]);
                    if (di1 == NULL)
                        mexErrMsgTxt("Error getting indexing array pointer.");

                    /* depending on data class ID */
                    switch (cid) {
                        case  mxDOUBLE_CLASS:
                            IA_EXP1D( mxDOUBLE_CLASS,          double, dbl, odbl, di1); break;
                        case  mxSINGLE_CLASS:
                            IA_EXP1D( mxSINGLE_CLASS,          float , sng, osng, di1); break;
                        case  mxUINT32_CLASS:
                            IA_EXP1D( mxUINT32_CLASS, unsigned int   , uin, ouin, di1); break;
                        case  mxUINT16_CLASS:
                            IA_EXP1D( mxUINT16_CLASS, unsigned short , ush, oush, di1); break;
                        case   mxUINT8_CLASS:
                            IA_EXP1D(  mxUINT8_CLASS, unsigned char  , uch, ouch, di1); break;
                        case   mxINT32_CLASS:
                            IA_EXP1D(  mxINT32_CLASS,   signed int   , sin, osin, di1); break;
                        case   mxINT16_CLASS:
                            IA_EXP1D(  mxINT16_CLASS,   signed short , ssh, ossh, di1); break;
                        case    mxINT8_CLASS:
                            IA_EXP1D(   mxINT8_CLASS,   signed char  , sch, osch, di1); break;
                        case mxLOGICAL_CLASS:
                            IA_EXP1D(mxLOGICAL_CLASS, unsigned char  , uch, ouch, di1); break;
                        default:
                            mexErrMsgTxt("Unsupported class.");
                            break;
                    }
                    break;

                /* UINT32 index */
                case mxUINT32_CLASS:

                    /* get indexing array pointers */
                    ui1 = (const unsigned int *) mxGetData(prhs[1]);
                    if (ui1 == NULL)
                        mexErrMsgTxt("Error getting indexing array pointer.");

                    /* depending on data class ID */
                    switch (cid) {
                        case  mxDOUBLE_CLASS:
                            IA_EXP1D( mxDOUBLE_CLASS,          double, dbl, odbl, ui1); break;
                        case  mxSINGLE_CLASS:
                            IA_EXP1D( mxSINGLE_CLASS,          float , sng, osng, ui1); break;
                        case  mxUINT32_CLASS:
                            IA_EXP1D( mxUINT32_CLASS, unsigned int   , uin, ouin, ui1); break;
                        case  mxUINT16_CLASS:
                            IA_EXP1D( mxUINT16_CLASS, unsigned short , ush, oush, ui1); break;
                        case   mxUINT8_CLASS:
                            IA_EXP1D(  mxUINT8_CLASS, unsigned char  , uch, ouch, ui1); break;
                        case   mxINT32_CLASS:
                            IA_EXP1D(  mxINT32_CLASS,   signed int   , sin, osin, ui1); break;
                        case   mxINT16_CLASS:
                            IA_EXP1D(  mxINT16_CLASS,   signed short , ssh, ossh, ui1); break;
                        case    mxINT8_CLASS:
                            IA_EXP1D(   mxINT8_CLASS,   signed char  , sch, osch, ui1); break;
                        case mxLOGICAL_CLASS:
                            IA_EXP1D(mxLOGICAL_CLASS, unsigned char  , uch, ouch, ui1); break;
                        default:
                            mexErrMsgTxt("Unsupported class.");
                            break;
                    }
                    break;
                default:
                    mexErrMsgTxt("Unsupported class.");
                    break;
            }

            break;

        /* two-dimensional input */
        case 2:

            /* check that number of input dimensions isn't too much */
            od = mxGetNumberOfDimensions(prhs[1]);
            if (od > 4)
                mexErrMsgTxt("Cannot expand beyond 4D.");

            /* get number of elements for first argument */
            c1 = mxGetNumberOfElements(prhs[1]);

            /* for now only supported if first argument is the expansion */
            idim = mxGetDimensions(prhs[1]);
            if ((c1 == *idim) ||
                (c1 == idim[1]))
                mexErrMsgTxt("Currently only expansion of first argument supported.");

            /* test whether second dimension is one! */
            if (idim[1] != 1)
                mexErrMsgTxt("Cannot expand over existing second dimension.");

            /* copy dimensions */
            for (c2 = 0; c2 < od; ++c2)
                odim[c2] = idim[c2];

            /* get counters for dimensions */
            s1 = *idim;
            s2 = mxGetNumberOfElements(prhs[2]);
            s3 = (int) (((double) c1) / ((double) *idim));

            /* copy second dimension from second argument */
            odim[1] = s2;

            /* re-set second dim size */
            n2 = mxGetN(*prhs);

            /* simply index into 1D array, depending on index class ID */
            switch (mxGetClassID(prhs[1])) {

                /* DOUBLE index class */
                case mxDOUBLE_CLASS:

                    /* get indexing array pointers */
                    di1 = (const double *) mxGetData(prhs[1]);
                    if (di1 == NULL)
                        mexErrMsgTxt("Error getting 1st indexing array pointer.");
                    di2 = (const double *) mxGetData(prhs[2]);
                    if (di2 == NULL)
                        mexErrMsgTxt("Error getting 2nd indexing array pointer.");

                    /* depending on data class ID */
                    switch (cid) {
                        case  mxDOUBLE_CLASS:
                            IA_EXP2D13( mxDOUBLE_CLASS,          double, dbl, dblc, odbl, di1, di2, di3, di4); break;
                        case  mxSINGLE_CLASS:
                            IA_EXP2D13( mxSINGLE_CLASS,          float , sng, sngc, osng, di1, di2, di3, di4); break;
                        case  mxUINT32_CLASS:
                            IA_EXP2D13( mxUINT32_CLASS, unsigned int   , uin, uinc, ouin, di1, di2, di3, di4); break;
                        case  mxUINT16_CLASS:
                            IA_EXP2D13( mxUINT16_CLASS, unsigned short , ush, ushc, oush, di1, di2, di3, di4); break;
                        case   mxUINT8_CLASS:
                            IA_EXP2D13(  mxUINT8_CLASS, unsigned char  , uch, uchc, ouch, di1, di2, di3, di4); break;
                        case   mxINT32_CLASS:
                            IA_EXP2D13(  mxINT32_CLASS,   signed int   , sin, sinc, osin, di1, di2, di3, di4); break;
                        case   mxINT16_CLASS:
                            IA_EXP2D13(  mxINT16_CLASS,   signed short , ssh, sshc, ossh, di1, di2, di3, di4); break;
                        case    mxINT8_CLASS:
                            IA_EXP2D13(   mxINT8_CLASS,   signed char  , sch, schc, osch, di1, di2, di3, di4); break;
                        case mxLOGICAL_CLASS:
                            IA_EXP2D13(mxLOGICAL_CLASS, unsigned char  , uch, uchc, ouch, di1, di2, di3, di4); break;
                        default:
                            mexErrMsgTxt("Unsupported class.");
                            break;
                    }
                    break;

                /* UINT32 index */
                case mxUINT32_CLASS:

                    /* get indexing array pointers */
                    ui1 = (const unsigned int *) mxGetData(prhs[1]);
                    if (ui1 == NULL)
                        mexErrMsgTxt("Error getting indexing array pointer.");

                    /* depending on data class ID */
                    switch (cid) {
                        case  mxDOUBLE_CLASS:
                            IA_EXP2D13( mxDOUBLE_CLASS,          double, dbl, dblc, odbl, ui1, ui2, ui3, ui4); break;
                        case  mxSINGLE_CLASS:
                            IA_EXP2D13( mxSINGLE_CLASS,          float , sng, sngc, osng, ui1, ui2, ui3, ui4); break;
                        case  mxUINT32_CLASS:
                            IA_EXP2D13( mxUINT32_CLASS, unsigned int   , uin, uinc, ouin, ui1, ui2, ui3, ui4); break;
                        case  mxUINT16_CLASS:
                            IA_EXP2D13( mxUINT16_CLASS, unsigned short , ush, ushc, oush, ui1, ui2, ui3, ui4); break;
                        case   mxUINT8_CLASS:
                            IA_EXP2D13(  mxUINT8_CLASS, unsigned char  , uch, uchc, ouch, ui1, ui2, ui3, ui4); break;
                        case   mxINT32_CLASS:
                            IA_EXP2D13(  mxINT32_CLASS,   signed int   , sin, sinc, osin, ui1, ui2, ui3, ui4); break;
                        case   mxINT16_CLASS:
                            IA_EXP2D13(  mxINT16_CLASS,   signed short , ssh, sshc, ossh, ui1, ui2, ui3, ui4); break;
                        case    mxINT8_CLASS:
                            IA_EXP2D13(   mxINT8_CLASS,   signed char  , sch, schc, osch, ui1, ui2, ui3, ui4); break;
                        case mxLOGICAL_CLASS:
                            IA_EXP2D13(mxLOGICAL_CLASS, unsigned char  , uch, uchc, ouch, ui1, ui2, ui3, ui4); break;
                        default:
                            mexErrMsgTxt("Unsupported class.");
                            break;
                    }
                    break;
                default:
                    mexErrMsgTxt("Unsupported class.");
                    break;
            }

            break;

        /* three-dimensional input */
        case 3:

            /* for now unsupported */
            mexErrMsgTxt("Not yet supported.");
            break;

        default:
            mexErrMsgTxt("Dimensional expansion only for up to 3D input.");
            break;
    }
}

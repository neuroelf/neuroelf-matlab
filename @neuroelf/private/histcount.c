/*

counts of elements with value (or monotonic range)

FORMAT:       hc = histcount(v, from, to, [step, [v2, from2, to2, step2]]);

Input fields:

      v, v2       value array
      from, from2 range definition (begin)
      to, to2     range definition (end)
      step, step2 range definition (stepsize, default: 1)

Output fields:

      hc          histogram count

Note: if v2 is set to a 1x1 double, it gives the dimension along which
      v is histogramed

Note: if v2 is set to a double array with the same size as v1 and from2
      is set to empty, v2 is treated as a weighting array; to combine
      dimension and weighting, give dimension first, then weighting as
      from2

% Version:  v0.9d
% Build:    14052616
% Date:     May-26 2014, 4:51 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2010 - 2014, Jochen Weber
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

#include "isinfnan.h"
#include "mex.h"
#include "math.h"

#define HIST1D(DATAP)                   \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) ((val1 - mn1) * f1); \
        ++(odbl[i1]);                   \
    }

#define HIST1DD(DATAP)                  \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) ((val1 - mn1) * f1); \
        ++(odbl[i1]);                   \
    }

#define HIST1DNOF(DATAP)                \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) (val1 - mn1);        \
        ++(odbl[i1]);                   \
    }

#define HIST1DNOFD(DATAP)               \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) (val1 - mn1);        \
        ++(odbl[i1]);                   \
    }

#define HIST1DNOO(DATAP)                \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        if (val1 < 0.0)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) (val1 * f1);         \
        ++(odbl[i1]);                   \
    }

#define HIST1DNOOD(DATAP)               \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        if (val1 < 0.0)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) (val1 * f1);         \
        ++(odbl[i1]);                   \
    }

#define HIST1DNOOF(DATAP)               \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        if (val1 < 0.0)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) val1;                \
        ++(odbl[i1]);                   \
    }

#define HIST1DNOOFD(DATAP)              \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        if (val1 < 0.0)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) val1;                \
        ++(odbl[i1]);                   \
    }

#define HIST2D(DATAP1, DATAP2)          \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP1++;      \
        val2 = (double) *DATAP2++;      \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        if (val2 < mn2)                 \
            continue;                   \
        else if(val2 > mx2)             \
            val2 = mx2;                 \
        i1 = (int) ((val1 - mn1) * f1); \
        i2 = (int) ((val2 - mn2) * f2); \
        ++(odbl[i1 + m * i2]);          \
    }

#define HIST2DD(DATAP1, DATAP2)         \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP1++;      \
        val2 = (double) *DATAP2++;      \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        if (val2 < mn2)                 \
            continue;                   \
        else if(val2 > mx2)             \
            val2 = mx2;                 \
        i1 = (int) ((val1 - mn1) * f1); \
        i2 = (int) ((val2 - mn2) * f2); \
        ++(odbl[i1 + m * i2]);          \
    }

#define HIST2DS(DATAP1)                 \
    switch (cid2) {                     \
        case mxDOUBLE_CLASS:                                       HIST2DD(DATAP1,dbl2); break; \
        case mxSINGLE_CLASS: sng2 = (const          float *) dbl2; HIST2DD(DATAP1,sng2); break; \
        case mxUINT32_CLASS:  ul2 = (const unsigned int   *) dbl2; HIST2D(DATAP1,  ul2); break; \
        case  mxINT32_CLASS:  sl2 = (const   signed int   *) dbl2; HIST2D(DATAP1,  sl2); break; \
        case mxUINT16_CLASS:  us2 = (const unsigned short *) dbl2; HIST2D(DATAP1,  us2); break; \
        case  mxINT16_CLASS:  ss2 = (const   signed short *) dbl2; HIST2D(DATAP1,  ss2); break; \
        case mxLOGICAL_CLASS:           \
        case  mxUINT8_CLASS:  uc2 = (const unsigned char  *) dbl2; HIST2D(DATAP1,  uc2); break; \
        case   mxINT8_CLASS:  sc2 = (const   signed char  *) dbl2; HIST2D(DATAP1,  sc2); break; \
        default:                        \
            mexErrMsgTxt("Invalid input class.");                       \
    }

#define HIST2DSD(DATAP1)                 \
    switch (cid2) {                     \
        case mxDOUBLE_CLASS:                                       HIST2DD(DATAP1,dbl2); break; \
        case mxSINGLE_CLASS: sng2 = (const          float *) dbl2; HIST2DD(DATAP1,sng2); break; \
        case mxUINT32_CLASS:  ul2 = (const unsigned int   *) dbl2; HIST2DD(DATAP1, ul2); break; \
        case  mxINT32_CLASS:  sl2 = (const   signed int   *) dbl2; HIST2DD(DATAP1, sl2); break; \
        case mxUINT16_CLASS:  us2 = (const unsigned short *) dbl2; HIST2DD(DATAP1, us2); break; \
        case  mxINT16_CLASS:  ss2 = (const   signed short *) dbl2; HIST2DD(DATAP1, ss2); break; \
        case mxLOGICAL_CLASS:           \
        case  mxUINT8_CLASS:  uc2 = (const unsigned char  *) dbl2; HIST2DD(DATAP1, uc2); break; \
        case   mxINT8_CLASS:  sc2 = (const   signed char  *) dbl2; HIST2DD(DATAP1, sc2); break; \
        default:                        \
            mexErrMsgTxt("Invalid input class.");                       \
    }

#define HISTND(DATAP)                               \
    i3 = i2 = i1 = 0;                               \
    if (m == 1) {                                   \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c) {              \
                val1 = (double) *DATAP++;           \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) ((val1 - mn1) * f1);     \
                ++(odbl[ii]);                       \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i2) {        \
                if (i2 >= n) {                      \
                    i2 = 0;                         \
                    i3 += o;                        \
                }                                   \
                val1 = (double) *DATAP++;           \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) ((val1 - mn1) * f1);     \
                ++(odbl[ii + i3]);                  \
            }                                       \
        }                                           \
    } else {                                        \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m)                        \
                    i1 = 0;                         \
                val1 = (double) *DATAP++;           \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) ((val1 - mn1) * f1);     \
                ++(odbl[i1 + ii * m]);              \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m) {                      \
                    i1 = 0;                         \
                    ++i2;                           \
                    if (i2 >= n) {                  \
                        i2 = 0;                     \
                        i3 += o;                    \
                    }                               \
                }                                   \
                val1 = (double) *DATAP++;           \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) ((val1 - mn1) * f1);     \
                ++(odbl[i1 + ii * m + i3]);         \
            }                                       \
        }                                           \
    }

#define HISTNDNOF(DATAP)                            \
    i3 = i2 = i1 = 0;                               \
    if (m == 1) {                                   \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c) {              \
                val1 = (double) *DATAP++;           \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 - mn1);            \
                ++(odbl[ii]);                       \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i2) {        \
                if (i2 >= n) {                      \
                    i2 = 0;                         \
                    i3 += o;                        \
                }                                   \
                val1 = (double) *DATAP++;           \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 - mn1);            \
                ++(odbl[ii + i3]);                  \
            }                                       \
        }                                           \
    } else {                                        \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m)                        \
                    i1 = 0;                         \
                val1 = (double) *DATAP++;           \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 - mn1);            \
                ++(odbl[i1 + ii * m]);              \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m) {                      \
                    i1 = 0;                         \
                    ++i2;                           \
                    if (i2 >= n) {                  \
                        i2 = 0;                     \
                        i3 += o;                    \
                    }                               \
                }                                   \
                val1 = (double) *DATAP++;           \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 - mn1);            \
                ++(odbl[i1 + ii * m + i3]);         \
            }                                       \
        }                                           \
    }

#define HISTNDNOO(DATAP)                            \
    i3 = i2 = i1 = 0;                               \
    if (m == 1) {                                   \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c) {              \
                val1 = (double) *DATAP++;           \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 * f1);             \
                ++(odbl[ii]);                       \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i2) {        \
                if (i2 >= n) {                      \
                    i2 = 0;                         \
                    i3 += o;                        \
                }                                   \
                val1 = (double) *DATAP++;           \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 * f1);             \
                ++(odbl[ii + i3]);                  \
            }                                       \
        }                                           \
    } else {                                        \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m)                        \
                    i1 = 0;                         \
                val1 = (double) *DATAP++;           \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 * f1);             \
                ++(odbl[i1 + ii * m]);              \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m) {                      \
                    i1 = 0;                         \
                    ++i2;                           \
                    if (i2 >= n) {                  \
                        i2 = 0;                     \
                        i3 += o;                    \
                    }                               \
                }                                   \
                val1 = (double) *DATAP++;           \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 * f1);             \
                ++(odbl[i1 + ii * m + i3]);         \
            }                                       \
        }                                           \
    }

#define HISTNDNOOF(DATAP)                           \
    i3 = i2 = i1 = 0;                               \
    if (m == 1) {                                   \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c) {              \
                val1 = (double) *DATAP++;           \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) val1;                    \
                ++(odbl[ii]);                       \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i2) {        \
                if (i2 >= n) {                      \
                    i2 = 0;                         \
                    i3 += o;                        \
                }                                   \
                val1 = (double) *DATAP++;           \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) val1;                    \
                ++(odbl[ii + i3]);                  \
            }                                       \
        }                                           \
    } else {                                        \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m)                        \
                    i1 = 0;                         \
                val1 = (double) *DATAP++;           \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) val1;                    \
                ++(odbl[i1 + ii * m]);              \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m) {                      \
                    i1 = 0;                         \
                    ++i2;                           \
                    if (i2 >= n) {                  \
                        i2 = 0;                     \
                        i3 += o;                    \
                    }                               \
                }                                   \
                val1 = (double) *DATAP++;           \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) val1;                    \
                ++(odbl[i1 + ii * m + i3]);         \
            }                                       \
        }                                           \
    }

#define WHIST1D(DATAP, WEIGHTP)         \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        val2 = *WEIGHTP++;              \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) ((val1 - mn1) * f1); \
        odbl[i1] += val2;               \
    }

#define WHIST1DD(DATAP, WEIGHTP)        \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        val2 = *WEIGHTP++;              \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) ((val1 - mn1) * f1); \
        odbl[i1] += val2;               \
    }

#define WHIST1DNOF(DATAP, WEIGHTP)      \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        val2 = *WEIGHTP++;              \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) (val1 - mn1);        \
        odbl[i1] += val2;               \
    }

#define WHIST1DNOFD(DATAP, WEIGHTP)     \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        val2 = *WEIGHTP++;              \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < mn1)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) (val1 - mn1);        \
        odbl[i1] += val2;               \
    }

#define WHIST1DNOO(DATAP, WEIGHTP)      \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        val2 = *WEIGHTP++;              \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < 0.0)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) (val1 * f1);         \
        odbl[i1] += val2;               \
    }

#define WHIST1DNOOD(DATAP, WEIGHTP)     \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        val2 = *WEIGHTP++;              \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < 0.0)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) (val1 * f1);         \
        odbl[i1] += val2;               \
    }

#define WHIST1DNOOF(DATAP, WEIGHTP)     \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        val2 = *WEIGHTP++;              \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < 0.0)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) val1;                \
        odbl[i1] += val2;               \
    }

#define WHIST1DNOOFD(DATAP, WEIGHTP)    \
    for (c = ne; c > 0; --c) {          \
        val1 = (double) *DATAP++;       \
        val2 = *WEIGHTP++;              \
        IF_IS_BAD_VAL(val1)             \
            continue;                   \
        IF_IS_BAD_VAL(val2)             \
            continue;                   \
        if (val1 < 0.0)                 \
            continue;                   \
        else if(val1 > mx1)             \
            val1 = mx1;                 \
        i1 = (int) val1;                \
        odbl[i1] += val2;               \
    }


#define WHISTND(DATAP, WEIGHTP)                     \
    i3 = i2 = i1 = 0;                               \
    if (m == 1) {                                   \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c) {              \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) ((val1 - mn1) * f1);     \
                odbl[ii] += val2;                   \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i2) {        \
                if (i2 >= n) {                      \
                    i2 = 0;                         \
                    i3 += o;                        \
                }                                   \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) ((val1 - mn1) * f1);     \
                odbl[ii + i3] += val2;              \
            }                                       \
        }                                           \
    } else {                                        \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m)                        \
                    i1 = 0;                         \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) ((val1 - mn1) * f1);     \
                odbl[i1 + ii * m] += val2;          \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m) {                      \
                    i1 = 0;                         \
                    ++i2;                           \
                    if (i2 >= n) {                  \
                        i2 = 0;                     \
                        i3 += o;                    \
                    }                               \
                }                                   \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) ((val1 - mn1) * f1);     \
                odbl[i1 + ii * m + i3] += val2;     \
            }                                       \
        }                                           \
    }

#define WHISTNDNOF(DATAP, WEIGHTP)                  \
    i3 = i2 = i1 = 0;                               \
    if (m == 1) {                                   \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c) {              \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 - mn1);            \
                odbl[ii] += val2;                   \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i2) {        \
                if (i2 >= n) {                      \
                    i2 = 0;                         \
                    i3 += o;                        \
                }                                   \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 - mn1);            \
                odbl[ii + i3] += val2;              \
            }                                       \
        }                                           \
    } else {                                        \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m)                        \
                    i1 = 0;                         \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 - mn1);            \
                odbl[i1 + ii * m] += val2;          \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m) {                      \
                    i1 = 0;                         \
                    ++i2;                           \
                    if (i2 >= n) {                  \
                        i2 = 0;                     \
                        i3 += o;                    \
                    }                               \
                }                                   \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < mn1)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 - mn1);            \
                odbl[i1 + ii * m + i3] += val2;     \
            }                                       \
        }                                           \
    }

#define WHISTNDNOO(DATAP, WEIGHTP)                  \
    i3 = i2 = i1 = 0;                               \
    if (m == 1) {                                   \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c) {              \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 * f1);             \
                odbl[ii] += val2;                   \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i2) {        \
                if (i2 >= n) {                      \
                    i2 = 0;                         \
                    i3 += o;                        \
                }                                   \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 * f1);             \
                odbl[ii + i3] += val2;              \
            }                                       \
        }                                           \
    } else {                                        \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m)                        \
                    i1 = 0;                         \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 * f1);             \
                odbl[i1 + ii * m] += val2;          \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m) {                      \
                    i1 = 0;                         \
                    ++i2;                           \
                    if (i2 >= n) {                  \
                        i2 = 0;                     \
                        i3 += o;                    \
                    }                               \
                }                                   \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) (val1 * f1);             \
                odbl[i1 + ii * m + i3] += val2;     \
            }                                       \
        }                                           \
    }

#define WHISTNDNOOF(DATAP, WEIGHTP)                 \
    i3 = i2 = i1 = 0;                               \
    if (m == 1) {                                   \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c) {              \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) val1;                    \
                odbl[ii] += val2;                   \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i2) {        \
                if (i2 >= n) {                      \
                    i2 = 0;                         \
                    i3 += o;                        \
                }                                   \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) val1;                    \
                odbl[ii + i3] += val2;              \
            }                                       \
        }                                           \
    } else {                                        \
        if (o == 1) {                               \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m)                        \
                    i1 = 0;                         \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) val1;                    \
                odbl[i1 + ii * m] += val2;          \
            }                                       \
        } else {                                    \
            for (c = 0; c < ne; ++c, ++i1) {        \
                if (i1 >= m) {                      \
                    i1 = 0;                         \
                    ++i2;                           \
                    if (i2 >= n) {                  \
                        i2 = 0;                     \
                        i3 += o;                    \
                    }                               \
                }                                   \
                val1 = (double) *DATAP++;           \
                val2 = (double) *WEIGHTP++;         \
                IF_IS_BAD_VAL(val2)                 \
                    continue;                       \
                if (val1 < 0.0)                     \
                    continue;                       \
                else if(val1 > mx1)                 \
                    val1 = mx1;                     \
                ii = (int) val1;                    \
                odbl[i1 + ii * m + i3] += val2;     \
            }                                       \
        }                                           \
    }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxClassID cid1, cid2 = mxDOUBLE_CLASS;
	const int *dim1, *dim2;
    const signed char *sc1, *sc2;
    const unsigned char *uc1, *uc2;
    const signed short *ss1, *ss2;
    const unsigned short *us1, *us2;
    const signed int *sl1, *sl2;
    const unsigned int *ul1, *ul2;
    const float *sng1, *sng2;
    const double *dbl1, *dbl2;
    double f1 = 1.0, f2 = 1.0, mn1 = 0.0, mn2 = 0.0, mx1 = 1.0, mx2 = 1.0,
           val1, val2;
    int c, i1, i2, i3, ii, m = 1, n = 1, nd, ne, o;
    int odim[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    double *odbl;
    /* char vstr[256]; */

    /* isinfnan.h variables */
    VARS_FOR_ISINFNAN

    /* test input arguments (first array configuration) */
	if (nrhs < 3 || nrhs > 8 || nlhs > 1)
		mexErrMsgTxt("Bad number of input/output arguments.");

    /* test array1 from:to range */
    if ((!mxIsDouble(prhs[1])) ||
        (mxGetNumberOfElements(prhs[1]) != 1) ||
        (!mxIsDouble(prhs[2])) ||
        (mxGetNumberOfElements(prhs[2]) != 1))
        mexErrMsgTxt("Range from and to must be 1x1 double scalars.");
    dbl1 = (const double *) mxGetData(prhs[1]);
    dbl2 = (const double *) mxGetData(prhs[2]);
    if ((dbl1 == NULL) ||
        (dbl2 == NULL))
        mexErrMsgTxt("Error getting range from/to data pointer.");
    if (mxIsInf(*dbl1) ||
        mxIsNaN(*dbl1) ||
        mxIsInf(*dbl2) ||
        mxIsNaN(*dbl2))
        mexErrMsgTxt("Bad range from/to value.");

    /* get array1 from:to range */
    mn1 = *dbl1;
    mx1 = *dbl2;
    if (mx1 < mn1)
        mx1 = mn1;

    /* test array1 range factor/stepsize */
    if ((nrhs > 3) &&
        ((!mxIsDouble(prhs[3])) ||
         (mxGetNumberOfElements(prhs[3]) != 1)))
        mexErrMsgTxt("Range stepsize must be a 1x1 double scalar.");

    /* get array1 range factor/stepsize */
    if (nrhs > 3) {
        dbl1 = (const double *) mxGetData(prhs[3]);
        if (dbl1 == NULL)
            mexErrMsgTxt("Error getting range stepsize data pointer.");
        if (mxIsInf(*dbl1) ||
            mxIsNaN(*dbl1) ||
            (*dbl1 <= 0.0))
            mexErrMsgTxt("Bad range stepsize value.");
        f1 = 1.0 / *dbl1;
    }

    /* initialize IS_* macros */
    INIT_INF_NAN_BAD_VAL()

    /* compute array1 required size (number of bins) */
    n = 1 + (int) ((mx1 - mn1) * f1);
    cid1 = mxGetClassID(*prhs);
    dim1 = mxGetDimensions(*prhs);
    nd = mxGetNumberOfDimensions(*prhs);
    ne = mxGetNumberOfElements(*prhs);

    /* additional arguments */
    if (nrhs > 4) {

        /* special case: array2 is a 1x1 double -> create ND histograms! */
        if (mxIsDouble(prhs[4]) &&
            (mxGetNumberOfElements(prhs[4]) == 1)) {

            /* not valid for more than 8D */
            if (nd > 8)
                mexErrMsgTxt("ND histogram only valid for up to 8D.");

            /* get dimension to use */
            dbl1 = (const double *) mxGetData(prhs[4]);
            if (mxIsInf(*dbl1) ||
                mxIsNaN(*dbl1) ||
                (*dbl1 < 1.0) ||
                (*dbl1 > ((double) nd)))
                mexErrMsgTxt("Bad dimension flag.");
            m = (int) *dbl1;

            /* copy input to output dimensions */
            for (c = 0; c < nd; ++c)
                odim[c] = dim1[c];

            /* patch output dimension at correct dim */
            i2 = m - 1;
            odim[i2] = n;

            /* create output */
            *plhs = mxCreateNumericArray(nd, odim, mxDOUBLE_CLASS, mxREAL);
            if (*plhs == NULL)
                mexErrMsgTxt("Error creating output array.");

            /* get output pointer */
            odbl = (double *) mxGetData(*plhs);
            if (odbl == NULL)
                mexErrMsgTxt("Error getting output data pointer.");

            /* compute access factors */
            m = 1;
            for (c = 0; c < i2; ++c)
                m *= odim[c];
            if ((i2 + 1) >= nd)
                o = 1;
            else
                o = m * n;

            /* reset n to original dim */
            n = dim1[i2];

            /* get input pointer (main arrays) */
            dbl1 = (const double *) mxGetData(*prhs);

            /* weighting array given */
            if ((nrhs > 5) &&
                (mxGetNumberOfElements(prhs[5]) == ne) &&
                mxIsDouble(prhs[5])) {

                /* get pointer */
                dbl2 = (const double *) mxGetData(prhs[5]);

                /* factor/offset needed */
                if ((mn1 != 0.0) &&
                    (f1 != 1.0)) {

                    /* switch over classes */
                    switch (cid1) {
                        case mxDOUBLE_CLASS:                                       WHISTND(dbl1, dbl2); break;
                        case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; WHISTND(sng1, dbl2); break;
                        case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; WHISTND( ul1, dbl2); break;
                        case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; WHISTND( sl1, dbl2); break;
                        case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; WHISTND( us1, dbl2); break;
                        case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; WHISTND( ss1, dbl2); break;
                        case mxLOGICAL_CLASS:
                        case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; WHISTND( uc1, dbl2); break;
                        case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; WHISTND( sc1, dbl2); break;
                        default:
                            mexErrMsgTxt("Invalid input class.");
                    }

                /* no factor needed */
                } else if (f1 == 1.0) {

                    /* no offset needed */
                    if (mn1 == 0.0) {

                        /* switch over classes */
                        switch (cid1) {
                            case mxDOUBLE_CLASS:                                       WHISTNDNOOF(dbl1, dbl2); break;
                            case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; WHISTNDNOOF(sng1, dbl2); break;
                            case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; WHISTNDNOOF( ul1, dbl2); break;
                            case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; WHISTNDNOOF( sl1, dbl2); break;
                            case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; WHISTNDNOOF( us1, dbl2); break;
                            case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; WHISTNDNOOF( ss1, dbl2); break;
                            case mxLOGICAL_CLASS:
                            case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; WHISTNDNOOF( uc1, dbl2); break;
                            case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; WHISTNDNOOF( sc1, dbl2); break;
                            default:
                                mexErrMsgTxt("Invalid input class.");
                        }

                    /* only offset needed */
                    } else {

                        /* switch over classes */
                        switch (cid1) {
                            case mxDOUBLE_CLASS:                                       WHISTNDNOF(dbl1, dbl2); break;
                            case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; WHISTNDNOF(sng1, dbl2); break;
                            case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; WHISTNDNOF( ul1, dbl2); break;
                            case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; WHISTNDNOF( sl1, dbl2); break;
                            case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; WHISTNDNOF( us1, dbl2); break;
                            case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; WHISTNDNOF( ss1, dbl2); break;
                            case mxLOGICAL_CLASS:
                            case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; WHISTNDNOF( uc1, dbl2); break;
                            case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; WHISTNDNOF( sc1, dbl2); break;
                            default:
                                mexErrMsgTxt("Invalid input class.");
                        }
                    }

                /* factor needed, but no offset! */
                } else {

                    /* switch over classes */
                    switch (cid1) {
                        case mxDOUBLE_CLASS:                                       WHISTNDNOO(dbl1, dbl2); break;
                        case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; WHISTNDNOO(sng1, dbl2); break;
                        case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; WHISTNDNOO( ul1, dbl2); break;
                        case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; WHISTNDNOO( sl1, dbl2); break;
                        case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; WHISTNDNOO( us1, dbl2); break;
                        case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; WHISTNDNOO( ss1, dbl2); break;
                        case mxLOGICAL_CLASS:
                        case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; WHISTNDNOO( uc1, dbl2); break;
                        case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; WHISTNDNOO( sc1, dbl2); break;
                        default:
                            mexErrMsgTxt("Invalid input class.");
                    }
                }

                /* return already */
                return;
            }

            /* factor/offset needed */
            if ((mn1 != 0.0) &&
                (f1 != 1.0)) {

                /* switch over classes */
                switch (cid1) {
                    case mxDOUBLE_CLASS:                                       HISTND(dbl1); break;
                    case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HISTND(sng1); break;
                    case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HISTND( ul1); break;
                    case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HISTND( sl1); break;
                    case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HISTND( us1); break;
                    case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HISTND( ss1); break;
                    case mxLOGICAL_CLASS:
                    case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HISTND( uc1); break;
                    case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HISTND( sc1); break;
                    default:
                        mexErrMsgTxt("Invalid input class.");
                }

            /* no factor needed */
            } else if (f1 == 1.0) {

                /* no offset needed */
                if (mn1 == 0.0) {

                    /* switch over classes */
                    switch (cid1) {
                        case mxDOUBLE_CLASS:                                       HISTNDNOOF(dbl1); break;
                        case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HISTNDNOOF(sng1); break;
                        case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HISTNDNOOF( ul1); break;
                        case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HISTNDNOOF( sl1); break;
                        case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HISTNDNOOF( us1); break;
                        case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HISTNDNOOF( ss1); break;
                        case mxLOGICAL_CLASS:
                        case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HISTNDNOOF( uc1); break;
                        case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HISTNDNOOF( sc1); break;
                        default:
                            mexErrMsgTxt("Invalid input class.");
                    }

                /* only offset needed */
                } else {

                    /* switch over classes */
                    switch (cid1) {
                        case mxDOUBLE_CLASS:                                       HISTNDNOF(dbl1); break;
                        case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HISTNDNOF(sng1); break;
                        case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HISTNDNOF( ul1); break;
                        case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HISTNDNOF( sl1); break;
                        case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HISTNDNOF( us1); break;
                        case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HISTNDNOF( ss1); break;
                        case mxLOGICAL_CLASS:
                        case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HISTNDNOF( uc1); break;
                        case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HISTNDNOF( sc1); break;
                        default:
                            mexErrMsgTxt("Invalid input class.");
                    }
                }

            /* factor needed, but no offset! */
            } else {

                /* switch over classes */
                switch (cid1) {
                    case mxDOUBLE_CLASS:                                       HISTNDNOO(dbl1); break;
                    case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HISTNDNOO(sng1); break;
                    case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HISTNDNOO( ul1); break;
                    case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HISTNDNOO( sl1); break;
                    case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HISTNDNOO( us1); break;
                    case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HISTNDNOO( ss1); break;
                    case mxLOGICAL_CLASS:
                    case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HISTNDNOO( uc1); break;
                    case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HISTNDNOO( sc1); break;
                    default:
                        mexErrMsgTxt("Invalid input class.");
                }
            }

            /* return already ! */
            return;
        }

        /* standard routine! -> create 2D histogram over combined array1/array2 scatter */
        m = n;
        if (nrhs > 6) {
            if ((!mxIsDouble(prhs[5])) ||
                (mxGetNumberOfElements(prhs[5]) != 1) ||
                (!mxIsDouble(prhs[6])) ||
                (mxGetNumberOfElements(prhs[6]) != 1))
                mexErrMsgTxt("Range from and to must be 1x1 double scalars.");
            dbl1 = (const double *) mxGetData(prhs[5]);
            dbl2 = (const double *) mxGetData(prhs[6]);
            if ((dbl1 == NULL) ||
                (dbl2 == NULL))
                mexErrMsgTxt("Error getting range from/to data pointer.");
            if (mxIsInf(*dbl1) ||
                mxIsNaN(*dbl1) ||
                mxIsInf(*dbl2) ||
                mxIsNaN(*dbl2))
                mexErrMsgTxt("Bad range from/to value.");
            mn2 = *dbl1;
            mx2 = *dbl2;
            if (mx2 < mn2)
                mx2 = mn2;
            if (nrhs > 7) {
                dbl1 = (const double *) mxGetData(prhs[7]);
                if (dbl1 == NULL)
                    mexErrMsgTxt("Error getting range stepsize data pointer.");
                if (mxIsInf(*dbl1) ||
                    mxIsNaN(*dbl1) ||
                    (*dbl1 <= 0.0))
                    mexErrMsgTxt("Bad range stepsize value.");
                f2 = 1.0 / *dbl1;
            }
            n = 1 + (int) ((mx2 - mn2) * f2);
        } else {
            mn2 = mn1;
            mx2 = mx1;
            f2 = f1;
        }
        cid2 = mxGetClassID(prhs[4]);
        i1 = mxGetNumberOfDimensions(prhs[4]);
        i2 = mxGetNumberOfElements(prhs[4]);
        if ((i1 != nd) ||
            (i2 != ne))
            mexErrMsgTxt("Two arrays must match in dimensions/elements.");
        dim2 = mxGetDimensions(prhs[4]);
        for (c = 0; c < nd; ++c)
            if (dim1[c] != dim2[c])
                mexErrMsgTxt("Two arrays must match in size.");
    }

    /* get input pointer (main arrays) */
    dbl1 = (const double *) mxGetData(*prhs);

    /* weighted histogram */
    if (nrhs == 6 &&
        cid2 == mxDOUBLE_CLASS &&
        mxGetNumberOfElements(prhs[5]) == 0)

        /* fix output size back to 1 */
        m = 1;

    /* create output array */
    *plhs = mxCreateDoubleMatrix(m, n, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error creating output array.");
    odbl = (double *) mxGetData(*plhs);
    if (odbl == NULL)
        mexErrMsgTxt("Error getting output data pointer.");

    /* two-dimensional case (or weighted histogram) */
    if (nrhs > 4) {

        /* get input pointer (2nd array) */
        dbl2 = (const double *) mxGetData(prhs[4]);

        /* weighted histogram */
        if (nrhs == 6 &&
            cid2 == mxDOUBLE_CLASS &&
            mxGetNumberOfElements(prhs[5]) == 0) {

            /* factor/offset needed */
            if ((mn1 != 0.0) &&
                (f1 != 1.0)) {

                /* switch over classes */
                switch (cid1) {
                    case mxDOUBLE_CLASS:                                       WHIST1DD(dbl1,dbl2); break;
                    case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; WHIST1DD(sng1,dbl2); break;
                    case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; WHIST1D( ul1, dbl2); break;
                    case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; WHIST1D( sl1, dbl2); break;
                    case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; WHIST1D( us1, dbl2); break;
                    case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; WHIST1D( ss1, dbl2); break;
                    case mxLOGICAL_CLASS:
                    case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; WHIST1D( uc1, dbl2); break;
                    case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; WHIST1D( sc1, dbl2); break;
                    default:
                        mexErrMsgTxt("Invalid input class.");
                }

            /* no factor needed */
            } else if (f1 == 1.0) {

                /* no offset needed */
                if (mn1 == 0.0) {

                    /* switch over classes */
                    switch (cid1) {
                        case mxDOUBLE_CLASS:                                       WHIST1DNOOFD(dbl1,dbl2); break;
                        case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; WHIST1DNOOFD(sng1,dbl2); break;
                        case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; WHIST1DNOOF( ul1, dbl2); break;
                        case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; WHIST1DNOOF( sl1, dbl2); break;
                        case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; WHIST1DNOOF( us1, dbl2); break;
                        case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; WHIST1DNOOF( ss1, dbl2); break;
                        case mxLOGICAL_CLASS:
                        case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; WHIST1DNOOF( uc1, dbl2); break;
                        case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; WHIST1DNOOF( sc1, dbl2); break;
                        default:
                            mexErrMsgTxt("Invalid input class.");
                    }

                /* only offset needed */
                } else {

                    /* switch over classes */
                    switch (cid1) {
                        case mxDOUBLE_CLASS:                                       WHIST1DNOFD(dbl1,dbl2); break;
                        case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; WHIST1DNOFD(sng1,dbl2); break;
                        case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; WHIST1DNOF( ul1, dbl2); break;
                        case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; WHIST1DNOF( sl1, dbl2); break;
                        case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; WHIST1DNOF( us1, dbl2); break;
                        case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; WHIST1DNOF( ss1, dbl2); break;
                        case mxLOGICAL_CLASS:
                        case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; WHIST1DNOF( uc1, dbl2); break;
                        case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; WHIST1DNOF( sc1, dbl2); break;
                        default:
                            mexErrMsgTxt("Invalid input class.");
                    }
                }

            /* factor needed, but no offset! */
            } else {

                /* switch over classes */
                switch (cid1) {
                    case mxDOUBLE_CLASS:                                       WHIST1DNOOD(dbl1,dbl2); break;
                    case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; WHIST1DNOOD(sng1,dbl2); break;
                    case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; WHIST1DNOO( ul1, dbl2); break;
                    case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; WHIST1DNOO( sl1, dbl2); break;
                    case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; WHIST1DNOO( us1, dbl2); break;
                    case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; WHIST1DNOO( ss1, dbl2); break;
                    case mxLOGICAL_CLASS:
                    case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; WHIST1DNOO( uc1, dbl2); break;
                    case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; WHIST1DNOO( sc1, dbl2); break;
                    default:
                        mexErrMsgTxt("Invalid input class.");
                }
            }

            /* return already */
            return;
        }

        /* switch over classes */
        switch (cid1) {
            case mxDOUBLE_CLASS:                                       HIST2DSD(dbl1);break;
            case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HIST2DSD(sng1);break;
            case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HIST2DS( ul1); break;
            case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HIST2DS( sl1); break;
            case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HIST2DS( us1); break;
            case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HIST2DS( ss1); break;
            case mxLOGICAL_CLASS:
            case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HIST2DS( uc1); break;
            case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HIST2DS( sc1); break;
            default:
                mexErrMsgTxt("Invalid input class.");
        }

    /* one-dimensional case */
    } else {

        /* factor/offset needed */
        if ((mn1 != 0.0) &&
            (f1 != 1.0)) {

            /* switch over classes */
            switch (cid1) {
                case mxDOUBLE_CLASS:                                       HIST1DD(dbl1);break;
                case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HIST1DD(sng1);break;
                case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HIST1D( ul1); break;
                case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HIST1D( sl1); break;
                case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HIST1D( us1); break;
                case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HIST1D( ss1); break;
                case mxLOGICAL_CLASS:
                case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HIST1D( uc1); break;
                case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HIST1D( sc1); break;
                default:
                    mexErrMsgTxt("Invalid input class.");
            }

        /* no factor needed */
        } else if (f1 == 1.0) {

            /* no offset needed */
            if (mn1 == 0.0) {

                /* switch over classes */
                switch (cid1) {
                    case mxDOUBLE_CLASS:                                       HIST1DNOOFD(dbl1);break;
                    case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HIST1DNOOFD(sng1);break;
                    case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HIST1DNOOF( ul1); break;
                    case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HIST1DNOOF( sl1); break;
                    case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HIST1DNOOF( us1); break;
                    case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HIST1DNOOF( ss1); break;
                    case mxLOGICAL_CLASS:
                    case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HIST1DNOOF( uc1); break;
                    case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HIST1DNOOF( sc1); break;
                    default:
                        mexErrMsgTxt("Invalid input class.");
                }

            /* only offset needed */
            } else {

                /* switch over classes */
                switch (cid1) {
                    case mxDOUBLE_CLASS:                                       HIST1DNOFD(dbl1);break;
                    case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HIST1DNOFD(sng1);break;
                    case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HIST1DNOF( ul1); break;
                    case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HIST1DNOF( sl1); break;
                    case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HIST1DNOF( us1); break;
                    case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HIST1DNOF( ss1); break;
                    case mxLOGICAL_CLASS:
                    case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HIST1DNOF( uc1); break;
                    case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HIST1DNOF( sc1); break;
                    default:
                        mexErrMsgTxt("Invalid input class.");
                }
            }

        /* factor needed, but no offset! */
        } else {

            /* switch over classes */
            switch (cid1) {
                case mxDOUBLE_CLASS:                                       HIST1DNOOD(dbl1);break;
                case mxSINGLE_CLASS: sng1 = (const          float *) dbl1; HIST1DNOOD(sng1);break;
                case mxUINT32_CLASS:  ul1 = (const unsigned int   *) dbl1; HIST1DNOO( ul1); break;
                case  mxINT32_CLASS:  sl1 = (const   signed int   *) dbl1; HIST1DNOO( sl1); break;
                case mxUINT16_CLASS:  us1 = (const unsigned short *) dbl1; HIST1DNOO( us1); break;
                case  mxINT16_CLASS:  ss1 = (const   signed short *) dbl1; HIST1DNOO( ss1); break;
                case mxLOGICAL_CLASS:
                case  mxUINT8_CLASS:  uc1 = (const unsigned char  *) dbl1; HIST1DNOO( uc1); break;
                case   mxINT8_CLASS:  sc1 = (const   signed char  *) dbl1; HIST1DNOO( sc1); break;
                default:
                    mexErrMsgTxt("Invalid input class.");
            }
        }
    }
}

/*
  quattrans.h

% Version:  v0.9b
% Build:    11050510
% Date:     Jul-24 2010, 1:01 PM EST
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


/* quaternion transformation with simple variables */


/* variables */
#define QUATVARS                                    \
double                                              \
    qt11 = 0.0, qt12 = 0.0, qt13 = 0.0, qt14 = 0.0, \
    qt21 = 0.0, qt22 = 0.0, qt23 = 0.0, qt24 = 0.0, \
    qt31 = 0.0, qt32 = 0.0, qt33 = 0.0, qt34 = 0.0, \
    qtp1 = 0.0, qtp2 = 0.0;

#define QUATVARS2                                   \
double                                              \
    q211 = 0.0, q212 = 0.0, q213 = 0.0, q214 = 0.0, \
    q221 = 0.0, q222 = 0.0, q223 = 0.0, q224 = 0.0, \
    q231 = 0.0, q232 = 0.0, q233 = 0.0, q234 = 0.0;

/* get vars from 4x4 input */
#define GET_QUATVARS(PTR)   \
    qt11 = (double) *PTR++; \
    qt21 = (double) *PTR++; \
    qt31 = (double) *PTR++; \
    ++PTR;                  \
    qt12 = (double) *PTR++; \
    qt22 = (double) *PTR++; \
    qt32 = (double) *PTR++; \
    ++PTR;                  \
    qt13 = (double) *PTR++; \
    qt23 = (double) *PTR++; \
    qt33 = (double) *PTR++; \
    ++PTR;                  \
    qt14 = (double) *PTR++; \
    qt24 = (double) *PTR++; \
    qt34 = (double) *PTR++; \
    ++PTR;

/* get vars from 4x4 input */
#define GET_QUATVARS2(PTR)  \
    q211 = (double) *PTR++; \
    q221 = (double) *PTR++; \
    q231 = (double) *PTR++; \
    ++PTR;                  \
    q212 = (double) *PTR++; \
    q222 = (double) *PTR++; \
    q232 = (double) *PTR++; \
    ++PTR;                  \
    q213 = (double) *PTR++; \
    q223 = (double) *PTR++; \
    q233 = (double) *PTR++; \
    ++PTR;                  \
    q214 = (double) *PTR++; \
    q224 = (double) *PTR++; \
    q234 = (double) *PTR++; \
    ++PTR;

/* multiplication */
#define QUAT_MULT(C1,C2,C3)                         \
    qtp1 = qt11 * C1 + qt12 * C2 + qt13 * C3 + qt14;\
    qtp2 = qt21 * C1 + qt22 * C2 + qt23 * C3 + qt24;\
      C3 = qt31 * C1 + qt32 * C2 + qt33 * C3 + qt34;\
      C1 = qtp1;                                    \
      C2 = qtp2;
#define QUAT_MULT_0B(C1,C2,C3)                              \
    qtp1 = qt11 * C1 + qt12 * C2 + qt13 * C3 + qt14 - 1.0;  \
    qtp2 = qt21 * C1 + qt22 * C2 + qt23 * C3 + qt24 - 1.0;  \
      C3 = qt31 * C1 + qt32 * C2 + qt33 * C3 + qt34 - 1.0;  \
      C1 = qtp1;                                            \
      C2 = qtp2;
#define QUAT_MULT2(C1,C2,C3)                        \
    qtp1 = q211 * C1 + q212 * C2 + q213 * C3 + q214;\
    qtp2 = q221 * C1 + q222 * C2 + q223 * C3 + q224;\
      C3 = q231 * C1 + q232 * C2 + q233 * C3 + q234;\
      C1 = qtp1;                                    \
      C2 = qtp2;
#define QUAT_MULT2_0B(C1,C2,C3)                             \
    qtp1 = q211 * C1 + q212 * C2 + q213 * C3 + q214 - 1.0;  \
    qtp2 = q221 * C1 + q222 * C2 + q223 * C3 + q224 - 1.0;  \
      C3 = q231 * C1 + q232 * C2 + q233 * C3 + q234 - 1.0;  \
      C1 = qtp1;                                            \
      C2 = qtp2;

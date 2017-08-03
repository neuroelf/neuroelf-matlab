/* flexinterpn

FORMAT:       idata = flexinterpn(data, coords [, k, ks [, h [, t]]])

Input fields:

      data        up to 4D data
      coords      CxD (number of coords -by- number of dims) coordinates
      k           Kx1 interpolation kernel (see flexinterpn_method.m)
      ks          kernel sampling
      h           optional "hold" (default) value, default: 0.0
      t           if a coordinate range is given (see Note), apply 4x4
                  quaternion multiplication to range coordinates

Output fields:

      idata       interpolated data

Note: instead of a CxD array, a special 4xD array can be given. if the
      for "row" are Inf, the following three rows are treated as
      if they were used in Matlab's colon operator (row2:row3:row4)
      this allows for very fast subsampling

% Version:  v0.9d
% Build:    14062016
% Date:     Jun-20 2014, 4:21 PM EST
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


/* Matlab default includes */
#include "mex.h"
#include "math.h"

/* faster isinf/isnan rejection */
#include "isinfnan.h"

/* simple quaternion multiplication */
#include "quattrans.h"

/* different kinds of jobs */
enum flexinterpn_jobkind {
    coordslist = 0,
    coordrange,
    transrange
};

/* constants:
   - default kernel weights pointer (linear interpolator)
   - range safety margin (approx sqrt(eps) to ensure upper limit is hit)
   - number of default weights
   - default kernel stepsize */
const double
    defkernelptr[3] = {0.0, 1.0, 0.0},
    positivemargin = 0.0000000149,
    negativemargin = -0.0000000149;
const int
    defnrofweights = 3,
    defstepsize = 1;


/* now come some preprocessor macros used to simplify code writing */


/* - if the sum of weights is not zero (the kernel was certainly used)
   -     store the (correctly weighted) sampled value
   - else if the sample value still is non-zero (we yet have hit data)
   -     store the sampled value anyway (but we cannot re-weigh it !)
   - else we didn't hit anything
   - 	 put in the default
   - make sure samplev and sumv are in good shape for the next loop
                                                    >> NO SCOPE CHANGE */
#define ADD_TO_SAMPLED                                                  \
    if (sumweights != 0.0)                                              \
        *sampled++ = sampleval / sumweights;                            \
    else if (sampleval != 0.0)                                          \
        *sampled++ = sampleval;                                         \
    else                                                                \
        *sampled++ = defsample;                                         \
    sampleval = sumweights = 0.0;

/* - add weighting to sum (for which we can also use the inner dim!)
   - and add product with voxval to sample          >> NO SCOPE CHANGE */
#define ADD_TO_SUM_AND_SAMPLE                                           \
    sumweights += weightval1;                                           \
    sampleval += weightval1 * voxval;

/* allocating the useweight/scoord/weights cache for all dims on coords */
#define ALLOCATE_COORD_CACHE                                            \
    weightptr1 = (double *) mxCalloc(nrofdims * kernelelements1, sizeof(double)); \
    scoordptr1 = (int *) mxCalloc(nrofdims * kernelelements1, sizeof(int)); \
    useweight1 = (unsigned char *) mxCalloc(nrofdims * kernelelements1, sizeof(unsigned char)); \
    if ((weightptr1 == NULL) ||                                         \
        (scoordptr1 == NULL) ||                                         \
        (useweight1 == NULL)) {                                         \
        mxDestroyArray(*plhs);                                          \
        if (weightptr1 != NULL) mxFree(weightptr1);                     \
        if (scoordptr1 != NULL) mxFree(scoordptr1);                     \
        if (useweight1 != NULL) mxFree(useweight1);                     \
        mexErrMsgTxt("Error allocating weight array(s).");              \
    }                                                                   \
    SHIFT_CACHE_POINTERS(1,kernelsize1)                                \
    if (nrofdims > 1) {                                                 \
        GET_CACHEPTRS(2,1)                                              \
        if (nrofdims > 2) {                                             \
            GET_CACHEPTRS(3,2)                                          \
            if (nrofdims > 3) {                                         \
                GET_CACHEPTRS(4,3)                                      \
            }                                                           \
        }                                                               \
    }

/* freeing the coord-based cache */
#define ALLOCATE_COORD_CACHE_FREE                                       \
    SHIFT_CACHE_POINTERS(1,-kernelsize1)                                \
    mxFree(useweight1);                                                 \
    mxFree(scoordptr1);                                                 \
    mxFree(weightptr1);

/* allocating the useweight/scoord/weights cache for one dim */
#define ALLOCATE_RANGE_CACHE(DIM)                                       \
    weightptr##DIM = (double *) mxCalloc(                               \
        kernelelements##DIM * outputdim##DIM, sizeof(double));          \
    scoordptr##DIM = (int *) mxCalloc(                                  \
        kernelelements##DIM * outputdim##DIM, sizeof(int));             \
    useweight##DIM = (unsigned char *) mxCalloc(                        \
        kernelelements##DIM * outputdim##DIM, sizeof(unsigned char));   \
    ALLOCATE_RANGE_CACHE_CHECK(DIM)

#define ALLOCATE_RANGE_CACHE_FREE(DIM)                                  \
    ALLOCATE_RANGE_CACHE_FREE_SAFE(weightptr##DIM)                      \
    ALLOCATE_RANGE_CACHE_FREE_SAFE(scoordptr##DIM)                      \
    ALLOCATE_RANGE_CACHE_FREE_SAFE(useweight##DIM)

#define ALLOCATE_RANGE_CACHE_FREEALL                                    \
    ALLOCATE_RANGE_CACHE_FREE(1)                                        \
    ALLOCATE_RANGE_CACHE_FREE(2)                                        \
    ALLOCATE_RANGE_CACHE_FREE(3)                                        \
    ALLOCATE_RANGE_CACHE_FREE(4)

#define ALLOCATE_RANGE_CACHE_FREE_SAFE(PTR)                             \
    if (PTR != NULL)                                                    \
        mxFree(PTR);

#define ALLOCATE_RANGE_CACHE_CHECK(DIM)                                 \
    if ((weightptr##DIM == NULL) ||                                     \
        (scoordptr##DIM == NULL) ||                                     \
        (useweight##DIM == NULL)) {                                     \
        ALLOCATE_RANGE_CACHE_FREEALL                                    \
        mexErrMsgTxt("Error allocating memory for the weight caches."); \
    }

/* single loop to cache all kernel values (for coordinate access)
                               >> ADDING ONE SCOPE TO FOLLOWING CODE ! */
#define CACHE_KERNEL_OPENLOOP(DIM)                                      \
    for (                                                               \
        kernelcount1 = -kernelsize##DIM;                                \
        kernelcount1 <= kernelsize##DIM;                                \
        ++kernelcount1) {

/* caching samplecoords, etc. for one kernel dim */
#define CACHE_KERNEL_SPEC(DIM)                                          \
    samplecoord##DIM = basecoord##DIM - kernelcount1;                   \
    kernelsubcoord##DIM = kernoffset##DIM - kernelcount1 * stepsize##DIM; \
    if ((samplecoord##DIM < 0) || (samplecoord##DIM >= inputdim##DIM))  \
        useweight##DIM[kernelcount1] = 0;                               \
    else {                                                              \
        kernelsubcoord##DIM = kernoffset##DIM - kernelcount1 * stepsize##DIM; \
        if ((kernelsubcoord##DIM < kernellimneg##DIM) ||                \
            (kernelsubcoord##DIM >= kernellimpos##DIM))                 \
            useweight##DIM[kernelcount1] = 0;                           \
        else {                                                          \
            useweight##DIM[kernelcount1] = 1;                           \
            scoordptr##DIM[kernelcount1] =                              \
                samplecoord##DIM * incindim##DIM;                       \
            weightptr##DIM[kernelcount1] =                              \
                (1.0 + coordval##DIM) * kernelptr##DIM[kernelsubcoord##DIM] \
                - coordval##DIM * kernelptr##DIM[kernelsubcoord##DIM + 1]; \
        }                                                               \
    }

#define CACHE_KERNEL_SPEC_UPTO2                                         \
    COMPUTE_KERNEL_OFFSET(1)                                            \
    COMPUTE_KERNEL_OFFSET(2)                                            \
    CACHE_KERNEL_OPENLOOP(1)                                            \
        CACHE_KERNEL_SPEC(1)                                            \
    CLOSE_LOOP                                                          \
    CACHE_KERNEL_OPENLOOP(2)                                            \
        CACHE_KERNEL_SPEC(2)                                            \
    CLOSE_LOOP

#define CACHE_KERNEL_SPEC_UPTO3                                         \
    COMPUTE_KERNEL_OFFSET(1)                                            \
    COMPUTE_KERNEL_OFFSET(2)                                            \
    COMPUTE_KERNEL_OFFSET(3)                                            \
    CACHE_KERNEL_OPENLOOP(1)                                            \
        CACHE_KERNEL_SPEC(1)                                            \
    CLOSE_LOOP                                                          \
    CACHE_KERNEL_OPENLOOP(2)                                            \
        CACHE_KERNEL_SPEC(2)                                            \
    CLOSE_LOOP                                                          \
    CACHE_KERNEL_OPENLOOP(3)                                            \
        CACHE_KERNEL_SPEC(3)                                            \
    CLOSE_LOOP

#define CACHE_KERNEL_SPEC_UPTO4                                         \
    COMPUTE_KERNEL_OFFSET(1)                                            \
    COMPUTE_KERNEL_OFFSET(2)                                            \
    COMPUTE_KERNEL_OFFSET(3)                                            \
    COMPUTE_KERNEL_OFFSET(4)                                            \
    CACHE_KERNEL_OPENLOOP(1)                                            \
        CACHE_KERNEL_SPEC(1)                                            \
    CLOSE_LOOP                                                          \
    CACHE_KERNEL_OPENLOOP(2)                                            \
        CACHE_KERNEL_SPEC(2)                                            \
    CLOSE_LOOP                                                          \
    CACHE_KERNEL_OPENLOOP(3)                                            \
        CACHE_KERNEL_SPEC(3)                                            \
    CLOSE_LOOP                                                          \
    CACHE_KERNEL_OPENLOOP(4)                                            \
        CACHE_KERNEL_SPEC(4)                                            \
    CLOSE_LOOP

#define CHECK_COORDDIM_STEP1(DIM)                                       \
    IF_IS_BAD_COORD_CONT(coordval##DIM)
#define CHECK_COORDDIM_STEP1_UPTO1 CHECK_COORDDIM_STEP1(1)
#define CHECK_COORDDIM_STEP1_UPTO2 CHECK_COORDDIM_STEP1_UPTO1 CHECK_COORDDIM_STEP1(2)
#define CHECK_COORDDIM_STEP1_UPTO3 CHECK_COORDDIM_STEP1_UPTO2 CHECK_COORDDIM_STEP1(3)
#define CHECK_COORDDIM_STEP1_UPTO4 CHECK_COORDDIM_STEP1_UPTO3 CHECK_COORDDIM_STEP1(4)

/* checking the range validity of a coordinate
   - compute the base coordinate (nearest full voxel)
   - if the base is outside the volume of data
   -     put default sample in output
         and continue (with next coordinate)        >> NO SCOPE CHANGE */
#define CHECK_COORDDIM_STEP2(DIM)                                       \
    basecoord##DIM = (int) floor(coordval##DIM + 0.5);                  \
    if ((basecoord##DIM < 0) ||                                         \
        (basecoord##DIM >= inputdim##DIM)) {                            \
        *sampled++ = defsample;                                         \
        continue;                                                       \
    }
#define CHECK_COORDDIM_STEP2_RANGE(DIM)                                 \
    basecoord##DIM = (int) floor(coordval##DIM + 0.5);                  \
    if ((basecoord##DIM < 0) ||                                         \
        (basecoord##DIM >= inputdim##DIM)) {                            \
        for (                                                           \
            largecount = incoutdim##DIM;                                \
            largecount > 0;                                             \
            --largecount)                                               \
            *sampled++ = defsample;                                     \
        continue;                                                       \
    }

/* shortcuts to check dims 1 ... 4                  >> NO SCOPE CHANGE */
#define CHECK_COORDDIM_STEP2_UPTO1 CHECK_COORDDIM_STEP2(1)
#define CHECK_COORDDIM_STEP2_UPTO2 CHECK_COORDDIM_STEP2_UPTO1 CHECK_COORDDIM_STEP2(2)
#define CHECK_COORDDIM_STEP2_UPTO3 CHECK_COORDDIM_STEP2_UPTO2 CHECK_COORDDIM_STEP2(3)
#define CHECK_COORDDIM_STEP2_UPTO4 CHECK_COORDDIM_STEP2_UPTO3 CHECK_COORDDIM_STEP2(4)

/* check the range of a dimension to see how many elements it takes;
   the from and to value are reduced by 1 to take care of the Matlab
   indexing, and also the security margin is added so that the range
   truly hits the "to" value; if the test is successful the range is
   once "walked" to determine the number of values in that dimension

   - is the first element Inf or NaN?
   -     get the from coordinate (-1 for 0-based indices)
   -     get the stepsize for this dim
   -     get the "to" coordinate and add (or subtract) the error margin
   -     if stepsize is 0 or the range would be empty
   -         set coordrange to false
   -     otherwise
   -         set output dim to zero
   -         iterate from "from" to "to" by "step"
   -             and increase dim accordingly
   - so the first element was not Inf or Nan?
   -     then set coordrange to false !
                                                    >> NO SCOPE CHANGE */
#define CHECK_GET_RANGE_OF_DIM(DIM)                                     \
    IF_IS_BAD_VAL(*coordptr##DIM) {                                     \
        coordfrom##DIM = coordptr##DIM[1] - 1.0;                        \
        coordstep##DIM = coordptr##DIM[2];                              \
        coordto##DIM = coordptr##DIM[3] - 1.0 +                         \
            (SIGN_OF_VALUE(coordstep##DIM) * positivemargin);           \
        if (((coordstep##DIM > 0.0) &&                                  \
             (coordto##DIM > coordfrom##DIM)) ||                        \
            ((coordstep##DIM < 0.0) &&                                  \
             (coordto##DIM < coordfrom##DIM))) {                        \
            outputdim##DIM = 0;                                         \
            dblstepsize##DIM = coordfrom##DIM;                          \
            if (coordstep##DIM > 0.0)                                   \
                for (;                                                  \
                    dblstepsize##DIM <= coordto##DIM;                   \
                    dblstepsize##DIM += coordstep##DIM)                 \
                        ++outputdim##DIM;                               \
            else                                                        \
                for (;                                                  \
                    dblstepsize##DIM >= coordto##DIM;                   \
                    dblstepsize##DIM += coordstep##DIM)                 \
                        ++outputdim##DIM;                               \
        }                                                               \
    } else                                                              \
        jobkind = coordslist;

/* closing a LOOP                     >> REDUCING SCOPE COUNT BY 1 !!! */
#define CLOSE_LOOP }

/* computing the offset into the kernel for that dim
   - reduce value by the base coordinate
   - compute the offset into the kernel             >> NO SCOPE CHANGE */
#define COMPUTE_KERNEL_OFFSET(DIM)                                      \
    coordval##DIM -= (double) basecoord##DIM;                           \
    kernoffset##DIM = -(int) ceil(dblstepsize##DIM * coordval##DIM);    \
    coordval##DIM += ((double) kernoffset##DIM) * dblstepsizequot##DIM; \
    coordval##DIM *= dblstepsize##DIM;                                  \
    if (coordval##DIM > 0.0) {                                          \
        --coordval##DIM;                                                \
        --kernoffset##DIM;                                              \
    }

/* fills the cache for a range */
#define FILL_RANGE_CACHE(DIM)                                           \
    largecount = 0;                                                     \
    for (coordcount##DIM = 0,                                           \
        coordpos##DIM = coordfrom##DIM;                                 \
        coordcount##DIM < outputdim##DIM;                               \
        coordpos##DIM += coordstep##DIM, ++coordcount##DIM) {           \
        coordval##DIM = coordpos##DIM;                                  \
        basecoord##DIM = (int) floor(coordval##DIM + 0.5);              \
        if ((basecoord##DIM >= 0) &&                                    \
            (basecoord##DIM < inputdim##DIM)) {                         \
            coordval##DIM -= (double) basecoord##DIM;                   \
            kernoffset##DIM = -(int) ceil(dblstepsize##DIM * coordval##DIM); \
            coordval##DIM += ((double) kernoffset##DIM) * dblstepsizequot##DIM; \
            coordval##DIM *= dblstepsize##DIM;                          \
            if (coordval##DIM > 0.0) {                                  \
                --coordval##DIM;                                        \
                --kernoffset##DIM;                                      \
            }                                                           \
            for (                                                       \
                kernelcount##DIM = -kernelsize##DIM;                    \
                kernelcount##DIM <= kernelsize##DIM;                    \
                ++kernelcount##DIM) {                                   \
                samplecoord##DIM = basecoord##DIM - kernelcount##DIM;   \
                if ((samplecoord##DIM < 0) ||                           \
                    (samplecoord##DIM >= inputdim##DIM)) {              \
                    useweight##DIM[kernelcount##DIM + largecount] = 0;  \
                    continue;                                           \
                }                                                       \
                kernelsubcoord##DIM =                                   \
                    kernoffset##DIM - kernelcount##DIM * stepsize##DIM; \
                if ((kernelsubcoord##DIM < kernellimneg##DIM) ||        \
                    (kernelsubcoord##DIM >= kernellimpos##DIM)) {       \
                    useweight##DIM[kernelcount##DIM + largecount] = 0;  \
                    continue;                                           \
                }                                                       \
                weightval##DIM =                                        \
                    (1.0 + coordval##DIM) * kernelptr##DIM[kernelsubcoord##DIM] \
                    - coordval##DIM * kernelptr##DIM[kernelsubcoord##DIM + 1]; \
                if ((weightval##DIM >= positivemargin) ||               \
                    (weightval##DIM <= negativemargin)) {               \
                    useweight##DIM[kernelcount##DIM + largecount] = 1;  \
                    scoordptr##DIM[kernelcount##DIM + largecount] =     \
                        samplecoord##DIM * incindim##DIM;               \
                    weightptr##DIM[kernelcount##DIM + largecount] =     \
                        weightval##DIM;                                 \
                } else                                                  \
                    useweight##DIM[kernelcount##DIM + largecount] = 0;  \
            }                                                           \
        }                                                               \
        largecount += kernelelements##DIM;                              \
    }

/* get cache pointers from previous dim */
#define GET_CACHEPTRS(DIM,FDIM)                                         \
    weightptr##DIM = &weightptr##FDIM[kernelelements1];                 \
    scoordptr##DIM = &scoordptr##FDIM[kernelelements1];                 \
    useweight##DIM = &useweight##FDIM[kernelelements1];

/* get coordinate value from array and
   continue if inf/nan
   then check the coordinate for data range         >> NO SCOPE CHANGE */
#define GET_COORDDIM(DIM)                                               \
    coordval##DIM = *coordptr##DIM++ - 1.0;                             \

/* shortcuts for number of dimensions               >> NO SCOPE CHANGE */
#define GET_COORD_UPTODIMS_1 GET_COORDDIM(1)
#define GET_COORD_UPTODIMS_2 GET_COORD_UPTODIMS_1 GET_COORDDIM(2)
#define GET_COORD_UPTODIMS_3 GET_COORD_UPTODIMS_2 GET_COORDDIM(3)
#define GET_COORD_UPTODIMS_4 GET_COORD_UPTODIMS_3 GET_COORDDIM(4)

/* get coord pointer and increments for dim */
#define GET_COORDPTR_AND_INC(DIM,FDIM)                                  \
    coordptr##DIM = &coordptr##FDIM [nrofcoords];                       \
    incindim##DIM = incindim##FDIM * inputdim##FDIM;                    \
    inputdim##DIM = inputdim[FDIM];

/* compute the weights for two adjacent kernel elements as
   the linear interpolation between the two */
#define GET_WEIGHT_ON_KDIM(DIM)                                         \
    weightval##DIM =                                                    \
        (1.0 + coordval##DIM) * kernelptr##DIM[kernelsubcoord##DIM]     \
        - coordval##DIM * kernelptr##DIM[kernelsubcoord##DIM + 1];

/* put default value into output and continue for bad coordinates */
#define IF_IS_BAD_COORD_CONT(VAL)                                       \
    IF_IS_BAD_VAL(VAL) {                                                \
        *sampled++ = defsample;                                         \
        continue;                                                       \
    }

/* (only) continue in the kernel sampling current dim ! */
#define IF_IS_BAD_VAL_CONT(VAL)                                         \
    IF_IS_BAD_VAL(VAL)                                                  \
        continue;

/* kernel sampling inner loop when at least two dims are given */
#define INNER_LOOP_STATEMENTS(DATAPTR)                                  \
    OPEN_LOOP_KDIM_LOWER(1,2)                                           \
        voxval = (double) DATAPTR[samplecoord1];                        \
        ADD_TO_SUM_AND_SAMPLE                                           \
    CLOSE_LOOP

/* kernel sampling inner loop for multiple dims, check floating points */
#define INNER_LOOP_STATEMENTS_FP(DATAPTR)                               \
    OPEN_LOOP_KDIM_LOWER(1,2)                                           \
        voxval = (double) DATAPTR[samplecoord1];                        \
        IF_IS_BAD_VAL_CONT(voxval)                                      \
        ADD_TO_SUM_AND_SAMPLE                                           \
    CLOSE_LOOP

/* kernel sampling inner loop when at least two dims are given */
#define INNER_LOOP_STATEMENTS_RANGE(DATAPTR)                            \
    OPEN_LOOP_KDIM_RANGE_LOWER(1,2)                                     \
        voxval = (double) DATAPTR[samplecoord1];                        \
        ADD_TO_SUM_AND_SAMPLE                                           \
    CLOSE_LOOP

/* kernel sampling inner loop when at least two dims are given, FP check */
#define INNER_LOOP_STATEMENTS_RANGE_FP(DATAPTR)                         \
    OPEN_LOOP_KDIM_RANGE_LOWER(1,2)                                     \
        voxval = (double) DATAPTR[samplecoord1];                        \
        IF_IS_BAD_VAL_CONT(voxval)                                      \
        ADD_TO_SUM_AND_SAMPLE                                           \
    CLOSE_LOOP

/* kernel sampling inner loop when only one dim is given */
#define INNER_LOOP_STATEMENTS_SINGLEDIM(DATAPTR)                        \
    OPEN_LOOP_KDIM_SINGLE                                               \
        voxval = (double) DATAPTR[samplecoord1];                        \
        GET_WEIGHT_ON_KDIM(1)                                           \
        ADD_TO_SUM_AND_SAMPLE                                           \
    CLOSE_LOOP

/* kernel sampling inner loop when only one dim is given, FP check */
#define INNER_LOOP_STATEMENTS_SINGLEDIM_FP(DATAPTR)                     \
    OPEN_LOOP_KDIM_SINGLE                                               \
        voxval = (double) DATAPTR[samplecoord1];                        \
        IF_IS_BAD_VAL_CONT(voxval)                                      \
        GET_WEIGHT_ON_KDIM(1)                                           \
        ADD_TO_SUM_AND_SAMPLE                                           \
    CLOSE_LOOP

/* open loop over coordinates */
#define OPEN_LOOP_COORDS                                                \
    for (                                                               \
        coordcount = nrofcoords;                                        \
        coordcount > 0;                                                 \
        --coordcount) {                                                 \

/* open loop to iterate over one dim of the kernel */
#define OPEN_LOOP_KDIM(DIM)                                             \
    for (                                                               \
        kernelcount##DIM = -kernelsize##DIM;                            \
        kernelcount##DIM <= kernelsize##DIM;                            \
        ++kernelcount##DIM) {                                           \
            if (useweight##DIM[kernelcount##DIM] == 0)                  \
                continue;                                               \
            samplecoord##DIM = scoordptr##DIM[kernelcount##DIM];        \
            weightval##DIM = weightptr##DIM[kernelcount##DIM];          \

/* open loop to iterate over one dim of the kernel */
#define OPEN_LOOP_KDIM_LOWER(DIM,UDIM)                                  \
    for (                                                               \
        kernelcount##DIM = -kernelsize##DIM;                            \
        kernelcount##DIM <= kernelsize##DIM;                            \
        ++kernelcount##DIM) {                                           \
            if (useweight##DIM[kernelcount##DIM] == 0)                  \
                continue;                                               \
            samplecoord##DIM =                                          \
                samplecoord##UDIM + scoordptr##DIM[kernelcount##DIM];   \
            weightval##DIM =                                            \
                weightval##UDIM * weightptr##DIM[kernelcount##DIM];

/* open loop to iterate over one dim of the kernel (cached weights) */
#define OPEN_LOOP_KDIM_RANGE(DIM)                                       \
    for (                                                               \
        kernelcount##DIM = -kernelsize##DIM;                            \
        kernelcount##DIM <= kernelsize##DIM;                            \
        ++kernelcount##DIM) {                                           \
            if (useweight##DIM[cachepos##DIM + kernelcount##DIM] == 0)  \
                continue;                                               \
            samplecoord##DIM =                                          \
                scoordptr##DIM[cachepos##DIM + kernelcount##DIM];       \
            weightval##DIM =                                            \
                weightptr##DIM[cachepos##DIM + kernelcount##DIM];

/* open loop to iterate over one dim of the kernel (cached weights) */
#define OPEN_LOOP_KDIM_RANGE_LOWER(DIM,UDIM)                            \
    for (                                                               \
        kernelcount##DIM = -kernelsize##DIM;                            \
        kernelcount##DIM <= kernelsize##DIM;                            \
        ++kernelcount##DIM) {                                           \
            if (useweight##DIM[cachepos##DIM + kernelcount##DIM] == 0)  \
                continue;                                               \
            samplecoord##DIM = samplecoord##UDIM +                      \
                scoordptr##DIM[cachepos##DIM + kernelcount##DIM];       \
            weightval##DIM = weightval##UDIM *                          \
                weightptr##DIM[cachepos##DIM + kernelcount##DIM];

/* open loop to iterate over single dim of the kernel */
#define OPEN_LOOP_KDIM_SINGLE                                           \
    for (                                                               \
        kernelcount1 = -kernelsize1;                                    \
        kernelcount1 <= kernelsize1;                                    \
        ++kernelcount1) {                                               \
            samplecoord1 =                                              \
                basecoord1 - kernelcount1;                              \
            if ((samplecoord1 < 0) ||                                   \
                (samplecoord1 >= inputdim1))                            \
                continue;                                               \
            kernelsubcoord1 =                                           \
                kernoffset1 - kernelcount1 * stepsize1;                 \
            if ((kernelsubcoord1 < kernellimneg1) ||                    \
                (kernelsubcoord1 >= kernellimpos1))                     \
                continue;

/* open a loop for one of the range dims */
#define OPEN_RANGE_LOOP_DIM(DIM)                                        \
    for (                                                               \
        cachepos##DIM = -kernelelements##DIM,                           \
            coordpos##DIM = coordfrom##DIM,                             \
            coordcount##DIM = outputdim##DIM;                           \
        coordcount##DIM > 0;                                            \
        --coordcount##DIM) {                                            \
        cachepos##DIM += kernelelements##DIM;                           \
        coordval##DIM = coordpos##DIM;                                  \
        coordpos##DIM += coordstep##DIM;                                \
        CHECK_COORDDIM_STEP2_RANGE(DIM)

#define OPEN_RANGE_LOOP_DIMS1 OPEN_RANGE_LOOP_SINGLEDIM
#define OPEN_RANGE_LOOP_DIMS2 OPEN_RANGE_LOOP_DIM(2) OPEN_RANGE_LOOP_DIM(1)
#define OPEN_RANGE_LOOP_DIMS3 OPEN_RANGE_LOOP_DIM(3) OPEN_RANGE_LOOP_DIMS2
#define OPEN_RANGE_LOOP_DIMS4 OPEN_RANGE_LOOP_DIM(4) OPEN_RANGE_LOOP_DIMS3

/* open a loop for one of the range dims */
#define OPEN_RANGE_LOOP_SINGLEDIM                                       \
    for (                                                               \
        coordpos1 = coordfrom1,                                         \
            coordcount1 = outputdim1;                                   \
        coordcount1 > 0;                                                \
        --coordcount1) {                                                 \
        coordval1 = coordpos1;                                          \
        coordpos1 += coordstep1;                                        \
        CHECK_COORDDIM_STEP2_RANGE(1)

/* open a loop for one of the range dims */
#define OPEN_TRANGE_LOOP_DIM(DIM)                                       \
    for (                                                               \
        coordpos##DIM = coordfrom##DIM,                                 \
            coordcount##DIM = outputdim##DIM;                           \
        coordcount##DIM > 0;                                            \
        --coordcount##DIM) {                                            \
        coordpos##DIM += coordstep##DIM;

/* open a loop for one of the range dims for caching */
#define PREPARE_COORDS_AND_KERNEL_FOR_DIM2                              \
    GET_COORD_UPTODIMS_2                                                \
    CHECK_COORDDIM_STEP1_UPTO2                                          \
    CHECK_COORDDIM_STEP2_UPTO2                                          \
    CACHE_KERNEL_SPEC_UPTO2
#define PREPARE_COORDS_AND_KERNEL_FOR_DIM3                              \
    GET_COORD_UPTODIMS_3                                                \
    CHECK_COORDDIM_STEP1_UPTO3                                          \
    CHECK_COORDDIM_STEP2_UPTO3                                          \
    CACHE_KERNEL_SPEC_UPTO3
#define PREPARE_COORDS_AND_KERNEL_FOR_DIM4                              \
    GET_COORD_UPTODIMS_4                                                \
    CHECK_COORDDIM_STEP1_UPTO4                                          \
    CHECK_COORDDIM_STEP2_UPTO4                                          \
    CACHE_KERNEL_SPEC_UPTO4

/* sign of value */
#define SIGN_OF_VALUE(VAL) ((VAL < 0.0) ? -1.0 : 1.0)

/* shift all cache pointers */
#define SHIFT_CACHE_POINTERS(DIM,POS)                                   \
    SHIFT_POINTER(weightptr##DIM,POS)                                   \
    SHIFT_POINTER(scoordptr##DIM,POS)                                   \
    SHIFT_POINTER(useweight##DIM,POS)

/* shift a pointer */
#define SHIFT_POINTER(PTR,SHIFT) PTR = &PTR[SHIFT];

#define SWITCH_COORDS_OVER_DIMS(DATAPTR,DATATYPE)                       \
    DATAPTR = (const DATATYPE *) inputvoid;                             \
    switch (nrofdims) {                                                 \
        case 1:                                                         \
            OPEN_LOOP_COORDS                                            \
                GET_COORD_UPTODIMS_1                                    \
                CHECK_COORDDIM_STEP1_UPTO1                              \
                CHECK_COORDDIM_STEP2_UPTO1                              \
                COMPUTE_KERNEL_OFFSET(1)                                \
                INNER_LOOP_STATEMENTS_SINGLEDIM(DATAPTR)                \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 2:                                                         \
            OPEN_LOOP_COORDS                                            \
                PREPARE_COORDS_AND_KERNEL_FOR_DIM2                      \
                OPEN_LOOP_KDIM(2)                                       \
                    INNER_LOOP_STATEMENTS(DATAPTR)                      \
                CLOSE_LOOP                                              \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 3:                                                         \
            OPEN_LOOP_COORDS                                            \
                PREPARE_COORDS_AND_KERNEL_FOR_DIM3                      \
                OPEN_LOOP_KDIM(3)                                       \
                    OPEN_LOOP_KDIM_LOWER(2,3)                           \
                        INNER_LOOP_STATEMENTS(DATAPTR)                  \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 4:                                                         \
            OPEN_LOOP_COORDS                                            \
                PREPARE_COORDS_AND_KERNEL_FOR_DIM4                      \
                OPEN_LOOP_KDIM(4)                                       \
                    OPEN_LOOP_KDIM_LOWER(3,4)                           \
                        OPEN_LOOP_KDIM_LOWER(2,3)                       \
                            INNER_LOOP_STATEMENTS(DATAPTR)              \
                        CLOSE_LOOP                                      \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
    }

#define SWITCH_COORDS_OVER_DIMS_FP(DATAPTR, DATATYPE)                   \
    DATAPTR = (const DATATYPE *) inputvoid;                             \
    switch (nrofdims) {                                                 \
        case 1:                                                         \
            OPEN_LOOP_COORDS                                            \
                GET_COORD_UPTODIMS_1                                    \
                CHECK_COORDDIM_STEP1_UPTO1                              \
                CHECK_COORDDIM_STEP2_UPTO1                              \
                COMPUTE_KERNEL_OFFSET(1)                                \
                INNER_LOOP_STATEMENTS_SINGLEDIM_FP(DATAPTR)             \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 2:                                                         \
            OPEN_LOOP_COORDS                                            \
                PREPARE_COORDS_AND_KERNEL_FOR_DIM2                      \
                OPEN_LOOP_KDIM(2)                                       \
                    INNER_LOOP_STATEMENTS_FP(DATAPTR)                   \
                CLOSE_LOOP                                              \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 3:                                                         \
            OPEN_LOOP_COORDS                                            \
                PREPARE_COORDS_AND_KERNEL_FOR_DIM3                      \
                OPEN_LOOP_KDIM(3)                                       \
                    OPEN_LOOP_KDIM_LOWER(2,3)                           \
                        INNER_LOOP_STATEMENTS_FP(DATAPTR)               \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 4:                                                         \
            OPEN_LOOP_COORDS                                            \
                PREPARE_COORDS_AND_KERNEL_FOR_DIM4                      \
                OPEN_LOOP_KDIM(4)                                       \
                    OPEN_LOOP_KDIM_LOWER(3,4)                           \
                        OPEN_LOOP_KDIM_LOWER(2,3)                       \
                            INNER_LOOP_STATEMENTS_FP(DATAPTR)           \
                        CLOSE_LOOP                                      \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
    }

/* switch over nrofdims for one datatype, coordrange */
#define SWITCH_CRANGE_OVER_DIMS(DATAPTR,DATATYPE)                       \
    DATAPTR = (const DATATYPE *) inputvoid;                             \
    switch (nrofdims) {                                                 \
        case 1:                                                         \
            OPEN_RANGE_LOOP_DIMS1                                       \
                COMPUTE_KERNEL_OFFSET(1)                                \
                INNER_LOOP_STATEMENTS_SINGLEDIM(DATAPTR)                \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 2:                                                         \
            OPEN_RANGE_LOOP_DIMS2                                       \
                    OPEN_LOOP_KDIM_RANGE(2)                             \
                        INNER_LOOP_STATEMENTS_RANGE(DATAPTR)            \
                    CLOSE_LOOP                                          \
                    ADD_TO_SAMPLED                                      \
                CLOSE_LOOP                                              \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 3:                                                         \
            OPEN_RANGE_LOOP_DIMS3                                       \
                        OPEN_LOOP_KDIM_RANGE(3)                         \
                            OPEN_LOOP_KDIM_RANGE_LOWER(2,3)             \
                                INNER_LOOP_STATEMENTS_RANGE(DATAPTR)    \
                            CLOSE_LOOP                                  \
                        CLOSE_LOOP                                      \
                        ADD_TO_SAMPLED                                  \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 4:                                                         \
            OPEN_RANGE_LOOP_DIMS4                                       \
                            OPEN_LOOP_KDIM_RANGE(4)                     \
                                OPEN_LOOP_KDIM_RANGE_LOWER(3,4)         \
                                    OPEN_LOOP_KDIM_RANGE_LOWER(2,3)     \
                                        INNER_LOOP_STATEMENTS_RANGE(DATAPTR) \
                                    CLOSE_LOOP                          \
                                CLOSE_LOOP                              \
                            CLOSE_LOOP                                  \
                            ADD_TO_SAMPLED                              \
                        CLOSE_LOOP                                      \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
            CLOSE_LOOP                                                  \
            break;                                                      \
    }

/* switch over nrofdims for one datatype, coordrange, FP check */
#define SWITCH_CRANGE_OVER_DIMS_FP(DATAPTR,DATATYPE)                    \
    DATAPTR = (const DATATYPE *) inputvoid;                             \
    switch (nrofdims) {                                                 \
        case 1:                                                         \
            OPEN_RANGE_LOOP_DIMS1                                       \
                COMPUTE_KERNEL_OFFSET(1)                                \
                INNER_LOOP_STATEMENTS_SINGLEDIM_FP(DATAPTR)             \
                ADD_TO_SAMPLED                                          \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 2:                                                         \
            OPEN_RANGE_LOOP_DIMS2                                       \
                    OPEN_LOOP_KDIM_RANGE(2)                             \
                        INNER_LOOP_STATEMENTS_RANGE_FP(DATAPTR)         \
                    CLOSE_LOOP                                          \
                    ADD_TO_SAMPLED                                      \
                CLOSE_LOOP                                              \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 3:                                                         \
            OPEN_RANGE_LOOP_DIMS3                                       \
                        OPEN_LOOP_KDIM_RANGE(3)                         \
                            OPEN_LOOP_KDIM_RANGE_LOWER(2,3)             \
                                INNER_LOOP_STATEMENTS_RANGE_FP(DATAPTR) \
                            CLOSE_LOOP                                  \
                        CLOSE_LOOP                                      \
                        ADD_TO_SAMPLED                                  \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
            CLOSE_LOOP                                                  \
            break;                                                      \
        case 4:                                                         \
            OPEN_RANGE_LOOP_DIMS4                                       \
                            OPEN_LOOP_KDIM_RANGE(4)                     \
                                OPEN_LOOP_KDIM_RANGE_LOWER(3,4)         \
                                    OPEN_LOOP_KDIM_RANGE_LOWER(2,3)     \
                                        INNER_LOOP_STATEMENTS_RANGE_FP(DATAPTR) \
                                    CLOSE_LOOP                          \
                                CLOSE_LOOP                              \
                            CLOSE_LOOP                                  \
                            ADD_TO_SAMPLED                              \
                        CLOSE_LOOP                                      \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
            CLOSE_LOOP                                                  \
            break;                                                      \
    }

/* a range of coordinates to be transformed, nrofdims == 3 ! */
#define SWITCH_TRANGE_OVER_DIMS(DATAPTR,DATATYPE)                       \
    DATAPTR = (const DATATYPE *) inputvoid;                             \
    coordfrom1 += (1.0 - coordstep1);                                   \
    coordfrom2 += (1.0 - coordstep2);                                   \
    coordfrom3 += (1.0 - coordstep3);                                   \
    OPEN_TRANGE_LOOP_DIM(3)                                             \
       OPEN_TRANGE_LOOP_DIM(2)                                          \
            OPEN_TRANGE_LOOP_DIM(1)                                     \
                coordval3 = coordpos3;                                  \
                coordval2 = coordpos2;                                  \
                coordval1 = coordpos1;                                  \
                QUAT_MULT_0B(coordval1, coordval2, coordval3)           \
                CHECK_COORDDIM_STEP2_UPTO3                              \
                CACHE_KERNEL_SPEC_UPTO3                                 \
                OPEN_LOOP_KDIM(3)                                       \
                    OPEN_LOOP_KDIM_LOWER(2,3)                           \
                        INNER_LOOP_STATEMENTS_FP(DATAPTR)               \
                    CLOSE_LOOP                                          \
                CLOSE_LOOP                                              \
            ADD_TO_SAMPLED                                              \
            CLOSE_LOOP                                                  \
        CLOSE_LOOP                                                      \
    CLOSE_LOOP

/* make name1, name2, name3, name4 (for variable lists) */
#define VARNAME_MAKE4(NAME)                                             \
    NAME##1 ,                                                           \
    NAME##2 ,                                                           \
    NAME##3 ,                                                           \
    NAME##4

/* make name1, name2, name3, name4 (for variable lists) */
#define VARNAME_MAKE4I(NAME,VAL)                                        \
    NAME##1 = VAL ,                                                     \
    NAME##2 = VAL ,                                                     \
    NAME##3 = VAL ,                                                     \
    NAME##4 = VAL

/* make *name1, *name2, *name3, *name4 (for pointer lists) */
#define VARNAME_MAKE4P(NAME)                                            \
    * NAME##1 = NULL ,                                                  \
    * NAME##2 = NULL ,                                                  \
    * NAME##3 = NULL ,                                                  \
    * NAME##4 = NULL



/* main function (mexFunction)
   - Number of Left Hand Side,
   - Pointer to Left Hand Side,
   - Number of Right Hand Side,
   - Point to Right Hand Side (const) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* various variable definitions */

    /* debugging ? */
    /* char vstr[256]; */

    /* - flag: coordinates list (CxD) or range (4xD with 1. row Inf/NaN) */
    enum flexinterpn_jobkind jobkind = coordslist;

    /* - data dim pointer,
       - coordinate dim pointer,
       - weight dim pointer,
       - output dim array/pointer and vars */
    const int
        *inputdim,
        *coorddim,
        *weightdim;

    /* - pointers to coordinates
       - pointer to weights
       - pointer to the input 1x1 stepsize double
       - default sample value pointer
       - quaternion pointer */
    const double
        VARNAME_MAKE4P(coordptr) ,
        VARNAME_MAKE4P(kernelptr) ,
        *stepsizeptr,
        *quotptr;

    /* - input data pointers */
    const           void    *inputvoid;
    const           double  *inputdbl;
    const           float   *inputflt;
    const unsigned  int     *inputul;
    const   signed  int     *inputsl;
    const unsigned  short   *inputus;
    const   signed  short   *inputss;
    const unsigned  char    *inputuc;
    const   signed  char    *inputsc;

    /* - cell array access */
    const mxArray *ckernelptr, *ckernelsize;

    /* cheaper isinf / isnan test helpers */
    VARS_FOR_ISINFNAN

    /* - coordinate counter,
       - number of data dims,
       - number of coordinates,
       - number of elements,
       - number of weights in kernel
       - input class ID
       - size of output array
       - base (nearest voxel) and sample (plus kernel steps) coordinates
       - position within one dims' cache
       - counter for the output range (outputdim)
       - step size in int
       - input dims (for fast comparison with base coords)
       - dimension increments for input and output
       - a "large" counter
       - negative and positive kernel "limits" (so as not to access beyond)
       - kernel offsets and subcoordinates */
    int
        coordcount,
        nrofdims,
        nrofcoords,
        nrofelems = 1,
        VARNAME_MAKE4I(nrofweights, 1) ,
        inputclass = -1,
        inputdimx[4] = {1, 1, 1, 1},
    	outputdim[4] = {1, 1, 1, 1},
        VARNAME_MAKE4(basecoord) ,
        VARNAME_MAKE4(samplecoord) ,
        VARNAME_MAKE4(cachepos) ,
        VARNAME_MAKE4(coordcount) ,
        VARNAME_MAKE4P(scoordptr) ,
        VARNAME_MAKE4I(stepsize, 1) ,
        VARNAME_MAKE4I(inputdim, 1) ,
        VARNAME_MAKE4I(incindim, 1) ,
        VARNAME_MAKE4I(outputdim, 1) ,
        VARNAME_MAKE4I(incoutdim, 1) ,
        largecount,
        VARNAME_MAKE4I(kernellimpos, 0) ,
        VARNAME_MAKE4I(kernellimneg, 0) ,
        VARNAME_MAKE4I(kernelelements, 1) ,
        VARNAME_MAKE4(kernoffset) ,
        VARNAME_MAKE4(kernelsubcoord) ;

    /* - double for value from input data (casted)
       - sample value (must be 0.0 to begin with!)
       - interpolation weight sum (must be 0.0 to begin with!)
       - sampled data pointer
       - default value
       - sample coordinate values and weight values
       - pointer to shortlist of weights for given coordinate
       - range (4 dims, from, to, step and pos for counter, 0-based)
       - kernel stepsize as double
       - and 1 / stepsize, for fast computation
       - testvalue for isinfnan cheap method */
    double
        voxval,
        sampleval = 0.0,
        sumweights = 0.0,
    	*sampled,
        defsample = 0.0,
        VARNAME_MAKE4(coordval) ,
        VARNAME_MAKE4(weightval) ,
        VARNAME_MAKE4P(weightptr) ,
        VARNAME_MAKE4I(coordfrom, 0.0) ,
        VARNAME_MAKE4I(coordto, 1.0) ,
        VARNAME_MAKE4I(coordstep, 1.0) ,
        VARNAME_MAKE4I(coordpos, 0.0) ,
        VARNAME_MAKE4(dblstepsize) ,
        VARNAME_MAKE4(dblstepsizequot) ;

    /* quaternion rotation variables */
    QUATVARS

    /* - shortlist flags, whether to use the weight or not */
    unsigned char
        VARNAME_MAKE4P(useweight) ;

    /* small int counter, kernel size (in coords), kernel index (voxels!) */
    signed char
        smallcount1 = 0,
        smallcount2 = 0;
    signed short
        VARNAME_MAKE4I(kernelsize, 1) ,
        VARNAME_MAKE4(kernelcount) ;

    /* check input */
    if (nrhs < 2 || nrhs > 6 || nlhs > 1)
        mexErrMsgTxt("Bad number of input/output arguments.");
    if (!mxIsNumeric(*prhs))
        mexErrMsgTxt("Input data must be of numeric type.");
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("Sample coordinates must be double.");
    if (mxGetNumberOfDimensions(prhs[1]) > 2)
        mexErrMsgTxt("Sample coordinates must be CxD.");

    nrofdims = mxGetNumberOfDimensions(*prhs);
    if (nrofdims > 4)
        mexErrMsgTxt("Only up to 4 dims supported for now.");
    inputdim = (const int *) mxGetDimensions(*prhs);
    if ((nrofdims == 2) && (inputdim[1] == 1))
        --nrofdims;
    for (smallcount1 = 0; smallcount1 < nrofdims; ++smallcount1) {
        nrofelems *= inputdim[smallcount1];
        inputdimx[smallcount1] = inputdim[smallcount1];
    }
    inputdim = inputdimx;
    if (nrofelems == 0)
        mexErrMsgTxt("Input must not be empty");
    inputclass = mxGetClassID(*prhs);
    inputvoid = (const void *) mxGetData(*prhs);

    /* check coordinates */
    coorddim = (const int *) mxGetDimensions(prhs[1]);
    if ((coorddim[1] < nrofdims) ||
        (coorddim[1] > 4))
        mexErrMsgTxt("Sample coordinates must have at least NrOfDataDims columns.");
    nrofdims = coorddim[1];
    nrofcoords = *coorddim;

    INIT_INF_NAN_BAD_VAL()

    /* get pointers to sample coordinates and compute increments */
    coordptr1 = (const double *) mxGetData(prhs[1]);
    incindim1 = 1;
    inputdim1 = *inputdim;
    if (nrofdims > 1) {
        GET_COORDPTR_AND_INC(2,1)
        if (nrofdims > 2) {
            GET_COORDPTR_AND_INC(3,2)
            if (nrofdims > 3) {
                GET_COORDPTR_AND_INC(4,3)
            }
        }
    }

    /* allow special case: range selection */
    if (nrofcoords == 4) {
        for (smallcount2 = 0; smallcount2 < nrofdims; ++smallcount2) {
            for (smallcount1 = 1; smallcount1 < 4; ++smallcount1) {
                IF_IS_BAD_VAL(coordptr1[smallcount1 + 4 * smallcount2])
                    break;
            }
            if (smallcount1 != 4)
                break;
        }
        if (smallcount2 == nrofdims) {
            jobkind = coordrange;
            incoutdim1 = 1;
            outputdim4 = outputdim3 = outputdim2 = outputdim1 = 1;
            CHECK_GET_RANGE_OF_DIM(1)
            if ((jobkind == coordrange) &&
                (nrofdims > 1)) {
                incoutdim2 = outputdim1;
                CHECK_GET_RANGE_OF_DIM(2)
            }
            if ((jobkind == coordrange) &&
                (nrofdims > 2)) {
                incoutdim3 = incoutdim2 * outputdim2;
                CHECK_GET_RANGE_OF_DIM(3)
            }
            if ((jobkind == coordrange) &&
                (nrofdims > 3)) {
                incoutdim4 = incoutdim3 * outputdim3;
                CHECK_GET_RANGE_OF_DIM(4)
            }
        }
    }

    /* check weight input */
    nrofweights1 = defnrofweights;
    kernelptr1 = defkernelptr;
    if ((nrhs > 3) &&
        mxIsDouble(prhs[2]) &&
        mxIsDouble(prhs[3]) &&
        (mxGetNumberOfElements(prhs[2]) > 2) &&
        (mxGetNumberOfElements(prhs[3]) == 1)) {
        stepsizeptr = mxGetData(prhs[3]);
        nrofweights1 = mxGetNumberOfElements(prhs[2]);
        weightdim = (const int *) mxGetDimensions(prhs[2]);
        if (!mxIsNaN(*stepsizeptr) &&
            !mxIsInf(*stepsizeptr) &&
            (*stepsizeptr >= 1) &&
            ((nrofweights1 == *weightdim) ||
             (nrofweights1 == weightdim[1]))) {
            stepsize1 = (int) *stepsizeptr;
            kernelptr1 = mxGetData(prhs[2]);
            for (largecount = nrofweights1 - 1; largecount >= 0; --largecount) {
                IF_IS_BAD_VAL(*kernelptr1++)
                    break;
            }
            if (largecount >= 0) {
                nrofweights1 = defnrofweights;
                stepsize1    = defstepsize;
                kernelptr1   = defkernelptr;
            } else
                kernelptr1 = mxGetData(prhs[2]);
        } else {
            nrofweights1 = defnrofweights;
            stepsize1    = defstepsize;
            kernelptr1   = defkernelptr;
        }
        nrofweights4 = nrofweights3 = nrofweights2 = nrofweights1;
        stepsize4 = stepsize3 = stepsize2 = stepsize1;
        kernelptr4 = kernelptr3 = kernelptr2 = kernelptr1;
    } else if ((jobkind == coordrange) &&
       (nrhs > 3) &&
        mxIsCell(prhs[2]) &&
        mxIsCell(prhs[3]) &&
       (mxGetNumberOfElements(prhs[2]) == nrofdims) &&
       (mxGetNumberOfElements(prhs[3]) == nrofdims)) {
        for (smallcount1 = 0; smallcount1 < nrofdims; ++smallcount1) {
            ckernelptr = (const mxArray*) mxGetCell(prhs[2], smallcount1);
            ckernelsize = (const mxArray*) mxGetCell(prhs[3], smallcount1);
            if (mxIsDouble(ckernelptr) &&
                mxIsDouble(ckernelsize) &&
                (mxGetNumberOfElements(ckernelptr) > 2) &&
                (mxGetNumberOfElements(ckernelsize) == 1)) {
                stepsizeptr = mxGetData(ckernelsize);
                nrofweights4 = mxGetNumberOfElements(ckernelptr);
                weightdim = (const int *) mxGetDimensions(ckernelptr);
                if (!mxIsNaN(*stepsizeptr) &&
                    !mxIsInf(*stepsizeptr) &&
                    (*stepsizeptr >= 1) &&
                    ((nrofweights4 == *weightdim) ||
                     (nrofweights4 == weightdim[1]))) {
                    stepsize4 = (int) *stepsizeptr;
                    kernelptr4 = mxGetData(ckernelptr);
                    for (largecount = nrofweights4 - 1; largecount >= 0; --largecount) {
                        IF_IS_BAD_VAL(*kernelptr4++)
                            break;
                    }
                    if (largecount >= 0) {
                        break;
                    } else
                        kernelptr4 = mxGetData(ckernelptr);
                } else
                    break;
            }
            if (smallcount1 == 0) {
                nrofweights1 = nrofweights4;
                stepsize1    = stepsize4;
                kernelptr1   = kernelptr4;
            } else if (smallcount1 == 1) {
                nrofweights2 = nrofweights4;
                stepsize2    = stepsize4;
                kernelptr2   = kernelptr4;
            } else if (smallcount1 == 2) {
                nrofweights3 = nrofweights4;
                stepsize3    = stepsize4;
                kernelptr3   = kernelptr4;
            }
        }
        if (smallcount1 < nrofdims) {
            nrofweights4 = nrofweights3 = nrofweights2 = nrofweights1 = defnrofweights;
            stepsize4 = stepsize3 = stepsize2 = stepsize1 = defstepsize;
            kernelptr4 = kernelptr3 = kernelptr2 = kernelptr1 = defkernelptr;
        } else {
        }
    } else {
        nrofweights4 = nrofweights3 = nrofweights2 = nrofweights1 = defnrofweights;
        stepsize4 = stepsize3 = stepsize2 = stepsize1 = defstepsize;
        kernelptr4 = kernelptr3 = kernelptr2 = kernelptr1 = defkernelptr;
    }
    dblstepsize1 = (double) stepsize1;
    dblstepsizequot1 = 1.0 / dblstepsize1;
    dblstepsize2 = (double) stepsize2;
    dblstepsizequot2 = 1.0 / dblstepsize2;
    dblstepsize3 = (double) stepsize3;
    dblstepsizequot3 = 1.0 / dblstepsize3;
    dblstepsize4 = (double) stepsize4;
    dblstepsizequot4 = 1.0 / dblstepsize4;

    /* default value given */
    if ((nrhs > 4) &&
        mxIsDouble(prhs[4]) &&
        (mxGetNumberOfElements(prhs[4]) == 1)) {
        defsample = *((const double *) mxGetData(prhs[4]));
    }

    /* default value given */
    if ((nrhs > 5) &&
        (nrofdims == 3) &&
        mxIsDouble(prhs[5]) &&
        (mxGetNumberOfElements(prhs[5]) == 16)) {
        quotptr = (const double *) mxGetData(prhs[5]);
        for (largecount = 0; largecount < 16; ++largecount) {
            IF_IS_BAD_VAL(quotptr[largecount])
                break;
        }
        if ((largecount == 16) &&
            (quotptr[3] == 0) &&
            (quotptr[7] == 0) &&
            (quotptr[11] == 0) &&
            (quotptr[15] == 1)) {
            GET_QUATVARS(quotptr);
            if (jobkind == coordrange)
                jobkind = transrange;
        }
    }

    /* compute kernel size, kernel midpoint and re-set weight pointer */
    kernelsize1 = (int) (((double) (nrofweights1 - 1) / ((double) (2 * stepsize1))));
    kernelelements1 = 2 * kernelsize1 + 1;
    kernellimpos1 = (int) (((double) (nrofweights1 - 1)) / 2.0);
    kernellimneg1 = -kernellimpos1;
    kernelptr1 = &kernelptr1[kernellimpos1];
    if (nrofdims > 1) {
        kernelsize2 = (int) (((double) (nrofweights2 - 1) / ((double) (2 * stepsize2))));
        kernelelements2 = 2 * kernelsize2 + 1;
        kernellimpos2 = (int) (((double) (nrofweights2 - 1)) / 2.0);
        kernellimneg2 = -kernellimpos2;
        kernelptr2 = &kernelptr2[kernellimpos2];
    }
    if (nrofdims > 2) {
        kernelsize3 = (int) (((double) (nrofweights3 - 1) / ((double) (2 * stepsize3))));
        kernelelements3 = 2 * kernelsize3 + 1;
        kernellimpos3 = (int) (((double) (nrofweights3 - 1)) / 2.0);
        kernellimneg3 = -kernellimpos3;
        kernelptr3 = &kernelptr3[kernellimpos3];
    }
    if (nrofdims > 3) {
        kernelsize4 = (int) (((double) (nrofweights4 - 1) / ((double) (2 * stepsize4))));
        kernelelements4 = 2 * kernelsize4 + 1;
        kernellimpos4 = (int) (((double) (nrofweights4 - 1)) / 2.0);
        kernellimneg4 = -kernellimpos4;
        kernelptr4 = &kernelptr4[kernellimpos4];
    }

    /* for a coordinate range */
    if (jobkind != coordslist) {

        /* create output as an array of input dims */
        *outputdim = outputdim1;
        outputdim[1] = outputdim2;
        outputdim[2] = outputdim3;
        outputdim[3] = outputdim4;
        *plhs = mxCreateNumericArray(nrofdims, outputdim, mxDOUBLE_CLASS, mxREAL);

        /* get pointer to output */
        sampled = (double *) mxGetData(*plhs);

        /* without transformation matrix */
        if (jobkind == coordrange) {

            /* then allocate weight arrays */
            if (nrofdims > 1) {
                ALLOCATE_RANGE_CACHE(1)
                ALLOCATE_RANGE_CACHE(2)
                if (nrofdims > 2) {
                    ALLOCATE_RANGE_CACHE(3)
                    if (nrofdims > 3) {
                        ALLOCATE_RANGE_CACHE(4)
                    }
                }

                /* shift pointers as necessary */
                SHIFT_CACHE_POINTERS(1,kernelsize1)
                FILL_RANGE_CACHE(1)
                SHIFT_CACHE_POINTERS(2,kernelsize2)
                FILL_RANGE_CACHE(2)
                if (nrofdims > 2) {
                    SHIFT_CACHE_POINTERS(3,kernelsize3)
                        FILL_RANGE_CACHE(3)
                    if (nrofdims > 3) {
                        SHIFT_CACHE_POINTERS(4,kernelsize4)
                        FILL_RANGE_CACHE(4)
                    }
                }
            }

            /* processing depends in input class */
            switch (inputclass) {
                case mxDOUBLE_CLASS: SWITCH_CRANGE_OVER_DIMS_FP(inputdbl,double) break;
                case mxSINGLE_CLASS: SWITCH_CRANGE_OVER_DIMS_FP(inputflt,float) break;
                case  mxINT32_CLASS: SWITCH_CRANGE_OVER_DIMS(inputsl,signed int) break;
                case mxUINT32_CLASS: SWITCH_CRANGE_OVER_DIMS(inputul,unsigned int) break;
                case  mxINT16_CLASS: SWITCH_CRANGE_OVER_DIMS(inputss,signed short) break;
                case mxUINT16_CLASS: SWITCH_CRANGE_OVER_DIMS(inputus,unsigned short) break;
                case   mxINT8_CLASS: SWITCH_CRANGE_OVER_DIMS(inputsc,signed char) break;
                case  mxUINT8_CLASS:
                case mxLOGICAL_CLASS:SWITCH_CRANGE_OVER_DIMS(inputuc,unsigned char) break;
                default: mexErrMsgTxt("Datatype not supported.");
            }

            /* finally de-allocate weight arrays */
            if (nrofdims > 1) {

                /* un-shift pointers as necessary */
                if (nrofdims > 2) {
                    if (nrofdims > 3) {
                        SHIFT_CACHE_POINTERS(4,-kernelsize4)
                        ALLOCATE_RANGE_CACHE_FREE(4)
                    }
                    SHIFT_CACHE_POINTERS(3,-kernelsize3)
                    ALLOCATE_RANGE_CACHE_FREE(3)
                }
                SHIFT_CACHE_POINTERS(2,-kernelsize2)
                ALLOCATE_RANGE_CACHE_FREE(2)
                SHIFT_CACHE_POINTERS(1,-kernelsize1)
                ALLOCATE_RANGE_CACHE_FREE(1)
            }

        /* transformed coord range (data/coords in 3-D!) */
        } else {

            /* allocate weight arrays (only small ones, coordinates one by one) */
            ALLOCATE_COORD_CACHE

            /* processing depends in input class */
            switch (inputclass) {
                case mxDOUBLE_CLASS: SWITCH_TRANGE_OVER_DIMS(inputdbl,double) break;
                case mxSINGLE_CLASS: SWITCH_TRANGE_OVER_DIMS(inputflt,float) break;
                case  mxINT32_CLASS: SWITCH_TRANGE_OVER_DIMS(inputsl,signed int) break;
                case mxUINT32_CLASS: SWITCH_TRANGE_OVER_DIMS(inputul,unsigned int) break;
                case  mxINT16_CLASS: SWITCH_TRANGE_OVER_DIMS(inputss,signed short) break;
                case mxUINT16_CLASS: SWITCH_TRANGE_OVER_DIMS(inputus,unsigned short) break;
                case   mxINT8_CLASS: SWITCH_TRANGE_OVER_DIMS(inputsc,signed char) break;
                case  mxUINT8_CLASS:
                case mxLOGICAL_CLASS:SWITCH_TRANGE_OVER_DIMS(inputuc,unsigned char) break;
                default: mexErrMsgTxt("Datatype not supported.");
            }

            /* freeing the cache */
            ALLOCATE_COORD_CACHE_FREE
        }

    /* list of specific coordinates given */
    } else {

        /* create output as a Cx1 vector */
        *plhs = mxCreateDoubleMatrix(nrofcoords, 1, mxREAL);

        /* get pointer to output */
        sampled = (double *) mxGetData(*plhs);

        /* allocate weight arrays (only small ones, coordinates one by one) */
        ALLOCATE_COORD_CACHE

        /* processing depends on input class */
        switch (inputclass) {
            case mxDOUBLE_CLASS: SWITCH_COORDS_OVER_DIMS_FP(inputdbl,double) break;
            case mxSINGLE_CLASS: SWITCH_COORDS_OVER_DIMS_FP(inputflt,float) break;
            case  mxINT32_CLASS: SWITCH_COORDS_OVER_DIMS(inputsl,signed int) break;
            case mxUINT32_CLASS: SWITCH_COORDS_OVER_DIMS(inputul,unsigned int) break;
            case  mxINT16_CLASS: SWITCH_COORDS_OVER_DIMS(inputss,signed short) break;
            case mxUINT16_CLASS: SWITCH_COORDS_OVER_DIMS(inputus,unsigned short) break;
            case   mxINT8_CLASS: SWITCH_COORDS_OVER_DIMS(inputsc,signed char) break;
            case  mxUINT8_CLASS:
            case mxLOGICAL_CLASS:SWITCH_COORDS_OVER_DIMS(inputuc,unsigned char) break;
            default: mexErrMsgTxt("Datatype not supported.");
        }

        /* freeing the cache */
        ALLOCATE_COORD_CACHE_FREE
    }
}

/*

mesh_normals  - computing normals from a list of vertices and faces

FORMAT:       nrm = mesh_normals(c, t)

Input fields:

      c           Cx3 coordinates
      t           Tx3 triangles (faces)

Output fields:

      nrm         Cx3 normals


% Version:  v0.9d
% Build:    14062714
% Date:     Jun-27 2014, 2:26 PM EST
% Author:   Dirk-Jan Kroon, University of Twente, Enschede, NL
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

Copyright (c) 2010, 2014, Dirk-Jan Kroon
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* All inputs */
    const double *FacesA, *FacesB, *FacesC, *VerticesX, *VerticesY, *VerticesZ;

    /* All outputs, Vertex Normals */
    double *NormalsX, *NormalsY, *NormalsZ;

    /* Temporary Face Normals */
    double *FaceNormalsX, *FaceNormalsY, *FaceNormalsZ;

    /* Temporary Face angles */
    double *AnglesA, *AnglesB, *AnglesC;

    /* Number of faces */
    const mwSize *FacesDims;
    int FacesN = 0;

    /* Number of vertices */
    const mwSize *VertexDims;
    int VertexN = 0;
    int VertexNA[1] = {0};

    /* 1D Index  */
    int index0, index1, index2;

    /* Edge coordinates and lenght */
    double e0x, e0y, e0z, e0l;
    double e1x, e1y, e1z, e1l;
    double e2x, e2y, e2z, e2l;

    /* Length of normal */
    double nl;

    /* Loop variable */
    int i;

    /* Check for proper number of arguments. */
    if (nrhs != 2)
        mexErrMsgTxt("2 inputs are required.");
    else if (nlhs > 2)
        mexErrMsgTxt("Up to 1 output supported.");

    /* Get number of VertexN */
    VertexDims = mxGetDimensions(*prhs);
    if ((VertexDims[1] != 3) ||
        (mxGetClassID(*prhs) != mxDOUBLE_CLASS))
        mexErrMsgTxt("Coordinates must be Cx3 double.");
    VertexN = *VertexDims;
    VerticesX = (const double *) mxGetData(*prhs);
    VerticesY = &VerticesX[VertexN];
    VerticesZ = &VerticesY[VertexN];

    /* Get number of FacesN */
    FacesDims = mxGetDimensions(prhs[1]);
    if ((FacesDims[1] != 3) ||
        (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS))
        mexErrMsgTxt("Triangles must be Tx3 double.");
    FacesN = *FacesDims;

    /* Read all inputs (faces and vertices) */
    FacesA = (const double *) mxGetData(prhs[1]);
    FacesB = &FacesA[FacesN];
    FacesC = &FacesB[FacesN];

    /* Create Output arrays for the Normal coordinates */
    *plhs = mxCreateDoubleMatrix(VertexN, 3, mxREAL);
    if (*plhs == NULL)
        mexErrMsgTxt("Error allocating output array.");
    NormalsX = (double *) mxGetData(*plhs);
    if (NormalsX == NULL)
        mexErrMsgTxt("Error retrieving output array pointer.");
    NormalsY = &NormalsX[VertexN];
    NormalsZ = &NormalsY[VertexN];

    /* temporary memory */
    if (nlhs == 1)
        FaceNormalsX = (double *) malloc(6 * FacesN * sizeof(double));
    else {
        plhs[1] = mxCreateDoubleMatrix(FacesN, 3, mxREAL);
        if (plhs[1] == NULL)
            mexErrMsgTxt("Error allocating second output array.");
        FaceNormalsX = (double *) mxGetData(plhs[1]);
    }
    if (FaceNormalsX == NULL)
        mexErrMsgTxt("Error allocating temporary memory.");
    FaceNormalsY = &FaceNormalsX[FacesN];
    FaceNormalsZ = &FaceNormalsY[FacesN];
    if (nlhs == 1)
        AnglesA = &FaceNormalsZ[FacesN];
    else {
        AnglesA = (double *) malloc(3 * FacesN * sizeof(double));
        if (AnglesA == NULL)
            mexErrMsgTxt("Error allocating temporary memory.");
    }        
    AnglesB = &AnglesA[FacesN];
    AnglesC = &AnglesB[FacesN];

    /* Calculate all face normals and angles */
    for (i = 0; i < FacesN; ++i) {

        /* Get indices of face vertices */
        index0 = (int) FacesA[i] - 1;
        index1 = (int) FacesB[i] - 1;
        index2 = (int) FacesC[i] - 1;

        /* Make edge vectors */
        e0x = VerticesX[index0] - VerticesX[index1];
        e0y = VerticesY[index0] - VerticesY[index1];
        e0z = VerticesZ[index0] - VerticesZ[index1];

        e1x = VerticesX[index1] - VerticesX[index2];
        e1y = VerticesY[index1] - VerticesY[index2];
        e1z = VerticesZ[index1] - VerticesZ[index2];

        e2x = VerticesX[index2] - VerticesX[index0];
        e2y = VerticesY[index2] - VerticesY[index0];
        e2z = VerticesZ[index2] - VerticesZ[index0];

        /* Normalize the edge vectors */
        e0l = sqrt(e0x * e0x + e0y * e0y + e0z * e0z) + 5e-16;
        e1l = sqrt(e1x * e1x + e1y * e1y + e1z * e1z) + 5e-16;
        e2l = sqrt(e2x * e2x + e2y * e2y + e2z * e2z) + 5e-16;
        e0x /= e0l; e0y /= e0l; e0z /= e0l;
        e1x /= e1l; e1y /= e1l; e1z /= e1l;
        e2x /= e2l; e2y /= e2l; e2z /= e2l;

        /* Calculate angles of face seen from vertices */
        AnglesA[i] = acos(e0x * (-e2x) + e0y * (-e2y) + e0z * (-e2z));
        AnglesB[i] = acos(e1x * (-e0x) + e1y * (-e0y) + e1z * (-e0z));
        AnglesC[i] = acos(e2x * (-e1x) + e2y * (-e1y) + e2z * (-e1z));

        /* Normal of the face */
        FaceNormalsX[i] = e0y * e2z - e0z * e2y;
        FaceNormalsY[i] = e0z * e2x - e0x * e2z;
        FaceNormalsZ[i] = e0x * e2y - e0y * e2x;
    }

    /* Calculate all vertex normals and angles */
    for (i = 0; i < FacesN; ++i) {

        index0 = (int) FacesA[i] - 1;
        index1 = (int) FacesB[i] - 1;
        index2 = (int) FacesC[i] - 1;
        NormalsX[index0] += FaceNormalsX[i] * AnglesA[i];
        NormalsY[index0] += FaceNormalsY[i] * AnglesA[i];
        NormalsZ[index0] += FaceNormalsZ[i] * AnglesA[i];
        NormalsX[index1] += FaceNormalsX[i] * AnglesB[i];
        NormalsY[index1] += FaceNormalsY[i] * AnglesB[i];
        NormalsZ[index1] += FaceNormalsZ[i] * AnglesB[i];
        NormalsX[index2] += FaceNormalsX[i] * AnglesC[i];
        NormalsY[index2] += FaceNormalsY[i] * AnglesC[i];
        NormalsZ[index2] += FaceNormalsZ[i] * AnglesC[i];
    }

    /* Normalize the Normals */
    for (i = 0; i < VertexN; ++i) {

        nl = sqrt(NormalsX[i] * NormalsX[i] + NormalsY[i] * NormalsY[i] + NormalsZ[i] * NormalsZ[i]) + 5e-16;
        NormalsX[i] /= nl;
        NormalsY[i] /= nl;
        NormalsZ[i] /= nl;
    }

    /* Free memory */
    if (nlhs == 1)
        free(FaceNormalsX);
    else {

        /* also normalize normals for faces */
        for (i = 0; i < FacesN; ++i) {
            e0x = FaceNormalsX[i];
            e0y = FaceNormalsY[i];
            e0z = FaceNormalsZ[i];
            e0l = sqrt(e0x * e0x + e0y * e0y + e0z * e0z);
            FaceNormalsX[i] /= e0l;
            FaceNormalsY[i] /= e0l;
            FaceNormalsZ[i] /= e0l;
        }

        /* then free memory */
        free(AnglesA);
    }
}

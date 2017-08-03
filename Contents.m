% NeuroElf
% Version 1.1        (6291) May2016
%
% Copyright (c) 2010 - 2016, Jochen Weber
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
%
% Matlab toolbox developed, written and maintained by Jochen Weber
% with the help and inspiration of
%
% - Hester Breman, Brain Innovation, B.V.
% - Federico Demartino, Brain Innovation, B.V.
% - Fabrizio Esposito, Brain Innovation, B.V.
% - Elia Formisano, Faculty of Psychology, University of Maastricht
% - Rainer Goebel, Faculty of Psychology, University of Maastricht
% - Armin Heinecke, Brain Innovation, B.V.
% - Pim Pullens, Brain Innovation, B.V.
% - Alard Roebroeck, Faculty of Psychology, University of Maastricht
% - Tor Wager, Department of Psychology, Columbia University
%
% additional thanks go to the SCAN Unit of Columbia University to allow
% me to continue working on this toolbox.
%
% core classes:
%
% @neuroelf     - class to "hide" functions from cluttering path
% @transimg     - class to handle stacking of (transparent) images (overlay)
% @transio      - transparent IO access into binary files (memory saving)
% @xff          - class to allow file I/O of diverse file formats
% @xfigure      - GUI handling
% @xini         - configuration file handling
% @xprogress    - simplified progress bar (without xfigure)
%
% main @neuroelf methods:
%
% acpc2tal        - convert ACPC based coordinates to TAL
% acsvread        - read tabular ASCII file into cell/struct array
% alphasim        - simulate smoothed maps for cluster size estimate
% analyzetype     - resolving Analyze datatypes
% any2ascii       - convert datatypes into ascii representation
% applybvtrf      - apply one of BV's 4x4 TRF matrices to coordinates
% applyfdr        - apply the FDR correction to threshold data
% asciiread       - read a ascii text file into a 1xN char array
% asciitab        - create custom ASCII file of 2D double array
% asciiwrite      - write a ascii text file from a 1xN char array
% autocorr        - compute auto-correlation coefficient
% averagevmrs     - build an average of a list of VMRs
% barycoord       - compute barycentric coordinates
% binread         - read a binary file into a 1xN uint8 array
% binwrite        - write a binary file from a 1xN char/uint8 array
% bitdump         - dump 8-bit bitfield of uint8 array
% bvcoordconv     - conversion chain between different systems
% bvtrf           - compute TRF matrix from 3D volume tools values
% bvxaddvartofile - adds a numeric variable to a BVX object file
% bvxcreatefile   - creates an empty BVX file container
% calcbetas       - compute beta values from design/time course data
% checkstruct     - check and populate structures
% checksyntax     - find basic errors in code
% clearxffobjects - clear the passed objects from the global storage
% clustercoords   - cluster a (binary) volume
% clustermeshmap  - cluster a surface map with threshold and options
% clustervol      - cluster a (data) volume with threshold and options
% colorpicker     - UI-based color selector
% conjval         - conjunction value
% conjvalp        - conjunction value (p-value logic)
% conv3d          - relatively fast 3D convolution
% convones        - convolute a vector with an all-1's vector
% correlinvtstat  - compute correlation value from t stat
% correlpvalue    - convert an r into a p value
% correltstat     - convert an r into a t value
% cov_nd          - calculate covariance on N-D matrices
% cpfile          - copy a file with IO
% crawford_abnorm - calculate abnormlity of single subject's score
% crawford_diff   - calculate difference of single subject's score
% crawford_diss   - calculate dissociation of single subject's score
% createfmr       - create an FMR object from input files
% createvmr       - create a VMR object from DICOM or PAR/REC files
% ddeblank        - deblank strings on both ends
% degclust        - degree of clustering
% depcorrt        - dependent correlations t (XrY > XrZ); (X1rY1 > X2rY2)
% dicom4tofmr     - convert a 4D DICOM file to FMR
% dilate3d        - 3d array dilation
% dispslice       - display 2-D data as image
% dispslicemovie  - display a 3-D array as series of images (movie)
% emptysrf        - create a truly empty SRF
% erode3d         - 3d array erosion
% exceltocoords   - get numeric coords from excel range string
% extcaller       - get name of file calling the caller
% fdr_thresholds  - compute FDR thresholds for MCP
% fileguessendian - guess endianness of file
% filesize        - get filesize
% fileswapendian  - swap endianness of file
% findfiledir     - find a named file or directory
% findfiles       - find files in a path (dir search)
% findfirst       - MEX function, stop after finding first
% fisherr2z       - apply Fisher z-transform on r values
% fitrobustbisquare - fit a model to data robustly with bisquare weights
% flexinterpn     - flexible interpolation (with arbitrary kernels)
% flexinterpn_method - interpolation applying standard kernels
% flexmask        - apply flexible data copying/masking
% floodfill3      - floodfill a (binary) volume from a seed
% fmriqasheet     - visualize fmriquality output
% fmriquality     - assess fMRI data quality
% gammapdf        - return gamma PDF for X, h, l
% glmtstat        - compute GLM t statistic
% gluetostring    - create a CSV list from a cell array
% gradsmooth      - apply gradient weighted smoothing (like sigma filter)
% grep            - Linux/Unix like grep command for files/strings
% hasfields       - checks a struct if it (only) has certain fields
% heartbeats      - detect heartbeats in ECG data
% hexdump         - hexdump of uint8/char/double content
% hrf             - calculate canonical HRF (single- or two-gamma)
% hxdouble        - internal double representation (hex)
% hxsingle        - internal single representation (hex)
% icc             - intra class correlation
% importbesa      - import BESA file formats
% importvmpfromspms - import VMP file from SPM maps/images
% initialalignment - compute IA matrix from FMR/VMR headers
% invsystem       - inverse system output arguments order
% isabsolute      - check paths for absolute content
% isrealvarname   - check whether a string is a varname (enhanced)
% lsqueeze        - make input Nx1 element vector
% ltriasc         - ASCII form of a lower triangle matrix (for Mx)
% mainver         - return MATLAB main version (for compatibility)
% makelabel       - make a label from a string
% maxpos          - only return position of maximum (uses max)
% minarray        - minimize array with thresholds
% minmaxbbox      - generate minimum and maximum bounding boxes
% minpos          - only return position of minimum (uses min)
% mkadir          - make all dir (incl. parents)
% mni2tal         - convert MNI to TAL space coordinates
% mstrrep         - apply strrep multiple times
% multiset        - set properties for multiple handles
% multislicetovol - convert a multislice raw image to an img/hdr vol
% newnatresvmp    - create empty NR-VMP
% ne_gzip         - pass gzip on to system command (if available)
% normvecs        - normalize vectors (to SS/length 1.0)
% orthimage       - return orthogonality image (scaled 0 .. 1)
% orthsdm         - orthogonalize an SDM object
% orthvec         - orthogonalize vectors
% orthvecs        - orthogonalize set of vectors
% ostype          - general hub for OS specific tests
% osz             - ones with the size of the input argument
% parseopts       - parse options (varargin)
% progresscount   - progress counter (console)
% psctrans        - percent signal transformation
% pwelch          - create Welch periodogram
% randm           - random numbers from different distributions
% readbesa        - read BESA file formats
% readpar         - Philips PAR/REC reader
% relfilename     - build relative filename
% renamefields    - rename fields in a struct
% renamefile      - lazy file renaming
% render_volume   - simple 3D rendering of 3-D array data
% robscatter      - robustly fit data vector Y to model X and scatter
% robustnsamplet  - perform robust N-sample t-test
% robustt         - compute t-scores on contrasts after robust fitting
% samplefmrspace  - sample BV's FMR space at given coordinates
% scaleimage      - display 2-dim data scaled automatically
% smoothdata3     - 3D data smoothing
% smoothest3d     - 3D smoothness estimation
% smoothkern      - gaussian kernel generation (arbitrary dims)
% spherecoords    - compute spherical from euclidean coordinates
% spherecoordsinv - compute euclidean from spherical coordinates
% spherevmr       - create pseudo VMR with "spherical" blob
% splitclustercoords - find sub-clusters in larger cluster (watershed)
% splittocell     - split a CSV list into a cell array
% spmmat2prt      - convert SPM2 design matrix into a PRT
% subget          - get a (sub) list of properties as a struct
% tal2acpc        - convert TAL coordinates to ACPC space
% tcplot          - fancy time course plotting
% tdclient        - neat interface to tdlocal2
% tdlocal2        - local TD client, based on v2.0
% tfmatrix        - create 4x4 transformation matrices
% unzerostring    - remove trailing zeros from string
% uunion          - unsorted union
% vtc_concat      - concatenate VTCs
% writestcfiles   - create STC files from a 4-D array
% zerodstring     - create zero padded string
% zsz             - zeros with the size of the input argument
% ztrans          - z transformation
%
% UI function:
%
% neuroelf_gui  - GUI for the toolbox
%
% UI sub-functions:
%
% screenshot    - create a screenshot (file) from a window or axes

% Version:  v1.1
% Build:    16053017
% Date:     May-30 2016, 5:54 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

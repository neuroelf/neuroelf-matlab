function cfr = hdr_CoordinateFrame(xo, vol)
% HDR::CoordinateFrame  - get coordinate frame from Analyze image
%
% FORMAT:       cframe = hdr.CoordinateFrame([vol])
%
% Input fields:
%
%       vol         1x1 volume number (default: last for 4-D files)
%
% Output fields:
%
%       cframe      struct with fields
%        .DimX/Y/Z  spatial dimension of image
%        .DimT      number of volumes for multi-volume image, default 1
%        .Dimensions combined list
%        .ResX/Y/Z  spatial resolution (in TAL/MNI/real world mm)
%        .Resolution combined list
%        .Slice1Center  center coordinate of first slice
%        .SliceNCenter  center coordinate of last slice
%        .SystemOrigin  coordinate system origin
%        .RowDir    equivalent to AffineTransX(1:3)
%        .ColDir    equivalent to AffineTransY(1:3)
%        .SlcDir    equivalent to AffineTransZ(1:3)
%        .IsRadiological  flag whether coordinate system is left-handed
%        .Trf       4x4 quaternion matrix, so that
%                   cframe.Trf * voxel := mm
%          - and -  inv(cframe.Trf) * mm := voxel
%
% Using: spmq2m, spmtrf.

% Version:  v1.1
% Build:    16021412
% Date:     Feb-14 2016, 12:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% neuroelf library and global config needed
global ne_methods xffsngl;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isa(vol, 'double') || numel(vol) ~= 1 || isinf(vol) || isnan(vol) || vol < 1
    vol = Inf;
else
    vol = floor(vol);
end
bc = xo.C;
cfr = struct;

% check for information from additional MAT file
mfilec = [];
mfilel = false;
try
    if isfield(bc.RunTimeVars, 'Mat44') && (isequal(size(bc.RunTimeVars.Mat44), [4, 4, 1]) || ...
        isequal(size(bc.RunTimeVars.Mat44), [4, 4, size(bc.VoxelData, 4)]))
        vol = min(vol, size(bc.RunTimeVars.Mat44, 3));
        mfilec = bc.RunTimeVars.Mat44(:, :, vol);
        if ~all(mfilec(:) == 0)
            mfilel = true;
        else
            mfilec = [];
        end
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end

% get dimensions
cfr.DimX = bc.ImgDim.Dim(2);
cfr.DimY = bc.ImgDim.Dim(3);
if bc.ImgDim.Dim(1) > 2
    cfr.DimZ = bc.ImgDim.Dim(4);
    resz = bc.ImgDim.PixSpacing(4);
else
    cfr.DimZ = 1;
    resz = 1;
end
if bc.ImgDim.Dim(1) > 3
    cfr.DimT = bc.ImgDim.Dim(5);
else
    cfr.DimT = 1;
end
cfr.Dimensions = [cfr.DimX, cfr.DimY, cfr.DimZ, cfr.DimT];
dimh = 0.5 + 0.5 * cfr.Dimensions(1:3);

% get resolution
cfr.ResX = bc.ImgDim.PixSpacing(2);
cfr.ResY = bc.ImgDim.PixSpacing(3);
cfr.ResZ = resz;

% set full resolution element
cfr.Resolution = [cfr.ResX, cfr.ResY, cfr.ResZ];

% compute 1-based transformation matrix (as used in SPM)
if ~mfilel

    % for real (full) NIFTI objects
    if bc.NIIFileType > 0

        % get NIftI data
        q = bc.DataHist.NIftI1;

        % first check SFormCode
        if q.SFormCode > 0

            % use AffineTransX/Y/Z
            mfilec = [q.AffineTransX; q.AffineTransY; q.AffineTransZ; 0, 0, 0, 1];
            if ~all(mfilec(1:3, 4) == 0)
                mfilec(1:3, 4) = mfilec(1:3, 4) - mfilec(1:3, 1:3) * [1; 1; 1];
                mfilel = true;
            end
        end

        % next check QFormCode
        if ~mfilel && q.QFormCode > 0

            % use that information instead
            if bc.ImgDim.PixSpacing(1) < 0
                mfilec = ne_methods.spmtrf([q.QuatOffsetX, q.QuatOffsetY, q.QuatOffsetZ]) * ...
                    ne_methods.spmq2m([q.QuaternionB, q.QuaternionC, q.QuaternionD]) * ...
                    ne_methods.spmtrf([0,0,0], [0,0,0], [1, 1, -1] .* bc.ImgDim.PixSpacing(2:4));
            else
                mfilec = ne_methods.spmtrf([q.QuatOffsetX, q.QuatOffsetY, q.QuatOffsetZ]) * ...
                    ne_methods.spmq2m([q.QuaternionB, q.QuaternionC, q.QuaternionD]) * ...
                    ne_methods.spmtrf([0,0,0], [0,0,0], bc.ImgDim.PixSpacing(2:4));
            end
            if ~all(mfilec(1:3, 4) == 0)
                mfilec(1:3, 4) = mfilec(1:3, 4) - mfilec(1:3, 1:3) * [1; 1; 1];
                mfilel = true;
            end
        end

    % let's try the old school stuff...
    else

        % an originator is given (SPM2)
        if ~all(bc.DataHist.OriginSPM(1:3) == 0)

            % use this information
            dho = bc.DataHist.OriginSPM(1:3);
            if all(bc.DataHist.OriginSPM(1:3) < 0)
                dho = cfr.Dimensions(1:3) + dho;
            end

            % then create a transformation matrix
            mfilec = [ ...
                cfr.ResX,     0   ,     0   , -cfr.ResX * dho(1); ...
                    0   , cfr.ResY,     0   , -cfr.ResY * dho(2); ...
                    0   ,     0   , cfr.ResZ, -cfr.ResZ * dho(3); ...
                    0   ,     0   ,     0   ,     1];
            mfilel = true;

        % for older images
        elseif bc.NIIFileType == 0

            % use the orientation field
            switch (double(bc.DataHist.Orientation))
                case 0
                    mfilec = [ ...
                        cfr.ResX,     0   ,     0   , -cfr.ResX * dimh(1); ...
                            0   , cfr.ResY,     0   , -cfr.ResY * dimh(2); ...
                            0   ,     0   , cfr.ResZ, -cfr.ResZ * dimh(3); ...
                            0   ,     0   ,     0   ,      1];
                case 1
                    mfilec = [ ...
                        cfr.ResX,     0   ,     0   , -cfr.ResX * dimh(1); ...
                            0   ,     0   , cfr.ResZ, -cfr.ResZ * dimh(3); ...
                            0   , cfr.ResY,     0   , -cfr.ResY * dimh(2); ...
                            0   ,     0   ,     0   ,      1];
                case 2
                    mfilec = [ ...
                            0   , cfr.ResY,     0   , -cfr.ResY * dimh(2); ...
                            0   ,     0   , cfr.ResZ, -cfr.ResZ * dimh(3); ...
                        cfr.ResX,     0   ,     0   , -cfr.ResX * dimh(1); ...
                            0   ,     0   ,     0   ,      1];
                case 3
                    mfilec = [ ...
                        cfr.ResX,     0   ,     0   , -cfr.ResX * dimh(1); ...
                            0   ,-cfr.ResY,     0   ,  cfr.ResY * dimh(2); ...
                            0   ,     0   , cfr.ResZ, -cfr.ResZ * dimh(3); ...
                            0   ,     0   ,     0   ,      1];
                case 4
                    mfilec = [ ...
                        cfr.ResX,     0   ,     0   , -cfr.ResX * dimh(1); ...
                            0   ,     0   ,-cfr.ResZ,  cfr.ResZ * dimh(3); ...
                            0   , cfr.ResY,     0   , -cfr.ResY * dimh(2); ...
                            0   ,     0   ,     0   ,      1];
                case 5
                    mfilec = [ ...
                            0   , cfr.ResY,     0   , -cfr.ResY * dimh(2); ...
                            0   ,     0   ,-cfr.ResZ,  cfr.ResZ * dimh(3); ...
                        cfr.ResX,     0   ,     0   , -cfr.ResX * dimh(1); ...
                            0   ,     0   ,     0   ,      1];
                otherwise
                    warning('neuroelf:xff:invalidObject', ...
                        'Unknown DataHist.Orientation value: %d.', double(bc.DataHist.Orientation));
            end
        end
    end
end

% last resort
if ~mfilel
    mfilec = [ ...
        cfr.ResX,     0   ,     0   , -cfr.ResX * dimh(1); ...
            0   , cfr.ResY,     0   , -cfr.ResY * dimh(2); ...
            0   ,     0   , cfr.ResZ, -cfr.ResZ * dimh(3); ...
            0   ,     0   ,     0   ,      1];
end

% convention guessing needed ?
if ~mfilel && (bc.ImgDim.PixSpacing(1) < 0 || ...
    (bc.ImgDim.PixSpacing(1) == 0 && bc.DataHist.Orientation == 0 && xffsngl.CONF.type.hdr.assumeflipped))
    if cfr.ResX < 0
        cfr.ResX = -cfr.ResX;
    end

    % perform x-flip to get to TAL
    mfilec(1, :) = -mfilec(1, :);
end

% compute slice 1 and N center coordinates (Trf is one-based!)
sl1c = mfilec * [0.5 + cfr.DimX / 2; 0.5 + cfr.DimY / 2; 0.5; 1];
slnc = mfilec * [0.5 + cfr.DimX / 2; 0.5 + cfr.DimY / 2; 0.5 + cfr.DimZ; 1];

% get direction components
dcp = -mfilec(1:3, 1:3);
dcp = dcp ./ repmat(sqrt(sum(dcp .* dcp, 1)), [3, 1]);
dcp(:, 4) = -cross(dcp(:, 1), dcp(:, 2));

% make settings
cfr.Slice1Center = sl1c(1:3)';
cfr.SliceNCenter = slnc(1:3)';
cfr.SystemOrigin = mfilec(1:3, 4)';
cfr.RowDir = dcp(:, 1)';
cfr.ColDir = dcp(:, 2)';
cfr.SlcDir = dcp(:, 3)';
cfr.IsRadiological = (sum(dcp(:, 3) .* dcp(:, 4)) > 0);
cfr.Trf = mfilec;

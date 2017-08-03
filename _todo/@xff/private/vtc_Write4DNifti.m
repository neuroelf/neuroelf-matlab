function xo = vtc_Write4DNifti(xo, outfile, range)
% VTC::Write4DNifti  - writes analyze vols for the (entire) TC
%
% FORMAT:       vtc.Write4DNifti(outfile, range)
%
% Input fields:
%
%       outfile     output filename (if empty, replace .vtc with .nii)
%       range       range of volumes (default: all)
%
% No output fields.
%
% Using: bvcoordconv.

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:03 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library
global ne_methods;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get number of volumes
bc = xo.C;
vr = bc.Resolution;
vs = size(bc.VTCData);
numvol = vs(1);
vs(1) = [];

% check further arguments
if nargin < 2 || ~ischar(outfile) || numel(outfile) < 5
    if isempty(xo.F)
        error('neuroelf:xff:filenameRequired', 'VTC unsaved, filename required.');
    end
    outfile = [xo.F(1:end-4) '.nii'];
else
    outfile = outfile(:)';
    if ~strcmpi(outfile(end-3:end), '.nii')
        outfile = [outfile '.nii'];
    end
end
if nargin < 3 || ~isa(range, 'double') || isempty(range) || numel(range) ~= max(size(range)) || ...
    any(isinf(range) | isnan(range) | range < 1 | range > numvol | range ~= fix(range))
    range = [];
else
    range = range(:)';
end

% big TRY/CATCH
nio = {[]};
try

    % create temporary NII object
    nio{1} = xff('new:nii');
    nii = nio{1}.C;

    % general settings
    nii.NIIFileType = 2;
    nii.FileMagic = 'n+1';

    % set VoxelData
    nii.VoxelData = permute(bc.VTCData(:, end:-1:1, end:-1:1, end:-1:1), [4, 2, 3, 1]);
    if ~isempty(range)
        nii.VoxelData = nii.VoxelData(:, :, :, range);
    end

    % set size, type, etc.
    nii.ImgDim.Dim(1:5) = [4, size(nii.VoxelData)];
    if strcmpi(class(nii.VoxelData), 'uint16')
        nii.ImgDim.DataType = 512;
        nii.ImgDim.BitsPerPixel = 16;
    else
        nii.ImgDim.DataType = 16;
        nii.ImgDim.BitsPerPixel = 32;
    end
    nii.ImgDim.PixSpacing(2:4) = vr;

    % generate transformation matrix rows
    ti = vr .* (vs - ne_methods.bvcoordconv([0, 0, 0], 'tal2bvc', aft_BoundingBox(xo)));
    nii.DataHist.NIftI1.AffineTransX = [vr,  0,  0, -ti(3)];
    nii.DataHist.NIftI1.AffineTransY = [ 0, vr,  0, -ti(1)];
    nii.DataHist.NIftI1.AffineTransZ = [ 0,  0, vr, -ti(2)];

    % set in object
    nio{1}.C = nii;

    % save file
    aft_SaveAs(nio{1}, outfile);
catch xfferror
    clearxffobjects(nio);
    rethrow(xfferror);
end

% clear object
clearxffobjects(nio);

function nh = newhdr(dim, mat, opts)
% newhdr  - create a new HDR/NII object with given properties
%
% FORMAT:       nh = newhdr(dim, mat [, opts])
%
% Input fields:
%
%       dim         1x3 or 1x4 size of data
%       mat         4x4 transformation from voxel to mm coordinates
%       opts        optional settings
%        .data      data to put into .VoxelData (must match dim!)
%        .dtype     datatype (either string or numeric code)
%        .scaling   1x2 scaling intercept (b0) and scaling slope (b1)
%
% Output fields:
%
%       nh          HDR object with requested properties

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:28 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% argument check
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
   ~any([3, 4] == numel(dim)) || ...
    numel(dim) ~= size(dim, 2) || ...
    any(isinf(dim) | isnan(dim) | dim < 1 | dim ~= fix(dim)) || ...
   ~isa(mat, 'double') || ...
   ~isequal(size(mat), [4, 4]) || ...
    any(isinf(mat(:)) | isnan(mat(:))) || ...
    any(mat(4, :) ~= [0, 0, 0, 1])
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin > 2 && ...
    isnumeric(opts) && ...
    any(ndims(opts) == [3, 4])
    opts = struct('data', opts);
    opts.dtype = class(opts.data);
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'data') || ...
   ~isnumeric(opts.data) || ...
   ~isequal(size(opts.data), dim)
    opts.data = 0;
end
if ~isempty(opts.data)
    opts.dtype = class(opts.data);
elseif ~isfield(opts, 'dtype') || ...
   ((~ischar(opts.dtype) || ...
     ~any(strcmpi(opts.dtype(:)', ...
        {'double', 'int16', 'int32', 'int8', 'single', 'uint16', 'uint32', 'uint8'}))) && ...
    (~isa(opts.dtype, 'double') || ...
      numel(opts.dtype) ~= 1 || ...
      isinf(opts.dtype) || ...
      isnan(opts.dtype) || ...
     ~any(dtype == [2, 4, 8, 16, 64, 130, 132, 136])))
    opts.dtype = 'single';
end
if ischar(opts.dtype)
    switch (lower(opts.dtype(:)'))
        case {'double'}
            opts.dtype = 64;
        case {'int16'}
            opts.dtype = 4;
        case {'int32'}
            opts.dtype = 8;
        case {'int8'}
            opts.dtype = 130;
        case {'single'}
            opts.dtype = 16;
        case {'uint8'}
            opts.dtype = 2;
        case {'uint16'}
            opts.dtype = 132;
        case {'uint32'}
            opts.dtype = 136;
    end
end
[edata, stype, bsize] = analyzetype(opts.dtype);
if ~isequal(size(opts.data), dim)
    opts.data = edata;
    if numel(dim) == 3
        opts.data(dim(1), dim(2), dim(3)) = 0;
    else
        opts.data(dim(1), dim(2), dim(3), dim(4)) = 0;
    end
end
if ~isfield(opts, 'scaling') || ...
   ~isa(opts.scaling, 'double') || ...
    numel(opts.scaling) ~= 2 || ...
    any(isinf(opts.scaling) | isnan(opts.scaling)) || ...
    opts.scaling(2) == 0
    opts.scaling = [0, 1];
end
spc = sqrt(sum(mat(1:3, 1:3) .^ 2));
org = round(inv(mat) * [0; 0; 0; 1])';
vx0 = mat * ones(4, 1);
mat(1:3, 4) = vx0(1:3);

% create output
nh = xff('new:hdr');

% make settings
nh.FileMagic = 'n+1';
nh.NIIFileType = 2;
nh.ImgDim.Dim(1:1+numel(dim)) = [numel(dim), dim];
nh.ImgDim.DataType = opts.dtype;
nh.ImgDim.BitsPerPixel = 8 * bsize;
nh.ImgDim.PixSpacing(2:4) = spc;
nh.ImgDim.ScalingSlope = opts.scaling(2);
nh.ImgDim.ScalingIntercept = opts.scaling(1);
nh.DataHist.OriginSPM(1:3) = org(1:3);
nh.DataHist.NIftI1.QFormCode = 2;
nh.DataHist.NIftI1.SFormCode = 2;
nh.DataHist.NIftI1.QuatOffsetX = vx0(1);
nh.DataHist.NIftI1.QuatOffsetY = vx0(2);
nh.DataHist.NIftI1.QuatOffsetZ = vx0(3);
nh.DataHist.NIftI1.AffineTransX = mat(1, :);
nh.DataHist.NIftI1.AffineTransY = mat(2, :);
nh.DataHist.NIftI1.AffineTransZ = mat(3, :);
nh.VoxelData = opts.data;

function xo = hdr_InhomogeneityCorrect(xo, opts)
% HDR::InhomogeneityCorrect  - attempt automatic inhomogeneity correction
%
% FORMAT:       [hdr = ] hdr.InhomogeneityCorrect([opts])
%
% Input fields:
%
%       opts        optional struct with settings
%        .mask      either 3D uint8/logical data or HDR object with preseg
%                   if omitted, try automatic mask detection
%        .model     either of 'log', {'mult'}
%        .numpasses number of passes, default 3 (valid: 1 through 5)
%        .order     polynomial order, default 3 (valid: 2 through 7)
%        .xmask     use mask in conjunction with autodetected mask
%
% Using: pmbfilter.

% Version:  v1.1
% Build:    16020515
% Date:     Feb-05 2016, 3:36 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2012, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if any(bc.ImgDim.DataType == [128, 1536, 2048, 2304])
    error('neuroelf:xff:badArgument', 'Not defined for complex datatypes.');
end
if ndims(bc.VoxelData) ~= 3
    error('neuroelf:xff:badArgument', 'Only implemented for 3D data.');
end
if nargin < 2 || numel(opts) ~= 1 || ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'mask') || isempty(opts.mask)
    opts.mask = [];
end
if numel(opts.mask) == 1 && xffisobject(opts.mask, true, 'hdr')
    mbc = opts.mask.C;
    opts.mask = mbc.VoxelData(:, :, :);
    if ~isequal(size(opts.mask), size(bc.VoxelData))
        opts.mask = [];
    end
    if islogical(opts.mask)
        opts.mask = uint8(240) .* uint8(opts.mask);
    end
    if ~isa(opts.mask, 'uint8')
        opts.mask = uint8([]);
    end
    if ~isempty(opts.mask)
        opts.mask(opts.mask < 226) = 0;
        opts.mask(opts.mask > 225) = opts.mask(opts.mask > 225) - 225;
        um = double(unique(opts.mask(:))) + 1;
        ur = uint8(1:max(um));
        ur(um) = 0:(numel(um)-1);
        opts.mask = uint8(ur(opts.mask + 1));
    end
end
if ~isfield(opts, 'model') || ~ischar(opts.model) || ...
   ~any(strcmpi(opts.model(:)', {'l', 'log', 'm', 'mult'}))
    opts.model = 'mult';
else
    opts.model = lower(opts.model(1));
    if opts.model == 'l'
        opts.model = 'log';
    else
        opts.model = 'mult';
    end
end
if ~isfield(opts, 'numpasses') || numel(opts.numpasses) ~= 1 || ...
   ~isa(opts.numpasses, 'double') || isnan(opts.numpasses) || ...
   ~any((1:5) == opts.numpasses)
    opts.numpasses = 3;
end
if ~isfield(opts, 'order') || numel(opts.order) ~= 1 || ...
   ~isa(opts.order, 'double') || isnan(opts.order) || ~any((2:7) == opts.order)
    opts.order = 3;
end
if ~isfield(opts, 'xmask') || numel(opts.xmask) ~= 1 || ~islogical(opts.xmask)
    opts.xmask = false;
end

% get data
vd = aft_GetVolume(xo, 1);

% apply correction (pre-filter)
for pc = 1:(opts.numpasses-1)
    vd = ne_methods.pmbfilter(vd, opts.order, opts.mask, ...
        struct('bcutoff', 1, 'xmask', opts.xmask));
end

% apply final pass
vd = ne_methods.pmbfilter(vd, opts.order, opts.mask, struct( ...
    'bcutoff', 0.5, 'robust', true, 'xmask', opts.xmask));
vd(vd < 0) = 0;

% set back
bc.ImgDim.DataType = 16;
bc.ImgDim.BitsPerPixel = 32;
bc.ImgDim.ScalingSlope = 1;
bc.ImgDim.ScalingIntercept = 0;
bc.VoxelData = single(vd);

% set in output
xo.C = bc;

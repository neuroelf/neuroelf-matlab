function vmr = importtifasvmr(tiffile, opts)
% importtifasvmr  - import multi-slice TIF image and store as VMR
%
% FORMAT:       vmr = importtifasvmr(tiffile [, opts])
%
% Input fields:
%
%       tiffile     filename of TIF file to be imported (blank: uigetfile)
%       opts        optional settings
%        .res       VMR resolution (default: auto-set)
%
% Output fields:
%
%       vmr         VMR object storing the TIF data

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:21 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, 2016, Jochen Weber
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

% filename not given
vmr = [];
if nargin < 1 || ...
   ~ischar(tiffile) || ...
    isempty(tiffile) || ...
    size(tiffile, 1) ~= 1 || ...
    exist(tiffile, 'file') ~= 2
    [tiffile, tifpath] = uigetfile({'*.tif;*.tiff', 'TIF files (*.tif, *.tiff)'}, ...
        'Please select the TIF file to import as VMR...');
    if isempty(tiffile) || ...
        isequal(tiffile, 0) || ...
        isequal(tifpath, 0)
        return;
    end
    if isempty(tifpath)
        tifpath = pwd;
    end
    tiffile = [tifpath filesep tiffile];
    if exist(tiffile, 'file') ~= 2
        return;
    end
end

% options
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'res') || ...
   ~isa(opts.res, 'double') || ...
    numel(opts.res) ~= 1 || ...
    isinf(opts.res) || ...
    isnan(opts.res) || ...
    opts.res > 1 || ...
    opts.res < 0
    opts.res = 0;
end

% get info on file
try
    tifinfo = imfinfo(tiffile);
    
    % pre-allocate memory
    m = uint8(0);
    m(tifinfo(1).Height, tifinfo(1).Width, numel(tifinfo)) = 0;
    
    % determine resolution
    if opts.res == 0
        opts.res = 2 ^ floor(log2(256 / max(size(m))));
    end
catch ne_eo;
    rethrow(ne_eo);
end

% read into memory
try
    for sc = 1:size(m, 3)
        m(:, :, sc) = imread(tiffile, sc);
    end
catch ne_eo;
    rethrow(ne_eo);
end

% create VMR
vmr = xff('new:vmr');

% set data
vmr.VMRData = m;
vmr.DimX = size(m, 1);
vmr.DimY = size(m, 2);
vmr.DimZ = size(m, 3);

% set resolution
vmr.VoxResX = opts.res;
vmr.VoxResY = opts.res;
vmr.VoxResZ = opts.res;
vmr.SliceThickness = opts.res;
vmr.FramingCube = 2 ^ ceil(log2(max(size(m))));
vmr.NRows = vmr.FramingCube;
vmr.NCols = vmr.FramingCube;

% compute offsets
vmr.OffsetX = floor(0.5 * (vmr.FramingCube - size(m, 1)));
vmr.OffsetY = floor(0.5 * (vmr.FramingCube - size(m, 2)));
vmr.OffsetZ = floor(0.5 * (vmr.FramingCube - size(m, 3)));

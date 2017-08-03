function amr = fmr_CreateAMR(xo, opts)
% FMR::CreateAMR  - generate an AMR file
%
% FORMAT:       amr = fmr.CreateAMR([opts]);
%
% Input fields:
%
%       opts        optional struct with settings
%        .autoctr   auto contrast, default: false
%        .invert    invert intensities, default: true
%        .thresh    threshold, default: 0
%        .vol       1x1 or 1xN array with volumes to sample, default: 1
%
% Note: as the FMR object is NOT returned, the AMR reference is NOT set!
%
% Using: flexinterpn_method, minmaxmean.

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:54 AM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'autoctr') || numel(opts.autoctr) ~= 1 || ~islogical(opts.autoctr)
    opts.autoctr = false;
end
if ~isfield(opts, 'dim') || numel(opts.dim) ~= 2 || ~isa(opts.dim, 'double') || ...
    any(isnan(opts.dim) | opts.dim < 64 | opts.dim > 512 | opts.dim ~= fix(opts.dim))
    opts.dim = [256, 256];
else
    opts.dim = 2 * round(opts.dim(:)' / 2);
end
if ~isfield(opts, 'invert') || numel(opts.invert) ~= 1 || ~islogical(opts.invert)
    opts.invert = true;
end
if ~isfield(opts, 'thresh') || numel(opts.thresh) ~= 1 || ~isa(opts.thresh, 'double') || ...
    isnan(opts.thresh) || opts.thresh < 0 || opts.thresh > 32767
    opts.thresh = uint16(0);
else
    opts.thresh = uint16(floor(opts.thresh));
end
if ~isfield(opts, 'vol') || isempty(opts.vol) || ~isa(opts.vol, 'double') || ...
    any(isnan(opts.vol(:)) | opts.vol(:) < 1 | opts.vol(:) > 512 | bc.NrOfVolumes | ...
        opts.vol(:) ~= fix(opts.vol(:)))
    opts.vol = 1;
else
    opts.vols = unique(opts.vol(:)');
end

% make sure data is loaded
try
    if isempty(bc.Slice) || ~isstruct(bc.Slice) || ~isfield(bc.Slice, 'STCData')
        fmr_LoadSTC(xo);
        bc = xo.C;
    end
catch xfferror
    error('neuroelf:xff:internalError', 'Error loading slice data: %s.', xfferror.message);
end

% get FMR resolution
dmx = bc.ResolutionX;
dmy = bc.ResolutionY;
dmz = bc.NrOfSlices;

% build output
amr = xff('new:amr');
amrc = amr.C;
slc = amrc.Slice(1);

% init slice object
slc.BITMAPFILEHEADER.bfSize = 1078 + prod(opts.dim);
slc.BITMAPINFOHEADER.bfWidth = opts.dim(2);
slc.BITMAPINFOHEADER.bfHeight = opts.dim(1);
slc.BITMAPINFOHEADER.biImageSize = prod(opts.dim);
slc(2:dmz) = slc(1);

% compute interpolation matrix
stps = [dmy, dmx] ./ opts.dim(1:2);
is12 = [Inf, Inf; 1, 1; stps; 1 + (opts.dim(1:2) - 1) .* stps];

% build data
amrd = zeros([opts.dim, dmz]);

% interpolate data
flexinterpn_method = ne_methods.flexinterpn_method;
for sc = 1:dmz

    % depends on file version
    if bc.FileVersion < 5 || bc.DataStorageFormat == 1
        amrd(:, :, sc) = flexinterpn_method(mean( ...
            bc.Slice(sc).STCData(:, :, opts.vol), 3), is12, 0, 'cubic');
    elseif bc.DataStorageFormat == 2
        amrd(:, :, sc) = flexinterpn_method(mean( ...
            bc.Slice.STCData(:, :, opts.vol, sc), 3), is12, 0, 'cubic');
    else
        delete(amr);
        error('neuroelf:xff:invalidObject', 'Unsupported DataStorageFormat.');
    end
end

% contrast
amm = ne_methods.minmaxmean(amrd);
amrd = amrd - amm(1);
amrd = amrd .* (225 ./ (amm(2) - amm(1)));

% thresholding
if opts.thresh > 0
    opts.thresh = opts.thresh * (225 / (amm(2) - amm(1)));
    amrd(amrd < opts.thresh) = 0;
end

% store in AMR slices
for sc = 1:dmz
    slc(sc).AMRData = uint8(round(amrd(:, end:-1:1, sc)));
end

% fill object
amrc.NrOfSlices = dmz;
amrc.Slice = slc;
amr.C = amrc;

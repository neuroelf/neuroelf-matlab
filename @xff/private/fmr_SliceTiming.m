function xo = fmr_SliceTiming(xo, opts)
% FMR::SliceTiming  - correct for slice timing differences
%
% FORMAT:       [fmr = ] fmr.SliceTiming([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .interp    interpolation method, one of {'cubic'}, 'lanczos3'
%        .order     slice order, either 1xS spatial scanning order or one
%                   of {'asc', 'aint1', 'aint2', 'des', 'dint1', 'dint2'}
%        .refslice  reference slice, default: 1
%        .tshift    if given, either a 1x1 or 1xS double in ms
%
% Output fields:
%
%       fmr         FMR with corrected data
%
% Using: flexinterpn_method.

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:15 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
flexinterpn_method = ne_methods.flexinterpn_method;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'interp') || ~ischar(opts.interp) || ~any(strcmpi(opts.interp(:)', ...
    {'linear', 'cubic', 'lanczos2', 'lanczos3', 'lanczos5'}))
    opts.interp = 'cubic';
else
    opts.interp = lower(opts.interp(:)');
end
if ~isfield(opts, 'order') || isempty(opts.order) || ~any(strcmpi(class(opts.order), {'char', 'double'}))
    opts.order = bc.SliceAcquisitionOrder;
end
cord = {'asc', 'aint1', 'aint2', 'des', 'dint1', 'dint2'};
if ischar(opts.order) && ~any(strcmpi(opts.order(:)', cord))
    opts.order = bc.SliceAcquisitionOrder;
else
    opts.order = find(strcmpi(opts.order(:)', cord)) - 1;
    if opts.order > 2
        opts.order = opts.order + 7;
    end
end
opts.order = round(opts.order(:)');
opts.order(isinf(opts.order) | isnan(opts.order) | opts.order < 1 | opts.order > bc.NrOfSlices) = [];
if ~any(numel(opts.order) == [1, bc.NrOfSlices])
    opts.order = bc.SliceAcquisitionOrder;
end
if numel(opts.order) == 1
    switch (opts.order)
        case 1
            opts.order = [1:2:bc.NrOfSlices, 2:2:bc.NrOfSlices];
        case 2
            opts.order = [2:2:bc.NrOfSlices, 1:2:bc.NrOfSlices];
        case 10
            opts.order = bc.NrOfSlices:-1:1;
        case 11
            opts.order = [bc.NrOfSlices:-2:1, (bc.NrOfSlices - 1):-2:1];
        case 12
            opts.order = [(bc.NrOfSlices - 1):-2:1, bc.NrOfSlices:-2:1];
        otherwise
            warning('neuroelf:xff:unknownOrder', 'Unknown slice order %d.', opts.order);
            opts.order = 1:bc.NrOfSlices;
    end
end
if ~isfield(opts, 'refslice') || ~isa(opts.refslice, 'double') || numel(opts.refslice) ~= 1 || ...
   ~any(opts.refslice == (1:bc.NrOfSlices))
    opts.refslice = 1;
end
if ~isfield(opts, 'tshift') || ~isa(opts.tshift, 'double') || ...
   ~any(numel(opts.tshift) == [1, bc.NrOfSlices]) || numel(opts.tshift) ~= max(size(opts.tshift)) || ...
    any(isinf(opts.tshift) | isnan(opts.tshift) | opts.tshift <= -1 | opts.tshift >= 1) || ...
   (numel(opts.tshift) == 1 && (opts.tshift <= 0 || opts.tshift > (bc.TR / bc.NrOfSlices)))
    opts.tshift = bc.TR / bc.NrOfSlices;
end
if numel(opts.tshift) == 1
    opts.tshift = (opts.tshift / bc.TR) * ((0:(bc.NrOfSlices - 1)) - (opts.refslice - 1));
else
    opts.tshift = (1 / bc.TR) * opts.tshift(:)';
end

% depending on file version
if bc.FileVersion > 4 && numel(bc.Slice) == 1

    % do work
    if istransio(bc.Slice.STCData)
        bc.Slice.STCData = resolve(bc.Slice.STCData);
    end
    stcd = bc.Slice.STCData;

% otherwise pack and unpack
else

    stcd = bc.Slice(1).STCData;
    if istransio(stcd)
        stcd = resolve(stcd);
    end
    stcd(1, 1, 1, numel(bc.Slice)) = 0;
    for sc = 2:numel(bc.Slice)
        stcd(:, :, :, sc) = bc.Slice(sc).STCData(:, :, :);
    end
end

% datatype?
if strcmpi(class(stcd), 'uint16')
    isui16 = true;
else
    isui16 = false;
end

% extend time course by one at the beginning and the end
szs = size(stcd);
stcd(1, 1, end + 2, 1) = 0;
stcd(:, :, 2:end-1, :) = stcd(:, :, 1:end-2, :);
stcx = stcd(:, :, 2:6, :);
for ic = 1:4
    stcx = flexinterpn_method(stcx, [Inf * ones(1, 4); ...
        1, 1, 0.75, 1; ones(1, 4); size(stcx)], 0, 'lanczos3');
end
stcd(:, :, 1, :) = stcx(:, :, 1, :);
stcx = stcd(:, :, end-5:end-1, :);
szs(3) = szs(3) + 1;
for ic = 1:4
    stcx = flexinterpn_method(stcx, [Inf * ones(1, 4); ...
         1, 1, 1.25, 1; ones(1, 4); size(stcx) + [0, 0, 1, 0]], 0, 'lanczos3');
end
stcd(:, :, end, :) = stcx(:, :, end, :);

% iterate over slices
for sc = 1:numel(opts.order)

    % don't do anything for a shift of 0!
    if opts.tshift(sc) == 0
        continue;
    end

    % interpolate new data
    scslc = flexinterpn_method(stcd(:, :, :, opts.order(sc)), ...
        [Inf * ones(1, 3); 1, 1, 2 - opts.tshift(sc); ones(1, 3); szs(1:3)], ...
        0, opts.interp);

    % put back into data
    if bc.FileVersion > 4 && numel(bc.Slice) == 1
        bc.Slice.STCData(:, :, :, opts.order(sc)) = max(0, scslc);
    else
        if isui16
            bc.Slice(opts.order(sc)).STCData = uint16(round(max(0, scslc)));
        else
            bc.Slice(opts.order(sc)).STCData = single(max(0, scslc));
        end
    end
end

% set into content
xo.C = bc;

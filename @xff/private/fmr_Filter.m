function xo = fmr_Filter(xo, opts)
% FMR::Filter  - filter the STC data in an FMR
%
% FORMAT:       [fmr =] fmr.Filter(opts)
%
% Input fields:
%
%       opts        mandatory struct but with optional fields
%        .spat      enable spatial filtering (default: false)
%        .spkern    smoothing kernel in mm (default: [6, 6, 6])
%        .temp      enable temporal filtering (default: false)
%        .tempdct   DCT-based filtering (min. wavelength, default: Inf)
%        .tempdt    detrend (default: true, is overriden by dct/sc)
%        .templp    temporal lowpass (smoothing) kernel in secs (def: 0)
%        .tempsc    sin/cos set of frequencies (number of pairs, def: 0)
%
% Output fields:
%
%       fmr         FMR with filtered data
%
% Using: tempfilter.

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:45 AM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr') || ...
   ~isstruct(opts) || numel(opts) ~= 1 || (~isfield(opts, 'spat') && ~isfield(opts, 'temp'))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if ~isfield(opts, 'spat') || (~islogical(opts.spat) && ~isnumeric(opts.spat)) || numel(opts.spat) ~= 1
    opts.spat = false;
else
    opts.spat = (opts.spat ~= 0);
end
if ~isfield(opts, 'temp') || (~islogical(opts.temp) && ~isnumeric(opts.temp)) || numel(opts.temp) ~= 1
    opts.temp = false;
else
    opts.temp = (opts.temp ~= 0);
end
if ~isfield(opts, 'spkern') || ~isa(opts.spkern, 'double') || numel(opts.spkern) ~= 3 || ...
    any(isinf(opts.spkern) | isnan(opts.spkern) | opts.spkern < 0 | opts.spkern > 20)
    opts.spkern = [6, 6, 6];
else
    opts.spkern = opts.spkern(:)';
end
opts.spkern = opts.spkern ./ [bc.InplaneResolutionX, bc.InplaneResolutionY, bc.SliceThickness + bc.GapThickness];
if all(opts.spkern < 0.5)
    opts.spat = false;
end
if opts.temp
    if ~isfield(opts, 'tempdct') || ~isa(opts.tempdct, 'double') || numel(opts.tempdct) ~= 1 || ...
        isinf(opts.tempdct) || isnan(opts.tempdct) || opts.tempdct < 30
        opts.tempdct = Inf;
    else
        opts.tempdct = opts.tempdct ./ (0.001 * bc.TR);
    end
    if ~isfield(opts, 'tempdt') || ~islogical(opts.tempdt) || numel(opts.tempdt) ~= 1
        opts.tempdt = true;
    end
    if ~isfield(opts, 'templp') || ~isa(opts.templp, 'double') || numel(opts.templp) ~= 1 || ...
        isinf(opts.templp) || isnan(opts.templp) || opts.templp < 0 || opts.templp > 20
        opts.templp = 0;
    else
        opts.templp = opts.templp ./ (0.001 * bc.TR);
    end
    if ~isfield(opts, 'tempsc') || ~isa(opts.tempsc, 'double') || numel(opts.tempsc) ~= 1 || ...
       ~any(opts.tempsc == (1:12))
        opts.tempsc = 0;
    end
    if isinf(opts.tempdct) && ~opts.tempdt && opts.templp == 0 && opts.tempsc == 0
        opts.temp = false;
    elseif opts.tempsc > 0
        opts.tempdct = Inf;
        opts.tempdt = false;
    end
end
if ~opts.spat && ~opts.temp
    return;
end

% patch opts
opts.tdim = 3;

% depending on file version
if bc.FileVersion > 4 && numel(bc.Slice) == 1

    % do work
    bc.Slice.STCData = ne_methods.tempfilter(bc.Slice.STCData, opts);

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

    % do work
    stcd = ne_methods.tempfilter(stcd, opts);

    % re-pack
    for sc = 1:numel(bc.Slice);
        bc.Slice(sc).STCData = stcd(:, :, :, sc);
    end
end

% but content into array
xo.C = bc;

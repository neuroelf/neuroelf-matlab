function sdm = vtc_MeanSRFTimeCourse(xo, srf, varargin)
% VTC::MeanSRFTimeCourse  - extract mean SRF-based time course
%
% FORMAT:       sdm = vtc.CreateMTC(srf [, srf2, ...] [ , opts])
%
% Input fields:
%
%       srf         required surface file
%       srf2, ...   optionally additional surface files
%       opts        1x1 struct with optional fields
%        .interp    method ('nearest', {'linear'}, 'cubic')
%        .ipfrom    interpolate from P + n * normal vector, default: -3
%        .ipstep    interpolation stepsize, default: 1 (in normal vectors)
%        .ipto      interpolate to P + n * normal vector, default: 1
%        .method    method to get value ({'mean'}, 'median')
%        .remmean   remove mean from extracts (default: true)
%
% Output fields:
%
%       sdm         SDM object with one regressor per SRF
%
% Note: results slightly differ from BV's results (sampling method)
%
% Using: flexinterpn_method.

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:53 PM EST
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
flexinterpn_method = ne_methods.flexinterpn_method;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
[null, vtcf] = fileparts(xo.F);
if isempty(vtcf)
    vtcf = 'unknown';
end
bc = xo.C;
nsrf = nargin - 1;
ncrd = zeros(1, nsrf);
srff = cell(1, nsrf);
srfs = cell(1, nsrf);
srfs{1} = srf;
[null, srff{1}] = fileparts(srf.F);
if isempty(srff{1})
    srff{1} = 'unknown1';
end
ncrd(1) = size(srf.C.VertexCoordinate, 1);
for ac = 1:nargin-2
    if numel(varargin{ac}) == 1 && xffisobject(varargin{ac}, true, 'srf')
        srfs{ac+1} = varargin{ac};
        [null, srff{ac+1}] = fileparts(srfs{ac+1}.F);
        if isempty(srff{ac+1})
            srff{ac+1} = sprintf('unknown%d', ac+1);
        end
        ncrd(ac+1) = size(srfs{ac+1}.C.VertexCoordinate, 1);
    end
end
crds = zeros(sum(ncrd), 3);
nrms = zeros(sum(ncrd), 3);
for ac = 1:numel(srfs)
    if ncrd(ac) > 0
        crds(1+sum(ncrd(1:ac-1)):sum(ncrd(1:ac)), :) = srfs{ac}.C.VertexCoordinate(:, :);
        nrms(1+sum(ncrd(1:ac-1)):sum(ncrd(1:ac)), :) = srfs{ac}.C.VertexNormal(:, :);
    end
end
usrf = (ncrd > 0);
if ~all(usrf)
    srff(~usrf) = [];
    ncrd(~usrf) = [];
end
if nargin < 3 || ~isstruct(varargin{end}) || numel(varargin{end}) ~= 1
    opts = struct;
else
    opts = varargin{end};
end
if ~isfield(opts, 'interp') || ~ischar(opts.interp) || ...
   ~any(strcmpi(opts.interp(:)', {'cubic', 'linear', 'nearest'}))
    opts.interp = 'linear';
else
    opts.interp = lower(opts.interp(:))';
end
if ~isfield(opts, 'ipfrom') || ~isa(opts.ipfrom, 'double') || numel(opts.ipfrom) ~= 1 || ...
    isnan(opts.ipfrom) || opts.ipfrom < -10 || opts.ipfrom > 10
    opts.ipfrom = -3;
end
if ~isfield(opts, 'ipstep') || ~isa(opts.ipstep, 'double') || numel(opts.ipstep) ~= 1 || ...
    isnan(opts.ipstep) || opts.ipstep <= 0 || opts.ipstep > 10
    opts.ipstep = 1;
end
if ~isfield(opts, 'ipto') || ~isa(opts.ipto, 'double') || numel(opts.ipto) ~= 1 || ...
    isnan(opts.ipto) || opts.ipto < opts.ipfrom || opts.ipto > 10
    opts.ipto = 1;
end
if ~isfield(opts, 'method') || ~ischar(opts.method) || ...
   ~any(strcmpi(opts.method(:)', {'mean', 'median'}))
    opts.method = 'mean';
else
    opts.method = lower(opts.method(:))';
end
if ~isfield(opts, 'remmean') || ~islogical(opts.remmean) || numel(opts.remmean) ~= 1
    opts.remmean = true;
end
vtcres  = (1 / bc.Resolution);
ipsamp  = vtcres .* (opts.ipfrom:opts.ipstep:opts.ipto);
ipsnum  = numel(ipsamp);
if istransio(bc.VTCData)
    bc.VTCData = resolve(bc.VTCData);
end
numvols = size(bc.VTCData, 1);

% create output
sdm = xff('new:sdm');
sdmc = sdm.C;
sdmc.NrOfPredictors = numel(ncrd);
sdmc.NrOfDataPoints = numvols;
sdmc.IncludesConstant = 0;
sdmc.FirstConfoundPredictor = numel(ncrd) + 1;
sdmc.PredictorColors = floor(255.999 .* rand(numel(ncrd), 3));
sdmc.PredictorNames = cell(1, numel(ncrd));
for ac = 1:numel(ncrd)
    sdmc.PredictorNames{ac} = sprintf('%s_of_%s', srff{ac}, vtcf);
end
sdmm = zeros(numvols, numel(ncrd));

% subtract XStart, YStart, ZStart
crds = 1 + vtcres .* [crds(:, 1) - bc.XStart, crds(:, 2) - bc.YStart, crds(:, 3) - bc.ZStart];

% prepare interpolation coordinates
ipcrds = crds(:, :, ones(1, ipsnum));
for ipc = 1:ipsnum
    ipcrds(:, :, ipc) = crds + ipsamp(ipc) .* nrms;
end
ipcrds = cat(2, ones(size(crds, 1), 1, ipsnum), ipcrds);
ipvals = zeros(size(crds, 1), ipsnum);

% iterate over maps
for vc = 1:numvols

    % interpolate according to method
    ipcrds(:, 1, :) = vc;
    for ipc = 1:ipsnum
        ipvals(:, ipc) = flexinterpn_method(bc.VTCData, ipcrds(:, :, ipc), 0, opts.interp);
    end

    % what summary method
    switch opts.method
        case 'mean'
            for ac = 1:numel(ncrd)
                sdmm(vc, ac) = mean(mean(ipvals(1+sum(ncrd(1:ac-1)):sum(ncrd(1:ac)), :)));
            end
        case 'median'
            for ac = 1:numel(ncrd)
                sdmm(vc, ac) = mean(median(ipvals(1+sum(ncrd(1:ac-1)):sum(ncrd(1:ac)), :), 2));
            end
    end
end

% remove mean
if opts.remmean
    sdmm = sdmm - ones(numvols, 1) * mean(sdmm);
end

% put content in new object
sdmc.SDMMatrix = sdmm;
sdmc.RTCMatrix = sdmm;
sdm.C = sdmc;

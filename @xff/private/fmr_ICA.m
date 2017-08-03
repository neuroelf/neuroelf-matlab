function ica = fmr_ICA(xo, opts)
% FMR::ICA  - compute ICA over FMR data
%
% FORMAT:       ica = fmr.ICA([opts])
%
% Input fields:
%
%       opts        options settings
%        .mask      binary mask (must match)
%        .mthresh   masking threshold (only used if mask is empty, 100)
%        .numics    number of components to estimate (default: 30)
%        .tfiltfrq  number of temporal filtering frequencies (default: 0)
%        .tfilttyp  filtering type, either of {'DCT'}, 'Fourier'
%        .trobust   perform filtering robustly
%                   !! additionally all fields of ne_fastica.m supported !!
%
% Output fields:
%
%       ica         ICA object (FMR-based)
%
% Using: autocorr, degclust, findfirst, kurt, linscale, ne_fastica, pwelch,
%        skew, spatent, tempfilter.

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2015, 2016, Jochen Weber
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
if bc.FileVersion < 4 || bc.DataStorageFormat < 2
    vd = permute(bc.Slice(1).STCData(:, :, :), [3, 1, 2]);
    vd(1, 1, 1, numel(bc.Slice)) = 0;
    for slc = 2:numel(bc.Slice)
        vd(:, :, :, slc) = permute(bc.Slice(slc).STCData(:, :, :), [3, 1, 2]);
    end
else
    vd = permute(bc.Slice.STCData(:, :, :, :), [3, 1, 2, 4]);
end
vsz = size(vd);
ntp = vsz(1);
vsz(1) = [];
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'mask') || ~islogical(opts.mask) || ~isequal(vsz, size(opts.mask))
    opts.mask = [];
end
if ~isfield(opts, 'mthresh') || numel(opts.mthresh) ~= 1 || ~isnumeric(opts.mthresh) || ...
    isinf(opts.mthresh) || isnan(opts.mthresh) || opts.mthresh <= 0
    opts.mthresh = 100;
end
if ~isfield(opts, 'numics') || ~isa(opts.numics, 'double') || numel(opts.numics) ~= 1 || ...
    isinf(opts.numics) || isnan(opts.numics) || opts.numics < 1 || opts.numics > ntp || ...
    opts.numics ~= fix(opts.numics)
    opts.numics = 30;
end
if ~isfield(opts, 'tfiltfrq') || numel(opts.tfiltfrq) ~= 1 || ~isa(opts.tfiltfrq, 'double') || ...
    isnan(opts.tfiltfrq) || opts.tfiltfrq < 0 || opts.tfiltfrq > 12
    opts.tfiltfrq = 0;
end
if ~isfield(opts, 'tfilttyp') || ~ischar(opts.tfilttyp) || ...
   ~any(strcmpi(opts.tfilttyp(:)', {'dct', 'fourier'}))
    opts.tfilttyp = 'dct';
else
    opts.tfilttyp = lower(opts.tfilttyp(:)');
end
if ~isfield(opts, 'trobust') || ~islogical(opts.trobust) || numel(opts.trobust) ~= 1
    opts.trobust = false;
end

% force data to double
if ~isa(vd, 'double')
    vd = double(vd);
end

% create impromptu mask
mask = opts.mask;
if isempty(mask)
    mask = squeeze(sum(vd > opts.mthresh) >= 0.5);
end
mask = mask(:);

% for small mask, really do use indices instead!
if sum(mask) < (0.25 * numel(mask))
    mask = find(mask);
end

% get masked time-course
vd = vd(:, mask);

% filter content
if opts.tfiltfrq > 0

    % prepare tempfilter options
    topts = opts;
    topts.spat = false;
    topts.tdim = 1;
    topts.temp = true;
    topts.robust = opts.trobust;
    if opts.tfilttyp(1) == 'd'
        topts.tempdct = ceil(size(vd, 1) * bc.TR / opts.tfiltfrq);
        topts.tempsc = 0;
    else
        topts.tempdct = Inf;
        topts.tempsc = opts.tfiltfrq;
    end

    % temp filter data of first object
    vd = ne_methods.tempfilter(vd, topts);
end

% compute ica
[unmixed, tcs, w] = ne_methods.ne_fastica(vd, opts);

% get some values
nc = size(unmixed, 1);

% create ica
ica = xff('new:ica');
icac = ica.C;
icac.FileVersion = 6;
icac.DocumentType = 0;
icac.NrOfMaps = nc;
icac.NrOfTimePoints = ntp;
icac.NrOfMapParameters = 16;
icac.ShowParamsRangeFrom = 0;
icac.ShowParamsRangeTo = 5;
icac.FingerprintParamsRangeFrom = 5;
icac.FingerprintParamsRangeTo = 16;
icac.NrOfColumns = vsz(1);
icac.NrOfRows = vsz(2);
icac.NrOfSlices = vsz(3);
icac.OriginatingXTC = xo.F;
icac.LinkedPRT = '<none>';
if ~isempty(bc.ProtocolFile)
    if ischar(bc.ProtocolFile)
        icac.LinkedPRT = bc.ProtocolFile(:)';
    elseif iscell(bc.ProtocolFile)
        icac.LinkedPRT = bc.ProtocolFile{1}(:)';
    end
end
icac.Map.Type = 12;
icac.Map.LowerThreshold = 2;
icac.Map.UpperThreshold = 10;
icac.Map.UseRGBColor = 0;
icac.Map.Name = 'IC 1';
icac.Map.TimePointData = zeros(ntp, 1);
icac.Map.CMPData = single(zeros(vsz));
icac.Map = icac.Map(1, ones(1, nc));
icac.MapParameter(1).Name = 'DM pred max';
icac.MapParameter(1).Values = zeros(1, nc);
icac.MapParameter = icac.MapParameter(1, ones(1, nc));
icac.MapParameter(2).Name = 'DM pred index';
icac.MapParameter(3).Name = 'RMS';
icac.MapParameter(4).Name = 'Spatial Template (max corr)';
icac.MapParameter(5).Name = 'Spatial Template (VOI index)';
icac.MapParameter(6).Name = 'Degree of Clustering';
icac.MapParameter(7).Name = 'Skewness';
icac.MapParameter(8).Name = 'Kurtosis';
icac.MapParameter(9).Name = 'Spatial Entropy';
icac.MapParameter(10).Name = '1-Lag Autocorr';
icac.MapParameter(11).Name = 'Temporal Entropy';
icac.MapParameter(12).Name = 'Power Spectrum Band 1';
icac.MapParameter(13).Name = 'Power Spectrum Band 2';
icac.MapParameter(14).Name = 'Power Spectrum Band 3';
icac.MapParameter(15).Name = 'Power Spectrum Band 4';
icac.MapParameter(16).Name = 'Power Spectrum Band 5';

% initialize fingerprint data
fps = zeros(nc, 16);
fps(:, 2) = 1;

% find fingerprint frequency entries
[rp, rw] = ne_methods.pwelch(randn(ntp, 1));
rw = (1000 / (2 * pi * bc.TR)) .* rw;
wf = [0, 0.008, 0.02, 0.05, 0.1, 0.25];
wi = cell(1, 5);
for wc = 1:5
    wif = ne_methods.findfirst(rw > wf(wc));
    if isempty(wif)
        wif = numel(rw);
    end
    wit = ne_methods.findfirst(rw >= wf(wc + 1));
    if isempty(wit)
        wit = numel(rw);
    end
    wi{wc} = wif:wit;
end

% settings
res = [bc.InplaneResolutionX, bc.InplaneResolutionY, bc.SliceThickness + bc.GapThickness];
clt = ceil(270 / prod(res));

% iterate over components
autocorr = ne_methods.autocorr;
degclust = ne_methods.degclust;
kurt     = ne_methods.kurt;
linscale = ne_methods.linscale;
pwelch   = ne_methods.pwelch;
skew     = ne_methods.skew;
spatent  = ne_methods.spatent;
for cc = 1:nc

    % component and time-course (BV-style)
    cp = unmixed(cc, :);
    tc = linscale(tcs(:, cc));

    % set to map
    icac.Map(cc).Name = sprintf('IC %d', cc);
    icac.Map(cc).CMPData(mask) = cp;
    icac.Map(cc).TimePointData = tc;

    % compute fingerprint parameters -> RMS
    fps(cc, 3) = sqrt(sum((tc - mean(tc)) .^ 2) / ntp);

    % degree of clustering, skewness, kurtosis, spatial entropy
    fps(cc, 6) = degclust(abs(icac.Map(cc).CMPData), clt, 2, 2.5);
    fps(cc, 7) = log(abs(skew(cp)));
    fps(cc, 8) = abs(log(abs(kurt(cp) - 3)));
    fps(cc, 9) = abs(log(spatent(cp)));

    % lag-1 autocorrelation, temporal entropy
    fps(cc, 10) = abs(autocorr(tc));
    fps(cc, 11) = exp(abs(spatent(tc - mean(tc), 20)));

    % compute Welch's periodogram power spectrum
    pw = pwelch(tc);
    pw = (1 / sum(pw)) .* pw;

    % get energies and store last five parameters
    for pc = 1:5
        fps(cc, pc + 11) = sum(pw(wi{pc}));
    end
end

% for power-spectrum values, first make sum per component = 1
fps(:, 12:16) = ((1 ./ sum(fps(:, 12:16), 2)) * ones(1, 5)) .* fps(:, 12:16);

% now scale a couple of the parameters (all but DOC / 1-lag auto-corr.)
fps(:, [7:9, 11:16]) = linscale(fps(:, [7:9, 11:16]), 1, 0, 1);

% set parameters
for pc = 1:16
    icac.MapParameter(pc).Values = fps(:, pc)';
end

% keep w in RunTimeVars for now
icac.RunTimeVars.W_matrix = w;

% set content
ica.C = icac;

% and also set a reference to the FMR
aft_SetHandle(ica, 'FMR', xo);

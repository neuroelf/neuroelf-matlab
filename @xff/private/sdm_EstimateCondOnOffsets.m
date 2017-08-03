function conds = sdm_EstimateCondOnOffsets(xo, opts)
% SDM::EstimateCondOnOffsets  - estimate on and offsets from design matrix
%
% FORMAT:       conds = sdm.EstimateCondOnOffsets([opts])
%
% Input fields:
%
%       opts        optional settings
%        .ndsp      dispersion of negative gamma PDF (default: 1)
%        .nttp      time to negative (undershoot) peak (default: 15 secs)
%        .ons       onset of the hrf (default: 0 secs, OK: [-5 .. 5])
%        .pdsp      dispersion of positive gamma PDF (default: 1)
%        .pnr       pos-to-neg ratio (default: 6, OK: [1 .. Inf])
%        .pttp      time to positive (response) peak (default: 5 secs)
%        .regspaced flag, assume onsets are regularly spaced (false)
%        .restcond  flag, add rest condition to PRT (false)
%        .s         sampling range (default: [0, ons + 2 * (nttp + 1)])
%        .sbins     number of slices (default: 30, only if srbin > 1)
%        .sbinta    acquisition time for all slices (default: prtr)
%        .sbinte    time of echo (duration to acquire one slice, prtr/sbins)
%        .shape     hrf general shape, one of 'boynton', {'twogamma'}
%        .srbin     reference bin (ref slice in slice-time correction: 1)
%        .tr        TR (temporal resolution) of SDM in ms (2000)
%
% Output fields:
%
%       conds       1xC conditions struct with OnOffsets (PRT format)
%
% Using: convones, cov_nd, hrf, maxpos.

% Version:  v1.1
% Build:    16021210
% Date:     Feb-12 2016, 10:07 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, 2016, Jochen Weber
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
convones = ne_methods.convones;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'sdm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'measure') || ~ischar(opts.measure) || isempty(opts.measure) || ...
   ~any(strcmpi(opts.measure(:)', {'corr', 'sdiff'}))
    opts.measure = 'corr';
else
    opts.measure = opts.measure(:)';
end
if ~isfield(opts, 'ndsp') || ~isa(opts.ndsp, 'double') || numel(opts.ndsp) ~= 1 || ...
    isinf(opts.ndsp) || isnan(opts.ndsp) || opts.ndsp <= 0
    opts.ndsp = 1;
end
if ~isfield(opts, 'nttp') || ~isa(opts.nttp, 'double') || numel(opts.nttp) ~= 1 || ...
    isinf(opts.nttp) || isnan(opts.nttp) || opts.nttp < 1
    opts.nttp = 15;
end
if ~isfield(opts, 'ons') || ~isa(opts.ons, 'double') || numel(opts.ons) ~= 1 || ...
    isinf(opts.ons) || isnan(opts.ons) || abs(opts.ons) > 5
    opts.ons = 0;
end
if ~isfield(opts, 'pdsp') || ~isa(opts.pdsp, 'double') || numel(opts.pdsp) ~= 1 || ...
    isinf(opts.pdsp) ||isnan(opts.pdsp) || opts.pdsp <= 0
    opts.pdsp = 1;
end
if ~isfield(opts, 'pnr') || ~isa(opts.pnr, 'double') || numel(opts.pnr) ~= 1 || ...
    isnan(opts.pnr) || opts.pnr < 1
    opts.pnr = 6;
end
if ~isfield(opts, 'pttp') || ~isa(opts.pttp, 'double') || numel(opts.pttp) ~= 1 || ...
    isinf(opts.pttp) || isnan(opts.pttp) || opts.pttp <= 0
    opts.pttp = 5;
end
if opts.pttp >= opts.nttp
    error('neuroelf:xff:badArgument', 'Bad time-to-peak values in call to %s.', mfilename);
end
if ~isfield(opts, 'regspaced') || ~islogical(opts.regspaced) || numel(opts.regspaced) ~= 1
    opts.regspaced = false;
end
if ~isfield(opts, 'restcond') || ~islogical(opts.restcond) || numel(opts.restcond) ~= 1
    opts.restcond = false;
end
if ~isfield(opts, 's') || ~isa(opts.s, 'double') || numel(opts.s) ~= 1 || ...
    isinf(opts.s) || isnan(opts.s) || opts.s <= 0
    opts.s = 2 * opts.nttp + opts.ons + 1;
end
if ~isfield(opts, 'shape') || ~ischar(opts.shape) || isempty(opts.shape) || ...
   ~any(strcmpi(opts.shape(:)', {'boynton', 'twogamma'}))
    opts.shape = 'twogamma';
else
    opts.shape = lower(opts.shape(:)');
end
if ~isfield(opts, 'tr') || ~isa(opts.tr, 'double') || numel(opts.tr) ~= 1 || ...
    isinf(opts.tr) || isnan(opts.tr) || opts.tr <= 0
    opts.tr = 2000;
end
if ~isfield(opts, 'sbins') || ~isa(opts.sbins, 'double') || numel(opts.sbins) ~= 1 || ...
    isinf(opts.sbins) || isnan(opts.sbins) || opts.sbins < 1
    opts.sbins = 30;
else
    opts.sbins = round(opts.sbins);
end
if ~isfield(opts, 'sbinta') || ~isa(opts.sbinta, 'double') || numel(opts.sbinta) ~= 1 || ...
    isinf(opts.sbinta) || isnan(opts.sbinta) || opts.sbinta < 1 || opts.sbinta > opts.tr
    opts.sbinta = opts.tr;
end
if ~isfield(opts, 'sbinte') || ~isa(opts.sbinte, 'double') || numel(opts.sbinte) ~= 1 || ...
    isinf(opts.sbinte) || isnan(opts.sbinte) || opts.sbinte < 1 || opts.sbinte > (opts.tr / opts.sbinta)
    opts.sbinte = max(1, floor(opts.tr / opts.sbins));
end
if ~isfield(opts, 'srbin') || ~isa(opts.srbin, 'double') || numel(opts.srbin) ~= 1 || ...
    isinf(opts.srbin) || isnan(opts.srbin) || opts.srbin < 1 || opts.srbin > opts.sbins
    opts.srbin = 1;
end

% get content
bc = xo.C;

% last predictor to transcode
lp = bc.FirstConfoundPredictor - 1;
lsdm = size(bc.SDMMatrix, 1);
maxo = (lsdm - 1) * opts.tr;

% generate output
conds = repmat(struct( ...
    'ConditionName', {{'Condition'}}, ...
    'NrOfOnOffsets', 0, ...
    'OnOffsets',     zeros(0, 2), ...
    'Weights',       [], ...
    'Color',         [0, 0, 0]), 1, lp + double(opts.restcond));

% rest condition
if opts.restcond
    conds(end).ConditionName{1} = 'Rest';
    conds(end).Color = [255, 255, 255];
end

% nothing else to do
if lp < 1
    if opts.restcond
        conds(end).NrOfOnOffsets = 1;
        conds(end).OnOffsets = [0, opts.tr * lsdm];
    end
    return;
end

% get design matrix right
sdm = bc.SDMMatrix(:, 1:lp);
pnames = bc.PredictorNames(1:lp);
pcols = bc.PredictorColors(1:lp, :);

% get HRF as specified
h = ne_methods.hrf(opts.shape, 0.001, opts.pttp, opts.nttp, opts.pnr, opts.ons, ...
    opts.pdsp, opts.ndsp, opts.s, 0);

% generate figure
f = figure;
set(f, 'Name', 'SDM to PRT onset estimation');
ax1 = subplot(2, 1, 1);
ax2 = subplot(2, 1, 2);
vp = plot(ax1, [0, 0]);
sp = plot(ax2, (0:opts.tr:maxo), sdm(:, [1, 1]));
xlabel(ax1, 'Iteration');
ylabel(ax1, '-log(1 - % variance explained)');
xlabel(ax2, 'Time (ms)');
ylabel(ax2, 'HRF convolution');

% cache
coc = cell(10000, 1);

% iterate over conditions
for cc = 1:lp

    % set name and color
    conds(cc).ConditionName = pnames(cc);
    conds(cc).Color = pcols(cc, :);

    % title figure
    title(ax1, [pnames{cc}, ' onset variance explained by search']);
    title(ax2, [pnames{cc}, ' original and estimated convolution']);

    % empty regressor
    svec = sdm(:, cc);
    set(sp(1), 'YData', svec);
    set(sp(2), 'YData', svec);
    if all(svec == svec(1))
        continue;
    end

    % find maximum peak (in ms)
    [mv, mp] = max(svec);

    % and continue, while any peaks found
    op = find(svec > 0.75 * mv);
    o = zeros(0, 2);
    tvec = svec;
    tvecs = sum(tvec .* tvec);
    while ~isempty(op)

        % compute onset
        omp = round(opts.tr * (mp - 1) - 1000 * (opts.pttp + tvec(mp, 1)));

        % and duration
        dur = 2500 .* (tvec(mp, 1) + log(1 + tvec(mp, 1)) .^ 4);
        ons = max(0, round([omp - dur, omp + dur]));
        o(end+1, :) = ons;

        % the convolve
        [coc, ohc] = convonsets(coc, ...
            h, ons, lsdm, opts.tr, opts.srbin, opts.sbins, opts.sbinta, opts.sbinte, convones);

        % and subtract from the curve
        tvec = tvec - ohc;
        tvec(ohc > (0.5 * max(ohc))) = 0;

        % copy over
        mp = ne_methods.maxpos(tvec);
        op = find(tvec > 0.75 * mv);
    end

    % assume this one onset, and then take it from there
    ot = repmat({o}, 257, 1);

    % apply some randomness
    ot(2:end, 1) = grandom(ot(2:end, 1), maxo, opts.regspaced, 0);

    % then convolve each of them
    coa = ot;
    for oc = 1:numel(ot)
        [coc, coa{oc}] = convonsets(coc, ...
            h, ot{oc}, lsdm, opts.tr, opts.srbin, opts.sbins, opts.sbinta, opts.sbinte, convones);
    end

    % now genetically run this
    maxcorr = 0;
    maxiter = 8;
    tmaxiter = 1000;
    maxcorrp = NaN .* zeros(tmaxiter, 1);
    if opts.measure(1) == 'c'
        svecm = repmat(svec', 257, 1);
    else
        svecm = repmat(svec, 1, 257);
    end
    while maxiter > 0 && tmaxiter > 0

        % compute correlations
        switch opts.measure
            case 'corr'
                [cv, cr] = ne_methods.cov_nd(svecm, cat(2, coa{:})');
            case 'sdiff'
                cr = 1 - lsqz(sum(abs(svecm - cat(2, coa{:})) .^ 2)) ./ tvecs;
        end
        maxlast = maxcorr;
        maxcorr = max(cr);
        maxcorrp(1001 - tmaxiter) = -log(1 - maxcorr * maxcorr);

        % keep top 100, multiply by 10
        [cr, cv] = sort(cr);
        ot = [repmat(ot(cv(end:-1:end-15)), 16, 1); ot(end)];
        coa = [repmat(coa(cv(end:-1:end-15)), 16, 1); ot(end)];
        if tmaxiter < 1000
            set(vp, 'XData', 1:(1001-tmaxiter), 'YData', maxcorrp(1:1001-tmaxiter));
            set(sp(2), 'YData', coa{1});
            drawnow;
        end

        % then reiterate
        ot(2:end-1, 1) = grandom(ot(2:end-1, 1), maxo, opts.regspaced, maxcorr);
        ose = true;
        for oc = 2:16
            if numel(ot{oc}) ~= numel(ot{1})
                ose = false;
                break;
            end
        end
        if ose
            ot{end} = round(mean(cat(3, ot{1:16}), 3));
        end
        for oc = 17:numel(ot)
            [coc, coa{oc}] = convonsets(coc, ...
                h, ot{oc}, lsdm, opts.tr, opts.srbin, opts.sbins, opts.sbinta, opts.sbinte, convones);
        end

        % and keep track of iteration
        if maxcorr > maxlast
            maxiter = 8;
        else
            maxiter = maxiter - 1;
            tmaxiter = tmaxiter + 1;
        end
        tmaxiter = tmaxiter - 1;
    end

    % set on-offsets
    conds(cc).NrOfOnOffsets = size(ot{1}, 1);
    conds(cc).OnOffsets = ot{1};
    conds(cc).Weights = zeros(size(ot{1}, 1), 0);
end

% delete figure
delete(f);
convonsets(0);


% sub-function to rapidly concolve on/offsets with hrf and resample
function [coc, rch] = convonsets(coc, h, o, l, prtr, srbin, sbins, sbinta, sbinte, convones)

% persistent zt
persistent nsdmsp sdmsp zt;

% neuroelf library
global ne_methods;

if isempty(zt)
    zt = zeros(prtr * l, 1);

    % sampling points
    sdmsp = round(1 + (srbin - 1) * (sbinta / sbins));
    sdmsp = sdmsp:(sdmsp + sbinte - 1);
    nsdmsp = numel(sdmsp);
    sdmsp = sdmsp(:) * ones(1, l) + ones(nsdmsp, 1) * (0:prtr:(l-1)*prtr);
elseif nargin < 2
    zt = [];
    o = [];
end

% nothing to do
if isempty(o)
    rch = zeros(l, 1);
    return;
end
zt = 0 .* zt;
nzt = numel(zt);

% generate train of zeros
o = o + 1;

% compute lengths
ol = max(1, diff(o, 1, 2));

% iterate
for oc = 1:size(o, 1)

    % convolve
    if ol(oc) <= 10000 && ~isempty(coc{ol(oc)})
        co = coc{ol(oc)};
    else
        co = ne_methods.convones(h, ol(oc));
        if ol(oc) <= 10000
            coc{ol(oc)} = co;
        end
    end

    % then store into array
    tz = o(oc, 1) + numel(co) - 1;
    if tz > nzt
        zt(o(oc, 1):end) = zt(o(oc, 1):end) + co(1:nzt-(o(oc, 1)-1));
    else
        zt(o(oc, 1):tz) = zt(o(oc, 1):tz) + co;
    end
end

% re-sample
rch = lsqz(mean(reshape(zt(sdmsp(:)), nsdmsp, l), 1));


% randomly jitter onsets
function ot = grandom(ot, maxo, regspaced, mc)

% randomization factor
rf = 500 * sqrt(1 - mc * mc);

% random values
if rf > 100
    rv = ceil((9 + double(regspaced)) .* rand(numel(ot), 1));
else
    rv = 3 + ceil((6 + double(regspaced)) .* rand(numel(ot), 1));
end

% iterate
for oc = 1:numel(ot)

    % current onset list
    o = ot{oc};

    % what to do
    switch rv(oc)

        % add another onset
        case 1
            o(end+1, :) = sort(ceil(maxo .* rand(2, 1)))';

        % duplicate one onset
        case 2
            so = ceil(size(o, 1) * rand(1, 1));
            o(end+1, :) = min(maxo, max(0, o(so, :) + round(20000 * randn(1, 1))));

        % remove random onset
        case 3
            if size(o, 1) > 1
                o(ceil(size(o, 1) * rand(1, 1)), :) = [];
            end

        % change onset time of one onset
        case 4
            so = ceil(size(o, 1) * rand(1, 1));
            o(so, 1) = max(0, min(o(so, 1) + round(rf * randn(1, 1)), o(so, 2) - 1));

        % change offset time of one onset
        case 5
            so = ceil(size(o, 1) * rand(1, 1));
            o(so, 2) = min(maxo, max(o(so, 1) + 1, o(so, 2) + round(rf * randn(1, 1))));

        % change on- and offset of one onset
        case 6
            so = ceil(size(o, 1) * rand(1, 1));
            o(so, :) = min(maxo, max(0, o(so, :) + round(rf * randn(1, 1))));

        % change onset time of all onsets
        case 7
            o(:, 1) = max(0, min(o(:, 1) + round(rf * randn(1, 1)), o(:, 2) - 1));

        % change offset time of all onsets
        case 8
            o(:, 2) = min(maxo, max(o(:, 1) + 1, o(:, 2) + round(rf * randn(1, 1))));

        % change on- and offset of all onset
        case 9
            o = min(maxo, max(0, o + round(rf * randn(1, 1))));

        % harmonize lengths
        case 10
            if size(o, 1) > 1
                o(:, 2) = o(:, 1) + round(mean(diff(o, 1, 2)) + randn(1, 1));
            end
    end

    % sort, store
    o = sortrows(o);

    % make sure nothing is double
    if size(o, 1) > 1
        o(o([2:end, end], 1) == o([1:end-1, 1], 1), :) = [];
    end

    % don't allow emptys
    if isempty(o)
        o = sort(ceil(maxo .* rand(2, 1)))';
    end

    % store
    ot{oc} = o;
end


% local implementation of lsqueeze
function v = lsqz(v)
v = v(:);

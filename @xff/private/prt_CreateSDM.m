function [xo2, hf] = prt_CreateSDM(xo, params)
% PRT::CreateSDM  - create an SDM from a PRT
%
% FORMAT:       [sdm, bfs] = prt.CreateSDM(params);
%
% Input fields:
%
%       params      1x1 struct with optional fields
%        .collapse  Cx2 or Cx3 cell array with collapsing arguments
%        .diffreg   Dx2 cell array with regressors to subtract
%        .erlen     duration of event-related (0-length) stimuli (100ms)
%        .hshape    HRF shape to use ({'twogamma'}, 'boynton')
%        .hpttp     HRF time to peak for positive response (5)
%        .hnttp     HRF time to peak for negative response (15)
%        .hpnr      HRF positive/negative ratio (6)
%        .hons      HRF onset (0)
%        .hpdsp     HRF positive response dispersion (1)
%        .hndsp     HRF negative response dispersion (1)
%        .ndcreg    deconvolution number of lags (regressors, 8)
%        .nderiv    add derivatives to HRF (1xN double, [])
%        .noderiv   Cx1 cell array with conditions without derivatives
%        .nvol      number of volumes (data points, 200)
%        .ortho     if given and true, orthogonalize derivatives
%        .params    1xN struct for parametrical designs (empty)
%          .cond    condition name or number
%          .name    parameter name
%          .opts    optional settings for parameter
%            .comp  compress parameter, 'log', 'sqrt'
%            .norm  z-normalization of parameter (otherwise as is)
%            .ortho normalize convolved regressor against design
%          .pval    parameter values (must match number of onsets)
%        .pnorm     normalize parameters from Cond(x).Weights (true)
%        .portho    orthogonalize params from Cond(x).Weights (false)
%        .ppicond   list of condition(s) to interact ppitc with
%        .ppitc     time-course from VOI/phys data (in TR resolution)
%        .ppitf     1x2 cell array with temporal filtering for .ppitc
%        .prtr      TR (in ms, 2000)
%        .rcond     Conditions to remove (rest, 1)
%        .regi      VxR regressors of interest or SDM (not orthogonalized)
%        .regni     VxR regressors of no interest or SDM (orthogonalized)
%        .rnorm     re-normalization of maximum (default: set plateau to 1)
%        .sbins     number of slices (default: 30, only if srbin > 1)
%        .sbinta    acquisition time for all slices (default: prtr)
%        .sbinte    time of echo (duration to acquire one slice, prtr/sbins)
%        .srbin     reference bin (ref slice in slice-time correction: 1)
%        .tshift    Temporally shift convoluted regressor (in ms, 0)
%        .type      HRF or FIR model ({'hrf'}, 'fir')
%
% Output fields:
%
%       sdm         SDM object (design matrix)
%       bfs         basis function set (tbin-by-functions)
%
% Note: parametrical weights in Cond(c).Weights will be *ADDED* to the
%       parameters after the evaluation of .param!
%
% Using: checkstruct, convones, ddeblank, emptystruct, findfirst,
%        flexinterpn_method, hrf, meannoinfnan, orthvec, splittocell,
%        tempfilter.

% Version:  v1.1
% Build:    16041209
% Date:     Apr-12 2016, 9:49 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010 - 2012, 2014, 2016, Jochen Weber
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

% import from neuroelf library into workspace
global ne_methods;
convones    = ne_methods.convones;
ddeblank    = ne_methods.ddeblank;
findfirst   = ne_methods.findfirst;
hrf         = ne_methods.hrf;
orthvec     = ne_methods.orthvec;
splittocell = ne_methods.splittocell;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
sbf = xo.F;
if isempty(sbf)
    sbf = 'interactive.prt';
else
    [sbfp, sbf, sbfe] = fileparts(sbf);
    sbf = [sbf sbfe];
end
bc = xo.C;
rtv = bc.RunTimeVars;

% default parameters
defpar = { ...
    'collapse','cell',   'nonempty', {}; ...          % collapsing conds
    'diffreg', 'cell',   'nonempty', {}; ...          % difference of regs
    'erlen',   'double', 'nonempty', 100; ...         % value if duration=0
    'hshape',  'char',   'nonempty', 'twogamma'; ...  % HRF type
    'hpttp',   'double', 'nonempty', 5; ...           % pos. time-to-peak
    'hnttp',   'double', 'nonempty', 15; ...          % neg. time-to-peak
    'hpnr',    'double', 'nonempty', 6; ...           % pos/neg ratio
    'hons',    'double', 'nonempty', 0; ...           % HRF onset
    'hpdsp',   'double', 'nonempty', 1; ...           % pos. dispersion
    'hndsp',   'double', 'nonempty', 1; ...           % neg. dispersion
    'moco',    'double', 'noinfnan', []; ...          % motion correction
    'ndcreg',  'double', 'nonempty', 8; ...           % # of deconv. lags
    'nderiv',  'double', 'noinfnan', []; ...          % add derivatives
    'noderiv', 'cell',   'nonempty', {}; ...          % noderiv conditions
    'nvol',    'double', 'nonempty', 200; ...         % number of volumes
    'ortho',   'logical','nonempty', false; ...       % orthogonalization
    'pnorm',   'logical','nonempty', true; ...        % param normalization
    'portho',  'logical','nonempty', false; ...       % param orthogon.
    'ppicond', 'cell',   'nonempty', {}; ...          % PPI condition list
    'ppitc',   'double', 'nonempty', []; ...          % PPI time course
    'ppitf',   'cell',   'nonempty', {}; ...          % PPI temporal filter
    'prtr',    'double', 'nonempty', 2000; ...        % TR (ms)
    'rcond',   'double', '',         1; ...           % remove which cond.
    'regi',    'double', 'nonempty', []; ...          % regs of interest
    'regni',   'double', 'nonempty', []; ...          % regs of no interest
    'rnorm',   'double', 'nonempty', 0; ...           % normalization
    'sbins',   'double', 'nonempty', 30; ...          % number of slices
    'sbinta',  'double', 'nonempty', []; ...          % time-of-acquisition
    'sbinte',  'double', 'nonempty', []; ...          % echo-time per slice
    'srbin',   'double', 'nonempty', 1; ...           % reference slice/bin
    'tshift',  'double', 'nonempty', 0; ...           % temporal shift
    'type',    'char',   'nonempty', 'hrf' ...        % model type hrf/fir
};

% check optional paramters
if nargin > 1 && isstruct(params) && numel(params) == 1
    params = ne_methods.checkstruct(params, defpar);
else
    params = ne_methods.checkstruct(struct, defpar);
end
if ~isempty(params.collapse) && any(size(params.collapse, 2) == [2, 3])
    for cc = 1:size(params.collapse, 1)
        try
            prt_Collapse(xo, params.collapse{cc, :});
        catch xfferror
            error('neuroelf:xff:badArgument', ...
                'Invalid condition collapsing argument %d.', cc);
        end
    end
    bc = xo.C;
end
if isinf(params.ndcreg(1)) || isnan(params.ndcreg(1)) || ...
    params.ndcreg(1) < 1 || params.ndcreg(1) > 120
    ndcreg = [];
else
    ndcreg = round(params.ndcreg(1));
end
nderiv = unique(params.nderiv(:)');
nderiv(nderiv ~= fix(nderiv) | nderiv < 1 | nderiv > 2) = [];
ndnum = numel(nderiv);
if ~isempty(params.noderiv)
    params.noderiv = params.noderiv(:);
    for cc = numel(params.noderiv):-1:1
        if ~ischar(params.noderiv{cc}) || isempty(params.noderiv{cc})
            params.noderiv(cc) = [];
        end
    end
end
if isinf(params.nvol(1)) || isnan(params.nvol(1)) || ...
    params.nvol(1) < 1 || params.nvol(1) > 10000
    nvol = 200;
else
    nvol = round(params.nvol(1));
end
if numel(params.pnorm) > 1
    params.pnorm = params.pnorm(1);
end
if numel(params.portho) > 1
    params.portho = params.portho(1);
end
if numel(params.ppitc) ~= nvol
    params.ppitc = [];
end
if numel(params.ppitf) ~= 2 || ~ischar(params.ppitf{1}) || ...
   ~any(strcmp(params.ppitf{1}, {'tempdct', 'tempsc'})) || ...
   ~isa(params.ppitf{2}, 'double') || numel(params.ppitf{2}) ~= 1
    params.ppitf = {};
end
ortho = params.ortho(1);
olist = struct;
if isinf(params.prtr(1)) || isnan(params.prtr(1)) || ...
    params.prtr(1) < 10 || params.prtr(1) > 30000
    if ~isinf(params.prtr(1)) && ~isnan(params.prtr(1)) && ...
        params.prtr(1) >= 0.25 && params.prtr(1) <= 5
        prtr = round(1000 * params.prtr(1));
    else
        warning('neuroelf:xff:replacingWithDefault', ...
            'Invalid or missing TR value; replacing with 2000ms.');
        prtr = 2000;
    end
else
    prtr = round(params.prtr(1));
end
if isempty(ndcreg)
    ndcreg = ceil(24000 / prtr);
end
if ~isempty(params.regi) && size(params.regi, 1) == nvol && ndims(params.regi) == 2
    xregi = params.regi;
    xregi(:, any(isinf(xregi) | isnan(xregi))) = [];
else
    xregi = [];
end
nxregi = size(xregi, 2);
if ~isempty(params.regni) && size(params.regni, 1) == nvol && ndims(params.regni) == 2
    xregni = params.regni;
    xregni(:, any(isinf(xregni) | isnan(xregni))) = [];
else
    xregni = [];
end
nxregni = size(xregni, 2);
ptype = params.type(:)';
if ~strcmpi(ptype, 'hrf')
    ptype = false;
else
    ptype = true;
end
rcond = params.rcond(:)';
rcond(isinf(rcond) | isnan(rcond) | rcond < 1 | rcond ~= fix(rcond)) = [];
if isinf(params.rnorm(1)) || isnan(params.rnorm(1))
    params.rnorm = 1;
else
    params.rnorm = params.rnorm(1);
end
if isinf(params.sbins(1)) || isnan(params.sbins(1)) || ...
    params.sbins(1) < 1 || params.sbins(1) > 256
    params.sbins = 30;
else
    params.sbins = round(params.sbins(1));
end
if isempty(params.sbinta) || isinf(params.sbinta(1)) || isnan(params.sbinta(1))
    params.sbinta = params.prtr;
else
    params.sbinta = min(params.prtr, max(params.sbins, round(params.sbinta(1))));
end
if isempty(params.sbinte) || isinf(params.sbinte(1)) || isnan(params.sbinte(1)) || ...
    params.sbinte(1) < 1 || params.sbinte(1) > (params.prtr / params.sbins)
    params.sbinte = round(params.prtr / params.sbins);
else
    params.sbinte = round(params.sbinte(1));
end
if isinf(params.srbin(1)) || isnan(params.srbin(1)) || ...
    params.srbin(1) < 1 || params.srbin(1) > params.sbins
    params.srbin = 1;
else
    params.srbin = round(params.srbin(1));
end
if isinf(params.tshift(1)) || isnan(params.tshift(1)) || (numel(params.tshift) > 1 && ...
    (isinf(params.tshift(2)) || isnan(params.tshift(2)) || params.tshift(2) == 0))
    tshift = [0, 1];
elseif numel(params.tshift) == 1
    tshift = [round(params.tshift(1)), 1];
else
    tshift = [round(params.tshift(1)), params.tshift(2)];
end

% create empty SDM
xo2 = xff('new:sdm');
bc2 = xo2.C;

% get shape of HRF (and derivatives)
if ptype
    if ~isempty(nderiv)
        mxd = max(nderiv);
    else
        mxd = 0;
    end
    hf = hrf(params.hshape, 1e-3, params.hpttp, params.hnttp, ...
        params.hpnr, params.hons, params.hpdsp, params.hndsp, ...
        -params.hons:1e-3:(2 * (params.hnttp + 1)), mxd);
else
    hf = eye(params.ndcreg);
end

% check protocol for volumes
usevol = 0;
if ~isempty(bc.ResolutionOfTime) && lower(bc.ResolutionOfTime(1)) == 'v'
    usevol = prtr;
end

% get conditions
cond = bc.Cond;
icond = 1:numel(cond);

% remove indices and conditions
rcond(rcond > numel(cond)) = [];
cond(rcond) = [];
icond(rcond) = [];
ncond = numel(cond);

% check condition parameters struct
if ~isfield(params, 'params') || ~isstruct(params.params) || isempty(params.params) || ...
   ~isfield(params.params, 'cond') || ~isfield(params.params, 'name') || ~isfield(params.params, 'pval')
    params.params = ne_methods.emptystruct({'cond', 'name', 'opts', 'pval'});
end
params.params = params.params(:)';

% add parameters for weights
if bc.ParametricWeights > 0

    % get names
    if isfield(rtv, 'ParameterNames') && iscell(rtv.ParameterNames) && ...
        numel(rtv.ParameterNames) == bc.ParametricWeights && ...
        all(cellfun(@ischar, rtv.ParameterNames(:))) && ...
       ~any(cellfun('isempty', rtv.ParameterNames(:)))
        pnames = rtv.ParameterNames(:)';
        for pcc = 1:bc.ParametricWeights
            pnames{pcc} = pnames{pcc}(:)';
        end
    else
        pnames = cell(1, bc.ParametricWeights);
        for pcc = 1:bc.ParametricWeights
            pnames{pcc} = sprintf('p%d', pcc);
        end
    end
    for cc = 1:numel(icond)
        cci = icond(cc);
        if size(bc.Cond(cci).Weights, 1) == size(bc.Cond(cci).OnOffsets, 1) && ...
           ~isempty(bc.Cond(cci).Weights)
            for pcc = 1:size(bc.Cond(cci).Weights, 2)
                if any(diff(bc.Cond(cci).Weights(:, pcc)) ~= 0)
                    params.params(end+1) = struct( ...
                        'cond', cci, 'name', pnames{pcc}, ...
                        'opts', struct('comp', 'no', ...
                            'norm', params.pnorm, 'ortho', params.portho), ...
                        'pval', bc.Cond(cci).Weights(:, pcc));
                end
            end
        end
    end
end

% for deconvolution, discard parameters
if ~ptype
    params.params(:) = [];
end

% check params settings
if ~isfield(params.params, 'opts')
    params.params(1).opts = [];
end
for pc = 1:numel(params.params)
    if ~isstruct(params.params(pc).opts)
        params.params(pc).opts = struct('comp', 'no', 'norm', true, 'ortho', false);
    end
    if ~isfield(params.params(pc).opts, 'comp') || ~ischar(params.params(pc).opts.comp) || ...
       ~any(strcmpi(params.params(pc).opts.comp(:)', {'no', 'sqrt', 'log'}))
        params.params(pc).opts.comp = 'no';
    else
        params.params(pc).opts.comp = lower(params.params(pc).opts.comp(:)');
    end
    if ~isfield(params.params(pc).opts, 'norm') || ~islogical(params.params(pc).opts.norm) || ...
        numel(params.params(pc).opts.norm) ~= 1
        params.params(pc).opts.norm = true;
    end
    if ~isfield(params.params(pc).opts, 'ortho') || ~islogical(params.params(pc).opts.ortho) || ...
        numel(params.params(pc).opts.ortho) ~= 1
        params.params(pc).opts.ortho = true;
    end
end
cnames = cell(numel(icond), 1);
for cc = 1:numel(icond)
    cnames{cc} = bc.Cond(icond(cc)).ConditionName{1};
end
for pc = numel(params.params):-1:1
    p = params.params(pc);
    if isempty(p.cond) || (~ischar(p.cond) && (~isa(p.cond, 'double') || ...
         numel(p.cond) ~= 1 || isinf(p.cond) || isnan(p.cond) || ...
         p.cond < 1 || p.cond > numel(bc.Cond) || p.cond ~= fix(p.cond))) || ...
        isempty(p.name) || ~ischar(p.name) || isempty(p.pval) || ...
       ~isa(p.pval, 'double') || numel(p.pval) < 2 || any(isinf(p.pval(:)))
        warning('neuroelf:xff:badArgument', 'Invalid parameter %d.', pc);
        params.params(pc) = [];
        continue;
    end
    if any(isnan(p.pval(:)))
        if ischar(p.cond)
            fprintf(['Warning: parameter %d (condition %s in file %s) has %d/%d ' ...
                'missing values, replacing with mean.\n'], pc, p.cond, sbf, ...
                sum(isnan(p.pval(:))), numel(p.pval));
        else
            fprintf(['Warning: parameter %d (condition %d in file %s) has %d/%d ' ...
                'missing values, replacing with mean.\n'], pc, p.cond, sbf, ...
                sum(isnan(p.pval(:))), numel(p.pval));
        end
        p.pval(isnan(p.pval(:))) = ne_methods.meannoinfnan(p.pval(:));
    end
    if numel(unique(p.pval(:))) < 2
        if ischar(p.cond)
            fprintf(['Warning: invalid parameter %d (condition %s in file %s) ' ...
                'with only one unique value.\n'], pc, p.cond, sbf);
        else
            fprintf(['Warning: invalid parameter %d (condition %d in file %s) ' ...
                'with only one unique value.\n'], pc, p.cond, sbf);
        end
        params.params(pc) = [];
        continue;
    end
    if ischar(p.cond)
        cf = findfirst(strcmpi(p.cond(:)', cnames));
        if isempty(cf)
            fprintf('Warning: condition %s for parameter %d not found.\n', p.cond, pc);
            params.params(pc) = [];
            continue;
        end
        p.cond = icond(cf(1));
    end
    if ~any(icond == p.cond)
        params.params(pc) = [];
        continue;
    end
    pvs = size(p.pval);
    npv = size(bc.Cond(p.cond).OnOffsets, 1);
    svf = find(pvs == npv);
    if isempty(svf)
        warning('neuroelf:xff:badArgument', ...
            'Invalid number of parameter values for condition %s.', cnames{p.cond});
        params.params(pc) = [];
    end
    if svf(1) > 1
        p.pval = reshape(permute( ...
            p.pval, [svf(1), setdiff(1:ndims(p.pval), svf(1))]), [npv, numel(p.pval) / npv]);
    end
    if size(p.pval, 2) > 3
        warning('neuroelf:xff:badArgument', ...
            'Only up to three sub-parameters with interaction supported.');
        params.params(pc) = [];
        continue;
    elseif size(p.pval, 2) == 3
        p.pval(:, 7) = p.pval(:, 1) .* p.pval(:, 2) .* p.pval(:, 3);
        p.pval(:, 4) = p.pval(:, 1) .* p.pval(:, 2);
        p.pval(:, 5) = p.pval(:, 1) .* p.pval(:, 3);
        p.pval(:, 6) = p.pval(:, 2) .* p.pval(:, 3);
    elseif size(p.pval, 2) == 2
        p.pval(:, 3) = p.pval(:, 1) .* p.pval(:, 2);
    end
    if ~strcmp(p.opts.comp, 'no')
        pvs = sign(p.pval);
        p.pval = abs(p.pval);
        if strcmp(p.opts.comp, 'log')
            if any(p.pval(:) < 1)
                p.pval = p.pval + 1;
            end
            p.pval = pvs .* log(p.pval);
        else
            p.pval = pvs .* sqrt(p.pval);
        end
    end
    if p.opts.norm
        p.pval = p.pval - repmat(mean(p.pval), [size(p.pval, 1), 1]);
        p.pval = p.pval ./ repmat(std(p.pval), [size(p.pval, 1), 1]);
    end
    if size(p.pval, 2) > 1
        rep = repmat(p, [1, size(p.pval, 2)]);
        for rpc = 1:numel(rep)
            rep(rpc).pval = rep.(rpc).pval(:, rpc);
            if numel(rep) == 3
                if rpc < 3
                    rep(rpc).name = sprintf('%s - p%d', rep(rpc).name, rpc);
                else
                    rep(rpc).name = sprintf('%s - p1Xp2', rep(rpc).name);
                end
            else
                switch (rpc)
                    case {1, 2, 3}
                        rep(rpc).name = sprintf('%s - p%d', rep(rpc).name, rpc);
                    case {4}
                        rep(rpc).name = sprintf('%s - p1Xp2', rep(rpc).name);
                    case {5}
                        rep(rpc).name = sprintf('%s - p1Xp3', rep(rpc).name);
                    case {6}
                        rep(rpc).name = sprintf('%s - p2Xp3', rep(rpc).name);
                    case {7}
                        rep(rpc).name = sprintf('%s - p1Xp2Xp3', rep(rpc).name);
                end
            end
        end
        p = rep;
    end
    if numel(p) == 1
        params.params(pc) = p;
    else
        params.params = [params.params(1:(pc - 1)), p, params.params((pc + 1):end)];
    end
end

% get correct number of parameters
numparams = numel(params.params);
parcond = zeros(size(params.params));
for pc = 1:numparams
    parcond(pc) = find(icond == params.params(pc).cond);
end

% generate one regressor per predictor in MS resolution
nsamp = nvol * prtr;
if ptype
    ncols = 1 + ndnum;
    nrows = nsamp;
else
    ncols = ndcreg;
    nrows = nvol;
end

% reserve memory
try
    sdm = zeros(nrows, (ncond + numparams) * ncols + 1);

% if fails...
catch xfferror

    % check if memory was the problem
    l = lower(xfferror.message);
    if ~isempty(strfind(l, 'memory'))

        % if so, create the SDM column by column
        xsnm = cell(1, (ncond + numparams) * ncols + 1 + nxregi + nxregni);
        xscl = zeros(numel(xsnm), 3);
        xsdm = zeros(nvol, numel(xsnm));
        nxtc = 1;
        nparams = params;
        for fbc = 1:numel(icond)
            nparams.rcond = 1:numel(bc.Cond);
            nparams.rcond(icond(fbc)) = [];
            nsdm = prt_CreateSDM(xo, nparams);
            nsdmc = nsdm.C;
            delete(nsdm);
            xsnm(nxtc:nxtc+size(nsdmc.SDMMatrix, 2)-2) = nsdmc.PredictorNames(1:end-1);
            xscl(nxtc:nxtc+size(nsdmc.SDMMatrix, 2)-2, :) = nsdmc.PredictorColors(1:end-1, :);
            xsdm(:, nxtc:nxtc+size(nsdmc.SDMMatrix, 2)-2) = nsdmc.SDMMatrix(:, 1:end-1);
            nxtc = nxtc + size(nsdmc.SDMMatrix, 2) - 1;
        end
        xsnm(end) = nsdmc.PredictorNames(end);
        xscl(end, :) = nsdmc.PredictorColors(end, :);
        xsdm(:, end) = nsdmc.SDMMatrix(:, end);

        % account for xregi/xregni
        if nxregi > 0
            xsdm(:, end-(nxregi+nxregni):end-(1+nxregni)) = xregi;
            xscl(end-(nxregi+nxregni):end-(1+nxregni), :) = floor(255.999 * rand(nxregi, 3));
            for xrc = 1:nxregi
                xsnm{end+xrc-(nxregi+nxregni)} = sprintf('regi_%d', xrc);
            end
        end
        if nxregni > 0
            xsdm(:, end-nxregni:end-1) = xregni;
            xscl(end-nxregni:end-1, :) = 255;
            for xrc = 1:nxregni
                xsnm{end+xrc-nxregni} = sprintf('regni_%d', xrc);
            end
        end

        % put into new object
        bc2.NrOfPredictors = size(xsdm, 2);
        bc2.NrOfDataPoints = nvol;
        bc2.IncludesConstant = 1;
        bc2.FirstConfoundPredictor = size(xsdm, 2) - nxregni;
        bc2.PredictorColors = xscl;
        bc2.PredictorNames = xsnm;
        bc2.SDMMatrix = xsdm;
        bc2.RTCMatrix = xsdm(:, 1:end-(1+nxregni));
        xo2.C = bc2;

        % return
        return;
    end

    % otherwise rethrow...
    rethrow(xfferror);
end

% go on
sdm(:, end) = 1;
sdmn = cell(1, (ncond + numparams) * ncols + 1);
sdmn{end} = 'Constant';
sdml = zeros(size(sdm, 2), 3);
sdml(end, :) = 255;

% compute first columns per condition / parameter
colcond = ones(1, ncond + 1);
for cc = 2:(ncond + 1)
    colcond(cc) = colcond(cc - 1) + (1 + sum(parcond == (cc - 1))) * ncols;
end

% check size
if colcond(end) ~= size(sdm, 2)
    error('neuroelf:xff:internalError', 'xff implementation error.');
end

% create list of indices for conditions and p

% init SDM
bc2.NrOfPredictors = size(sdm, 2);
bc2.NrOfDataPoints = nvol;
bc2.IncludesConstant = 1;

% create empty cache
ccache = struct;

% iterate over conditions
for cc = 1:ncond

    % get first target column
    tcc = colcond(cc);
    sdml(tcc, :) = cond(cc).Color;

    % any parameters
    pcf = params.params(parcond == cc);

    % get condition on-/offsets and possible parameters
    oo = cond(cc).OnOffsets;
    if ptype
        for dc = [0, nderiv]
            if dc == 0
                sdmn{tcc} = cond(cc).ConditionName{1};
            else
                sdml(tcc, :) = min(255, max(0, ...
                    round((dc * [128, 128, 128] + cond(cc).Color) ./ (dc + 1))));
                sdmn{tcc} = sprintf('%s - %d.deriv', cond(cc).ConditionName{1}, dc);
            end
            tcc = tcc + 1;
        end
        for pc = 1:numel(pcf)
            pcf(pc).sdmcol1 = tcc;
            for dc = [0, nderiv]
                if dc == 0
                    sdml(tcc, :) = min(255, max(0, ...
                        floor(64 .* rand(1, 3) + 0.75 .* cond(cc).Color)));
                    sdmn{tcc} = sprintf('%s x %s', cond(cc).ConditionName{1}, pcf(pc).name);
                else
                    sdml(tcc, :) = min(255, max(0, round(64 * rand(1, 3) + ...
                        (dc * [128, 128, 128] + cond(cc).Color) ./ (dc + 1))));
                    sdmn{tcc} = sprintf('%s x %s - %d.deriv', ...
                        cond(cc).ConditionName{1}, pcf(pc).name, dc);
                end
                tcc = tcc + 1;
            end
        end
    else
        namepart = cond(cc).ConditionName{1};
        for nc = 0:(ndcreg - 1)
            sdmn{nc + tcc} = sprintf('%s_D%d', namepart, nc);
        end
        tcc = tcc + ndcreg;
        for pc = 1:numel(pcf)
            for nc = 0:(ndcreg - 1)
                sdmn{nc + tcc} = sprintf('%sX%s_D%d', namepart, pcf(pc).name, nc);
            end
            tcc = tcc + ndcreg;
        end
    end

    % iterate over onsets
    for oc = 1:size(oo, 1)

        % get first target column again
        tcc = colcond(cc);

        % for volume-PRTs
        if usevol > 0
            oon = round(tshift(2) * (oo(oc, 1) - 1) * usevol + 1 + tshift(1));
            oof = round(tshift(2) * oo(oc, 2) * usevol + tshift(1));

        % for millisecond-PRTs
        else
            oon = round(1 + tshift(2) * oo(oc, 1) + tshift(1));
            oof = round(tshift(2) * oo(oc, 2) + tshift(1));

            % for pure event-related design (duration = 0ms), assume the
            % shortest possible "stick function" length to be of this
            % length (in resolution of SDM creation, which is ms for now)
            if oof == (oon - 1)
                oof = oof + params.erlen;
            end
        end

        % continue if out of range
        if oon > nsamp || oof < oon
            fprintf('Warning: onset %d (%d:%d) out of range in condition %d (%s), in file ''%s''.\n', ...
                oc, oo(oc, 1), oo(oc, 2), cc, cond(cc).ConditionName{1}, sbf);
            continue;
        end

        % for HRF
        if ptype

            % generate vector of ones
            ov = 1 + oof - oon;
            ccn = sprintf('c%d_0', ov);

            % already conv
            if ~isfield(ccache, ccn)
                for ndc = 0:mxd
                    cv = convones(hf(:, ndc + 1), ov);
                    ccache.(sprintf('c%d_%d', ov, ndc)) = cv(:);
                end
            end
            lcv = numel(ccache.(ccn));

            % get indices
            sif = 1;
            sit = lcv;
            if oon < 1
                sif = 2 - oon;
            end
            if (oon + sit - 1) > nsamp
                sit = nsamp + 1 - oon;
            end
            if oon < 1
                sif = sif + 1 - oon;
                oon = 1;
            end

            % add to predictor at onset
            if sit > 0

                % iterate over derivatives
                fcc = colcond(cc);
                for dc = [0, nderiv]

                    % get column value
                    cv = ccache.(sprintf('c%d_%d', ov, dc));

                    % put regressor into SDM
                    sdm(oon:(oon+sit-sif), tcc) = sdm(oon:(oon+sit-sif), tcc) + cv(sif:sit);

                    % increase target column
                    tcc = tcc + 1;
                end

                % parameters
                for pc = 1:numel(pcf)

                    % iterate over derivatives
                    for dc = [0, nderiv]

                        % no derivative
                        if dc == 0

                            % put HRF into SDM
                            sdm(oon:(oon+sit-sif), tcc) = ...
                                sdm(oon:(oon+sit-sif), tcc) + pcf(pc).pval(oc) * cv(sif:sit);

                            % keep track for orthogonalization
                            if pcf(pc).opts.ortho
                                olist.(sprintf('r%d', tcc)) = [tcc, (fcc:tcc-1)];
                            end

                        % derivative
                        else

                            % build derivative and scale to 1000 (match HRF)
                            dcv = [zeros(dc, 1); diff(cv, dc)];
                            dcv = 1000 * dcv;
                            sdm(oon:(oon+sit-sif), tcc) = ...
                                sdm(oon:(oon+sit-sif), tcc) + pcf(pc).pval(oc) * dcv(sif:sit);

                            % keep track for orthogonalization
                            if ortho || pcf(pc).opts.ortho
                                olist.(sprintf('r%d', tcc)) = [tcc, (fcc:tcc-1)];
                            end
                        end

                        % increase target column
                        tcc = tcc + 1;
                    end
                end
            end

        % deconvolution design (possible improvements to be made!)
        else
            for nc = 1:ndcreg
                don = round(oon / prtr) + nc;
                if don > nvol
                    break;
                end
                sdm(don, tcc) = 1;
                tcc = tcc + 1;
            end
            for pc = 1:numel(pcf)
                for nc = 1:ndcreg
                    don = round(oon / prtr) + nc;
                    if don > nvol
                        break;
                    end
                    sdm(don, tcc) = pcf(pc).pval(oc);
                    tcc = tcc + 1;
                end
            end
        end
    end
end

% PPI
if ~isempty(params.ppicond) && ~isempty(params.ppitc)

    % resample signal
    ppitc = ne_methods.flexinterpn_method( ...
        params.ppitc(:), [Inf;1;1/prtr;1+(size(sdm, 1) - 0.25)/prtr], 'cubic');

    % initialize marker
    mark = cell2struct(cell(numel(params.ppicond), 1, 2), {'pos', 'neg'}, 3);

    % initialize vector with conditions to sum
    sumcond = false(1, numel(sdmn));

    % iterate over condition list
    for cc = numel(mark):-1:1

        % set marker
        mark(cc).pos = false(1, numel(sdmn));
        mark(cc).neg = false(1, numel(sdmn));

        % get name
        cname = params.ppicond{cc};

        % contrast?
        if any(cname == '>')

            % split
            parts = splittocell(cname, '>');
            parts{1} = ddeblank(parts{1});
            parts{2} = ddeblank(parts{2});
        else
            parts = {cname};
        end

        % split parts even further
        if numel(parts{1}) > 2 && parts{1}(1) == '(' && parts{1}(end) == ')'
            parts{1} = parts{1}(2:end-1);
        end
        parts{1} = splittocell(parts{1}, '+');
        passed = true;
        for cpc = 1:numel(parts{1})
            parts{1}{cpc} = ddeblank(parts{1}{cpc});
            condi = strcmpi(parts{1}{cpc}, sdmn);
            if ~any(condi)
                passed = false;
                break;
            end
            mark(cc).pos = mark(cc).pos | condi(:)';
            sumcond = sumcond | condi(:)';
        end
        if ~passed
            mark(cc) = [];
            params.ppicond(cc) = [];
            continue;
        end
        if numel(parts) > 1
            if numel(parts{2}) > 2 && parts{2}(1) == '(' && parts{2}(end) == ')'
                parts{2} = parts{2}(2:end-1);
            end
            parts{2} = splittocell(parts{2}, '+');
            for cpc = 1:numel(parts{2})
                parts{2}{cpc} = ddeblank(parts{2}{cpc});
                condi = strcmpi(parts{2}{cpc}, sdmn);
                if ~any(condi)
                    passed = false;
                    break;
                end
                mark(cc).neg = mark(cc).neg | condi(:)';
                sumcond = sumcond | condi(:)';
            end
            if ~passed
                mark(cc) = [];
                params.ppicond(cc) = [];
            end
        end
    end

    % no more PPI
    if isempty(mark)
        params.ppicond = {};
        params.ppitc = [];

    % still valid
    else

        % regress explained variance out of ppitc
        ppitc = ne_methods.tempfilter(ppitc, struct('nuisreg', sdm, 'temp', true, ...
            params.ppitf{:}));

        % patch SDM (with CHEAP check)
        sumcond = find(sumcond);
        if numel(sumcond) < numel(mark)
            error('neuroelf:xff:badArgument', ...
                'PPI conditions/contrasts not independent.');
        end

        % increase size
        ppitc(numel(ppitc), 1 + numel(mark)) = 0;

        % ppi names
        ppin = cell(1, size(ppitc, 2));
        ppin{1} = 'VOI';
        ppil = floor(255.99 .* rand(1, 3));
        ppil(1 + numel(mark), 3) = 0;

        % for each condition/contrast
        for cc = 1:numel(mark)

            % get pos/neg columns in smaller matrix
            mpos = find(mark(cc).pos);
            mneg = find(mark(cc).neg);

            % build new regressor
            if ~isempty(mneg)
                sdmsum = (1 / numel(mpos)) .* sum(sdm(:, mpos), 2) - ...
                    (1 / numel(mneg)) .* sum(sdm(:, mneg), 2);
                sdmscl = (1 / numel(mpos)) .* sum(sdml(mpos, :), 1) - ...
                    (1 / numel(mneg)) .* sum(sdml(mneg, :), 1);
            elseif numel(mpos) > 1
                sdmsum = (1 / numel(mpos)) .* sum(sdm(:, mpos), 2);
                sdmscl = (1 / numel(mpos)) .* sum(sdml(mpos, :), 1);
            else
                sdmsum = sdm(:, mpos);
                sdmscl = sdml(mpos, :);
            end

            % set new name
            ppin{cc + 1} = sprintf('VOI-x-%s', params.ppicond{cc});

            % and interact
            ppitc(:, cc + 1) = ppitc(:, 1) .* sdmsum;
            ppil(cc + 1, :) = max(0, min(255, round(0.5 .* (sdmscl + ppil(1, :)))));
        end

        % add interacted timecourses
        sdm(:, end:end+size(ppitc, 2)) = [ppitc, ones(size(ppitc, 1), 1)];
        sdmn(end:end+numel(ppin)) = [ppin, {'Constant'}];
        sdml(end:end+numel(ppin), :) = [ppil; sdml(end, :)];
    end
end

% subsample / normalize SDM (for HRF)
if ptype

    % build new SDM
    nsdm = zeros(nvol, size(sdm, 2));

    % for milliseconds PRTs
    if usevol == 0

        % get sample points
        sdmsp = round(1 + (params.srbin - 1) * (params.sbinta / params.sbins));
        sdmsp = sdmsp:(sdmsp + params.sbinte - 1);

        % sample
        for vc = 1:nvol
            nsdm(vc, :) = mean(sdm(sdmsp + (vc - 1) * params.prtr, :), 1);
        end

    % for volume-based PRTs
    else

        % iterate over volumes
        for vc = 1:nvol
            nsdm(vc, :) = mean(sdm(1 + (vc - 1) * usevol:vc * usevol, :));
        end
    end

    % put back
    sdm = nsdm;

    % normalize predictors
    if params.rnorm ~= 0
        sdm = params.rnorm * sdm / max(sdm(:));

        % but keep constant at 1
        sdm(:, end) = 1;
    end

    % orthogonalization
    if ortho

        % iterate over regressors
        oregs = fieldnames(olist);
        for rc = 1:length(oregs)

            % get columns to orthogonalize
            sdmc = olist.(oregs{rc});

            % orthoganlize sdmc(1) against the others
            oc = sdmc(1);
            for occ = sdmc(2:end)
                sdm(:, oc) = orthvec(sdm(:, oc), sdm(:, occ));
            end
        end
    end
end

% work in regi/regni
if nxregi > 0
    regin = cell(1, nxregi);
    for xrc = 1:nxregi
        regin{xrc} = sprintf('regi_%d', xrc);
    end
    sdml = [sdml(1:end-1, :); floor(255.999 .* randn(nxregi, 3)); sdml(end, :)];
    sdmn = [sdmn(1:end-1), regin, sdmn(end)];
    sdm = [sdm(:, 1:end-1), xregi, sdm(:, end)];
end
if nxregni > 0
    regnin = cell(1, nxregni);
    for xrc = 1:nxregni
        regnin{xrc} = sprintf('regni_%d', xrc);
    end
    sdml = [sdml(1:end-1, :); 255 .* ones(nxregni, 3); sdml(end, :)];
    sdmn = [sdmn(1:end-1), regnin, sdmn(end)];
    sdm = [sdm(:, 1:end-1), xregni, sdm(:, end)];
end

% remove derivatives for some regressors
for dc = 1:numel(params.noderiv)
    dci = find(~cellfun('isempty', regexp(sdmn, ['^' params.noderiv{dc} ' - \d\.deriv$'])));
    if ~isempty(dci)
        sdml(dci, :) = [];
        sdmn(dci) = [];
        sdm(:, dci) = [];
    end
end

% difference regressors
if size(params.diffreg, 2) == 2
    for dc = 1:size(params.diffreg, 1)
        if ~ischar(params.diffreg{dc, 1}) || isempty(params.diffreg{dc, 1}) || ...
           ~ischar(params.diffreg{dc, 2}) || isempty(params.diffreg{dc, 2}) || ...
            strcmpi(params.diffreg{dc, 1}(:)', params.diffreg{dc, 2}(:)')
            continue;
        end
        df1 = findfirst(strcmpi(sdmn, params.diffreg{dc, 1}(:)'));
        df2 = findfirst(strcmpi(sdmn, params.diffreg{dc, 2}(:)'));
        if isempty(df1) || df1 > (size(sdm, 2) - (nxregni + 1)) || ...
            isempty(df2) || df2 > (size(sdm, 2) - (nxregni + 1))
            continue;
        end
        sdml(df1, :) = max(0, sdml(df1, :) - sdml(df2, :));
        sdml(df2, :) = [];
        sdmn{df1} = [sdmn{df1} '-' sdmn{df2}];
        sdmn(df2) = [];
        sdm(:, df1) = sdm(:, df1) - sdm(:, df2);
        sdm(:, df2) = [];
    end
end

% update/put into output
bc2.NrOfPredictors = size(sdm, 2);
bc2.PredictorColors = sdml;
bc2.PredictorNames = sdmn;
bc2.FirstConfoundPredictor = size(sdm, 2) - nxregni;
bc2.SDMMatrix = sdm;
bc2.RTCMatrix = sdm(:, 1:end-(1+nxregni));
bc2.RunTimeVars.CreatedFromPRT = sbf;
xo2.C = bc2;

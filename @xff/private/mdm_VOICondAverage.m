function [mtc, mtcse, btc, tc, tcf, tr, onsets, sdms, tcv] = mdm_VOICondAverage(xo, voi, conds, opts, tc, tcf, tr, onsets, sdms)
% MDM::VOICondAverage  - average condition timecourses
%
% FORMAT:       [mtc, mtcse, btc, tc] = mdm.VOICondAverage(voi, conds [, opts])
%
% Input fields:
%
%       voi         VOI object
%       conds       1xC cell array of condition names
%       opts        optional 1x1 struct with fields
%        .avgtype   either of 'conv', 'deconv', or {'meanx'}
%        .avgwin    length of averaging window in ms (default: 20000)
%        .baseline  either of 'none' or {'trial'} (only meanx)
%        .basewin   baseline window in ms (default: -4000:2000:0)
%        .collapse  Cx2 or Cx3 cell array with PRT::Collapse arguments
%        .group     either of 'ffx', {'off'}, or 'rfx'
%        .interp    interpolation method, (see flexinterpn_method, 'cubic')
%        .motpars   motion parameters (Sx1 cell array with sdm/txt files)
%        .motparsd  also add diff of motion parameters (default: false)
%        .motparsq  also add squared motion parameters (default: false)
%        .ndcreg    number of deconvolution time bins (default: 12)
%        .pbar      1x1 progress bar (xfigure::progress or xprogress)
%        .remnuis   remove nuisance (only for 'meanx', default: true)
%        .restcond  remove rest condition (rest cond. name, default: '')
%        .robust    perform robust instead of OLS regression
%        .samptr    upsample TR (milliseconds, if not given use first data)
%        .sdse      either or {'SD'} or 'SE'
%        .subsel    cell array with subject IDs to work on
%        .subvois   subject specific VOIs, either 'sub_', '_sub', {'voi}
%        .tfilter   add filter regressors (cut-off in secs)
%        .tfilttype temporal filter type (either 'dct' or {'fourier'})
%        .trans     either of {'none'}, 'psc', 'z'
%        .unique    flag, only extract unique functional voxels (false)
%        .voisel    cell array with sub-VOI selection to use
%        .xconfound just as motpars, but without restriction on number
%
% Output fields:
%
%       mtc         TxCxS numeric array (mean Time-x-Condition-x-Subject)
%       mtcse       TxCxS numeric array (SD/SE Time-x-Condition-x-Subject)
%       btc         TxCxS numeric array (baseline data)
%       tc          optional Study-by-1 cell array with TxV raw time courses
%
% Note: this call is only valid if all model files in the MDM are PRTs!
%
%       for conv and deconv, a GLM regression is performed (not yet
%       implemented!!)
%
%       the baseline window (basewin) must be evenly spaced!
%
% Using: fitrobustbisquare_img, flexinterpn, flexinterpn_method,
%        meannoinfnan, multimatch, varc, ztrans.

% Version:  v1.1
% Build:    16030309
% Date:     Mar-03 2016, 9:43 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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
using(neuroelf, {'fitrobustbisquare_img', 'flexinterpn', 'flexinterpn_method', ...
    'meannoinfnan', 'multimatch', 'varc', 'ztrans'});

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm') || ...
    numel(voi) ~= 1 || ~xffisobject(voi, true, {'poi', 'voi'}) || ...
   ~iscell(conds) || isempty(conds)
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
cfs = struct('autofind', true, 'silent', true);
try
    if nargin < 9 || ~iscell(tcf) || ~iscell(tc) || isempty(tc) || ...
       ~isa(tr, 'double') || ~iscell(onsets) || isempty(onsets) || ...
       ~iscell(sdms) || numel(tcf) ~= numel(tc) || numel(tcf) ~= size(onsets, 1) || ...
        numel(tcf) ~= numel(tr) || numel(tcf) ~= numel(sdms) || ...
        any(isinf(tr(:)) | isnan(tr(:)) | tr(:) < 0) || numel(conds) ~= size(onsets, 2)
        mdm_CheckFiles(xo, cfs);
    end
catch xfferror
    rethrow(xfferror);
end
conds = conds(:);
dimc = numel(conds);
bc = xo.C;
hc = xo.H;
rtv = bc.RunTimeVars;
if ~any(strcmpi(bc.TypeOfFunctionalData, {'mtc'; 'vtc'})) || ...
   ~any(size(bc.XTC_RTC, 2) == [2, 3])
    error('neuroelf:xff:badArgument', 'MDM must be MTC-/VTC-based.');
end
if any(cellfun('isempty', regexpi(bc.XTC_RTC(:, end), '\.prt$')))
    error('neuroelf:xff:badArgument', 'MDM must be PRT-based.');
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'avgbegin') || ~isa(opts.avgbegin, 'double') || numel(opts.avgbegin) ~= 1 || ...
    isinf(opts.avgbegin) || isnan(opts.avgbegin) || opts.avgbegin < -30000 || opts.avgbegin > 30000
    opts.avgbegin = 0;
end
if ~isfield(opts, 'avgtype') || ~ischar(opts.avgtype) || ...
   ~any(strcmpi(opts.avgtype(:)', {'conv', 'deconv', 'meanx'}))
    opts.avgtype = 'meanx';
else
    opts.avgtype = lower(opts.avgtype(:)');
end
if ~isfield(opts, 'avgwin') || ~isa(opts.avgwin, 'double') || numel(opts.avgwin) ~= 1 || ...
    isinf(opts.avgwin) || isnan(opts.avgwin) || opts.avgwin <= 0
    opts.avgwin = 20000;
end
if ~isfield(opts, 'baseline') || ~ischar(opts.baseline) || ...
   ~any(strcmpi(opts.baseline(:)', {'none', 'trial'}))
    opts.baseline = 'trial';
else
    opts.baseline = lower(opts.baseline(:)');
end
if ~isfield(opts, 'basewin') || ~isa(opts.basewin, 'double') || isempty(opts.basewin) || ...
    any(isinf(opts.basewin(:)) | isnan(opts.basewin(:)))
    opts.basewin = -4000:2000:0;
else
    opts.basewin = opts.basewin(:)';
end
bwinmin = min(opts.basewin);
bwinmax = max(opts.basewin);
if bwinmax > bwinmin
    bwinstep = (bwinmax - bwinmin) / (numel(opts.basewin) - 1);
    if bwinstep == 0
        error('neuroelf:xff:badArgument', 'Invalid baseline window.');
    end
else
    bwinstep = 1;
end
if any(abs(opts.basewin - (bwinmin:bwinstep:bwinmax)) > (0.01 * bwinstep))
    error('neuroelf:xff:badArgument', 'Invalid baseline window spacing.');
end
bwinmax = bwinmax + 0.1 * bwinstep;
if ~isfield(opts, 'collapse') || ~iscell(opts.collapse) || ...
   ~any([2, 3] == size(opts.collapse, 2)) || ndims(opts.collapse) ~= 2
    opts.collapse = cell(0, 2);
end
if ~isfield(opts, 'globsigs') && isfield(rtv, 'GlobalSignals') && iscell(rtv.GlobalSignals)
    opts.globsigs = rtv.GlobalSignals;
end
if ~isfield(opts, 'group') || ~ischar(opts.group) || ...
   ~any(strcmpi(opts.group(:)', {'ffx', 'off', 'rfx'}))
    opts.group = 'off';
else
    opts.group = lower(opts.group(:)');
end
if ~isfield(opts, 'interp') || ~ischar(opts.interp) || ...
    isempty(regexpi(opts.interp(:)', '^(cubic|lanczos\d|linear|nearest)$'))
    opts.interp = 'cubic';
else
    opts.interp = lower(opts.interp(:)');
end
if isfield(opts, 'motpars') && islogical(opts.motpars) && numel(opts.motpars) == 1 && ...
    opts.motpars && isfield(rtv, 'MotionParameters')
    opts.motpars = rtv.MotionParameters;
end
if ~isfield(opts, 'motpars') || ~iscell(opts.motpars)
    opts.motpars = {};
else
    opts.motpars = opts.motpars(:);
end
if ~isfield(opts, 'motparsd') || ~islogical(opts.motparsd) || numel(opts.motparsd) ~= 1
    opts.motparsd = false;
end
if ~isfield(opts, 'motparsq') || ~islogical(opts.motparsq) || numel(opts.motparsq) ~= 1
    opts.motparsq = false;
end
if ~isfield(opts, 'ndcreg') || ~isa(opts.ndcreg, 'double') || numel(opts.ndcreg) ~= 1 || ...
    isinf(opts.ndcreg) || isnan(opts.ndcreg) || opts.ndcreg < 0 || opts.ndcreg > 120
    opts.ndcreg = 0;
else
    opts.ndcreg = floor(opts.ndcreg);
end
if ~isfield(opts, 'pbar') || numel(opts.pbar) ~= 1 || ...
   (~isxfigure(opts.pbar) && ~isa(opts.pbar, 'xprogress'))
    pbar = [];
else
    pbar = opts.pbar;
    pbarvis = pbar.Visible;
end
if ~isfield(opts, 'remnuis') || ~islogical(opts.remnuis) || numel(opts.remnuis) ~= 1
    opts.remnuis = true;
end
if ~isfield(opts, 'restcond') || (~ischar(opts.restcond) && ~iscell(opts.restcond)) || ...
    isempty(opts.restcond)
    opts.restcond = {};
elseif ischar(opts.restcond)
    opts.restcond = {lower(opts.restcond(:)')};
else
    for rcc = numel(opts.restcond):-1:1
        if ~ischar(opts.restcond{rcc}) || isempty(opts.restcond{rcc})
            opts.restcond(rcc) = [];
        else
            opts.restcond{rcc} = lower(opts.restcond{rcc}(:)');
        end
    end
end
if ~isfield(opts, 'robust') || numel(opts.robust) ~= 1 || ~islogical(opts.robust)
    opts.robust = false;
end
if ~isfield(opts, 'samptr') || ~isa(opts.samptr, 'double') || numel(opts.samptr) ~= 1 || ...
    isinf(opts.samptr) || isnan(opts.samptr) || opts.samptr <= 0
    opts.samptr = 0;
else
    opts.samptr = max(1, round(opts.samptr));
end
if ~isfield(opts, 'sdse') || ~ischar(opts.sdse) || ~any(strcmpi(opts.sdse(:)', {'sd', 'se'}))
    opts.sdse = 'sd';
else
    opts.sdse = lower(opts.sdse(:)');
end
if ~isfield(opts, 'subsel') || ~iscell(opts.subsel) || isempty(opts.subsel)
    opts.subsel = mdm_Subjects(xo);
else
    opts.subsel = opts.subsel(:);
    try
        ssm = multimatch(opts.subsel, mdm_Subjects(xo));
    catch xfferror
        rethrow(xfferror);
    end
    if any(ssm < 1)
        error('neuroelf:xff:badArgument', 'Invalid subject ID in selection.');
    end
end
dims = numel(opts.subsel);
if ~isfield(opts, 'tfilter') || ~isa(opts.tfilter, 'double') || numel(opts.tfilter) ~= 1 || ...
    isnan(opts.tfilter) || opts.tfilter < 60
    opts.tfilter = Inf;
end
if ~isfield(opts, 'tfilttype') || ~ischar(opts.tfilttype) || isempty(opts.tfilttype) || ...
    lower(opts.tfilttype(1)) ~= 'd'
    opts.tfilttype = 'fourier';
else
    opts.tfilttype = 'dct';
end
if ~isfield(opts, 'trans') || ~ischar(opts.trans) || isempty(opts.trans) || ...
   ~any('npz' == lower(opts.trans(1)))
    opts.trans = 'n';
    if bc.PSCTransformation > 0
        opts.trans = 'p';
    end
    if bc.zTransformation > 0
        opts.trans = 'z';
    end
end
if ~isfield(opts, 'xconfound') || ~iscell(opts.xconfound)
    opts.xconfound = {};
else
    opts.xconfound = opts.xconfound(:);
end

% required to extract timecourses
extc = true;
try
    if nargin > 8 && iscell(tcf) && iscell(tc) && ~isempty(tc) && isa(tr, 'double') && ...
        iscell(onsets) && ~isempty(onsets) && iscell(sdms) && ...
        numel(tcf) == numel(tc) && numel(tcf) == size(onsets, 1) && ...
        numel(tcf) == numel(tr) && numel(tcf) == numel(sdms) && ...
       ~any(isinf(tr(:)) | isnan(tr(:)) | tr(:) < 0) && numel(conds) == size(onsets, 2)
        tcf = tcf(:);
        tc = tc(:);
        tr = tr(:);
        tcm = multimatch(tcf, bc.XTC_RTC(:, end-1));
        extcp = true;
        numv = size(tc{1}, 2);
        if all(tcm > 0) && numel(unique(tcm(tcm > 0))) == sum(tcm > 0)
            for sc = 1:numel(tc)
                if (~isa(tc{sc}, 'double') && ~isa(tc{sc}, 'single')) || ...
                    ndims(tc{sc}) ~= 2 || isempty(tc{sc}) || size(tc{sc}, 2) ~= numv || ...
                   ~isstruct(sdms{sc}) || numel(sdms{sc}) ~= 1 || ...
                   ~isfield(sdms{sc}, 'FirstConfoundPredictor') || ...
                   ~isa(sdms{sc}.FirstConfoundPredictor, 'double') || ...
                    numel(sdms{sc}.FirstConfoundPredictor) ~= 1 || ...
                   ~isfield(sdms{sc}, 'SDMMatrix') || ~isa(sdms{sc}.SDMMatrix, 'double') || ...
                    size(sdms{sc}.SDMMatrix, 1) ~= size(tc{sc}, 1)
                    extcp = false;
                    break;
                end
            end
            extc = ~extcp;
            tcv = repmat({(1:numv)'}, numel(tcf), 1);
        end
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end

% progress
if ~isempty(pbar)
    pbar.Progress(0, 'Extracting time courses');
    pbar.Visible = 'on';
    drawnow;
end

% required to extract data?
if extc

    % try to extract timecourses
    try
        opts.subsel = unique(mdm_Subjects(xo, true));
        if strcmpi(bc.TypeOfFunctionalData, 'vtc')
            [tc, tcf, tcv, tr] = mdm_VOITimeCourses(xo, voi, opts);
        else
            [tc, tcf, tcv, tr] = mdm_POITimeCourses(xo, voi, opts);
        end
        tco = multimatch(bc.XTC_RTC(:, end-1), tcf);
        if sum(tco > 0) ~= numel(tcf) || numel(unique(tco(tco > 0))) ~= numel(tcf)
            error('neuroelf:xff:internalError', ...
                'Error resolving time course file reference.');
        end
        tc = tc(tco(tco > 0));
        tcf = tcf(tco(tco > 0));
        tcv = tcv(tco(tco > 0));
        tr = tr(tco(tco > 0));
    catch xfferror
        if ~isempty(pbar)
            pbar.Visible = pbarvis;
            drawnow;
        end
        rethrow(xfferror);
    end

    % set NrOfVolumes and TR
    nvol = zeros(numel(tc), 1);
    for stc = 1:numel(nvol)
        nvol(stc) = size(tc{stc}, 1);
    end
    hc.NrOfVolumes = nvol;
    hc.TR = tr;
    xo.H = hc;

    % get onsets
    try
        if ~isempty(pbar)
            pbar.Progress(1, 'Reading onsets...');
            drawnow;
        end
        opts.prtr = [];
        opts.store = true;
        onsets = mdm_ConditionOnsets(xo, conds, opts);
    catch xfferror
        if ~isempty(pbar)
            pbar.Visible = pbarvis;
            drawnow;
        end
        rethrow(xfferror);
    end

    % requires models
    if ~strcmp(opts.avgtype, 'meanx') || opts.remnuis

        % get models
        opts.asstruct = true;
        if ~isempty(pbar)
            pbar.Progress(1, 'Creating models for nuisance regression...');
            drawnow;
        end
        sdms = mdm_SDMs(xo, opts);
    else
        if nargout > 7
            sdms = cell(1, numel(tc));
            for sc = 1:numel(tc)
                sdms{sc} = struct('FirstConfoundPredictor', 1, ...
                    'SDMMatrix', ones(size(tc{sc}, 1), 1));
            end
        else
            sdms = {};
        end
    end

    % subject selection
    subids = mdm_Subjects(xo, true);

% with data delivered
else

    % get number of volumes
    nvol = zeros(numel(tc), 1);
    for stc = 1:numel(nvol)
        nvol(stc) = size(tc{stc}, 1);
    end

    % get subids from filenames
    subids = tcf;
    for sc = 1:numel(subids)
        [null, subids{sc}] = fileparts(subids{sc});
    end
    subids = regexprep(subids, '[_\.].*$', '');
end

% generate selection
subseli = multimatch(subids, opts.subsel);
subseld = (subseli > 0);
if ~any(subseld)
    if ~isempty(pbar)
        pbar.Visible = pbarvis;
        drawnow;
    end
    error('neuroelf:xff:badArgument', 'No subjects selected.');
end

% apply selection
if ~all(subseld)
    subids = subids(subseld);
    tc = tc(subseld);
    tr = tr(subseld);
    onsets = onsets(subseld, :);
    if ~isempty(sdms)
        sdms = sdms(subseld);
    end
    subseli = multimatch(subids, opts.subsel);
end
numstudy = numel(tc);

% get maximum number of onsets
numo = zeros(size(onsets));
for cc = 1:numel(onsets)
    numo(cc) = size(onsets{cc}, 1);
end
dimo = zeros(dims, 1);
for cc = 1:dims
    dimo(cc) = max(sum(numo(subseli == cc, :)));
end
dimo = max(dimo);

% progress
if ~isempty(pbar)
    pbar.Progress(1, 'Averaging...');
    drawnow;
end

% create data
if opts.samptr == 0
    opts.samptr = tr(1);
end
dimtm = ceil(opts.avgwin / opts.samptr);
dimt = dimtm + 1;
dimv = size(tc{1}, 2);
dimtb = numel(bwinmin:bwinstep:bwinmax);
btc = single(NaN .* zeros(dimtb, dimc, dimv, dims, dimo));
mtc = single(NaN .* zeros(dimt, dimc, dimv, dims, dimo));
mtcse = NaN .* zeros(dimt, dimc, dimv, dims);

% initialize condition-subject counters
csi = ones(dimc, dims);

% depending on type of averaging -> averaging
if strcmp(opts.avgtype, 'meanx')

    % get interpolation kernel
    [idata, ik] = flexinterpn_method(zeros(5, 1), [Inf; 1; 0.5; 5], 'cubic');

    % removing nuisance variance first
    if opts.remnuis

        % iterate over remaining studies
        for sc = 1:numstudy

            % regress out nuisance variance
            nsdm = sdms{sc};
            fcp = nsdm.FirstConfoundPredictor;
            if fcp > 0 && fcp <= size(nsdm.SDMMatrix, 2)
                nsdm = nsdm.SDMMatrix;
                nsdm = ztrans(nsdm);
                rsdm = find(sum(abs(nsdm)) == 0);
                if ~isempty(rsdm)
                    nsdm(:, rsdm) = [];
                    fcp = fcp - sum(rsdm < fcp);
                end
                if isempty(nsdm)
                    continue;
                end
                nsdm(:, end+1) = 1;
                if opts.robust
                    [b, w] = fitrobustbisquare_img(nsdm, tc{sc}');
                    bm = b;
                    b(:, end) = 0;
                    if fcp > 1
                        b(:, 1:(fcp-1)) = 0;
                    end
                    bm(:, fcp:end-1) = 0;
                    tc{sc} = w' .* (tc{sc} - nsdm * b') + (1 - w') .* (nsdm * bm');
                else
                    b = pinv(nsdm' * nsdm) * nsdm' * tc{sc};
                    b(end, :) = 0;
                    if fcp > 1
                        b(1:(fcp-1), :) = 0;
                    end
                    tc{sc} = tc{sc} - nsdm * b;
                end
            end
        end
    end

    % create samping window
    bwin = [Inf, Inf; 0, 1; 1, 1; 1, dimv];
    swin = bwin;
    bwinmsm = [bwinmin; bwinstep; bwinmax];

    % iterate over remaining studies
    for sc = 1:numstudy

        % get subject number
        s = subseli(sc);

        % get timecourse and TR
        stc = tc{sc};
        ntc = size(stc, 1);
        st1 = opts.avgbegin / tr(sc);
        str = opts.samptr / tr(sc);

        % iterate over conditions
        for cc = 1:dimc

            % get onsets (in sampling points resolution)
            ons = onsets{sc, cc};
            if isempty(ons)
                continue;
            end
            ons = 1 + ons(:, 1) ./ tr(sc);

            % iterate over onsets
            for oc = 1:numel(ons)

                % create baseline and sampling window
                bwin(2:4, 1) = ons(oc) * [1; 0; 1] + bwinmsm ./ tr(sc);
                swin(2:4, 1) = ons(oc) * [1; 0; 1] + [st1; str; st1 + str * (dimt - 0.5)];

                % skip if either window is undefined
                if all(bwin(4, 1) < 1 | bwin(2, 1) > ntc) || ...
                    all(swin(4, 1) < 1 | swin(2, 1) > ntc)
                    continue;
                end

                % sample values
                bdata = flexinterpn(stc, bwin, ik{:});
                sdata = flexinterpn(stc, swin, ik{:});

                % store
                bwinv = bwin(2):bwin(3):bwin(4);
                bwinv = (bwinv >= 1 & bwinv <= ntc);
                btc(bwinv, cc, :, s, csi(cc, s)) = bdata(bwinv, :);
                swinv = swin(2):swin(3):swin(4);
                swinv = (swinv >= 1 & swinv <= ntc);
                mtc(swinv, cc, :, s, csi(cc, s)) = sdata(swinv, :);

                % increase counter
                csi(cc, s) = csi(cc, s) + 1;
            end
        end
    end

    % baseline computation (remove trial-based influence)
    if strcmp(opts.baseline, 'trial')

        % average across time
        btc = meannoinfnan(btc, 1);

        % set all-zero baselines to bad values
        btc(btc == 0) = NaN;

        % subtract
        mtc = mtc - btc(ones(1, dimt), :, :, :, :);
    end

% -> regression
else


end

% FFX
if strcmp(opts.group, 'ffx')

    % remove subject as a factor
    mtc = reshape(mtc, [dimt, dimc, dimv, dims * dimo]);

    % remove onsets without data
    owd = squeeze(all(all(all(isnan(mtc), 1), 2), 3));
    mtc(:, :, :, owd) = [];
end

% unless rfx
if ~strcmp(opts.group, 'rfx')

    % compute SE measure
    mtcse = sqrt(varc(mtc, ndims(mtc), true));
end

% unless raw
if ~strcmp(opts.group, 'off')

    % build averages
    [mtc, ge, ges] = meannoinfnan(mtc, ndims(mtc));
end

% RFX grouping
if strcmp(opts.group, 'rfx')

    % all-0 responses are NaN!
    mtc(:, all(mtc(:, :) == 0, 1)) = NaN;

    % build averages along subject dimension
    mtcse = sqrt(varc(mtc, 4, true));
    [mtc, ge, ges] = meannoinfnan(mtc, 4);
end

% type of error
if strcmp(opts.sdse, 'se') && ~strcmp(opts.group, 'off')
    mtcse = (1 ./ sqrt(ges)) .* mtcse;
end

% progress bar visibility
if ~isempty(pbar)
    pbar.Visible = pbarvis;
end

% raw time courses
if nargout > 3
    for sc = 1:numel(tc)
        tc{sc} = single(tc{sc});
    end
end

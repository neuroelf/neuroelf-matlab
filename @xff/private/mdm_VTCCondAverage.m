function [vtc, rvtc, wvtc] = mdm_VTCCondAverage(xo, conds, opts)
% MDM::VTCCondAverage  - average condition timecourses
%
% FORMAT:       vtc = mdm.VTCCondAverage(conds [, opts])
%
% Input fields:
%
%       conds       1xC cell array of condition names (or {'listdlg'})
%       opts        optional 1x1 struct with fields
%        .avgwin    length of averaging window in ms (default: 20000)
%        .basewin   baseline window in ms (default: -4000:2000:0)
%        .collapse  Cx2 or Cx3 cell array with PRT::Collapse arguments
%        .condcols  Cx3 or Cx6 or Cx12 double with [0..255] RGB color codes
%        .ffx       compute fixed-effects mean/error (no weighting, true)
%        .gpthresh  global-signal f-thresh for inclusion (default: 0.05)
%        .ithresh   intensity threshold for masking, default: auto-detect
%        .mask      mask object in which to compute the average (reused)
%        .naive     remove only confound regressors from models (false)
%        .remgsig   remove variance from global signal (default: false)
%        .remgsmod  use events to determine GS voxels (default: true)
%        .rfx       compute random-effects mean/error (no weighting, false)
%        .robtune   robust tuning parameter (default: 4.685)
%        .robust    perform 2-pass robust instead of OLS regression (false)
%        .rsngtrial regress out all but the trials being extracted (false)
%        .samptr    sampling TR in ms (default: 400)
%        .smooth    temporal smoothing in ms (after nuisance removal, 0)
%        .subsel    cell array with subject IDs to work on (default: all)
%        .trans     either of 'none', 'psc', 'z' (default: from MDM)
%        .wrfx      compute weighted random-effects mean/error (false)
%
% Output fields:
%
%       vtc         extended VTC object(s) (RunTimeVars)
%
% Note: this call is only valid if all model files in the MDM are PRTs! to
%       compute the nuisance regressors, a call to MDM::SDMs is performed,
%       and all options to this call are carried forward (e.g. motpars).
%
%       the baseline window (basewin) must be evenly spaced!
%
%       the intensity threshold is only used to create an implicit mask!
%
%       the 'st' stats time course type is t-scores (for each time point)
%
%       if multiple subjects are selected, the VTC will contain both a
%       fixed-effects and a (weighted) random-effects error for the mean
%
% Using: calcbetas, findfirst, flexinterpn, flexinterpn_method, histcount,
%        makelabel, meannoinfnan, modelcomp, multimatch, psctrans, sdist,
%        smoothkern, ztrans.

% Version:  v1.1
% Build:    16052509
% Date:     May-25 2016, 9:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013 - 2016, Jochen Weber
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
ddeblank = ne_methods.ddeblank;
calcbetas = ne_methods.calcbetas;
findfirst = ne_methods.findfirst;
flexinterpn = ne_methods.flexinterpn;
meannoinfnan = ne_methods.meannoinfnan;
modelcomp = ne_methods.modelcomp;
multimatch = ne_methods.multimatch;
psctrans = ne_methods.psctrans;
sdist = ne_methods.sdist;
wvartally = ne_methods.wvartally;
ztrans = ne_methods.ztrans;

% preset outputs
vtc = [];
rvtc = [];
wvtc = [];

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm') || ...
   ~iscell(conds) || isempty(conds)
    error('neuroelf:xff:BadArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if ~strcmpi(bc.TypeOfFunctionalData, 'vtc') || size(bc.XTC_RTC, 2) ~= 2 || ...
    any(cellfun('isempty', regexpi(bc.XTC_RTC(:, 1), '\.vtc$')))
    error('neuroelf:xff:badArgument', 'MDM must be VTC-based.');
end
if any(cellfun('isempty', regexpi(bc.XTC_RTC(:, 2), '\.prt$')))
    error('neuroelf:xff:badArgument', 'MDM must be PRT-based.');
end

% options (main argument check)
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end

% list-based condition selection
if numel(conds) == 1 && ischar(conds{1}) && strcmpi(conds{1}(:)', 'listdlg') && ...
    size(xo.C.XTC_RTC, 2) == 2 && ~isempty(regexpi(xo.C.XTC_RTC{1, 2}, '\.prt$'))

    % allow for errors
    try
        prto = {[]};

        % load first protocol
        prto{1} = xff(xo.C.XTC_RTC{1, 2});
        if ~xffisobject(prto{1}, true, 'prt')
            error('neuroelf:xff:badObject', 'Condition selection requires PRT.');
        end

        % apply or initialize collapse argument if necessary
        if isfield(opts, 'collapse') && iscell(opts.collapse) && size(opts.collapse, 2) == 2
            for cc = 1:size(opts.collapse, 1)
                prt_Collapse(prto{1}, opts.collapse{cc, 1}, opts.collapse{cc, 2});
            end
        else
            opts.collapse = cell(0, 2);
        end

        % get current names
        prtc = prt_ConditionNames(prto{1});

        % too many conditions for reasonable PRT (single-trial guess)
        if numel(prtc) > 24
            vsel = 1;

            % repeat until nothing is selected
            while ~isempty(vsel)

                % select from conditions
                [vsel, vsok] = listdlg('ListString', prtc, ...
                    'Name', 'Condition collapsing for averaging...', ...
                    'PromptString', 'Please select conditions to collapse (if any):');
                if isequal(vsok, 0)
                    vsel = [];
                end

                % ask for collapsed name
                if ~isempty(vsel)
                    cname = inputdlg({'Collapsed name'}, 'Collapsed name', 1);
                    if isempty(cname) || isempty(ddeblank(cname{1}))
                        vsel = [];
                    end
                end

                % create collapse string from ALL names (not numeric!)
                if ~isempty(vsel)
                    rstring = sprintf('%s|', prtc{vsel});
                    rstring = ['^(' rstring(1:end-1) ')$'];

                    % perform collapse (error check)
                    prt_Collapse(prto{1}, rstring, ddeblank(cname{1}));

                    % keep track
                    opts.collapse(end+1, :) = {rstring, ddeblank(cname{1})};

                    % and re-get condition names
                    prtc = prt_ConditionNames(prto{1});
                end
            end
        end

        % clear temp object
        clearxffobjects(prto);

        % then get list of desired (collapsed) conditions
        [prtci, selok] = listdlg('ListString', prtc, 'SelectionMode', 'multiple', ...
            'ListSize', [max(420, 12 * size(char(prtc), 2)), 16 * numel(prtc)], ...
            'InitialValue', (1:numel(prtc))', 'Name', 'Condition Selection for averaging...');
        if ~isequal(selok, 1)
            return;
        end
        conds = prtc(prtci);
    catch xfferror
        clearxffobjects(prto);
        rethrow(xfferror);
    end
end
conds = conds(:);

% number of conditions as "dimc"
dimc = numel(conds);

% check condition names
for cc = 1:dimc
    if ~ischar(conds{cc}) || isempty(conds{cc})
        error('neuroelf:xff:badArgument', 'Invalid condition name (%d).', cc);
    end
    conds{cc} = conds{cc}(:)';
end
if numel(unique(conds)) ~= dimc
    error('neuroelf:xff:badArgument', 'Duplicate condition names.');
end

% allow contrast configurations
nconds = regexprep(conds, '^.*(\s*[\>\-]\s*)(\w+)$', '$2');
nconds(strcmp(nconds, conds)) = {''};
conds = regexprep(conds, '^(.*)(\s*[\>\-]\s*)(\w+)$', '$1');

% check referenced files
try
    mdm_CheckFiles(xo, struct('autofind', true, 'silent', true));
catch xfferror
    rethrow(xfferror);
end

% averaging window (default: 20000ms)
if ~isfield(opts, 'avgwin') || ~isa(opts.avgwin, 'double') || numel(opts.avgwin) ~= 1 || ...
    isinf(opts.avgwin) || isnan(opts.avgwin) || opts.avgwin <= 0
    opts.avgwin = 20000;
end

% baseline window
if ~isfield(opts, 'basewin') || ~isa(opts.basewin, 'double') || isempty(opts.basewin) || ...
    any(isinf(opts.basewin(:)) | isnan(opts.basewin(:)))
    basewin = -4000:2000:0;
else
    basewin = opts.basewin(:)';
end
bwinmin = min(basewin);
bwinmax = max(basewin);
if bwinmax > bwinmin
    bwinstep = ceil((bwinmax - bwinmin) / (numel(basewin) - 1));
else
    bwinstep = 1;
end
if any(abs(basewin - (bwinmin:bwinstep:bwinmax)) > (0.01 * bwinstep))
    error('neuroelf:xff:badArgument', 'Invalid baseline window spacing.');
end
bwinmax = bwinmax + 0.1 * bwinstep;
basewin = bwinmin:bwinstep:bwinmax;

% collapsing (if not yet defined from UI phase)
if ~isfield(opts, 'collapse') || ~iscell(opts.collapse) || ...
   ~any([2, 3] == size(opts.collapse, 2)) || ndims(opts.collapse) ~= 2
    opts.collapse = cell(0, 2);
end
% ocollapse = opts.collapse;
% for cc = 1:size(ocollapse, 1)
%     if ~ischar(ocollapse{cc, 1}) || isempty(ocollapse{cc, 1})
%         error('neuroelf:xff:badArgument', 'Invalid collapsing argument.');
%     end
% end

% condition colors
if ~isfield(opts, 'condcols') || ~isa(opts.condcols, 'double') || ...
    size(opts.condcols, 1) ~= dimc || size(opts.condcols, 2) < 3 || ...
    any(isinf(opts.condcols(:)) | isnan(opts.condcols(:)) | opts.condcols(:) < 0 | opts.condcols(:) > 255)
    opts.condcols = zeros(dimc, 0);
else
    opts.condcols(:, 17:end) = [];
    if all(opts.condcols(:) <= 1)
        opts.condcols = round(255 .* opts.condcols);
    end
    opts.condcols = uint8(round(opts.condcols));
end

% compute FFX (keep FFX tally)
if ~isfield(opts, 'ffx') || numel(opts.ffx) ~= 1 || ~islogical(opts.ffx)
    opts.ffx = true;
end

% global p (F) thresh for GS inclusion
if ~isfield(opts, 'gpthresh') || ~isa(opts.gpthresh, 'double') || numel(opts.gpthresh) ~= 1 || ...
    isinf(opts.gpthresh) || isnan(opts.gpthresh) || opts.gpthresh <= 0 || opts.gpthresh >= 0.5
    opts.gpthresh = 0.05;
else
    opts.gpthresh = min(0.25, opts.gpthresh);
end

% intensity threshold
if ~isfield(opts, 'ithresh') || ~isa(opts.ithresh, 'double') || numel(opts.ithresh) ~= 1 || ...
    isinf(opts.ithresh) || isnan(opts.ithresh) || opts.ithresh < 0
    opts.ithresh = -1;
end

% mask
if ~isfield(opts, 'mask') || ~xffisobject(opts.mask, true, 'msk')
    if ~isfield(opts, 'mask') || ~islogical(opts.mask) || ndims(opts.mask) ~= 3
        opts.mask = [];
    end
else
    opts.mask = opts.mask.C;
    opts.mask = (opts.mask.Mask(:, :, :) > 0);
end

% naive mode (remove all *other* regressors from nuisance model?)
if ~isfield(opts, 'naive') || numel(opts.naive) ~= 1 || ~islogical(opts.naive)
    opts.naive = false;
end

% TR (default: auto-detect)
if ~isfield(opts, 'prtr') || ~isa(opts.prtr, 'double') || numel(opts.prtr) ~= 1 || ...
    isinf(opts.prtr) || isnan(opts.prtr) || opts.prtr <= 0
    opts.prtr = [];
end

% remove global signal
if ~isfield(opts, 'remgsig') || numel(opts.remgsig) ~= 1 || ~islogical(opts.remgsig)
    opts.remgsig = false;
end
if ~isfield(opts, 'remgsmod') || numel(opts.remgsmod) ~= 1 || ~islogical(opts.remgsmod)
    opts.remgsmod = true;
end

% compute plain RFX (unweighted group tally)
if ~isfield(opts, 'rfx') || numel(opts.rfx) ~= 1 || ~islogical(opts.rfx)
    opts.rfx = false;
end

% robust tuning parameter
if ~isfield(opts, 'robtune') || ~isa(opts.robtune, 'double') || numel(opts.robtune) ~= 1 || ...
    isinf(opts.robtune) || isnan(opts.robtune) || opts.robtune <= 0
    opts.robtune = 4.685;
end

% use robust estimate
if ~isfield(opts, 'robust') || numel(opts.robust) ~= 1 || ~islogical(opts.robust)
    opts.robust = false;
end

% regress out all but single-trial (different nuisance model for each trial)
if ~isfield(opts, 'rsngtrial') || ~islogical(opts.rsngtrial) || numel(opts.rsngtrial) ~= 1
    opts.rsngtrial = false;
end

% sampling (virtual) TR
if ~isfield(opts, 'samptr') || ~isa(opts.samptr, 'double') || numel(opts.samptr) ~= 1 || ...
    isinf(opts.samptr) || isnan(opts.samptr) || opts.samptr <= 0
    samptr = 400;
else
    samptr = max(10, round(opts.samptr));
end

% temporal smoothing
if ~isfield(opts, 'smooth') || ~isa(opts.smooth, 'double') || numel(opts.smooth) ~= 1 || ...
    isinf(opts.smooth) || isnan(opts.smooth) || opts.smooth < 0
    opts.smooth = 0;
else
    opts.smooth = min(10000, opts.smooth);
end

% subject selection (default: all)
mdmsubs = mdm_Subjects(xo);
if ~isfield(opts, 'subsel') || ~iscell(opts.subsel) || isempty(opts.subsel)
    opts.subsel = mdmsubs;
else

    % ensure subjects are valid strings
    opts.subsel = opts.subsel(:);
    for sc = numel(opts.subsel):-1:1
        if ~ischar(opts.subsel{sc}) || isempty(opts.subsel{sc})
            opts.subsel(sc) = [];
        else
            opts.subsel{sc} = opts.subsel{sc}(:)';
        end
    end
    opts.subsel = unique(opts.subsel);

    % match against actual subjects
    try
        ssm = multimatch(opts.subsel, mdmsubs);
    catch xfferror
        rethrow(xfferror);
    end
    
    % any missing
    if any(ssm < 1)
        error('neuroelf:xff:badArgument', 'Invalid subject ID in selection.');
    end
end
numsubs = numel(opts.subsel);

% time-course transformation
transv = 0;
if ~isfield(opts, 'trans') || ~ischar(opts.trans) || isempty(opts.trans) || ...
   ~any('npz' == lower(opts.trans(1)))
    opts.trans = 'n';
    if bc.PSCTransformation > 0
        opts.trans = 'p';
        transv = 3;
    end
    if bc.zTransformation > 0
        opts.trans = 'z';
        transv = 1;
    end
elseif lower(opts.trans(1)) == 'p'
    opts.trans = 'p';
    transv = 3;
elseif lower(opts.trans(1)) == 'z'
    opts.trans = 'z';
    transv = 1;
end
tlim = 5;
if transv == 1
    tlim = 2;
elseif transv == 3
    tlim = 3;
end

% compute weighted RFX (variance-weighted group tally)
if ~isfield(opts, 'wrfx') || numel(opts.wrfx) ~= 1 || ~islogical(opts.wrfx)
    opts.wrfx = false;
end

% at least ONE computation
if ~opts.ffx && ~opts.rfx && ~opts.wrfx
    error('neuroelf:xff:badOption', 'At least one computation must be enabled.');
end

% objects and subject IDs
vtcs = bc.XTC_RTC(:, 1);
prts = bc.XTC_RTC(:, 2);
subids = mdm_Subjects(xo, true);
usevtc = (multimatch(subids, opts.subsel) > 0);
numvtcs = sum(usevtc);
if numvtcs < 1
    error('neuroelf:xff:badArgument', 'Invalid arguments or options.');
end
if ~all(usevtc)
    subids = subids(usevtc);
    vtcs = vtcs(usevtc);
    prts = prts(usevtc);
end

% length of averaging window
avgwin = 0:samptr:opts.avgwin;
if avgwin(end) < opts.avgwin
    avgwin(end+1) = avgwin(end) + avgwin(2);
end
avgwinmin = avgwin(1);
avgwinmax = avgwin(end) + 0.1 * samptr;
dima = numel(avgwin);

% how many time courses
dimt = dima * 2 * dimc;

% get onsets
try
    opts.store = true;
    onsets = mdm_ConditionOnsets(xo, conds, opts);
    nonsets = cell(size(onsets));
    nci = find(~cellfun('isempty', nconds));
    if ~isempty(nci)
        nonsets(:, nci) = mdm_ConditionOnsets(xo, nconds(nci), opts);
    end
catch xfferror
    rethrow(xfferror);
end

% get models
opts.asstruct = true;
try
    opts.sngtrial = false;
    sdms = mdm_SDMs(xo, opts);

    % every condition must have at least one match
    cmatch = false(dimc, 1);
    for sc = 1:numel(sdms)
        cmatch = cmatch | (multimatch(conds, sdms{sc}.PredictorNames(:)) > 0);
        if all(cmatch)
            break;
        end
    end
    if ~all(cmatch)
        error('neuroelf:xff:badArgument', 'Not all conditions covered.');
    end

    % for single-trial regression
    if opts.rsngtrial

        % get those SDMs as well
        %opts.collapse = {};
        opts.sngtrial = true;
        ssdms = mdm_SDMs(xo, opts);
%         for sc = 1:numel(ssdms)
%             pnames = ssdms{sc}.PredictorNames(1:ssdms{sc}.FirstConfoundPredictor-1);
%             tnums = regexprep(pnames, '^.*(_T\d+)$', '$1');
%             pnames = regexprep(pnames, '^(.*)_T\d+$', '$1');
%             for cc = 1:size(ocollapse, 1)
%                 pnames = regexprep(pnames, ocollapse{cc, 1}, ocollapse{cc, 2});
%             end
%             for cc = 1:numel(pnames)
%                 pnames{cc} = [pnames{cc} tnums{cc}];
%             end
%             ssdms{sc}.PredictorNames(1:ssdms{sc}.FirstConfoundPredictor-1) = pnames;
%         end
    end
catch xfferror
    rethrow(xfferror);
end

% get colors
if size(opts.condcols, 2) == 3
    ccol = opts.condcols;
    dcol = uint8(round(0.333 .* double(ccol)));
    oo = ones(size(ccol, 1), 1);
    opts.condcols = [dcol, uint8(64 .* oo), ccol, uint8(255 .* oo), ...
            uint8(255 .* oo) - dcol, uint8(64 .* oo), uint8(255 .* oo) - ccol, uint8(255 .* oo)];
elseif size(opts.condcols, 2) == 6
    ccol1 = opts.condcols(:, 1:3);
    dcol1 = uint8(round(0.333 .* double(ccol1)));
    ccol2 = opts.condcols(:, 4:6);
    dcol2 = uint8(round(0.333 .* double(ccol2)));
    oo = ones(size(ccol1, 1), 1);
    opts.condcols = [dcol1, uint8(64 .* oo), ccol1, uint8(255 .* oo), ...
            dcol2, uint8(64 .* oo), ccol2, uint8(255 .* oo)];
elseif size(opts.condcols, 2) == 12
    oo = ones(size(opts.condcols, 1), 1);
    opts.condcols = [opts.condcols(:, 1:3), uint8(64 .* oo), opts.condcols(:, 4:6), uint8(255 .* oo), ...
        opts.condcols(:, 7:9), uint8(64 .* oo), opts.condcols(:, 10:12), uint8(255 .* oo)];
elseif size(opts.condcols,2 ) ~= 16
    prednames = sdms{1}.PredictorNames;
    predcols = sdms{1}.PredictorColors;
    opts.condcols(end, 16) = 0;
    for cc = 1:dimc
        predi = findfirst(strcmpi(prednames(:), conds{cc}));
        if isempty(predi)
            ccol = uint8(floor(255.9999 .* rand(1, 3)));
        else
            ccol = uint8(round(predcols(predi, :)));
        end
        dcol = uint8(round(0.333 .* double(ccol)));
        opts.condcols(cc, :) = [dcol, uint8(64), ccol, uint8(255), ...
            uint8(255) - dcol, uint8(64), uint8(255) - ccol, uint8(255)];
    end
end

% apply selection
if ~all(usevtc)
    onsets = onsets(usevtc, :);
    if ~isempty(sdms)
        sdms = sdms(usevtc);
    end
end
numstudy = numel(vtcs);

% load VTC objects (transio access)
vtco = vtcs;
onsrem = zeros(numstudy, dimc);
try
    for sc = 1:numstudy
        vtco{sc} = xff(vtcs{sc}, 't');
        if ~xffisobject(vtco{sc}, true, 'vtc')
            error('neuroelf:xff:badObject', 'Not a valid VTC object (%s).', vtcs{sc});
        end

        % for first VTC
        tvtcc = vtco{sc}.C;
        tlay = aft_Layout(vtco{sc});
        tlnv = tlay(4);
        tlay(4) = [];
        if sc == 1

            % copy object (for output)
            vtc = aft_CopyObject(vtco{sc});

            % store layout
            vlay = tlay;
            vsz = vlay(1:3);
            numvox = prod(vsz);

        % for others
        else

            % check layout
            if ~isequal(vlay, tlay)
                error('neuroelf:xff:badObject', ...
                    'VTCs must match in dims and offsets (mismatch: %s).', vtcs{sc});
            end
        end

        % go over onsets
        lastonstime = tvtcc.TR * tlnv - 2500;
        for cc = 1:size(onsets, 2)
            if isempty(onsets{sc, cc})
                continue;
            end
            remons = find(onsets{sc, cc}(:, 1) >= lastonstime);
            if ~isempty(remons)
                onsets{sc, cc}(remons, :) = [];
                onsrem(sc, cc) = numel(remons);
            end
        end
    end
catch xfferror
    clearxffobjects(vtco);
    delete(vtc);
    rethrow(xfferror);
end

% total number of onsets
try
    nonsets = size(cat(1, onsets{:}), 1);
catch xfferror
    neuroelf_lasterr(xfferror);
    for cc = 1:numel(onsets)
        onsets{cc} = onsets{cc}(:, 1:2);
    end
    nonsets = size(cat(1, onsets{:}), 1);
end
consets = 0;

% no or invalid mask given
rtv = xo.C.RunTimeVars;
if ~islogical(opts.mask) || ~isequal(size(opts.mask), vsz)

    % check for averaging mask in handles
    if isfield(rtv, 'VTCAveragingMask') && iscell(rtv.VTCAveragingMask) && ...
        numel(rtv.VTCAveragingMask) == 2 && isequal(rtv.VTCAveragingMask{1}, vtcs)

        % re-use
        opts.mask = rtv.VTCAveragingMask{2};

    % no mask yet
    else

        % initialize progress bar
        try
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', 'Full-VTC averaging: mask generation...');
            xprogress(pbar, 0, sprintf('Averaging VTC 1/%d...', numel(vtcs)), 'visible', 0, numel(vtcs));
        catch xfferror
            neuroelf_lasterr(xfferror);
            pbar = [];
        end

        % create mask
        opts.mask = zeros(vsz);
        for rc = 1:numel(vtcs)
            vtccont = vtco{rc}.C;
            if ~isempty(pbar)
                xprogress(pbar, rc - 1, sprintf('Averaging VTC %d/%d...', rc, numel(vtcs)));
            end
            opts.mask = opts.mask + squeeze(mean(vtccont.VTCData(:, :, :, :), 1));
        end
        opts.mask = (1 / numel(vtcs)) .* opts.mask;

        % auto-detect threshold
        if ~isempty(pbar)
            xprogress(pbar, rc, 'Augmenting mask...');
        end
        if opts.ithresh < 0
            ithreshm = double(ceil(max(opts.mask(:))));
            ithreshc = ne_methods.histcount(opts.mask, 0, ithreshm, ithreshm / 511.000001);
            ithreshc = flexinterpn(ithreshc(:), [Inf; 1; 1; numel(ithreshc)], ne_methods.smoothkern(5), 1);
            opts.ithresh = 0.25 * ithreshm * (findfirst(diff(ithreshc) > 0, -1) / 511);
        end

        % intersect with colin brain if available
        colinbrain = [neuroelf_path('colin') '/colin_brain.vmr'];
        colinbicbm = [neuroelf_path('colin') '/colin_brain_ICBMnorm.vmr'];
        colinbtaln = [neuroelf_path('colin') '/colin_brain_TALnorm.vmr'];
        if exist(colinbrain, 'file') > 0
            colin = {[]};
            try
                colin{1} = xff(colinbrain);
                sc = aft_SampleBVBox(colin{1}, aft_BoundingBox(vtc), 1) > 0;
                delete(colin{1});
            catch xfferror
                clearxffobjects(colin);
                neuroelf_lasterr(xfferror);
                sc = true(size(opts.mask));
            end
        else
            sc = true(size(opts.mask));
        end
        if exist(colinbicbm, 'file') > 0
            colin = {[]};
            try
                colin{1} = xff(colinbicbm);
                rc = aft_SampleBVBox(colin{1}, aft_BoundingBox(vtc), 1) > 0;
                delete(colin{1});
            catch xfferror
                clearxffobjects(colin);
                neuroelf_lasterr(xfferror);
                rc = true(size(opts.mask));
            end
        else
            rc = true(size(opts.mask));
        end
        if exist(colinbtaln, 'file') > 0
            colin = {[]};
            try
                colin{1} = xff(colinbtaln);
                cc = aft_SampleBVBox(colin{1}, aft_BoundingBox(vtc), 1) > 0;
                delete(colin{1});
            catch xfferror
                clearxffobjects(colin);
                neuroelf_lasterr(xfferror);
                cc = true(size(opts.mask));
            end
        else
            cc = true(size(opts.mask));
        end

        % intersect
        opts.mask = ((opts.mask >= opts.ithresh) & (sc | rc | cc)) | ...
            ((opts.mask >= (0.25 .* opts.ithresh)) & sc & rc & cc);

        % save
        xo.C.RunTimeVars.AutoSave = true;
        xo.C.RunTimeVars.VTCAveragingMask = {vtcs, opts.mask};
        if ~isempty(xo.F)
            try
                aft_SaveRunTimeVars(xo);
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
        end

        % close bar
        if ~isempty(pbar)
            closebar(pbar);
        end
    end
end
mask = opts.mask;
smask = sum(mask(:));
osmask = ones(1, smask);

% initialize progress bar
try
    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 640, 36]);
    xprogress(pbar, 'settitle', sprintf('Full-VTC averaging %d onsets...', nonsets));
    xprogress(pbar, 0, 'Preparation...', 'visible', 0, nonsets);
catch xfferror
    neuroelf_lasterr(xfferror);
    pbar = [];
end

% get interpolation kernel
[tw, ipk] = ne_methods.flexinterpn_method([0; 1; 0], [Inf; 0; 1; 3], 'cubic');

% VTC settings
vtc.C.NameOfSourceFMR = xo.F;
vtc.C.NrOfCurrentPRT = 1;
vtc.C.DataType = 2;
vtc.C.NrOfVolumes = dimt;
vtc.C.TR = samptr;
vtc.C.VTCData = single(zeros([0, vlay(1:3)]));
vtc.C.VTCData(dimt, 1, 1, 1) = 0;
vtc.C.RunTimeVars.AutoSave = true;
vtc.C.RunTimeVars.AvgVTC = true;
vtc.C.RunTimeVars.AvgGlobSigRemoved = opts.remgsig;
vtc.C.RunTimeVars.AvgRobust = opts.robust;
vtc.C.RunTimeVars.AvgRobustTune = opts.robtune;
vtc.C.RunTimeVars.AvgSmooth = opts.smooth;
vtc.C.RunTimeVars.AvgTransformationType = transv;
vtc.C.RunTimeVars.AvgWindowFrom = 0;
vtc.C.RunTimeVars.AvgWindowStep = avgwin(2);
vtc.C.RunTimeVars.AvgWindowTo = avgwin(end);
vtc.C.RunTimeVars.BaseWindowFrom = basewin(1);
vtc.C.RunTimeVars.BaseWindowStep = bwinstep;
vtc.C.RunTimeVars.BaseWindowTo = basewin(end);
vtc.C.RunTimeVars.Discard = [];
vtc.C.RunTimeVars.DVARS = [];
vtc.C.RunTimeVars.FFXWeights = struct;
vtc.C.RunTimeVars.GlobSigs = cell(0, 2);
vtc.C.RunTimeVars.Map = repmat(struct( ...
    'Type', 30, 'LowerThreshold', 0, 'UpperThreshold', 0.5, ...
    'Name', '', 'RGBLowerThreshPos', [255, 0, 0], 'RGBUpperThreshPos', [255, 255, 0], ...
    'RGBLowerThreshNeg', [255, 0, 255], 'RGBUpperThreshNeg', [0, 0, 255], ...
    'UseRGBColor', 1, 'LUTName', '<default>', 'TransColorFactor', 1, ...
    'NrOfLags', 0, 'MinLag', 0, 'MaxLag', 0, 'CCOverlay', 0, ...
    'ClusterSize', 4, 'EnableClusterCheck', 0, 'UseValuesAboveThresh', 1, ...
    'DF1', 0, 'DF2', 0, 'ShowPositiveNegativeFlag', 1, ...
    'BonferroniValue', 0, 'NrOfFDRThresholds', 0, 'FDRThresholds', zeros(0, 3), ...
    'OverlayColors', []), 1, dimc);
vtc.C.RunTimeVars.MapSelection = {{}, []};
vtc.C.RunTimeVars.MotionParameters = zeros(0, 6);
vtc.C.RunTimeVars.MPFD = [];
vtc.C.RunTimeVars.NrOfConditions = dimc;
vtc.C.RunTimeVars.NrOfConditionOnsets = zeros(1, dimc);
vtc.C.RunTimeVars.NrOfTCsPerCondition = 2;
vtc.C.RunTimeVars.NrOfVolumesPerTC = dima;
vtc.C.RunTimeVars.NrOfSourceVTCs = numel(usevtc);
vtc.C.RunTimeVars.NrOfSubjects = numsubs;
vtc.C.RunTimeVars.ConditionColors = opts.condcols;
vtc.C.RunTimeVars.ConditionNames = conds;
vtc.C.RunTimeVars.ConditionOnsets = onsets;
vtc.C.RunTimeVars.ConditionThresholds = repmat(reshape([0, 2, 0.5, 5], [1, 2, 2]), dimc, 1);
vtc.C.RunTimeVars.SubjectNames = opts.subsel;
vtc.C.RunTimeVars.TCNames = {'mean', 'z'};
vtc.C.RunTimeVars.ScalingWindow = [-2, 2];
vtc.C.RunTimeVars.ScalingWindowLim = [0.25, 1];
vtc.C.RunTimeVars.SourcePRTs = prts;
vtc.C.RunTimeVars.SourceVTCs = vtcs;
vtc.C.RunTimeVars.SubMapVol = 1;
vtc.C.VTCData(:) = 0;
if opts.rfx
    rvtc = aft_CopyObject(vtc);
end
if opts.wrfx
    wvtc = aft_CopyObject(vtc);
end

% prepare tally for group
ftally = cell(1, dimc);
rtally = cell(1, dimc);
wtally = cell(1, dimc);
onscount = zeros(1, dimc);
osubids = subids;
osubid = '_NOTASUBJECT_';

% iterate over subjects
for sc = 1:numel(opts.subsel)

    % find VTCs
    usevtc = find(strcmp(subids, opts.subsel{sc}));

    % reset data
    tally = cell(1, dimc);

    % for all runs for this subject
    for rc = 1:numel(usevtc)

        % get data
        vtci = usevtc(rc);
        ons = onsets(vtci, :);
        if all(cellfun('isempty', ons))
            continue;
        end
        vtccont = vtco{vtci}.C;
        if isempty(opts.prtr)
            vtctr = vtccont.TR;
        else
            vtctr = opts.prtr;
        end
        obwinstep = bwinstep / vtctr;
        awinstep = samptr / vtctr;
        if ~isempty(pbar)
            [vtcfolder, vtcshort] = fileparts(vtcs{vtci});
            xprogress(pbar, consets, sprintf('Reading data from %s...', vtcshort));
        end
        vtcdata = vtccont.VTCData(:, :, :, :);
        nvtctp = size(vtcdata, 1);
        sdmmatrix = sdms{vtci}.SDMMatrix;
        sdmconds = sdms{vtci}.PredictorNames(:);
        sdmconf1 = sdms{vtci}.FirstConfoundPredictor;

        % for single-trial designs, get some additional info
        if opts.rsngtrial
            ssdmmatrix = ssdms{vtci}.SDMMatrix;
            ssdmconds = ssdms{vtci}.PredictorNames(:);

            % new onset counters
            if ~strcmp(osubid, osubids{vtci})
                sonscount = zeros(1, dimc);
                osubid = osubids{vtci};
            end
        end

        % smoothing
        if opts.smooth > 0
            smk = ne_methods.smoothkern(opts.smooth / vtctr);
        else
            smk = 1;
        end

        % remove empty regressors
        removebets = all(sdmmatrix == 0);
        sdmmatrix(:, removebets) = [];
        sdmconds(removebets) = [];

        % reshape
        if ~isempty(pbar)
            xprogress(pbar, consets, sprintf('Masking data from %s...', vtcshort));
        end
        vtcdata = double(reshape(vtcdata, nvtctp, numvox));
        vtcdata = vtcdata(:, mask);
        vtccsz = size(vtcdata);

        % progress
        if ~isempty(pbar)
            xprogress(pbar, consets, sprintf('Regressing out nuisance from %s...', vtcshort));
        end

        % for robust computation
        if opts.robust

            % compute weights (first pass)
            if opts.naive
                [betas, iXX, cc] = calcbetas(sdmmatrix(:, sdmconf1:end), vtcdata, 1);
            else
                [betas, iXX, cc] = calcbetas(sdmmatrix, vtcdata, 1);
            end
            w = vtcdata - cc;
            ws = opts.robtune .* max(std(w), sqrt(eps));
            w = repmat((1 ./ ws), nvtctp, 1) .* w;
            wo = (abs(w) < 1) .* (1 - w .^ 2) .^ 2;

            % compute global signal
            if opts.remgsig

                % full model
                if opts.remgsmod

                    % find regressors to drop (conditions)
                    gpidx = setdiff((1:size(sdmmatrix, 2))', multimatch(conds, sdmconds));

                    % compute F-stat
                    [gsigf, gdf1, gdf2] = modelcomp(sdmmatrix, sdmmatrix(:, gpidx), vtcdata, 1);

                    % compute threshold
                    gsfthresh = sdist('finv', opts.gpthresh, gdf1, gdf2, true);

                    % mask
                    gsig = sum(wo(:, gsigf <= gsfthresh) .* vtcdata(:, gsigf <= gsfthresh), 2) ./ ...
                        sum(wo(:, gsigf <= gsfthresh), 2);

                    % and add GS from mask
                    if ~any(isnan(gsig))
                        sdmmatrix = [sdmmatrix(:, 1:end-1), ztrans(gsig), sdmmatrix(:, end)];
                    end

                % only nuisance
                else

                    % compute F-stat of nuisance over mean-only
                    [gsigf, gdf1, gdf2] = ...
                        modelcomp(sdmmatrix(:, sdmconf1:end), sdmmatrix(:, end), vtcdata, 1);

                    % divide, then subtract (mask)
                    gsigf = min(2, gsigf ./ mean(gsigf));
                    gsigf(isinf(gsigf) | isnan(gsigf)) = 0;
                    gsigf = gsigf ./ (sum(gsigf) / numel(gsigf));

                    % add global signal
                    gsig = sum(wo .* vtcdata .* (ones(nvtctp, 1) * gsigf(:)'), 2) ./ ...
                        sum(wo, 2);
                    sdmmatrix = [sdmmatrix(:, 1:end-1), ztrans(gsig), sdmmatrix(:, end)];
                end
            end

            % then adapt data
            vtcdata = wo .* vtcdata + (1 - wo) .* cc;

            % then second pass
            if opts.naive
                [betas, iXX, cc] = calcbetas(sdmmatrix(:, sdmconf1:end), vtcdata, 1);
            else
                [betas, iXX, cc] = calcbetas(sdmmatrix, vtcdata, 1);
            end
            w = vtcdata - cc;
            ws = opts.robtune .* max(std(w), sqrt(eps));
            w = repmat((1 ./ ws), nvtctp, 1) .* w;
            w = (abs(w) < 1) .* (1 - w .^ 2) .^ 2;
            vtcdata = w .* vtcdata + (1 - w) .* cc;

        % remove global signal
        elseif opts.remgsig

            % run F-test
            if opts.remgsmod
                gpidx = setdiff((1:size(sdmmatrix, 2))', multimatch(conds, sdmconds));
                [gsigf, gdf1, gdf2] = modelcomp(sdmmatrix, sdmmatrix(:, gpidx), vtcdata, 1);
                gsfthresh = sdist('finv', opts.gpthresh, gdf1, gdf2, true);
                gsig = mean(vtcdata(:, gsigf <= gsfthresh), 2);
                if ~any(isnan(gsig))
                    sdmmatrix = [sdmmatrix(:, 1:end-1), ztrans(gsig), sdmmatrix(:, end)];
                end
            else
                [gsigf, gdf1, gdf2] = ...
                    modelcomp(sdmmatrix(:, sdmconf1:end), sdmmatrix(:, end), vtcdata, 1);
                gsfthresh = sdist('finv', opts.gpthresh, gdf1, gdf2, true);
                gsigf = max(0, (sqrt(gsigf) ./ sqrt(gsfthresh)) - 1);
                gsigf(isinf(gsigf) | isnan(gsigf)) = 0;
                if opts.robust
                    gsig = sum(wo .* vtcdata .* (ones(nvtctp, 1) * gsigf(:)'), 2) ./ ...
                        (sum(gsigf(:)) .* sum(wo, 2));
                else
                    gsig = sum(vtcdata .* (ones(nvtctp, 1) * gsigf(:)'), 2) ./ ...
                        (sum(gsigf(:)));
                end
                sdmmatrix = [sdmmatrix(:, 1:end-1), ztrans(gsig), sdmmatrix(:, end)];
            end
        end

        % transformation
        if transv == 1
            vtcdata = ztrans(vtcdata);
        elseif transv == 3
            vtcdata = psctrans(vtcdata);
        end

        % final (or first) pass
        betas = calcbetas(sdmmatrix, vtcdata, 1)';

        % for PSC trans
        if transv == 3
            vtcdata = repmat((100 ./ betas(end, :)), nvtctp, 1) .* vtcdata;
            if ~opts.rsngtrial
                betas = calcbetas(sdmmatrix, vtcdata, 1)';
            end
        end

        % regular weights (to account for sampling beyond data limits)
        tw = ones(nvtctp, 1);

        % iterate over conditions
        for cc = 1:dimc

            % no onsets, continue
            if isempty(ons{cc})
                continue;
            end
            o = 1 + ons{cc}(:, 1) ./ vtctr;
            no = size(o, 1);
            if no == 0
                continue;
            end

            % which predictor
            pidx = find(strcmp(sdmconds, conds{cc}));
            if isempty(pidx)
                consets = consets + no;
                continue;
            end
            onscount(cc) = onscount(cc) + no;

            % process condition
            if ~isempty(pbar)
                xprogress(pbar, consets, sprintf('Processing condition %s...', conds{cc}));
            end

            % for condition-wide processing
            if ~opts.rsngtrial

                % remove residual variance
                usebets = 1:size(betas, 1);
                usebets([pidx, end]) = [];
                if opts.naive
                    vtccorr = vtcdata;
                else
                    vtccorr = vtcdata - sdmmatrix(:, usebets) * betas(usebets, :);
                end

                % smoothing
                if numel(smk) ~= 1
                    vtccorr = flexinterpn(vtccorr, [Inf, Inf; 1, 1; 1, 1; vtccsz], ...
                        {smk, [0; 1; 0]}, {1, 1});
                end
            end

            % for each onset
            for oc = 1:no

                % for single-trial processing
                if opts.rsngtrial && ~opts.naive

                    % find trial number in single-design matrix
                    sonscount(cc) = sonscount(cc) + 1;
                    tidx = find(~cellfun('isempty', regexpi(ssdmconds, ...
                        sprintf('^%s_T0*%d$', conds{cc}, sonscount(cc)))));
                    if isempty(tidx) && no == 1
                        tidx = find(strcmpi(ssdmconds, conds{cc}));
                    end
                    if isempty(tidx)
                        if ~isempty(pbar)
                            closebar(pbar);
                        end
                        clearxffobjects(vtco);
                        delete(vtc);
                        delete(rvtc);
                        delete(wvtc);
                        error('neuroelf:xff:internalError', ...
                            'Error locating single-trial regressor for %s:%d.', conds{cc}, oc);
                    end

                    % first, take the original design matrix
                    rsdmmatrix = [sdmmatrix, zeros(size(sdmmatrix, 1), 1)];

                    % then remove influence of single trial
                    rsdmmatrix(:, pidx) = rsdmmatrix(:, pidx) - ssdmmatrix(:, tidx);

                    % and add this as additional regressor as last
                    rsdmmatrix(:, end) = ssdmmatrix(:, tidx);

                    % only one onset -> remove original regressor
                    if no == 1
                        rsdmmatrix(:, pidx) = [];
                    end

                    % still remove (almost) empty regressors
                    rsdmmatrixs = (sum(abs(rsdmmatrix), 1) < 1e-6);
                    if rsdmmatrixs(end)
                        continue;
                    end
                    rsdmmatrix(:, rsdmmatrixs) = [];

                    % then regress
                    betas = calcbetas(rsdmmatrix, vtcdata, 1)';
                    vtccorr = vtcdata - rsdmmatrix(:, 1:end-1) * betas(1:end-1, :);

                    % smoothing
                    if numel(smk) ~= 1
                        vtccorr = flexinterpn(vtccorr, [Inf, Inf; 1, 1; 1, 1; vtccsz], ...
                            {smk, [0; 1; 0]}, {1, 1});
                    end

                % fully-naive
                elseif opts.naive

                    % just copy data
                    vtccorr = vtcdata;

                    % smoothing
                    if numel(smk) ~= 1
                        vtccorr = flexinterpn(vtccorr, [Inf, Inf; 1, 1; 1, 1; vtccsz], ...
                            {smk, [0; 1; 0]}, {1, 1});
                    end
                end

                % sample baseline window
                obwinmin = o(oc) + bwinmin / vtctr;
                obwinmax = o(oc) + bwinmax / vtctr;
                basewin = [Inf, Inf; obwinmin, 1; obwinstep, 1; obwinmax, vtccsz(2)];
                obwin = meannoinfnan(flexinterpn(vtccorr, basewin, ipk{:}), 1, true);

                % sample onset window
                awinmin = o(oc) + avgwinmin / vtctr;
                awinmax = o(oc) + avgwinmax / vtctr;
                awin = [Inf, Inf; awinmin, 1; awinstep, 1; awinmax, vtccsz(2)];
                avgwin = flexinterpn(vtccorr, awin, ipk{:}) - repmat(obwin, dima, 1);
                tcwin = max(0, flexinterpn(tw, awin(:, 1), ipk{:}));

                % add to VTCData
                if opts.robust
                    wdata = max(0, flexinterpn(wo, awin, ipk{:}));
                else
                    wdata = repmat(tcwin, 1, vtccsz(2));
                end
                if opts.ffx
                    ftally{cc} = wvartally(ftally{cc}, avgwin, wdata);
                end
                if opts.rfx || opts.wrfx
                    tally{cc} = wvartally(tally{cc}, avgwin, wdata);
                end

                % update counter
                if ~isempty(pbar)
                    consets = consets + 1;
                    xprogress(pbar, consets, sprintf( ...
                        'Added data from onset %d of %d (condition %s, volumes %.2f++).', ...
                        consets, nonsets, conds{cc}, o(oc)));
                end
            end % onsets

            % keep track of removed onsets for single-trial designs
            if opts.rsngtrial
                sonscount(cc) = sonscount(cc) + onsrem(vtci, cc);
            end

        end % conditions
    end % runs
    
    % tally
    if opts.rfx || opts.wrfx
        for cc = 1:dimc
            if opts.wrfx
                [gm, gv] = wvartally(tally{cc});
                gv = 1 ./ gv;
                gv(isinf(gv) | isnan(gv)) = 0;
            else
                gm = wvartally(tally{cc});
            end
            if numel(gm) == 1 && isnan(gm)
                continue;
            end
            if opts.rfx
                rtally{cc} = wvartally(rtally{cc}, gm);
            end
            if opts.wrfx
                wtally{cc} = wvartally(wtally{cc}, gm, 1 ./ gv);
            end
        end
    end
end

% compute required summary measures
for cc = 1:dimc

    % first index
    tci = 1 + 2 * dima * (cc - 1);
    tce = tci + dima - 1;
    sci = tci + dima;
    sce = tce + dima;
    
    % evaluate tallies (patching DF to nominal!)
    if opts.ffx
        [tmean, tvar, terror] = wvartally(ftally{cc});
        tvar = (isnan(tmean) | isinf(terror) | isnan(terror));
        tmean(tvar) = 0;
        terror(tvar) = 1;
        tvar = abs(tmean ./ terror);
        tvar(tvar == 0) = 0;
        vtc.C.VTCData(tci:tce, mask) = tmean;
        vtc.C.VTCData(sci:sce, mask) = tvar;
        vtc.C.RunTimeVars.NrOfConditionOnsets = onscount;
        vtc.C.RunTimeVars.Map(cc).Name = [vtc.C.RunTimeVars.ConditionNames{cc} ' (FFX)'];
        vtc.C.RunTimeVars.Map(cc).LowerThreshold = vtc.C.RunTimeVars.ConditionThresholds(cc, 1, 1);
        vtc.C.RunTimeVars.Map(cc).UpperThreshold = vtc.C.RunTimeVars.ConditionThresholds(cc, 1, 2);
        vtc.C.RunTimeVars.Map(cc).RGBLowerThreshPos = double(vtc.C.RunTimeVars.ConditionColors(cc, 1:3));
        vtc.C.RunTimeVars.Map(cc).RGBUpperThreshPos = double(vtc.C.RunTimeVars.ConditionColors(cc, 5:7));
        vtc.C.RunTimeVars.Map(cc).RGBLowerThreshNeg = double(vtc.C.RunTimeVars.ConditionColors(cc, 9:11));
        vtc.C.RunTimeVars.Map(cc).RGBUpperThreshNeg = double(vtc.C.RunTimeVars.ConditionColors(cc, 13:15));
        vtc.C.RunTimeVars.Map(cc).UseRGBColor = 1;
        vtc.C.RunTimeVars.Map(cc).DF1 = onscount(cc) - 1;
    end
    if opts.rfx
        [tmean, tvar, terror] = wvartally(rtally{cc});
        tvar = (isnan(tmean) | isinf(terror) | isnan(terror));
        tmean(tvar) = 0;
        terror(tvar) = 1;
        tvar = abs(tmean ./ terror);
        rvtc.C.VTCData(tci:tce, mask) = tmean;
        rvtc.C.VTCData(sci:sce, mask) = tvar;
        rvtc.C.RunTimeVars.Map(cc).Name = [vtc.C.RunTimeVars.ConditionNames{cc} ' (RFX)'];
        rvtc.C.RunTimeVars.Map(cc).LowerThreshold = vtc.C.RunTimeVars.ConditionThresholds(cc, 1, 1);
        rvtc.C.RunTimeVars.Map(cc).UpperThreshold = vtc.C.RunTimeVars.ConditionThresholds(cc, 1, 2);
        rvtc.C.RunTimeVars.Map(cc).RGBLowerThreshPos = double(vtc.C.RunTimeVars.ConditionColors(cc, 1:3));
        rvtc.C.RunTimeVars.Map(cc).RGBUpperThreshPos = double(vtc.C.RunTimeVars.ConditionColors(cc, 5:7));
        rvtc.C.RunTimeVars.Map(cc).RGBLowerThreshNeg = double(vtc.C.RunTimeVars.ConditionColors(cc, 9:11));
        rvtc.C.RunTimeVars.Map(cc).RGBUpperThreshNeg = double(vtc.C.RunTimeVars.ConditionColors(cc, 13:15));
        rvtc.C.RunTimeVars.Map(cc).UseRGBColor = 1;
        rvtc.C.RunTimeVars.Map(cc).DF1 = numsubs - 1;
    end
    if opts.wrfx
        [tmean, tvar, terror, tdf] = wvartally(wtally{cc});
        tvar = (isnan(tmean) | isinf(terror) | isnan(terror));
        tmean(tvar) = 0;
        terror(tvar) = 1;
        tvar = abs(sdist('tinv', sdist('tcdf', -abs(tmean ./ terror), max(sqrt(eps), tdf - 1)), numsubs - 1));
        tvar(tvar == 0) = 0;
        wvtc.C.VTCData(tci:tce, mask) = tmean;
        wvtc.C.VTCData(sci:sce, mask) = tvar;
        wvtc.C.RunTimeVars.Map(cc).Name = [vtc.C.RunTimeVars.ConditionNames{cc} ' (weighted RFX)'];
        wvtc.C.RunTimeVars.Map(cc).LowerThreshold = vtc.C.RunTimeVars.ConditionThresholds(cc, 1, 1);
        wvtc.C.RunTimeVars.Map(cc).UpperThreshold = vtc.C.RunTimeVars.ConditionThresholds(cc, 1, 2);
        wvtc.C.RunTimeVars.Map(cc).RGBLowerThreshPos = double(vtc.C.RunTimeVars.ConditionColors(cc, 1:3));
        wvtc.C.RunTimeVars.Map(cc).RGBUpperThreshPos = double(vtc.C.RunTimeVars.ConditionColors(cc, 5:7));
        wvtc.C.RunTimeVars.Map(cc).RGBLowerThreshNeg = double(vtc.C.RunTimeVars.ConditionColors(cc, 9:11));
        wvtc.C.RunTimeVars.Map(cc).RGBUpperThreshNeg = double(vtc.C.RunTimeVars.ConditionColors(cc, 13:15));
        wvtc.C.RunTimeVars.Map(cc).UseRGBColor = 1;
        wvtc.C.RunTimeVars.Map(cc).DF1 = numsubs - 1;
    end
end

% update scaling window
if opts.wrfx
    aft_SetScalingWindow(wvtc, [-tlim, tlim], true);
end
if opts.rfx
    aft_SetScalingWindow(rvtc, [-tlim, tlim], true);
else
    if opts.wrfx
        rvtc = wvtc;
    end
end
if opts.ffx
    aft_SetScalingWindow(vtc, [-tlim, tlim], true);
else
    aft_ClearObject(vtc);
    if opts.rfx
        vtc = rvtc;
        if opts.wrfx
            rvtc = wvtc;
        end
    else
        vtc = wvtc;
    end
end

% close bar
if ~isempty(pbar)
    closebar(pbar);
end

% clear interim data
clearxffobjects(vtco);

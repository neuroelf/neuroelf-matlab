function [map, midx, midxdist] = glm_Searchlight(xo, opts)
% GLM::Searchlight  - run a searchlight analysis
%
% FORMAT:       [map, midx, midxdist] = glm.Searchlight(opts);
%
% Input fields:
%
%       opts        mandatory (and optional) settings
%        .condrep   Rx2 replacements for regexprep (e.g. to collapse, etc)
%                   default: {'_T\d+$', ''}
%        .condsel   condition selection (regexp, default: if condrep is
%                   given, use second column, otherwise all multi-trial)
%        .debug     flag, if true, print debug messages
%        .dist      distance (in mm) for which searchlight is created (6)
%        .dweight   distance weight function; either vector or handle, if
%                   a vector is given, the values will be interpolated
%                   using dweight(1 + 10 * dist), default []
%        .gfunc     cell array with first cell set to function handle
%                   being evaluated to compute a group statistic
%                   along subjects, string argument replacements:
%                   '$Cov'     -> covariate Cov values (for subjects)
%                   'data'     -> subjects-by-slstored array (see slstore)
%                   'data\d+'  -> subjects-by-1 array
%                   'sldata\d' -> subjects-by-slreturn array
%                   'subjects' -> cell array with subject names
%                   'subgroup' -> if RunTimeVars.Groups is defined, a
%                                 double array with group assignments
%                   default: {@neuroelf.slttest, 'data1'}
%        .gmdf      group map DF (default: number of selected subjects - 1)
%        .gmtype    group map type (default: 1)
%        .gonly     directly pass data from conditions into gfunc (false)
%        .gstore    output storage (group level), cell array with index
%                   into outputs, within output (numeric) index, and name
%                   default: {1, 1, 'Group: condition-difference Searchlight'}
%        .maponly   only create the SMP/VMP (with all maps) and return
%        .mask      mask object (or matching array)
%        .maxtime   maximum time (in seconds) before VMP object is stored,
%                   set this to a value a few minutes below any hard cutoff
%                   in case partial results are preferable to no results
%        .midx      cell array with target -> searchlight mappings
%        .midxdist  cell array with distances of mapped indices (required
%                   if dist-weighting is enabled)
%        .midxonly  create searchlight mappings and return as only outputs
%        .minnfeat  minimum number of features (default: 6)
%        .mlevel    multi-level (pass values from all subjects into func,
%                   and skip gfunc, directly store into group map)
%        .nullcorr  correct samples by median of NULL (default: false)
%        .nulld     number of null-distribution samples, default: 0
%        .nullt     use time courses+models to generate NULL samples by
%                   applying phase+label scrambling (instead of randn, false) 
%        .nullz     recode first output into Z value from NULL (default:
%                   true only if nulld > 500)
%        .progress  either {false} or a 1x1 xfigure::progress or xprogress
%        .recode    recode condition numbers from 1 to selected (true)
%        .remcmean  for trial-based data, remove per-condition mean (false)
%        .slfnout   number of output arguments to expect (default: 1)
%        .slfunc    cell array with function handle evaluated for each
%                   searchlight, in which the following string (1xN char)
%                   arguments will be replaced as follows:
%                   'cond'    -> double array with unique condition number
%                   'data'    -> condition-by-features array
%                   'datat'   -> features-by-condition array
%                   'opts'    -> pass-through of opts from this call
%                   'run'     -> column with run number of condition
%                   'subject' -> double array with subject index (mlevel)
%                   default: {@neuroelf.slsvmclassify, 'data', 'cond'}
%        .slstore   output storage (lower level), cell array with index
%                   into outputs, within output (numeric) index, and name
%                   default: {1, 1, 'Condition-difference Searchlight'}
%        .sortcond  sort conditions by name (for coding, true)
%        .srf       surface object (required if midx/midxdist not given!)
%        .subpart   1x2 double indicating which sub-partition of the space
%                   is to be processed (subpart(1) of subpart(2)), [1, 1]
%        .subsel    subject selection (default: all)
%        .t2sl      cell array, pass results from timefunc into slfunc,
%                   e.g. {'data', 'cond'} means tfunc must produce two
%                   outputs (data and cond) which will be passed to slfunc
%                   default {} (means: do not pass on to slfunc!)
%        .targetmap filename of target map (if given, will be loaded using
%                   transio, and results will be written into this file!)
%                   or existing SMP/VMP (if mapnames don't exist, add)
%        .tcclean   remove nuisance from time courses (default: true), if
%                   a char array, combination from letters
%                   'd' - derivatives of nuisance (motion) regressors
%                   'f' - temporal filters
%                   'g' - global signal regressors
%                   'm' - motion parameters
%                   'o' - other nuisance regressors
%        .tctrans   override time-course transform: 'none', 'psc', 'z'
%        .tfunc     function handle evaluated for extracted time courses
%                   from files, replacements:
%                   'cond\d+'  -> prt.Cond(\d+).OnOffsets
%                   'data'     -> time-by-features array
%                   'datat'    -> features-by-time array
%                   'prtcont'  -> full protocol contents
%                   'run'      -> Tx1 double vector, run number
%                   'subject'  -> Tx1 double vector, subject number
%                   'tr'       -> TR
%        .timeml    if true, extract time courses from all subjects
%        .transio   if true, do not resolve transio betamaps (slow, false)
%
% Output fields:
%
%       map         output object with map(s) (SMP or VMP)
%       midx        cell array with target -> searchlight mappings
%       midxdist    cell array with target -> searchlight mapping distances
%
% Using: flexinterpn, multimatch, sltargets (and others).

% Version:  v1.1
% Build:    17050418
% Date:     May-04 2017, 6:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, 2017, Jochen Weber
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
flexinterpn  = ne_methods.flexinterpn;

% check object argument
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end

% get object contents
bc = xo.C;
rtv = bc.RunTimeVars;

% check options (must be 1x1 struct)
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end

% and first, check whether debug messages are requested
if ~isfield(opts, 'debug') || ~islogical(opts.debug) || numel(opts.debug) ~= 1
    dodebug = false;
else
    dodebug = opts.debug;
    if dodebug
        dbnext = now + 1 / 8640;
    end
end

% get data from GLM
if dodebug
    fprintf('Parsing options...\n');
end

% get layout (dims, etc.), subject IDs, and predictor names
glmlay = aft_Layout(xo);
subids = glm_Subjects(xo);
stsubids = glm_Subjects(xo, true);
nsubs = numel(subids);
subpreds = glm_SubjectPredictors(xo);
studies = bc.Study;
stfiles = {studies.NameOfAnalyzedFile};

% covariates
if isfield(rtv, 'CovariatesData') && isfield(rtv, 'CovariatesNames') && ...
    isa(rtv.CovariatesData, 'double') && ndims(rtv.CovariatesData) == 2 && ...
    iscell(rtv.CovariatesNames) && numel(rtv.CovariatesNames) == size(rtv.CovariatesData, 2) && ...
    size(rtv.CovariatesData, 1) == nsubs
    covnames = rtv.CovariatesNames(:);
    covdata = rtv.CovariatesData;
else
    covnames = {};
    covdata = zeros(nsubs, 0);
end

% initialize time course objects array (for clear calls)
stobjs = {};

% get map size depending on RFX
isrfx = (bc.ProjectTypeRFX > 0);
if isrfx
    msize = size(bc.GLMData.RFXGlobalMap);
else
    msize = size(bc.GLMData.TimeCourseMean);
end

% number of map elements, and dims, as well as multi-map dim
mnumel = prod(msize);
mdim = sum(msize > 1);
cdim = mdim + 1;

% subscript access
sra = [repmat({':'}, 1, mdim), {[]}];

% now to options... -> by default, replace _T### with _T
if ~isfield(opts, 'condrep') || ~iscell(opts.condrep) || size(opts.condrep, 2) ~= 2
    opts.condrep = {'_T\d+$', '_T'};
end

% replace condition names (in a copy)
reppreds = subpreds;
for rc = 1:size(opts.condrep, 1)
    if ischar(opts.condrep{rc, 1}) && ischar(opts.condrep{rc, 2}) && ~isempty(opts.condrep{rc, 1})
        reppreds = regexprep(reppreds, opts.condrep{rc, 1}(:)', opts.condrep{rc, 2}(:)');
    end
end

% keep the copy (to locate the indices later)
oreppreds = reppreds;

% and find the unique names (to match against condsel)
[uconds, ucondi] = unique(reppreds);

% no selection was made
if ~isfield(opts, 'condsel') || (~ischar(opts.condsel) && ~iscell(opts.condsel)) || isempty(opts.condsel)

    % remove the unique elements from copy of condition names
    oreppreds(ucondi) = [];

    % get the unique of the remaining ones (those are multiple trials!)
    oreppreds = unique(oreppreds);

    % if not empty, the condition selection is than this selection
    if ~isempty(oreppreds)
        opts.condsel = oreppreds;

    % otherwise we use a catch all
    else
        opts.condsel = '^.*$';
    end
end

% for a character selection
if ischar(opts.condsel)

    % find all matching conditions
    selconds = find(~cellfun('isempty', regexpi(reppreds, opts.condsel)));

    % and their labels
    if ~isfield(opts, 'sortcond') || ~islogical(opts.sortcond) || ...
        numel(opts.sortcond) ~= 1 || opts.sortcond
        selcondi = ne_methods.multimatch(reppreds(selconds), uconds);
    else
        selcondi = unique(ne_methods.multimatch(reppreds(selconds), reppreds));
    end

% for a cell array (multiple conditions selected)
else

    % copy the selection
    selconds = opts.condsel;
    selcondi = cell(size(selconds));

    % iterate
    for rc = 1:numel(selconds)

        % find matches
        try
            selconds{rc} = find(~cellfun('isempty', regexpi(reppreds, selconds{rc})));
            if isempty(selconds{rc})
                error('neuroelf:xff:badConditionSelection', ...
                    'No conditions found for ''%s''.', opts.condsel{rc});
            end
        catch xfferror
            rethrow(xfferror);
        end
        selcondi{rc} = rc .* ones(numel(selconds{rc}), 1);
    end

    % catenate
    selconds = ne_methods.lsqueezecells(selconds);
    selconds = cat(1, selconds{:});
    if numel(selconds) ~= numel(unique(selconds))
        error('neuroelf:xff:badConditionSelection', 'Conditions not selected uniquely.');
    end
    selcondi = cat(1, selcondi{:});
end

% re-code from 1 to number of selected conditions
if ~isfield(opts, 'recode') || ~islogical(opts.recode) || numel(opts.recode) ~= 1 || opts.recode
    rselc = zeros(max(selcondi), 1);
    rselc(unique(selcondi)) = 1:numel(unique(selcondi));
    selcondi = rselc(selcondi);
end

% total number of selected "conditions" (maps per subject to sample from)
numconds = numel(selconds);
selpreds = subpreds(selconds);

% distance from searchlight center (in mm)
if ~isfield(opts, 'dist') || ~isa(opts.dist, 'double') || numel(opts.dist) ~= 1 || ...
    isinf(opts.dist) || isnan(opts.dist) || opts.dist <= 0
    opts.dist = 6;

% restrict to 30mm
else
    opts.dist = min(30, opts.dist);
end

% distance weight (array or function), default []
if ~isfield(opts, 'dweight') || isempty(opts.dweight) || (~isa(opts.dweight, 'double') && ~isa(opts.dweight, 'function_handle'))
    dweight = [];
    doweight = 0;

% array given
elseif isa(opts.dweight, 'double')

    % linearize
    dweight = opts.dweight(:);

    % use weights flag (1 = values)
    doweight = 1;

    % but only if at least 3 valid values given
    if numel(dweight) < 3 || any(isinf(dweight) | isnan(dweight))
        dweight = [];
        doweight = 0;
    end

    % fill up with one 0 value (for correct interpolation at the end)
    dweight(end+1) = 0;

% function handle
else

    % store and flag
    dweight = opts.dweight;
    doweight = 2;

    % then test
    try
        if dodebug
            fprintf('Testing weighting function...\n');
        end
        wsampl = feval(dweight, 3 .* abs(randn(30, 1)));
        if ~isa(wsample, 'double') || numel(wsampl) ~= 30 || any(isinf(wsampl(:)) | isnan(wsampl(:)))
            error('Bad return value.');
        end
    catch xfferror
        error('neuroelf:xff:fevalError', ['Error sampling weighting function: ' xfferror.message]);
    end
end

% group-function outputs (default: 1)
if ~isfield(opts, 'gfnout') || ~isa(opts.gfnout, 'double') || numel(opts.gfnout) ~= 1 || ...
    isinf(opts.gfnout) || isnan(opts.gfnout) || opts.gfnout < 1
    gfnout = 1;
else
    gfnout = ceil(opts.gfnout);
end

% group-based function
if ~isfield(opts, 'gfunc') || ~iscell(opts.gfunc) || isempty(opts.gfunc) || ...
   ~isa(opts.gfunc{1}, 'function_handle')

    % default: one-sample t-test of first searchlight output
    gfunc = {ne_methods.slttest, 'data1'};

% if given, linearize
else
    gfunc = opts.gfunc(:)';
end

% find which cells to replace with what arguments
gfrep = find(~cellfun('isempty', gfunc) & cellfun(@ischar, gfunc));
gfrep(cellfun('isempty', regexpi(gfunc(gfrep), ...
    '^(\$\(?[a-z0-9_\.\*\+\?\|]+\)?|data(\d+[\,0-9]*)?|opts|sldata\d|subjects|subgroup)$'))) = [];

% covariates
gfrcovs = gfrep(~cellfun('isempty', regexpi(gfunc(gfrep), '^\$\(?[a-z0-9_\.\*\+\?\|]+\)?$'))) - 1;
agfrcovs = any(gfrcovs);
if agfrcovs

    % error out of no covariates available
    if isempty(covnames)
        error('neuroelf:xff:noCovariates', 'No Covariates configured in GLM.');
    end

    % find indices
    for rc = 1:numel(gfrcovs)

        % try full match first
        gfrcovit = find(strcmpi(covnames, gfunc{gfrcovs(rc)+1}(2:end)));

        % more than one match?
        if numel(gfrcovit) > 1

            % invalid selection then (need to recode in GLM!)
            error('neuroelf:xff:covariateNotFound', ...
                'Covariate %s not available.', gfunc{gfrcovs(rc)+1}(2:end));

        % not found
        elseif isempty(gfrcovit)

            % try regexp instead
            gfrcpart = gfunc{gfrcovs(rc)+1}(2:end);
            gfrcovit = find(~cellfun('isempty', regexpi(covnames, gfrcpart)));

            % still not found
            if isempty(gfrcovit)
                error('neuroelf:xff:covariateNotFound', ...
                    'Covariate %s not available.', gfunc{gfrcovs(rc)+1}(2:end));
            end

            % but then only take the first if a regular string was used
            if ~any(gfrcpart == '*' | gfrcpart == '+' | gfrcpart == '?' | gfrcpart == '.' | gfrcpart == '|')
                gfrcovit = gfrcovit(1);
            end
        end

        % set in array (fix replace!)
        gfunc{gfrcovs(rc)+1} = covdata(:, gfrcovit);
    end
end

% all data (group function replacement)
gfrdata = gfrep(strcmpi(gfunc(gfrep), 'data')) - 1;
agfrdata = any(gfrdata);

% specific data (data\d+)
gfrdatan = gfrep(~cellfun('isempty', regexpi(gfunc(gfrep), '^data\d+'))) - 1;
agfrdatan = any(gfrdatan);
if agfrdatan
    gfrdatai = cellfun(ne_methods.u8str2double, strrep(gfunc(gfrdatan + 1), 'data', ''), ...
        'UniformOutput', false);
end

% options, subjects, and groups
gfropts = gfrep(strcmpi(gfunc(gfrep), 'opts')) - 1;
agfropts = any(gfropts);
gfrsubg = gfrep(strcmpi(gfunc(gfrep), 'subgroup')) - 1;
agfrsubg = any(gfrsubg);
gfrsubs = gfrep(strcmpi(gfunc(gfrep), 'subjects')) - 1;
agfrsubs = any(gfrsubs);

% full searchlight data (sldata\d)
gfrsldatan = gfrep(~cellfun('isempty', regexpi(gfunc(gfrep), '^sldata\d+'))) - 1;
agfrsldatan = any(gfrdatan);
if agfrsldatan
    gfrsldatai = str2double(strrep(gfunc(gfrsldatan + 1), 'sldata', ''));
end

% what to store from group function (default, to also allow no storage!)
if ~isfield(opts, 'gstore') || ~iscell(opts.gstore) || size(opts.gstore, 2) ~= 3 || isempty(opts.sltore)
    gstore = {1, 1, 'Group: Condition difference Searchlight'};

% if defined, limit to 2D
else
    gstore = opts.gstore(:, :, 1);
end

% test if all information matches
for fc = 1:size(gstore, 1)
    if ~isa(gstore{fc, 1}, 'double') || numel(gstore{fc, 1}) ~= 1 || ...
        isinf(gstore{fc, 1}) || isnan(gstore{fc, 1}) || ~any(gstore{fc, 1} == 1:gfnout) || ...
       ~isa(gstore{fc, 2}, 'double') || isempty(gstore{fc, 2}) || ...
        any(isinf(gstore{fc, 2}(:)) | isnan(gstore{fc, 2}(:)) | gstore{fc, 2}(:) < 1) || ...
       ~ischar(gstore{fc, 3}) || isempty(gstore{fc, 3})
        error('neuroelf:xff:badArgument', 'Bad gstore option.');
    else
        gstore{fc, 2} = ceil(gstore{fc, 2}(:)');
        gstore{fc, 3} = gstore{fc, 3}(:)';
        if numel(gstore{fc, 2}) > 1 && isempty(strfind(gstore{fc, 3}, '%d'))
            error('neuroelf:xff:badArgument', 'Bad gstore option.');
        end
    end
end

% concatenate and determine source indices (arguments, etc.)
gnstore = cat(2, gstore{:, 2});
gnstoren = zeros(size(gnstore));
tgnstore = 0;
for fc = 1:size(gstore, 1)
    gnstoren(tgnstore+1:tgnstore+numel(gstore{fc, 2})) = gstore{fc, 1};
    tgnstore = tgnstore + numel(gstore{fc, 2});
end

% group map type
if ~isfield(opts, 'gmtype') || ~isa(opts.gmtype, 'double') || ...
   ~any(numel(opts.gmtype) == [1, tgnstore]) || ...
    any(isinf(opts.gmtype(:)) | isnan(opts.gmtype(:)) | opts.gmtype(:) < 1)
    gmtype = ones(1, tgnstore);
elseif numel(opts.gmtype) == 1
    gmtype = opts.gmtype .* ones(1, tgnstore);
else
    gmtype = opts.gmtype(:);
end
gmtype = round(gmtype);

% VMP only requested (default false)
if ~isfield(opts, 'maponly') || ~islogical(opts.maponly) || numel(opts.maponly) ~= 1
    opts.maponly = false;
end

% no mask given
if ~isfield(opts, 'mask') || ...
    ((numel(opts.mask) ~= 1 || ~xffisobject(opts.mask, true, 'msk')) && ...
     (~isequal(size(opts.mask), msize) || ~isnumeric(opts.mask)))

    % already set from previous run
    if isfield(rtv, 'SearchlightMask') && islogical(rtv.SearchlightMask) && isequal(size(rtv.SearchlightMask), msize)
        opts.mask = rtv.SearchlightMask;
    else

        % set to []
        opts.mask = [];
    end

% MSK object given
elseif numel(opts.mask) == 1

    % get mask
    opts.mask = (opts.mask.C.Mask > 0);

    % and check size
    if ~isequal(size(opts.mask), msize)

        % error out
        error('neuroelf:xff:sizeMismatch', 'MSK object mismatches in size.');
    end
end

% adhoc mask (to avoid running across ALL voxels
if isempty(opts.mask)
    if dodebug
        fprintf('Creating ad-hoc mask...\n');
    end

    % for surface-based stats
    if bc.ProjectType == 2

        % simply use all data
        opts.mask = true(msize);

    % for RFX (non-surface stats
    elseif bc.ProjectTypeRFX > 0

        % start with a all-0 mask
        opts.mask = zeros(msize);

        % then iterate over subjects (all subjects for mask!)
        for sc = 1:numel(bc.GLMData.Subject);

            % and count all ~(Inf | Nan | 0) voxels (in any map)
            bvals = bc.GLMData.Subject(sc).BetaMaps;
            if istransio(bvals)
                bvals = resolve(bvals);
            end
            opts.mask = opts.mask + ~any(isinf(bvals) | isnan(bvals) | bvals == 0, cdim);
        end

        % finally threshold
        opts.mask = (opts.mask > (0.8 * numel(bc.GLMData.Subject)));

    % FFX mask
    else

        % start with a all-0 mask
        opts.mask = zeros(msize);

        % then iterate over subjects (all subjects for mask!)
        for sc = 1:bc.NrOfSubjects

            % and count all ~(Inf | Nan | 0) voxels (in any map)
            sra{end} = (1+(sc-1)*bc.NrOfSubjectsPredictors):(sc*bc.NrOfSubjectPredictors);
            bvals = bc.GLMData.BetaMaps(sra{:});
            opts.mask = opts.mask + ~any(isinf(bvals) | isnan(bvals) | bvals == 0, cdim);
        end

        % finally threshold
        opts.mask = (opts.mask > (0.8 * bc.NrOfSubjects));
    end

    % set (but don't save to disk)
    if any(opts.mask(:))
        xo.C.RunTimeVars.SearchlightMask = opts.mask;
    end
end

% number of mask candidates
nmask = numel(opts.mask);
if sum(opts.mask(:)) < 27
    error('neuroelf:xff:tooFewVoxels', 'Too few voxels in mask.');
end

% maximum run-time (mainly useful for running on a cluster)
if ~isfield(opts, 'maxtime') || ~isa(opts.maxtime, 'double') || numel(opts.maxtime) ~= 1

    % default: no restriction
    maxtime = Inf;
else

    % at least allow 2 minutes
    maxtime = max(120, opts.maxtime);

    % and inform
    if dodebug
        fprintf('Limiting runtime to %02d:%02d.%03d minutes.\n', ...
            floor(maxtime / 60), floor(mod(maxtime, 60)), maxtime - floor(maxtime));
    end

    % determine endtime
    maxtime = now + maxtime / 86400;
end

% mask-based indices and distances
midxdist = [];
if ~isfield(opts, 'midx') || ~iscell(opts.midx) || numel(opts.midx) ~= nmask

    % if not given
    if dodebug
        fprintf('Looking up searchlight targets...\n');
    end

    % for now, use either sltargets (@neuroelf) or for srf, require SRF
    if mdim == 3
        [midx, midxdist] = ne_methods.sltargets(opts.mask, struct( ...
            'dist', opts.dist, 'res', bc.Resolution .* [1, 1, 1]));
    elseif ~isfield(opts, 'srf') || numel(opts.srf) ~= 1 || ~xffisobject(opts.srf, true, 'srf') || ...
        size(opts.srf.C.Neighbors, 1) ~= bc.NrOfVertices
        error('neuroelf:xff:missingOption', 'Surface-based Searchlight requires SRF.');

    % for SRF
    else

        % use NeighborsNDegree
        [midx, fc, midxdist] = srf_NeighborsNDegree(opts.srf, floor(opts.dist / sqrt(srf_Area(opts.srf) / (pi * opts.srf.C.NrOfVertices))));

        % compute which are to be kept
        fc = cellfun(@le, midxdist, repmat({opts.dist}, size(midxdist)), 'UniformOutput', false);

        % and keep only those
        for fci = 1:numel(fc)
            midx{fci} = midx{fci}(fc{fci});
            midxdist{fci} = midxdist{fci}(fc{fci});
        end
    end

    % midx and midxdist only?
    if isfield(opts, 'midxonly') && islogical(opts.midxonly) && numel(opts.midxonly) == 1 && opts.midxonly
        map = midx;
        midx = midxdist;
        return;
    end

% linearize
else
    midx = opts.midx(:);
end

% number of sources per searchlight target
numidx = cellfun('prodofsize', midx);
maxnumidx = max(numidx);

% remove from mask
opts.mask(numidx == 0) = false;

% (how many) voxels within mask
xmask = find(opts.mask(:));
smask = numel(xmask);
if smask == 0
    error('neuroelf:xff:emptyMask', 'Nothing selected.');
end

% distance metric for searchlight sources
if isfield(opts, 'midxdist') && iscell(opts.midxdist) && isequal(size(opts.midxdist), size(midx)) && ...
    all(cellfun('prodofsize', opts.midxdist(:)) == cellfun('prodofsize', midx))
    midxdist = opts.midxdist;

% complain if not given (or invalid), but required for weighting
elseif doweight
    error('neuroelf:xff:missingOption', 'Option midxdist must be given for distance weighting.');
end

% minimum number of features
if ~isfield(opts, 'minnfeat') || ~isa(opts.minnfeat, 'double') || numel(opts.minnfeat) ~= 1 || ...
    isinf(opts.minnfeat) || isnan(opts.minnfeat) || opts.minnfeat < 1
    minnfeat = min(6, maxnumidx);
else
    minnfeat = min(round(opts.minnfeat), maxnumidx);
end

% multi-level
if ~isfield(opts, 'mlevel') || ~islogical(opts.mlevel) || numel(opts.mlevel) ~= 1
    mlevel = false;
else
    mlevel = opts.mlevel;
end

% correct for, and how many NULL samples
if isfield(opts, 'nullcorr') && islogical(opts.nullcorr) && numel(opts.nullcorr) == 1
    nullcorr = opts.nullcorr;
else
    nullcorr = false;
end
if isfield(opts, 'nulld') && isa(opts.nulld, 'double') && numel(opts.nulld) == 1 && ...
   ~isinf(opts.nulld) && ~isnan(opts.nulld) && opts.nulld >= 0
    nulld = min(100000, round(opts.nulld));
    if nulld < 100
        nullcorr = false;
    end
else
    nulld = 0;
    nullcorr = false;
end

% progress bar
pbar = [];
if ~isfield(opts, 'progress') || numel(opts.progress) ~= 1 || ...
   (~islogical(opts.progress) && ~isa(opts.progress, 'xprogress') && ~isxfigure(opts.progress, true))
    opts.progress = false;
    progress = false;
elseif ~islogical(opts.progress)
    progress = true;
    pbar = opts.progress;
else
    progress = opts.progress;
end

% remove condition mean?
if ~isfield(opts, 'remcmean') || ~islogical(opts.remcmean) || numel(opts.remcmean) ~= 1
    opts.remcmean = false;
end

% number of outputs for searchlight function (calls regardless of storage!)
if ~isfield(opts, 'slfnout') || ~isa(opts.slfnout, 'double') || numel(opts.slfnout) ~= 1 || ...
    isinf(opts.slfnout) || isnan(opts.slfnout) || opts.slfnout < 1
    slfnout = 1;
else
    slfnout = ceil(opts.slfnout);
end

% regular (subject-level) searchlight function
if ~isfield(opts, 'slfunc') || ~iscell(opts.slfunc) || isempty(opts.slfunc) || ...
   ~isa(opts.slfunc{1}, 'function_handle')

    % default: SVM classifier with data and condition labels
    slfunc = {ne_methods.slsvmclassify, 'data', 'cond'};

% otherwise, use function as given
else
    slfunc = opts.slfunc(:)';
end

% find indices for replacement
slfrep = find(~cellfun('isempty', slfunc) & cellfun(@ischar, slfunc));
slfrep(cellfun('isempty', regexpi(slfunc(slfrep), '^(cond|datat?|opts|run|subjects)$'))) = [];
slfrcond = slfrep(strcmpi(slfunc(slfrep), 'cond')) - 1;
aslfrcond = any(slfrcond);
slfrdata = slfrep(strcmpi(slfunc(slfrep), 'data')) - 1;
aslfrdata = any(slfrdata);
slfrdatat = slfrep(strcmpi(slfunc(slfrep), 'datat')) - 1;
aslfrdatat = any(slfrdatat);
slfropts = slfrep(strcmpi(slfunc(slfrep), 'opts')) - 1;
aslfropts = any(slfropts);
slfrrun = slfrep(strcmpi(slfunc(slfrep), 'run')) - 1;
aslfrrun = any(slfrrun);
slfrsubj = slfrep(strcmpi(slfunc(slfrep), 'subject')) - 1;
aslfrsubj = any(slfrsubj);

% what to store from the output of the searchlight function
if ~isfield(opts, 'slstore') || ~iscell(opts.slstore) || ...
    size(opts.slstore, 2) ~= 3 || isempty(opts.slstore)
    slstore = {1, 1, 'Condition difference Searchlight'};
else
    slstore = opts.slstore(:, :, 1);
end

% test storage requests
for fc = 1:size(slstore, 1)
    if ~isa(slstore{fc, 1}, 'double') || numel(slstore{fc, 1}) ~= 1 || ...
        isinf(slstore{fc, 1}) || isnan(slstore{fc, 1}) || ~any(slstore{fc, 1} == 1:slfnout) || ...
       ~isa(slstore{fc, 2}, 'double') || isempty(slstore{fc, 2}) || ...
        any(isinf(slstore{fc, 2}(:)) | isnan(slstore{fc, 2}(:)) | slstore{fc, 2}(:) < 1) || ...
       ~ischar(slstore{fc, 3}) || isempty(slstore{fc, 3})
        error('neuroelf:xff:badArgument', 'Bad slstore option.');
    else
        slstore{fc, 2} = ceil(slstore{fc, 2}(:)');
        slstore{fc, 3} = slstore{fc, 3}(:)';
        if numel(slstore{fc, 2}) > 1 && isempty(strfind(slstore{fc, 3}, '%d'))
            error('neuroelf:xff:badArgument', 'Bad slstore option.');
        end
    end
end
nstore = cat(2, slstore{:, 2});
nstoren = zeros(size(nstore));
tnstore = 0;
for fc = 1:size(slstore, 1)
    nstoren(tnstore+1:tnstore+numel(slstore{fc, 2})) = slstore{fc, 1};
    tnstore = tnstore + numel(slstore{fc, 2});
end

% which partition of the problem to solve in this run?
if ~isfield(opts, 'subpart') || ~isa(opts.subpart, 'double') || numel(opts.subpart) ~= 2 || ...
    any(isinf(opts.subpart) | isnan(opts.subpart) | opts.subpart < 1 | opts.subpart > nmask) || ...
    opts.subpart(1) > opts.subpart(2)
    opts.subpart = [1, 1];
else
    opts.subpart = ceil(opts.subpart(:)');
end

% splice the data (across the entire mask to avoid different runtimes)
spi = opts.subpart(1):opts.subpart(2):smask;
nspi = numel(spi);

% subject selection
if isfield(opts, 'subsel') && (iscell(opts.subsel) || isa(opts.subsel, 'double'))

    % cell array
    if iscell(opts.subsel)

        % linearize
        opts.subsel = opts.subsel(:);

        % invalid selection?
        if any(~cellfun(@ischar, opts.subsel) | cellfun('isempty', opts.subsel))
            error('neuroelf:xff:badArgument', 'Invalid subject selection.');
        end

        % match against subject IDs
        opts.subsel = ne_methods.multimatch(opts.subsel, subids);

        % any missing IDs
        if any(opts.subsel < 1)
            error('neuroelf:xff:badArgument', 'Invalid subject selection.');
        end

    % invalid selection (not a double or outside range)
    elseif any(isinf(opts.subsel(:)) | isnan(opts.subsel(:)) | opts.subsel(:) < 1 | opts.subsel(:) > nsubs)
        error('neuroelf:xff:badArgument', 'Invalid subject selection.');

    % linearize and round
    else
        opts.subsel = round(opts.subsel(:));
    end

    % force to unique values
    subsel = unique(opts.subsel);

% by default, use all
else
    subsel = (1:nsubs)';
end

% number of subjects and their IDs (for map names)
nssubs = numel(subsel);
selsubs = subids(subsel);

% apply to covariates
if ~isempty(gfrcovs)
    for fc = 1:numel(gfrcovs)
        gfunc{gfrcovs(rc)+1} = gfunc{gfrcovs(rc)+1}(subsel, :);
    end
end

% clean time courses
if ~isfield(opts, 'tcclean') || ~ischar(opts.tcclean) || isempty(opts.tcclean) || ...
   ~any(lower(opts.tcclean(1) == 'adfgmno'))
    opts.tcclean = 'n';
else
    opts.tcclean = lower(opts.tcclean(:)');
end

% time-course transformation (default: from GLM.TransformationType)
if ~isfield(opts, 'tctrans') || ~ischar(opts.tctrans) || isempty(opts.tctrans) || ...
   ~any(lower(opts.tctrans(1) == 'npz'))
    if bc.TransformationType == 1
        opts.tctrans = 'z';
    elseif bc.TransformationType == 3
        opts.tctrans = 'p';
    else
        opts.tctrans = 'n';
    end
else
    opts.tctrans = lower(opts.tctrans(1));
end

% function to work on time courses (default: none)
if ~isfield(opts, 'tfunc') || ~iscell(opts.tfunc) || isempty(opts.tfunc) || ...
   ~isa(opts.tfunc{1}, 'function_handle')
    tfunc = [];
    usetfunc = false;
else
    tfunc = opts.tfunc(:)';
    usetfunc = true;
end

% if using time course function
if usetfunc
    tfrep = find(~cellfun('isempty', tfunc) & cellfun(@ischar, tfunc));
    tfrep(cellfun('isempty', regexpi(gfunc(tfrep), '^(cond\d*|datat?|files|prtcont|run|opts|subject|tr)$'))) = [];
    tfrcond = gfrep(strcmpi(tfunc(tfrep), 'cond')) - 1;
    tfrcondn = gfrep(~cellfun('isempty', regexpi(tfunc(tfrep), '^cond\d+'))) - 1;
    if ~isempty(tfrcondn)
        tfrcondi = str2double(strrep(tfunc(tfrcondn + 1), 'cond', ''));
    end
    tfrdata = gfrep(strcmpi(tfunc(tfrep), 'data')) - 1;
    tfrdatat = gfrep(strcmpi(tfunc(tfrep), 'datat')) - 1;
    tfrprtc = gfrep(strcmpi(tfunc(tfrep), 'prtcont')) - 1;
    tfrrun = gfrep(strcmpi(tfunc(tfrep), 'run')) - 1;
    tfrsubj = gfrep(strcmpi(tfunc(tfrep), 'subject')) - 1;
    tfrtr = gfrep(strcmpi(tfunc(tfrep), 'tr')) - 1;
    if ~isfield(opts, 't2sl') || ~iscell(opts.t2sl) || isempty(opts.t2sl) || ...
       ~all(cellfun(@ischar, opts.t2sl(:)) & ~cellfun('isempty', opts.t2sl(:)))
        t2sl = {};
    else
        t2sl = opts.t2sl(:)';
        t2slcond = [];
        t2sldata = [];
        t2slsubj = [];
        for fc = 1:numel(t2sl)
            if strcmp(t2sl{fc}, 'cond')
                t2slcond = fc;
            elseif strcmp(t2sl{fc}, 'data')
                t2sldata = fc;
            elseif strcmp(t2sl{fc}, 'subject')
                t2slsubj = fc;
            end
        end
    end
    if ~isfield(opts, 'timeml') || ~islogical(opts.timeml) || numel(opts.timeml) ~= 1
        opts.timeml = false;
    end
end
if ~isfield(opts, 'transio') || ~islogical(opts.transio) || numel(opts.transio) ~= 1
    usetransio = false;
else
    usetransio = opts.transio;
end
copts = {opts};

% we need to figure out condition -> run indexing (per subject)
if aslfrrun || aslfrsubj
    slfsrun = cell(nssubs, 1);
    for fc = 1:nssubs
        stidx = find(strcmp(stsubids, selsubs{fc}));
        stcond = NaN(numconds, 1);
        for sfc = 1:numel(stidx)
            study = studies(stidx(sfc));
            if ~isfield(study, 'RunTimeVars') || ~isstruct(study.RunTimeVars) || ...
                numel(study.RunTimeVars) ~= 1 || ~isfield(study.RunTimeVars, 'Predictors') || ...
               ~isfield(study.RunTimeVars, 'SDMMatrix')
                error('neuroelf:xff:missingInformation', ...
                    'Study %d for subject %s doesn''t contain necessary information.', ...
                    sfc, selsubs{fc});
            end
            stpreds = study.RunTimeVars.Predictors;
            stsdm = study.RunTimeVars.SDMMatrix;
            for spc = 1:numel(stpreds)
                if all(stsdm(:, spc) == 0)
                    continue;
                end
                stpmatch = strcmpi(selpreds, stpreds{spc});
                if ~any(stpmatch)
                    continue;
                end
                stcond(stpmatch & ~isnan(stcond)) = 0;
                stcond(stpmatch & isnan(stcond)) = sfc;
            end
        end
        slfsrun{fc} = stcond;
    end
    if aslfrsubj
        csubj = cellfun('prodofsize', slfsrun);
        csubja = zeros(sum(csubj), 1);
        sfc = 1;
        for fc = 1:nssubs
            csubja(sfc:(sfc+csubj(fc)-1)) = fc;
            sfc = sfc + csubj(fc);
        end
        csubj = {csubja};
    end
end

% do we need to figure out grouping
if agfrsubg
    if ~isfield(rtv, 'Groups') || ~iscell(rtv.Groups) || size(rtv.Groups, 2) ~= 2 || ...
        isempty(rtv.Groups) || size(rtv.Groups, 1) < 2
        error('neuroelf:xff:badSetting', 'Bad grouping in RunTimeVars.');
    end
    subgroup = zeros(nssubs, 1);
    for fc = 1:size(rtv.Groups, 1)
        subgroup(ne_methods.multimatch(selsubs, subids(rtv.Groups{fc, 2})) > 0) = fc;
    end
    if any(subgroup == 0)
        error('neuroelf:xff:badSetting', 'Some subjects are not in a selected group.');
    end
    csubg = {subgroup};
end

% total number of ALL maps
tnmaps = tgnstore + tnstore * nssubs;

% copy argument list
gfargs = gfunc(2:end);
gfunc = gfunc{1};
gfout = cell(1, gfnout);
slfargs = slfunc(2:end);
slfunc = slfunc{1};
slfout = cell(nssubs, slfnout);

% progress bar
try
    if progress
        if isempty(pbar)
            pbar = xprogress;
            xprogress(pbar, 'settitle', 'GLM-based searchlight...');
            xprogress(pbar, 0, 'Checking function handles...', 'visible', 0, 1);
        else
            pbar.Progress(0, 'Checking function handles...');
        end
        pbint = 1 / 86400;
        pbnext = now;
    end
catch xfferror
    pbar = [];
    progress = false;
    neuroelf_lasterr(xfferror);
end

% time to figure out what to do -> time domain
if usetfunc

    % get all studys
    if any(~cellfun(@ischar, stfiles) | cellfun('isempty', stfiles))
        error('neuroelf:xff:fileNotAvailable', 'Study:AnalyzedFiles not all available.');
    end
    if any(cellfun(@exist, stfiles) ~= 2)
        error('neuroelf:xff:fileNotAvailable', 'Study:AnalyzedFiles not all available.');
    end
    stsubs = glm_Subjects(xo, true);
    stobjs = cell(numel(stfiles), 1);
    try
        if dodebug
            fprintf('Accessing time-course files (transio)...\n');
        end
        for fc = 1:numel(stobjs)
            if any(strcmp(selsubs, stsubs{fc}))
                stobjs{fc} = xff(stfiles{fc}, 't');
            end
        end
        selsubsst = selsubs;
        selsubsstd = selsubs;
        selsubsstn = selsubs;
        for fc = 1:nssubs
            selsubsst{fc} = stobjs(strcmp(stsubs, selsubs{fc}));
            selsubsstn(fc) = numel(selsubsst{fc});
            if selsubsstn(fc) == 0
                error('neuroelf:xff:noData', 'No time-course data for subjects %s.', selsubs{fc});
            end
            selsubsstd{fc} = cell(selsubsstn(fc), 1);
            for sc = 1:selsubsstn(fc)
                vtclay = aft_Layout(selsubsst{fc}{sc});
                if any(vtclay([1:3, 4:13]) ~= glmlay([1:3, 4:13]))
                    error('neuroelf:xff:nadData', 'Bad time-course data for subjects %s.', selsubs{fc});
                end
                selsubsst{fc}{sc} = selsubsst{fc}{sc}.C.VTCData;
            end
        end
    catch xfferror
        clearxffobjects(stobjs);
        rethrow(xfferror);
    end

    % try to access time course data (pick central voxel of slab)
    tvox = xmask(spi(1 + floor(0.5 * nspi)));
    for sc = 1:nssubs
    end
    
% use trial data
else

    % multilevel
    if mlevel
        error('neuroelf:xff:notYetImplemented', 'Multi-level Searchlight not yet implemented.');

    % single subject data
    else
        
        % data would be conditions-by-features (maxnumidx)
        data = randn(numconds, maxnumidx);

        % replace
        if aslfrcond
            slfargs(slfrcond) = {selcondi};
        end
        if aslfrdata
            slfargs(slfrdata) = {data};
        end
        if aslfrdatat
            slfargs(slfrdatat) = {data'};
        end
        if aslfropts
            slfargs(slfropts) = copts;
        end
        if aslfrrun
            slfargs(slfrrun) = slfsrun(1);
        end
        if aslfrsubj
            slfargs(slfrsubj) = csubj;
        end

        % evaluate function
        try
            if dodebug
                fprintf('Testing searchlight function (slfunc)...\n');
            end
            [slfout{1, :}] = feval(slfunc, slfargs{:});
        catch xfferror
            rethrow(xfferror);
        end

        % test output
        mnames = cell(tnstore, 1);
        mni = 1;
        for fc = 1:size(slstore, 1)
            if numel(slfout{1, slstore{fc, 1}}) < max(slstore{fc, 2})
                error('neuroelf:xff:tooFewOutputs', 'Searchlight function produces too few outputs.');
            end
            for sfc = 1:numel(slstore{fc, 2})
                if numel(slstore{fc, 2}) > 1 || ~isempty(strfind(slstore{fc, 3}, '%d'))
                    mnames{mni} = sprintf(slstore{fc, 3}, slstore{fc, 2}(sfc));
                else
                    mnames{mni} = slstore{fc, 3};
                end
                mni = mni + 1;
            end
        end

        % null-distribution
        nulldata = zeros(nulld, tnstore);
        if aslfrdata
            if dodebug && nulld > 0
                fprintf('Null-distribution (%d samples) of searchlight function (slfunc)...\n', nulld);
            end
            for fc = 1:nulld
                if progress && now > pbnext
                    pbar.Progress(fc / nulld, ...
                        sprintf('Null-distribution (sample %d/%d)...', fc, nulld));
                    pbnext = now + pbint;        
                end
                slfargs(slfrdata) = {randn(numconds, maxnumidx)};
                if aslfrrun
                    slfargs(slfrrun) = slfsrun(ceil(nssubs * rand(1, 1)));
                end
                [slfout{1, :}] = feval(slfunc, slfargs{:});
                for mc = 1:tnstore
                    nulldata(fc, mc) = slfout{1, nstoren(mc)}(nstore(mc));
                end
            end
        end

        % compute null-correction (chance/median value)
        if nullcorr
            ncdata = median(nulldata, 1);
        end

        % test group function
        if nulld >= nssubs
            gdata = nulldata(1:nssubs, :);
        else
            gdata = rand(nssubs, tnstore) - 0.5;
        end
        
        % replace simple arguments
        if agfrdata
            gfargs(gfrdata) = {gdata};
        end
        if agfropts
            gfargs(gfropts) = copts;
        end
        if agfrsubg
            gfargs(gfrsubg) = csubg;
        end
        if agfrsubs
            gfargs(gfrsubs) = selsub;
        end
        if agfrsldatan
            slfout = repmat(slfout(1, :), nssubs, 1);
        end
        
        % evaluate function
        try
            for fc = 1:numel(gfrdatan)
                gfargs{gfrdatan(fc)} = gdata(:, gfrdatai{fc});
            end
            for fc = 1:numel(gfrsldatan)
                gfargs{gfrsldatan(fc)} = cat(1, slfout{:, gfrsldatai});
            end
            if dodebug
                fprintf('Testing searchlight function (gfunc)...\n');
            end
            [gfout{:}] = feval(gfunc, gfargs{:});
        catch xfferror
            rethrow(xfferror);
        end

    end
end

% what output
if bc.ProjectType == 1
    otype = 'vmp';
    omap = 'VMPData';
elseif bc.ProjectType == 2
    otype = 'smp';
    omap = 'SMPData';
end

% prepare output
reusemap = false;
if isfield(opts, 'targetmap') && numel(opts.targetmap) == 1 && ...
    xffisobject(opts.targetmap, true, otype) && numel(opts.targetmap.C.Map) == tnmaps
    map = opts.targetmap;
    glmlay = aft_Layout(xo);
    maplay = aft_Layout(map);
    if ~isequal(glmlay([1:3, 5:10]), maplay([1:3, 5:10]))
        error('neuroelf:xff:badOption', 'Invalid targetmap specified.');
    end
    reusemap = true;
elseif isfield(opts, 'targetmap') && ischar(opts.targetmap) && ~isempty(opts.targetmap) && ...
   ~isempty(regexpi(opts.targetmap(:)', ['\.' otype '$']))
    try
        map = [];
        map = xff(opts.targetmap(:)', 't');
        if ~xffisobject(map, true, otype) || numel(map.C.Map) ~= tnmaps
            error('neuroelf:xff:badOption', 'Invalid targetmap specified.');
        end
        maplay = aft_Layout(map);
        if ~isequal(glmlay([1:3, 5:10]), maplay([1:3, 5:10]))
            error('neuroelf:xff:badOption', 'Invalid targetmap specified.');
        end
        reusemap = true;
    catch xfferror
        delete(map);
        rethrow(xfferror);
    end
else
    map = xff(['new:' otype]);

    % make setting
end

% access output
mapc = map.C;

% store NULL data
mapc.RunTimeVars.SLNull = nulldata;

% preset first map
if ~reusemap
    if bc.ProjectType == 1
        mapc.Resolution = bc.Resolution;
        mapc.XStart = bc.XStart;
        mapc.XEnd = bc.XEnd;
        mapc.YStart = bc.YStart;
        mapc.YEnd = bc.YEnd;
        mapc.ZStart = bc.ZStart;
        mapc.ZEnd = bc.ZEnd;
    else
        mapc.NrOfVertices = bc.NrOfVertices;
    end
    mapc.Map = mapc.Map(1);
    mapc.Map.Type = 93;
    mapc.Map.(omap) = single(zeros(msize));
    mapc.Map = mapc.Map(ones(1, tnmaps));
end

% generate maps and store indices (important later to collect data!)
smapi = repmat((1:tnstore)', 1, nssubs) + tnstore .* repmat(0:(nssubs-1), tnstore, 1);
mapi = numel(smapi);
gmapi = mapi + (1:tgnstore);
if ~reusemap
    for sc = 1:nssubs
        for fc = 1:tnstore
            mapc.Map(smapi(fc, sc)).Name = sprintf('Subject %s: %s', ...
                subids{subsel(sc)}, mnames{fc});
        end
    end
    for fc = 1:tgnstore
        mapc.Map(gmapi(fc)).Name = gstore{fc, 3};
    end
end

% return
if opts.maponly
    clearxffobjects(stobjs);
    if progress && islogical(opts.progress) && ~isempty(pbar)
        closebar(pbar);
    end
    return;
end

% collect data required for extraction
if ~usetfunc
    tmaps = cell(nssubs, 1);
    if dodebug
        fprintf('Accessing BetaMaps (for slfunc)...\n');
    end
    for fc = 1:nssubs
        tmaps{fc} = bc.GLMData.Subject(subsel(fc)).BetaMaps;
        if ~usetransio && istransio(tmaps{fc})
            tmaps{fc} = resolve(tmaps{fc});
        end
        tmaps{fc} = reshape(tmaps{fc}, mnumel, size(tmaps{fc}, cdim));
    end
end

% iterate over part to do
for ic = 1:nspi

    % voxel index
    sic = xmask(spi(ic));

    % debug output
    if dodebug && now > dbnext
        fprintf('Searchlight (voxel %d/%d; part %d/%d -- total voxel %d/%d)...\n', ...
            ic, nspi, opts.subpart(1), opts.subpart(2), spi(ic), smask);
        dbnext = now + 1 / 720;
    end

    % progress
    if progress && now > pbnext
        pbar.Progress(ic / nspi, sprintf('Searchlight (voxel %d/%d)...', ic, nspi));
        pbnext = now + pbint;
        map.C = mapc;
    end

    % searchlight indices
    slidx = midx{sic};
    if numel(slidx) < minnfeat
        continue;
    end
    
    % weight
    if doweight == 1
        wsampl = spdiag(flexinterpn(dweight, 1 + 10 .* midxdist{sic}(:)));
    elseif doweight == 2
        wsampl = dweight(midxdist{sic});
    end
    if doweight
        nws = numel(wsampl);
        nw = 1:nws;
        wsampl = sparse(nw, nw, wsampl, nws, nws, nws);
    end

    % multilevel
    if mlevel
        error('neuroelf:xff:notYetImplemented', 'Not yet implemented.');

    % lower level first
    else

        % prepare group data
        gdata(:) = NaN;

        % iterate over subjects
        for sc = 1:nssubs

            % time first
            if usetfunc

            % trials
            else

                % access trials
                data = double(tmaps{sc}(slidx, selconds)');

                % which samples are bad
                goodsamples = ~all(isinf(data) | isnan(data) | data == 0, 2);
                if all(goodsamples)
                    cond = selcondi;
                else
                    ngoodsamples = sum(goodsamples);
                    cond = selcondi(goodsamples);
                    data(~goodsamples, :) = [];
                end
                
                % weighting
                if doweight
                    data = data * wsampl;
                end

                % bad features
                goodfeatures = ~all(isinf(data) | isnan(data) | data == 0, 1);
                if ~all(goodfeatures)
                    ngoodfeatures = sum(goodfeatures);
                    if ngoodfeatures < minnfeat
                        continue;
                    end
                    data = data(:, goodfeatures);
                end

                % prepare arguments
                if aslfrcond
                    slfargs(slfrcond) = {cond};
                end
                if aslfrdata
                    slfargs(slfrdata) = {data};
                end
                if aslfrdatat
                    slfargs(slfrdatat) = {data'};
                end
                if aslfrrun
                    slfargs(slfrrun) = {slfsrun{sc}(goodsamples)};
                end

                % evaluate
                try
                    [slfout{sc, :}] = feval(slfunc, slfargs{:});
                catch xfferror
                    rethrow(xfferror);
                end

                % store values
                if nullcorr
                    for mc = 1:tnstore
                        gdata(sc, mc) = slfout{sc, nstoren(mc)}(nstore(mc)) - ncdata(mc);
                        mapc.Map(smapi(mc, sc)).(omap)(sic) = gdata(sc, mc);
                    end
                else
                    for mc = 1:tnstore
                        gdata(sc, mc) = slfout{sc, nstoren(mc)}(nstore(mc));
                        mapc.Map(smapi(mc, sc)).(omap)(sic) = gdata(sc, mc);
                    end
                end

                % return early
                if ~isinf(maxtime) && now > maxtime
                    if dodebug
                        fprintf('Maximum time reached -- quitting...\n');
                    end
                    if progress && islogical(opts.progress) && ~isempty(pbar)
                        closebar(pbar);
                    end
                    map.C = mapc;
                    return;
                end
            end
        end

        % group function replacements
        if agfrdata
            gfargs(gfrdata) = {gdata};
        end
        for fc = 1:numel(gfrdatan)
            gfargs{gfrdatan(fc)} = gdata(:, gfrdatai{fc});
        end
        for fc = 1:numel(gfrsldatan)
            gfargs{gfrsldatan(fc)} = cat(1, slfout{:, gfrsldatai});
        end

        % evaluation
        [gfout{:}] = feval(gfunc, gfargs{:});

        % and data storage
        for mc = 1:tgnstore
            mapc.Map(gmapi(mc)).(omap)(sic) = gfout{gnstoren(mc)}(gnstore(mc));
        end
    end
end

% clean up progress
if progress && islogical(opts.progress) && ~isempty(pbar)
    closebar(pbar);
end

% set back in VMP object
map.C = mapc;

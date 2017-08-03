function ava = glm_Anova(xo, factors, betas, subjects, group, covariate, am)
% GLM::Anova  - create an ANOVA object from a GLM
%
% FORMAT:       ava = glm.Anova(factors, betas)
%
% Input fields:
%
%       factors     1xF struct with fields:
%        .fname     factor name
%        .ftype     'between' or 'within'
%        .levels    1xL cell array with level labels
%       betas       N-dimensional array with unique beta assignments,
%                   number of dimensions = number of within factors
%       subjects    list of subjects to include (by default: all)
%       group       SxB double array with group assignment(s)
%                   (default: empty, no between factor)
%       covariate   CxS covariate data
%       addmean     addmean flag (default: false)
%
% Note: currently, only the models used by BrainVoyager QX can be created:
%       - 1 within factor
%       - 2 within factors
%       - 3 within factors
%       - 1 between factor X 1 within factor
%
%       for covariates-only models, please use the GLM::RFX_rMap
%
% Using: emptystruct, mrmanova.

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:06 AM EST
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
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
   ~isstruct(factors) || isempty(factors) || numel(factors) > 3 || ...
   ~isfield(factors, 'fname') || ~isfield(factors, 'ftype') || ~isfield(factors, 'levels') || ...
    isempty(betas) || ~isa(betas, 'double') || ...
    any(isinf(betas(:)) | isnan(betas(:)) | betas(:) < 1 | betas(:) ~= fix(betas(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if bc.ProjectTypeRFX ~= 1
    error('neuroelf:xff:invalidObject', 'Only defined for RFX GLM.');
end
betas = squeeze(betas);
if numel(betas) == size(betas, 2)
    betas = betas(:);
end
if any(betas(:) > (size(bc.GLMData.Subject(1).BetaMaps, ndims(bc.GLMData.Subject(1).BetaMaps)) - 1)) || ...
    numel(betas) ~= numel(unique(betas(:)))
    error('neuroelf:xff:invalidObject', 'Bad beta assignment.');
end
bfc = 0;
wfc = 0;
wfl = ones(1, 4);
for fc = 1:numel(factors)
    if ~ischar(factors(fc).fname) || isempty(factors(fc).fname) || ~ischar(factors(fc).ftype) || ...
       ~any(strcmpi(factors(fc).ftype(:)', {'b', 'between', 'w', 'within'})) || ...
       ~iscell(factors(fc).levels) || numel(factors(fc).levels) < 2
        error('neuroelf:xff:badArgument', 'Bad factors argument.');
    end
    factors(fc).ftype = lower(factors(fc).ftype(1));
    if factors(fc).ftype == 'b'
        bfc = bfc + 1;
    else
        wfc = wfc + 1;
        wfl(wfc) = numel(factors(fc).levels);
    end
    factors(fc).levels = factors(fc).levels(:);
    for lc = 1:numel(factors(fc).levels)
        if ~ischar(factors(fc).levels{lc}) || isempty(factors(fc).levels{lc})
            error('neuroelf:xff:badArgument', 'Bad factors argument.');
        end
        factors(fc).levels{lc} = factors(fc).levels{lc}(:)';
    end
end
if sum(size(betas) > 1) ~= wfc || ~isequal(size(betas), wfl(1:max(2, wfc)))
    error('neuroelf:xff:badArgument', 'Bad beta levels assignment.');
end

% for now only support certain models
if bfc > 1 || (bfc > 0 && wfc > 1) || wfc > 3 || (bfc == 0 && wfc == 0)
    error('neuroelf:xff:unsupported', ...
        'AVA files with this factor combination not yet supported.');
end

% check other arguments
if nargin < 4 || ~isa(subjects, 'double') || isempty(subjects) || ...
    numel(subjects) ~= length(subjects) || ndims(subjects) > 2 || ...
    any(isnan(subjects) | subjects < 1 | subjects > bc.NrOfSubjects | subjects ~= fix(subjects)) || ...
    numel(unique(subjects)) ~= numel(subjects)
    subjects = 1:bc.NrOfSubjects;
else
    subjects = subjects(:)';
end
numsubs = numel(subjects);
if nargin < 5 || ~isa(group, 'double') || size(group, 1) ~= numsubs || ...
    any(isnan(group(:)) | group(:) < 0 | group(:) >= numsubs | group(:) ~= fix(group(:)))
    group = zeros(numsubs, 0);
end
if size(group, 2) ~= bfc
    error('neuroelf:xff:badArgument', 'Group assignment incorrect.');
end
if nargin < 6 || ~isa(covariate, 'double') || size(covariate, 1) ~= numsubs || ...
    any(isinf(covariate(:)) | isnan(covariate(:)))
    covariate = zeros(numsubs, 0);
end
if nargin < 7 || ~islogical(am) || numel(am) ~= 1
    am = false;
end

% get some GLM settings
prjtype = bc.ProjectType + 1;
numsprd = bc.NrOfSubjectPredictors;
numstud = numel(bc.Study);
stfiles = cell(numstud, 1);
for stc = 1:numstud
    stfiles{stc} = bc.Study(stc).NameOfAnalyzedFile;
end

% generate predictor structs
bfs = ne_methods.emptystruct({'Type', 'Name', 'NrOfLevels', 'LevelNames'});
wfs = ne_methods.emptystruct({'Type', 'Name', 'NrOfLevels', 'LevelNames'});

% further parse factors
numfac = numel(factors);
bfc = 0;
wfc = 0;
blevs = [0, 0, 0];
wlevs = [0, 0, 0];
for fc = 1:numfac
    flevels = factors(fc).levels;
    switch (lower(factors(fc).ftype))
        case 'b'
            bfc = bfc + 1;
            blevs(bfc) = numel(flevels);
            if numel(unique(group(:, bfc))) ~= blevs(bfc)
                error('neuroelf:xff:badArgument', ...
                    'Invalid between-subjects factor structure.');
            end
            bfs(bfc).Type = 1;
            bfs(bfc).Name = factors(fc).fname;
            bfs(bfc).NrOfLevels = blevs(bfc);
            bfs(bfc).LevelNames = flevels;
        case 'w'
            wfc = wfc + 1;
            wlevs(wfc) = numel(flevels);
            wfs(wfc).Type = 1;
            wfs(wfc).Name = factors(fc).fname;
            wfs(wfc).NrOfLevels = wlevs(wfc);
            wfs(wfc).LevelNames = flevels(:);
    end
end

% determine model type
mtype = 0;
switch (wfc)
    case 0
        if bfc == 0
            mtype = 4;
        else
            mtype = 6;
        end
    case 1
        if bfc == 0
            mtype = 1;
        else
            mtype = 3;
        end
    case 2
        if bfc == 0
            mtype = 2;
        else
            mtype = 9;
        end
    case 3
        mtype = 5;
end
if mtype == 0
    error('neuroelf:xff:unsupported', ...
        'AVA object with this combination of factors unsupported.');
end

% generate arguments for mrmanova call
bsz = size(bc.GLMData.Subject(1).BetaMaps);
sd = numel(bsz);
bsi = repmat({':'}, 1, sd);
bse = repmat({':'}, 1, sd);
bse{sd} = bsz(end);
nd = sd + 1;
bti = repmat({':'}, 1, nd);
bsz(end) = numsubs;
bsz(end+1) = numel(betas);
bm = zeros(bsz);
for subc = 1:numsubs
    bti{sd} = subc;
    for betac = 1:numel(betas)
        bsi{sd} = betas(betac);
        bti{nd} = betac;
        if am
            bm(bti{:}) = bc.GLMData.Subject(subc).BetaMaps(bsi{:}) + ...
                bc.GLMData.Subject(subc).BetaMaps(bse{:});
        else
            bm(bti{:}) = bc.GLMData.Subject(subc).BetaMaps(bsi{:});
        end
    end
end
[w{1:sum(size(betas) > 1)}] = ind2sub(size(betas), (1:numel(betas))');
for wc = 1:numel(w)
    w{wc} = w{wc}(:)';
end
w = cat(1, w{:});

% compute anova output
[msmaps, msmapn, cmmaps] = ne_methods.mrmanova(bm, w, group, struct('cov', covariate));
cmmaps = reshape(permute(cmmaps, [1:(sd-1), sd + 1, sd]), ...
    [bsz(1:sd-1), size(cmmaps, sd) * size(cmmaps, sd + 1)]);
bsi(end) = [];

% generate new AVA object
ava = xff('new:ava');
avac = ava.C;
avac.ProjectType = prjtype;
avac.ModelType = mtype;
avac.NrOfBetweenFactors = bfc;
avac.NrOfWithinFactors = wfc;
avac.NrOfCovariates = 0;
if isfield(bc.RunTimeVars, 'TrfPlus')
    avac.RunTimeVars.TrfPlus = bc.RunTimeVars.TrfPlus;
end
if bfc > 0
    avac.NrOfBetweenFactorLevels = blevs(1:bfc);
else
    avac.NrOfBetweenFactorLevels = zeros(1, 0);
end
if wfc > 0
    avac.NrOfWithinFactorLevels = wlevs(1:wfc);
else
    avac.NrOfWithinFactorLevels = zeros(1, 0);
end
avac.NrOfSubjects = numsubs;
avac.NrOfPredictors = numsprd - 1;
avac.BetweenFactor = bfs;
avac.WithinFactor = wfs;
avac.ContrastName = '';
avac.Covariates = {};
avac.OriginatingGLM = xo.F;
avac.NrOfStudies = numstud;
avac.TimeCourseFiles = stfiles;
avac.Resolution = bc.Resolution;
avac.Maps = struct;

% further calculation is somewhat based on project type
switch (prjtype)
    case 2
        avac.XStart = bc.XStart;
        avac.XEnd   = bc.XEnd;
        avac.YStart = bc.YStart;
        avac.YEnd   = bc.YEnd;
        avac.ZStart = bc.ZStart;
        avac.ZEnd   = bc.ZEnd;
    case 3
        avac.NrOfVertices = bc.NrOfVertices;
    otherwise
        error('neuroelf:xff:badObject', 'Only VTC and MTC GLMs supported for now.');
end

% depending on model type
switch (mtype)

    % one within factor
    case {1}
        avac.Maps.MS_A = msmaps(bsi{:}, 2);
        avac.Maps.MS_ErrorA = msmaps(bsi{:}, 3);

    % two within factors
    case {2}

    % one within factor, one between factor
    case {3}
        avac.Maps.MS_A = msmaps(bsi{:}, 3);
        avac.Maps.MS_ErrorB = msmaps(bsi{:}, 2);
        avac.Maps.MS_B = msmaps(bsi{:}, 1);
        avac.Maps.MS_AxB = msmaps(bsi{:}, 4);
        avac.Maps.MS_ErrorAxB = msmaps(bsi{:}, 5);

    % three within factors
    case {5}

    % one between factor (only)
    case {6}

    % two within factors, one between factor
    case {9}
end

% add cell means to maps
avac.Maps.CellMeans = cmmaps;

% put back in storage
ava.C = avac;

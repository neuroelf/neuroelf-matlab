% PUBLIC FUNCTION ne_mkda_setplp: set one PLP as current object
function varargout = ne_mkda_setplp(varargin)

% Version:  v1.1
% Build:    16040614
% Date:     Apr-06 2016, 2:45 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2016, Jochen Weber
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

% global variable
global ne_gcfg;
hFig = ne_gcfg.h.MKDA.MKDAFig;
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get content of PLP dropdown
plps = ch.PLPs;
plpud = plps.UserData;
plpst = plps.String;
if ~iscell(plpst)
    plpst = cellstr(plpst);
end

% inputs
if nargin < 3 || ...
    numel(varargin{3}) ~= 1 || ...
   (~isa(varargin{3}, 'double') && ...
    ~isxff(varargin{3}, 'plp'))
    plpidx = ch.PLPs.Value;
elseif isa(varargin{3}, 'double')
    plpidx = varargin{3};
    if isinf(plpidx) || ...
        isnan(plpidx) || ...
        plpidx < 1 || ...
        plpidx > numel(plpst) || ...
        plpidx ~= fix(plpidx)
        return;
    end
else
    plpf = false;
    for plpidx = 1:size(plpud, 1)
        if numel(plpud{plpidx, 3}) == 1 && ...
            isxff(plpud{plpidx, 3}, 'plp') && ...
            plpud{plpidx, 3} == varargin{3}
            plpf = true;
            break;
        end
    end
    if ~plpf
        return;
    end
end
plps.Value = plpidx;

% disable PLP menu
hFig.SetGroupEnabled('PLPOK', 'off');
hFig.SetGroupEnabled('Gauss', 'off');
hFig.SetGroupEnabled('HasAnas', 'off');
ne_gcfg.fcfg.plp = [];
ch.Studies.Value = [];
ch.Studies.String = {};
ch.Columns.Value = [];
ch.Columns.String = {};
ch.Analyses.Value = 1;
ch.Analyses.String = {'<no analyses defined>'};
ch.Contrast.String = '<not yet specified>';
ch.ContColumn.Value = 1;
ch.ContColumn.String = {'<no columns>'};
ch.CndParts.Value = [];
ch.CndParts.String = {};
ch.CndColumn.Value = 1;
ch.CndColumn.String = {'<no columns>'};
ch.CndOperator.Value = 1;
ch.Weights.String = '1';
ch.PointsLabel.String = 'Selected points: (selecting has no effect!)';
ch.Points.Value = [];
ch.Points.String = {'none'};

% get plp
if isempty(plpud)
    return;
end
plp = plpud{plpidx, 3};
if numel(plp) ~= 1 || ...
   ~isxff(plp, 'plp')
    return;
end
rtv = plp.RunTimeVars;

% no study column
cn = plp.ColumnNames(:);
data = plp.Points(:, :);
if size(data, 1) > 1
    mcn = (any(diff(data, 1, 1) ~= 0, 1));
    mcn = mcn(:);
else
    mcn = true(size(data, 2), 1);
end
if ~any(strcmp(cn, 'Study'))
    uiwait(errordlg('PLP object must contain Study column.', 'NeuroElf - error', 'modal'));
    return;
end

% set as globally accessible object
ne_gcfg.fcfg.plp = plp;

% make sure object has settings
if ~isfield(rtv, 'Config') || ...
   ~isstruct(rtv.Config)
    plp.RunTimeVars.Config = ne_gcfg.fcfg.MKDA;
    rtv = plp.RunTimeVars;
end
sfld = fieldnames(ne_gcfg.fcfg.MKDA);
for sfc = 1:numel(sfld)
    if ~isfield(rtv.Config, sfld{sfc})
        plp.RunTimeVars.Config.(sfld{sfc}) = ne_gcfg.fcfg.MKDA.(sfld{sfc});
    end
end
rtv = plp.RunTimeVars;

% update settings
mcfg = rtv.Config;
ne_mkda_setoption(0, 0, 'ApplyMask', mcfg.ApplyMask);
ne_mkda_setoption(0, 0, 'ContrastComp', mcfg.ContrastComp);
ne_mkda_setoption(0, 0, 'ContrastCompExclWeight', mcfg.ContrastCompExclWeight);
ne_mkda_setoption(0, 0, 'GroupMapComp', mcfg.GroupMapComp);
ne_mkda_setoption(0, 0, 'JoinBlobComp', mcfg.JoinBlobComp);
ne_mkda_setoption(0, 0, 'LabelColumn', mcfg.LabelColumn);
ne_mkda_setoption(0, 0, 'KeepIndivMaps', mcfg.KeepIndivMaps);
ne_mkda_setoption(0, 0, 'PPSWeighting', mcfg.PPSWeighting);
ne_mkda_setoption(0, 0, 'SpatialNull', mcfg.SpatialNull);
ne_mkda_setoption(0, 0, 'SummaryVMP', mcfg.SummaryVMP);
ne_mkda_setoption(0, 0, 'UniqueUnitPoints', mcfg.UniqueUnitPoints);

% make sure object has analyses info
if ~isfield(rtv, 'MKDAAnalyses') || ...
   ~iscell(rtv.MKDAAnalyses) || ...
    isempty(rtv.MKDAAnalyses) || ...
    size(rtv.MKDAAnalyses, 2) ~= 2 || ...
    ndims(rtv.MKDAAnalyses) > 2
    plp.RunTimeVars.MKDAAnalyses = {'<no analyses defined>', struct( ...
        'CndParts',    {{}}, ...
        'ContColumn',  '', ...
        'Contrast',    '', ...
        'Iterations',  5000, ...
        'Mask',        'colin_brain_ICBMnorm.vmr', ...
        'NullDist',    'spatial', ...
        'Resolution',  3, ...
        'Scaling',     'gauss', ...
        'SphereSize',  12, ...
        'SphereTaper', 0, ...
        'StudyColumn', '', ...
        'Weights',     '1')};
    rtv = plp.RunTimeVars;
end

% set UICs
studies = plp.Study;
if min(studies) >= 1 && max(studies) <= numel(plp.Labels)
    [studies, ustudies] = unique(plp.Labels(studies));
else
    studies = splittocellc(sprintf('%d|', studies), '|');
    [studies, ustudies] = unique(studies(:));
end
ch.Studies.String = studies;
ch.Studies.Value = (1:numel(studies))';
ch.Studies.UserData = struct('ustudies', ustudies);
if any(strcmpi(cn, 'contrast'))
    clnames = uunion({'Study'; 'Contrast'}, cn(mcn));
else
    clnames = uunion({'Study'}, cn(mcn));
end
ch.Columns.String = [clnames(:); setdiff(cn(:), clnames(:))];
ch.Columns.Value = lsqueeze(find(multimatch(lower(clnames), ...
    {'x', 'y', 'z', 'study', 'contrast', 'n', 'fixedrandom'}) > 0));
if any(strcmpi('fixedrandom', clnames)) && ...
    all(diff(plp.Points(:, findfirst(strcmpi('fixedrandom', cn)))) == 0)
    ch.Columns.Value(ch.Columns.Value == findfirst(strcmpi('fixedrandom', clnames))) = [];
end
ch.Analyses.Value = 1;
ch.Analyses.String = plp.RunTimeVars.MKDAAnalyses(:, 1);
ch.CndColumn.String = clnames;
if any(strcmpi('n', cn))
    if any(strcmpi('fixedrandom', cn))
        if any(diff(plp.Points(:, findfirst(strcmpi('fixedrandom', cn))) ~= 0))
            ch.Weights.String = 'sqrt($N) .* (0.75 + 0.25 * ($FixedRandom > 0))';
        else
            ch.Weights.String = 'sqrt($N)';
        end
    else
        ch.Weights.String = 'sqrt($N)';
    end
elseif any(strcmpi('fixedrandom', cn))
    ch.Weights.String = '0.75 + 0.25 * ($FixedRandom > 0)';
end

% try to determine which columns *might* be statistical units
if ~isfield(rtv, 'StudyColumns') || ...
   ~iscell(rtv.StudyColumns) || ...
    isempty(rtv.StudyColumns)

    % copy column names
    sclnames = clnames;

    % remove clear no-shows
    sclnames(multimatch(lower(sclnames), ...
        {'x'; 'y'; 'z'; 'coordsys'; 'coordsystem'; 'usepoint'; 'value'; ...
         'df1'; 'df2'; 'clustersize'; 'vthresh'; 'kthresh'; 'label'; ...
         'software'; 'rfx'; 'correction'; 'group'; 'groupsize'; ...
         'symbol'; 'size'; 'color'; 'duration'; 'n'; 'fixedrandom'}) > 0) = [];

    % get actual study column index and data, as well as unique entries
    stcol = findfirst(strcmpi('study', cn));
    stdat = data(:, stcol);
    studat = unique(stdat);

    % for each unique entry
    stui = cell(1, numel(studat));
    stun = zeros(1, numel(stui));
    for cc = 1:numel(stun)

        % determine indices and number of indices
        stui{cc} = find(stdat == studat(cc));
        stun(cc) = numel(stui{cc});
    end

    % create test data for as much space as we need
    stumn = max(stun);
    stumf = ceil(stumn ./ stun);
    stutd = randn(numel(stun), stumn);

    % for each remaining column
    sclmatch = multimatch(lower(sclnames(:)), lower(cn));
    for cc = numel(sclnames):-1:1

        % fill test data
        for sc = 1:numel(stun)
            stutr = repmat(data(stui{sc}, sclmatch(cc))', 1, stumf(sc));
            stutd(sc, :) = stutr(1, 1:stumn);
        end

        % values must be unique across studies
        stutd = sort(stutd, 1);
        if any(any(diff(stutd, 1, 1) == 0))
            sclnames(cc) = [];
        end
    end
    plp.RunTimeVars.StudyColumns = sclnames;
    rtv = plp.RunTimeVars;
else
    sclnames = rtv.StudyColumns(:);
end
ch.StudyColumn.String = sclnames;

% remove non-sensical entries (for contrast column)
clnames(multimatch(lower(clnames), {'n'; 'study'; 'x'; 'y'; 'z'}) > 0) = [];
ch.ContColumn.String = clnames;
if any(strcmpi(clnames, 'contrast'))
    ch.ContColumn.Value = find(strcmpi(clnames, 'contrast'));
end

% has no analyses
if size(rtv.MKDAAnalyses, 1) == 1 && ...
    strcmp(rtv.MKDAAnalyses{1, 1}, '<no analyses defined>')
    ne_mkda_setcontcol(0, 0, false);

    % figure out which unit column has the most promising spread
    sclnumu = zeros(1, numel(sclnames));
    for sc = 1:numel(sclnumu)
        sclnumu(sc) = numel(unique(data(:, findfirst(strcmpi(sclnames{sc}, cn)))));
    end
    sclnumu = maxpos(sclnumu);
    ch.StudyColumn.Value = sclnumu;

    % now set in default analysis
    if isempty(rtv.MKDAAnalyses{1, 2}.ContColumn)
        plp.RunTimeVars.MKDAAnalyses{1, 2}.ContColumn = ...
            ch.ContColumn.String{ch.ContColumn.Value};
        plp.RunTimeVars.MKDAAnalyses{1, 2}.Contrast = ch.Contrast.String;
    end
    if isempty(rtv.MKDAAnalyses{1, 2}.StudyColumn)
        plp.RunTimeVars.MKDAAnalyses{1, 2}.StudyColumn = sclnames{sclnumu};
    end
    if strcmp(rtv.MKDAAnalyses{1, 2}.Weights, '1')
        plp.RunTimeVars.MKDAAnalyses{1, 2}.Weights = ch.Weights.String;
    end

% otherwise
else

    % enable group already
    hFig.SetGroupEnabled('HasAnas', 'on');
end

% make sure condition values are set
ne_mkda_setcondcol;

% configuration of studies exists
if isfield(rtv, 'StudySelection') && ...
    iscell(rtv.StudySelection)
    ssel = multimatch(lower(rtv.StudySelection), lower(studies));
    ch.Studies.Value = ssel(ssel > 0);
else
    plp.RunTimeVars.StudySelection = studies;
end

% simply read from config
ne_mkda_setana;

% re-enable PLP menu
hFig.SetGroupEnabled('PLPOK', 'on');

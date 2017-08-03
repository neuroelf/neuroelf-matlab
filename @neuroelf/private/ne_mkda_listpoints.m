% PUBLIC FUNCTION ne_mkda_listpoints: list matching points
function varargout = ne_mkda_listpoints(varargin)

% Version:  v0.9c
% Build:    13040223
% Date:     Nov-08 2011, 1:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only parse one at a time
if any(strcmp(ne_gcfg.c.blockcb, 'mkda_listpoints'))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'mkda_listpoints';

% pre-set points list value
ch.Points.Value = [];
ch.Points.String = {'<getting points list>'};

% get content of PLP dropdown
plps = ch.PLPs;
plpud = plps.UserData;
plpid = plps.Value;
try
    plp = plpud{plpid, 3};
    if numel(plp) ~= 1 || ...
       ~isxff(plp, 'plp')
        error( ...
            'neuroelf:GUI:BadPLP', ...
            'Bad PLP object.' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_listpoints')) = [];
    return;
end
rtv = plp.RunTimeVars;

% get condition particles
condstr = ch.CndParts.String;
if ~iscell(condstr)
    if isempty(condstr)
        condstr = {};
    else
        condstr = cellstr(condstr);
    end
end

% form string
if isempty(condstr)
    condstr = '';
else
    condstr = gluetostringc(condstr, ' ');
end

% extend by contrast particles
cparts = ch.Contrast.UserData;
if iscell(cparts{1})
    cparts = [cparts{1}(:); cparts{2}(:)];
else
    cparts = cparts(:);
end
ccol = ch.ContColumn.String{ch.ContColumn.Value};
ccoltext = ch.ContColumn.UserData;
if isempty(ccoltext)
    if isfield(rtv, 'ColumnIsText') && ...
        isfield(rtv.ColumnIsText, ccol)
        ccoltext = rtv.ColumnIsText.(ccol);
    end
end
for cc = 1:numel(cparts)
    if ccoltext || ...
       (isempty(regexpi(ddeblank(cparts{cc}), '^[+\-]?\d+(\.\d*)?(e[+\i]?\d+)?$')) && ...
        ~any(strcmpi(cparts{cc}, {'inf', '-inf', 'nan'})))
        cparts{cc} = sprintf(' $%s == ''%s'' ', ccol, cparts{cc});
    else
        cparts{cc} = sprintf(' $%s == %s ', ccol, cparts{cc});
    end
end
cparts = gluetostringc(cparts, '|');
if isempty(condstr)
    condstr = cparts;
else
    condstr = ['(' condstr ') & (' cparts ')'];
end

% extend by study selection
allstudies = plp.Study;
studies = unique(allstudies);
selstudies = studies(ch.Studies.Value);
if isempty(selstudies)
    condstr = 'false';
elseif numel(selstudies) < numel(studies)
    ustudies = ch.Studies.UserData.ustudies(:)';
    if numel(selstudies) > (0.5 * numel(studies))
        nonselstudies = setdiff(studies(:), selstudies(:));
        spart = sprintf('$Study ~= %d & ', allstudies(ustudies(nonselstudies(:)')));
    else
        spart = sprintf('$Study == %d | ', allstudies(ustudies(selstudies(:)')));
    end
    if condstr(1) == '('
        condstr = [condstr ' & (' spart(1:end-3) ')'];
    else
        condstr = ['(' condstr ') & (' spart(1:end-3) ')'];
    end
end

% get points selection
mskfile = ch.Mask.String{1};
if mskfile(1) == 'c'
    mskfile = [neuroelf_path('colin'), filesep, mskfile];
end
sel = plp.Select(condstr, mskfile);
onumsel = numel(sel);

% get column selection
colnames = ch.Columns.String(ch.Columns.Value);
colmatch = multimatch(lower(colnames(:)), lower(plp.ColumnNames(:)));
if any(colmatch < 1)
    ch.Points.String = {'<error processing points>'};
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_listpoints')) = [];
    return;
end

% get data
cn = plp.ColumnNames(:);
data = plp.Points(:, :);
datasel = (1:size(data, 1))';
labels = plp.Labels(:);

% how many units
stcol = findfirst(strcmpi(ch.StudyColumn.String{ch.StudyColumn.Value}, cn));
numunits = numel(unique(data(:, stcol)));
selunits = numel(unique(data(sel, stcol)));

% enforce selection
if ch.SelectedPts.Value > 0
    data = data(sel, :);
    datasel = datasel(sel);
    sel = sel(sel);
end

% get weights
w = lower(ddeblank(ch.Weights.String));

% no weighting
if isempty(w) || ...
    strcmp(w, '1')

    % set to 1
    data(:, end+1) = 1;

% work to do
else

    % compute weights
    cnl = zeros(numel(cn), 1);
    for cnc = 1:numel(cn)
        cnl(cnc) = numel(cn{cnc});
    end
    [cnl, cni] = sort(cnl, 'descend');
    cns = cn(cni);
    for cnc = 1:numel(cn)
        w = regexprep(w, ['\$' lower(cns{cnc}) '\s*~=\s*''([^'']*)'''], ...
            sprintf('~strcmpi(labels(data(:,%d)),''$1'')', cni(cnc)));
        w = regexprep(w, ['\$' lower(cns{cnc}) '\s*==\s*''([^'']*)'''], ...
            sprintf('strcmpi(labels(data(:,%d)),''$1'')', cni(cnc)));
        w = strrep(w, ['$' lower(cns{cnc})], sprintf('data(:,%d)', cni(cnc)));
    end

    % apply weights
    try
        data(:, end+1) = eval(w);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(warndlg('Invalid weighting function.', 'NeuroElf - error', 'modal'));
        data(:, end+1) = 1;
        ch.Weights.String = '1';
    end
end

% get data for selected columns
headers = cn(colmatch);
data = data(:, [colmatch(:)', end]);
headers{end+1} = 'Weight';
textcols = false(1, numel(headers));
if isfield(rtv, 'ColumnIsText') && ...
    isstruct(rtv.ColumnIsText) && ...
    numel(rtv.ColumnIsText) == 1
    for cc = 1:numel(headers)
        if isfield(rtv.ColumnIsText, headers{cc})
            textcols(cc) = rtv.ColumnIsText.(headers{cc});
        end
    end
end
table = cell(size(data, 1), numel(headers));
nlabels = numel(labels);
formats = cell(1, numel(headers));
csizes = zeros(1, numel(headers));
for cc = 1:numel(headers)
    if textcols(cc)
        for rc = 1:size(data, 1)
            if data(rc, cc) >= 1 && ...
                data(rc, cc) <= nlabels
                table{rc, cc} = labels{data(rc, cc)};
            else
                table{rc, cc} = sprintf('%.3g', data(rc, cc));
            end
        end
    else
        for rc = 1:size(data, 1)
            table{rc, cc} = sprintf('%.3g', data(rc, cc));
        end
    end
    csizes(cc) = max(numel(headers{cc}), size(char(table(:, cc)), 2));
    formats{cc} = sprintf('%%%ds', csizes(cc));
end
tformat = gluetostringc(formats, ' | ');

% print together
points = cell(size(table, 1) + 2, 1);
points{1} = sprintf(tformat, headers{:});
for cc = 1:numel(headers)
    headers{cc} = repmat('-', 1, csizes(cc));
end
points{2} = sprintf(tformat, headers{:});
for rc = 1:size(data, 1)
    points{rc+2} = sprintf(tformat, table{rc, :});
end

% set to points list
ch.Points.Value = [];
ch.Points.String = points;
if ch.SelectedPts.Value <= 0 && ...
    any(~sel)
    ch.Points.Value = find(sel) + 2;
    ch.Points.ListboxTop = findfirst(sel) + 1;
else
    ch.Points.ListboxTop = 1;
end
ch.Points.UserData = {condstr, datasel};
ch.PointsLabel.String = sprintf( ...
    'Selected points: (%d of %d points, %d of %d units selected; selecting has no effect!)', ...
    sum(sel), onumsel, selunits, numunits);

% update study selection
studies = ch.Studies.String;
if ~iscell(studies)
    studies = cellstr(studies);
end
studies = studies(ch.Studies.Value);
plp.RunTimeVars.StudySelection = studies(:);

% allow next update
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_listpoints')) = [];

% read out analysis
ne_mkda_updana;

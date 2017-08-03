function [t, tlist] = plp_CPointsAtCoord(xo, crd, opts)
% PLP::CPointsAtCoord  - create list of units with points at coordinate/s
%
% FORMAT:       [text, tlist] = plp.CPointsAtCoord(crd [, opts])
%
% Input fields:
%
%       crd         Cx3 coordinate(s) (in same notation as .X, .Y, .Z)
%       opts        optional settings
%        .columns   columns to output (default: {'Study', 'Contrast'})
%        .columntxt flag which columns are text (default: [true, false])
%        .cond      conditional statement, e.g.
%                   '$Study >= 1 & $Study <= 3 & $Type == 2'
%                   whereas $names are replaced by their column data
%                   and syntax as in '$Column == ''content''' is translated
%                   into a call to strcmpi (or regexpi for patterns)
%        .fringe    minimal distance of any point/coordinate (default: 5)
%        .gridres   coordinate grid (source map) resolution (default: 3)
%        .unitcol   column which represents units, default: 'Contrast'
%        .unitsel   unit selection (default: all unique values)
%        .uvmp      VMP object to extract data from
%        .uvmpmaps  maps in VMP that correspond to unitsel (default: first)
%        .uvmptype  type of data, either of {'max'}, 'mean', or 'min'
%
% Output fields:
%
%       text        text output
%       tlist       cell array with output arguments
%
% Using: bvcoordconv, findfirst, gluetostringc, minmaxmean, multimatch,
%        psetdists.

% Version:  v1.1
% Build:    16021210
% Date:     Feb-12 2016, 10:46 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
findfirst  = ne_methods.findfirst;
minmaxmean = ne_methods.minmaxmean;
psetdists  = ne_methods.psetdists;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'plp') || ...
   ~isa(crd, 'double') || isempty(crd) || ndims(crd) > 2 || ...
    size(crd, 2) ~= 3 || any(isinf(crd(:)) | isnan(crd(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
cn = bc.ColumnNames(:);
labels = bc.Labels(:);
p = bc.Points;
mcrd = mean(crd, 1);
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'columns') || ~iscell(opts.columns) || ...
   ~all(cellfun('ischar', opts.columns(:))) || any(cellfun('isempty', opts.columns(:)))
    opts.columns = {'Study'; 'Contrast'};
else
    opts.columns = opts.columns(:);
end
if ~isfield(opts, 'columntxt') || ~islogical(opts.columntxt) || ...
    numel(opts.columntxt) ~= numel(opts.columns)
    if numel(opts.columns) == 2 && all(strcmpi(opts.columns, {'study'; 'contrast'}))
        opts.columntxt = [true; false];
    else
        opts.columntxt = false(numel(opts.columns), 1);
    end
else
    opts.columntxt = opts.columntxt(:);
end
cnm = (ne_methods.multimatch(opts.columns, cn) < 1);
opts.columns(cnm) = [];
opts.columntxt(cnm) = [];
if isempty(opts.columns)
    error('neuroelf:xff:badObject', ...
        'PLP object requires at least one of the named columns.');
end
if ~isfield(opts, 'cond') || ~ischar(opts.cond)
    opts.cond = '';
    ocond = '';
else
    opts.cond = ['(' lower(opts.cond(:)') ')'];
    ocond = opts.cond;
    ocond(ocond == ' ') = [];
    while ~isempty(ocond) && ocond(1) == '(' && ocond(end) == ')'
        ocond = ocond(2:end-1);
    end
    cregx = regexp(opts.cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*==\s*''[^'']+'')', 'tokens');
    while ~isempty(cregx) && ~isempty(cregx{1})
        cregp = regexp(cregx{1}{1}, '^\$([a-zA-Z][a-zA-Z_0-9]*)\s*==\s*''([^'']+)''$', 'tokens');
        if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
            error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
        end
        if any(cregp{1}{2} == '*')
            opts.cond = strrep(opts.cond, cregx{1}{1}, sprintf( ...
                '~cellfun(''isempty'', regexpi(labels($%s), ''%s''))', cregp{1}{:}));
        else
            opts.cond = strrep(opts.cond, cregx{1}{1}, sprintf( ...
                'strcmpi(labels($%s), ''%s'')', cregp{1}{:}));
        end
        cregx = regexp(opts.cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*==\s*''[^'']+'')', 'tokens');
    end
    cregx = regexp(opts.cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*~=\s*''[^'']+'')', 'tokens');
    while ~isempty(cregx) && ~isempty(cregx{1})
        cregp = regexp(cregx{1}{1}, '^\$([a-zA-Z][a-zA-Z_0-9]*)\s*~=\s*''([^'']+)''$', 'tokens');
        if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
            error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
        end
        if any(cregp{1}{2} == '*')
            opts.cond = strrep(opts.cond, cregx{1}{1}, sprintf( ...
                'cellfun(''isempty'', regexpi(labels($%s), ''%s''))', cregp{1}{:}));
        else
            opts.cond = strrep(opts.cond, cregx{1}{1}, sprintf( ...
                '~strcmpi(labels($%s), ''%s'')', cregp{1}{:}));
        end
        cregx = regexp(opts.cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*~=\s*''[^'']+'')', 'tokens');
    end
end
cnl = zeros(numel(cn), 1);
for cnc = 1:numel(cn)
    cnl(cnc) = numel(cn{cnc});
end
[cnl, cni] = sort(cnl, 'descend');
cns = cn(cni);
for cnc = 1:numel(cn)
    opts.cond = strrep(opts.cond, ['$' lower(cns{cnc})], sprintf('p(:,%d)', cni(cnc)));
end
if isempty(opts.cond)
    opts.cond = '(true)';
end
if ~strcmp(opts.cond, '(true)')
    try
        pl = eval(opts.cond);
    catch xfferror
        error('neuroelf:xff:BadArgument', 'Bad condition given: %s.', xfferror.message);
    end
    p = p(pl, :);
end
if ~isfield(opts, 'fringe') || ~isa(opts.fringe, 'double') || numel(opts.fringe) ~= 1 || ...
    isinf(opts.fringe) || isnan(opts.fringe) || opts.fringe < 2
    opts.fringe = 2;
end
if ~isfield(opts, 'gridres') || ~isa(opts.gridres, 'double') || numel(opts.gridres) ~= 1 || ...
    isinf(opts.gridres) || isnan(opts.gridres) || opts.gridres < 0
    opts.gridres = 3;
end
gridres = sqrt(opts.gridres) + sqrt(eps);
if ~isfield(opts, 'unitcol') || ~ischar(opts.unitcol) || ...
    isempty(opts.unitcol) || ~any(strcmpi(cn, opts.unitcol(:)'))
    if any(strcmpi(cn, 'contrast'))
        opts.unitcol = 'Contrast';
    elseif any(strcmpi(cn, 'study'))
        opts.unitcol = 'Study';
    else
        error('neuroelf:xff:badObject', ...
            'PLP object requires either Contrast or Study column.');
    end
else
    opts.unitcol = opts.unitcol(:)';
end
ucol = findfirst(strcmpi(cn, opts.unitcol));
uselgiven = true;
if ~isfield(opts, 'unitsel') || ~isa(opts.unitsel, 'double') || ...
   ~any(size(opts.unitsel) == numel(opts.unitsel)) || ...
    any(isinf(opts.unitsel) | isnan(opts.unitsel))
    opts.unitsel = unique(p(:, ucol));
    uselgiven = false;
end
if ~isfield(opts, 'uvmp') || numel(opts.uvmp) ~= 1 || ~xffisobject(opts.uvmp, true, 'vmp')
    vc = struct('Map', struct('Name', '', 'Type', []));
    vc.Map(:) = [];
else
    vcx = unique(ne_methods.bvcoordconv(crd, 'tal2bvx', aft_BoundingBox(opts.uvmp)));
    vc = opts.uvmp.C;
end
if ~isfield(opts, 'uvmpmaps') || ~isa(opts.uvmpmaps, 'double') || ...
    numel(opts.uvmpmaps) ~= numel(opts.unitsel) || ...
    any(isinf(opts.uvmpmaps(:)) | isnan(opts.uvmpmaps(:)) | ...
        opts.uvmpmaps(:) ~= fix(opts.uvmpmaps(:)) | opts.uvmpmaps(:) < 1) || ...
    numel(unique(opts.uvmpmaps(:))) ~= numel(opts.unitsel) || ...
    any(opts.uvmpmaps(:) > numel(vc.Map))
    if uselgiven
        opts.uvmpmaps = (1:numel(opts.unitsel))';
    elseif ~isempty(vc.Map)
        mnum = (1:numel(vc.Map))';
        mtype = cat(1, vc.Map.Type);
        mnum(mtype ~= 16) = [];
        mrtv = {vc.Map(mnum).RunTimeVars};
        muse = zeros(numel(opts.unitsel), 1);
        for cc = 1:numel(mrtv)
            mcond = mrtv{cc}.Condition(mrtv{cc}.Condition ~= ' ');
            while ~isempty(mcond) && mcond(1) == '(' && mcond(end) == ')'
                mcond = mcond(2:end-1);
            end
            if mrtv{cc}.UnitColumn == ucol && strcmpi(ocond, mcond) && ...
                any(opts.unitsel == mrtv{cc}.UnitID)
                muse(opts.unitsel == mrtv{cc}.UnitID) = mnum(cc);
            end
        end
        if ~all(muse) > 0
            opts.uvmpmaps = [];
            vc.Map(:) = [];
        else
            opts.uvmpmaps = muse;
        end
    else
        opts.uvmpmaps = [];
    end
    if ~isempty(opts.uvmpmaps) && opts.uvmpmaps(end) > numel(vc.Map)
        vc.Map(:) = [];
        opts.uvmpmaps = [];
    end
end
if ~isfield(opts, 'uvmptype') || ~ischar(opts.uvmptype) || isempty(opts.uvmptype) || ...
   ~any(strcmpi(opts.uvmptype(:)', {'max', 'mean', 'min'}))
    uvtyp = 'a';
else
    uvtyp = lower(opts.uvmptype(2));
end
[unitsel, uia] = intersect(opts.unitsel(:), unique(p(:, ucol)));
if ~isempty(vc.Map)
    opts.uvmpmaps = opts.uvmpmaps(uia);
end

% get column numbers
tcol = zeros(1, numel(opts.columns));
tco2 = numel(tcol) + 1;
tco3 = numel(tcol) + 2;
for cc = 1:numel(tcol)
    tcol(cc) = findfirst(strcmpi(cn, opts.columns{cc}));
end
xyz = [findfirst(strcmpi(cn, 'x')), findfirst(strcmpi(cn, 'y')), findfirst(strcmpi(cn, 'z'))];

% create output
tlist = cell(numel(unitsel), numel(tcol) + 4 + double(~isempty(vc.Map)));

% iterate over each unit
cmdist = Inf * ones(numel(unitsel), 1);
cmtake = false(size(cmdist));
for cc = 1:numel(unitsel)

    % find rows that match
    ri = find(p(:, ucol) == unitsel(cc));

    % no points, nothing to do
    if isempty(ri)
        continue;
    end

    % compute minimal distance between sets of points
    pdist = min(psetdists(crd, p(ri, xyz)), [], 1);
    pdist(pdist <= gridres) = 0;
    if sum(pdist == 0) > 1
        [pdist, pdi] = min(psetdists(mcrd, p(ri, xyz)));
        pdist = 0;
    else
        [pdist, pdi] = min(pdist);
    end
    cmdist(cc) = pdist;

    % use this
    if cmdist(cc) <= opts.fringe
        cmtake(cc) = true;

        % fill data
        tlist{cc, tco2} = cmdist(cc);
        tlist(cc, end-2:end) = {p(ri(pdi), xyz(1)), p(ri(pdi), xyz(2)), p(ri(pdi), xyz(3))};
        for c = 1:numel(tcol)
            if opts.columntxt(c)
                try
                    tlist{cc, c} = labels{p(ri(pdi), tcol(c))};
                catch xfferror
                    neuroelf_lasterr(xfferror);
                    tlist{cc, c} = 'n/a';
                end
            else
                tlist{cc, c} = p(ri(1), tcol(c));
            end
        end

        % extract data
        if ~isempty(vc.Map)
            vmd = minmaxmean(vc.Map(opts.uvmpmaps(cc)).VMPData(vcx), 4);
            switch (uvtyp)
                case 'a'
                    tlist{cc, tco3} = vmd(2);
                case 'e'
                    tlist{cc, tco3} = vmd(3);
                case 'i'
                    tlist{cc, tco3} = vmd(1);
            end
        end
    end
end

% shrink list
tlist = tlist(cmtake, :);

% create text
text = cell(size(tlist, 1), 1);
patt = cell(1, size(tlist, 2));
for c = 1:numel(tcol)
    if opts.columntxt(c)
        patt{c} = sprintf('%%-%ds', size(char(tlist(:, c)), 2) + 1);
    else
        cdata = cat(1, tlist{:, c});
        if all(cdata == fix(cdata))
            patt{c} = sprintf('%%%dd', ceil(log10(max(abs(cdata)) + 1)));
        else
            patt{c} = '%g';
        end
    end
end
patt{tco2} = '(%6.3f)';
if ~isempty(vc.Map)
    patt{tco2} = '(%6.3f,';
    patt{tco3} = '%6.4f)';
end
patt(end-2:end) = {'[%3.0f,', '%3.0f,', '%3.0f]'};
patt = ne_methods.gluetostringc(patt, ' ');
for cc = 1:numel(text)
    text{cc} = sprintf(patt, tlist{cc, :});
end
t = ne_methods.gluetostringc(text, char(10));

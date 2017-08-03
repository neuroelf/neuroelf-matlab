function [ph, th] = plp_PlotOnSurface(xo, opts)
% PLP::PlotOnSurface  - plot points on a surface
%
% FORMAT:       [ph, tx] = plp.PlotOnSurface(opts)
%
% Input fields:
%
%       opts        options
%        .addrand   add random displacement (useful for clusters)
%        .axes      axes object, if not given, create new from GUI
%        .bcolor    optional 1x3 RGB color for border of dots (if filled)
%        .bgcolor   background RGB color
%        .bothhemi  force both hemispheres to be drawn
%        .color     either column name/number or 1x3 RGB color
%        .cond      conditional statement, e.g.
%                   '$Study >= 1 & $Study <= 3 & $Type == 2'
%                   whereas names are replaced by their column data
%        .enum      enumerate points (can be a double as start value)
%        .filled    logical value (default: true)
%        .hover     plot points N mm "hovering in front" of coordinates
%        .label     either column name/number or fixed char
%        .labcolor  either column name/number or 1x3 RGB color
%        .labfont   label font name (empty for default)
%        .labfontsz label font size (empty for default)
%        .labshift  either mm value (radius from center) or X/Y screen vals
%        .srf       single filename or cell array with filenames
%        .symbol    either symbol number or 1x1 char symbol
%        .symsize   either column name/number or negative for fixed size
%        .view      either of (mandatory field)
%                   'll' - left lateral
%                   'lm' - left medial
%                   'rl' - right lateral
%                   'rm' - right medial
%                   'fr' - frontal
%                   'oc' - occipetal
%                   'in' - inferior
%                   'su' - superior
%        .xrange    1x2 double, restrict X to between values
%        .yrange    1x2 double, restrict Y to between values
%        .zrange    1x2 double, restrict Z to between values
%
% Output fields:
%
%       ph          plot handles
%       tx          text handles (if any)
%
% Using: findfirst, spmtrf.

% Version:  v1.1
% Build:    16020112
% Date:     Feb-01 2016, 12:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library
global ne_methods;
findfirst = ne_methods.findfirst;
spmtrf    = ne_methods.spmtrf;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'plp') || ...
    numel(opts) ~= 1 || ~isstruct(opts) || ~isfield(opts, 'view') || ...
   ~ischar(opts.view) || ~any(strcmpi(opts.view(:)', ...
        {'f', 'fr', 'i', 'in', 'll', 'lm', 'o', 'oc', 'rl', 'rm', 's', 'su'}))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
rtv = bc.RunTimeVars;
if ~isfield(opts, 'addrand') || ~isa(opts.addrand, 'double') || numel(opts.addrand) ~= 1 || ...
    isinf(opts.addrand) || isnan(opts.addrand) || opts.addrand <= 0
    opts.addrand = 0;
end
if ~isfield(opts, 'axes') || (~isa(opts.axes, 'double') && ...
    ~isa(opts.axes, 'matlab.graphics.axis.Axes')) || ...
    numel(opts.axes) ~= 1 || ~ishandle(opts.axes) || ...
   ~strcmpi(get(opts.axes, 'Type'), 'axes')
    opts.axes = [];
end
if ~isfield(opts, 'bcolor') || ~isa(opts.bcolor, 'double') || numel(opts.bcolor) ~= 3 || ...
    any(isinf(opts.bcolor) | isnan(opts.bcolor) | opts.bcolor < 0 | opts.bcolor > 255)
    opts.bcolor = [];
elseif any(opts.bcolor > 1)
    opts.bcolor = (1 / 255) .* opts.bcolor(:)';
else
    opts.bcolor = opts.bcolor(:)';
end
if ~isfield(opts, 'bgcolor') || ~isa(opts.bgcolor, 'double') || numel(opts.bgcolor) ~= 3 || ...
    any(isinf(opts.bgcolor) | isnan(opts.bgcolor) | opts.bgcolor < 0 | opts.bgcolor > 255)
    opts.bgcolor = [];
else
    opts.bgcolor = (1 / 255) .* opts.bgcolor(:)';
end
if ~isfield(opts, 'bothhemi') || ~islogical(opts.bothhemi) || numel(opts.bothhemi) ~= 1
    opts.bothhemi = false;
end
if ~isfield(opts, 'color') || (~isa(opts.color, 'double') && ~ischar(opts.color))
    opts.color = findfirst(strcmpi(bc.ColumnNames, 'color'));
    if isempty(opts.color) || opts.color > size(bc.Points, 2)
        opts.color = [255, 0, 0];
    end
end
if ischar(opts.color)
    opts.color = findfirst(strcmpi(bc.ColumnNames, opts.color(:)'));
    if isempty(opts.color)
        opts.color = [255, 0, 0];
    end
elseif isa(opts.color, 'double') && (~any(numel(opts.color) == [1, 3]) || ...
    any(isinf(opts.color) | isnan(opts.color) | opts.color < 0))
    opts.color = [255, 0, 0];
elseif isa(opts.color, 'double') && numel(opts.color) == 1 && ...
   ~any(1:size(bc.Points, 2) == opts.color)
    opts.color = [255, 0, 0];
elseif isa(opts.color, 'double') && numel(opts.color) == 3
    opts.color = min(255, round(opts.color(:)'));
end
if ~isfield(opts, 'cond') || ~ischar(opts.cond)
    opts.cond = '';
else
    opts.cond = ['(' lower(opts.cond(:)') ')'];
end
if ~isfield(opts, 'enum') || numel(opts.enum) ~= 1 || (~islogical(opts.enum) && ~isa(opts.enum, 'double'))
    opts.enum = 0;
elseif islogical(opts.enum)
    opts.enum = double(opts.enum);
elseif isinf(opts.enum) || isnan(opts.enum) || opts.enum < 1
    opts.enum = 0;
else
    opts.enum = floor(opts.enum);
end
if opts.enum > 0
    opts.label = [];
end
if ~isfield(opts, 'filled') || ~islogical(opts.filled) || numel(opts.filled) ~= 1
    opts.filled = {'filled'};
else
    if opts.filled
        opts.filled = {'filled'};
    else
        opts.filled = {};
    end
end
if ~isfield(opts, 'hover') || ~isa(opts.hover, 'double') || numel(opts.hover) ~= 1 || ...
    isinf(opts.hover) || isnan(opts.hover)
    opts.hover = 5;
end
if ~isfield(opts, 'label') || (~isa(opts.label, 'double') && ~ischar(opts.label))
    opts.label = [];
end
opts.labeltxt = false;
if ischar(opts.label)
    labeli = findfirst(strcmpi(bc.ColumnNames, opts.label(:)'));
    if ~isempty(labeli) && labeli <= size(bc.Points, 2) && ~isempty(opts.label)
        opts.label = labeli;
        try
            opts.labeltxt = rtv.ColumnIsText.(bc.ColumnNames{labeli});
        catch xfferror
            neuroelf_lasterr(xfferror);
        end
    elseif isempty(opts.label)
        opts.label = [];
    end
elseif isa(opts.label, 'double') && (numel(opts.label) ~= 1 || isinf(opts.label) || isnan(opts.label) || ...
    ~any((1:size(bc.Points, 2)) == opts.label))
    opts.label = [];
elseif isa(opts.label, 'double')
    try
        opts.labeltxt = rtv.ColumnIsText.(bc.ColumnNames{opts.label});
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end
if ~isfield(opts, 'labcolor') || (~isa(opts.labcolor, 'double') && ~ischar(opts.labcolor))
    opts.labcolor = [];
end
if ischar(opts.labcolor)
    opts.labcolor = findfirst(strcmpi(bc.ColumnNames, opts.labcolor(:)'));
    if isempty(opts.labcolor) || opts.labcolor > size(bc.Points, 2)
        opts.labcolor = [];
    end
elseif isa(opts.labcolor, 'double') && (~any(numel(opts.labcolor) == [1, 3]) || ...
    any(isinf(opts.labcolor) | isnan(opts.labcolor) | opts.labcolor < 0))
    opts.labcolor = [];
elseif isa(opts.labcolor, 'double') && numel(opts.labcolor) == 1 && ~any(1:size(bc.Points, 2) == opts.labcolor)
    opts.labcolor = [];
elseif isa(opts.labcolor, 'double') && numel(opts.labcolor) == 3
    opts.labcolor = min(255, round(opts.labcolor(:)'));
end
if ~isfield(opts, 'labfont') || ~ischar(opts.labfont) || isempty(opts.labfont)
    opts.labfont = '';
end
if ~isfield(opts, 'labfontsz') || ~isa(opts.labfontsz, 'double') || numel(opts.labfontsz) ~= 1 || ...
    isinf(opts.labfontsz) || isnan(opts.labfontsz) || opts.labfontsz <= 0
    opts.labfontsz = [];
end
if ~isfield(opts, 'labshift') || ~isa(opts.labshift, 'double') || ~any([1, 2] == numel(opts.labshift)) || ...
    any(isinf(opts.labshift) | isnan(opts.labshift))
    opts.labshift = [4, -4];
elseif numel(opts.labshift) == 1
    opts.labshift = abs(opts.labshift);
end
if ~isfield(opts, 'srf') || isempty(opts.srf) || (~ischar(opts.srf) && ~iscell(opts.srf))
    opts.srf = {};
elseif ischar(opts.srf)
    opts.srf = {opts.srf(:)'};
end
for srfc = numel(opts.srf):-1:1
    if ~ischar(opts.srf{srfc}) || isempty(opts.srf{srfc}) || exist(opts.srf{srfc}(:)', 'file') ~= 2
        opts.srf(srfc) = [];
    else
        [isabs, opts.srf{srfc}] = isabsolute(opts.srf{srfc}(:)');
        opts.srf{srfc} = strrep(opts.srf{srfc}, '\', '/');
    end
end
if ~isfield(opts, 'symbol') || numel(opts.symbol) ~= 1 || ...
   (~isa(opts.symbol, 'double') && ~ischar(opts.symbol)) || ...
   (isa(opts.symbol, 'double') && (isinf(opts.symbol) || isnan(opts.symbol) || ...
    opts.symbol < 1 || opts.symbol > numel(bc.Symbols) || opts.symbol ~= fix(opts.symbol)))
    opts.symbol = 'o';
elseif isa(opts.symbol, 'double')
    opts.symbol = bc.Symbols{opts.symbol};
end
if ~isfield(opts, 'symsize') || (~isa(opts.symsize, 'double') && ~ischar(opts.symsize))
    opts.symsize = 'Size';
end
if ischar(opts.symsize)
    opts.symsize = findfirst(strcmpi(bc.ColumnNames, opts.symsize(:)'));
    if isempty(opts.symsize) || opts.symsize > size(bc.Points, 2)
        opts.symsize = [];
    end
elseif isa(opts.symsize, 'double') && (numel(opts.symsize) > 1 || ...
    any(isinf(opts.symsize) | isnan(opts.symsize) | opts.symsize == 0))
    opts.symsize = [];
end
if isa(opts.symsize, 'double') && ~isempty(opts.symsize)
    if opts.symsize > 0 && ~any((1:size(bc.Points, 2)) == opts.symsize)
        opts.symsize = [];
    end
end
opts.view = lower(opts.view(:)');
if ~isfield(opts, 'xrange')
    if opts.view(1) == 'l'
        opts.xrange = [-Inf, round(0.5 * opts.hover)];
    elseif opts.view(1) == 'r'
        opts.xrange = [-round(0.5 * opts.hover), Inf];
    end
end
if ~isfield(opts, 'xrange') || ~isa(opts.xrange, 'double') || numel(opts.xrange) ~= 2 || ...
    any(isnan(opts.xrange)) || opts.xrange(1) > opts.xrange(2)
    opts.xrange = [];
end
if ~isfield(opts, 'yrange') || ~isa(opts.yrange, 'double') || numel(opts.yrange) ~= 2 || ...
    any(isnan(opts.yrange)) || opts.yrange(1) > opts.yrange(2)
    opts.yrange = [];
end
if ~isfield(opts, 'zrange') || ~isa(opts.zrange, 'double') || numel(opts.zrange) ~= 2 || ...
    any(isnan(opts.zrange)) || opts.zrange(1) > opts.zrange(2)
    opts.zrange = [];
end

% get points
p = bc.Points;

% enumerate
if opts.enum > 0
    p(:, end + 1) = opts.enum:(opts.enum + size(p, 1) - 1);
end

% parse condition
labels = bc.Labels(:);
cond = opts.cond;
cregx = regexp(cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*==\s*''[^'']+'')', 'tokens');
while ~isempty(cregx) && ~isempty(cregx{1})
    cregp = regexp(cregx{1}{1}, '^\$([a-zA-Z][a-zA-Z_0-9]*)\s*==\s*''([^'']+)''$', 'tokens');
    if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
        error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
    end
    if any(cregp{1}{2} == '*')
        cond = strrep(cond, cregx{1}{1}, sprintf( ...
            '~cellfun(''isempty'', regexpi(labels($%s), ''%s''))', cregp{1}{:}));
    else
        cond = strrep(cond, cregx{1}{1}, sprintf( ...
            'strcmpi(labels($%s), ''%s'')', cregp{1}{:}));
    end
    cregx = regexp(cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*==\s*''[^'']+'')', 'tokens');
end
cregx = regexp(cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*~=\s*''[^'']+'')', 'tokens');
while ~isempty(cregx) && ~isempty(cregx{1})
    cregp = regexp(cregx{1}{1}, '^\$([a-zA-Z][a-zA-Z_0-9]*)\s*~=\s*''([^'']+)''$', 'tokens');
    if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
        error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
    end
    if any(cregp{1}{2} == '*')
        cond = strrep(cond, cregx{1}{1}, sprintf( ...
            'cellfun(''isempty'', regexpi(labels($%s), ''%s''))', cregp{1}{:}));
    else
        cond = strrep(cond, cregx{1}{1}, sprintf( ...
            '~strcmpi(labels($%s), ''%s'')', cregp{1}{:}));
    end
    cregx = regexp(cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*~=\s*''[^'']+'')', 'tokens');
end

% replace in condition
cns = bc.ColumnNames;
cnl = zeros(numel(cns), 1);
for cn = 1:numel(cns)
    cnl(cn) = numel(cns{cn});
end
[cnl, cni] = sort(cnl, 'descend');
cns = cns(cni);
for cn = 1:numel(cns)
    cond = strrep(cond, ['$' lower(cns{cn})], sprintf('p(:,%d)', cni(cn)));
end
if isempty(cond)
    cond = '(true)';
end

% add range to condition if needed
xc = findfirst(strcmpi(bc.ColumnNames, 'x'));
if isempty(xc)
    xc = 1;
end
if ~isempty(opts.xrange)
    cond = [cond sprintf(' & p(:,%d)>=%f & p(:,%d)<=%f', ...
        xc, opts.xrange(1), xc, opts.xrange(2))];
end
yc = findfirst(strcmpi(bc.ColumnNames, 'y'));
if isempty(yc)
    yc = 2;
end
if ~isempty(opts.yrange)
    cond = [cond sprintf(' & p(:,%d)>=%f & p(:,%d)<=%f', ...
        yc, opts.yrange(1), yc, opts.yrange(2))];
end
zc = findfirst(strcmpi(bc.ColumnNames, 'z'));
if isempty(zc)
    zc = 3;
end
if ~isempty(opts.zrange)
    cond = [cond sprintf(' & p(:,%d)>=%f & p(:,%d)<=%f', ...
        zc, opts.zrange(1), zc, opts.zrange(2))];
end

% try to parse condition
if strcmp(cond, '(true)')
    pl = 1:size(p, 1);
else
    try
        pl = find(eval(cond));
    catch xfferror
        error('neuroelf:xff:badArgument', 'Bad condition given: %s.', xfferror.message);
    end
end
p = p(pl, :);
np = size(p, 1);
opts.pcrd = [p(:, xc), p(:, yc), p(:, zc), ones(np, 1)];

% add random value?
if opts.addrand > 0
    opts.pcrd = opts.pcrd + (opts.addrand .* rand(size(opts.pcrd)) - 0.5 * opts.addrand);
end

% try to get colors
if numel(opts.color) == 1
    try
        opts.color = bc.Colors(p(:, opts.color), :);
    catch xfferror
        neuroelf_lasterr(xfferror);
        if ~isempty(bc.Colors)
            opts.color = bc.Color(ones(np, 1), :);
        else
            opts.color = 255 .* [ones(np, 1), zeros(np, 2)];
        end
    end
else
    opts.color = opts.color(ones(np, 1), :);
end
opts.color = (1 / 255) .* opts.color;

% labels
if ~isempty(opts.label)

    % prepare labels
    if ischar(opts.label)
        opts.label = repmat({opts.label}, np, 1);
    else
        try
            if opts.labeltxt
                opts.label = bc.Labels(p(:, opts.label));
            else
                lcci = opts.label;
                opts.label = cell(size(p, 1), 1);
                for lc = 1:numel(opts.label)
                    opts.label{lc} = sprintf('%g', p(lc, lcci));
                end
            end
        catch xfferror
            neuroelf_lasterr(xfferror);
            opts.label;
        end
    end
end

% try to get label colors
if numel(opts.labcolor) == 1
    try
        opts.labcolor = bc.Colors(p(:, opts.labcolor), :);
    catch xfferror
        neuroelf_lasterr(xfferror);
        opts.labcolor = [];
    end
elseif numel(opts.labcolor) == 3
    opts.labcolor = opts.labcolor(ones(np, 1), :);
end
opts.labcolor = (1 / 255) .* opts.labcolor;

% symbol size
if ~isempty(opts.symsize) && opts.symsize > 0
    try
        opts.symsize = p(:, opts.symsize);
        opts.symsize(opts.symsize < 0) = bc.SymbolSizes(-opts.symsize(opts.symsize < 0));
    catch xfferror
        neuroelf_lasterr(xfferror);
        opts.symsize = [];
    end
elseif numel(opts.symsize) == 1
    opts.symsize = repmat(abs(opts.symsize), size(p, 1), 1);
else
    opts.symsize = repmat(10, size(p, 1), 1);
end
opts.symsize(opts.symsize <= 0) = 40;

% make sure the required files are loaded
try
    [scn, scs] = neuroelf_gui('scenerylist');
catch xfferror
    rethrow(xfferror);
end
for lc = 1:size(scn, 1)
    scn{lc, 3} = aft_FilenameOnDisk(scn{lc, 4});
end
if ~isempty(scn)
    scn(:, 3) = strrep(scn(:, 3), '\', '/');
end

% for empty list of requested files
if isempty(opts.srf)

    % use the LH and RH colin brain hemispheres
    lh = strrep([neuroelf_path('colin') '/colin_LH_SPH.srf'], '\', '/');
    rh = strrep([neuroelf_path('colin') '/colin_RH_SPH.srf'], '\', '/');
    if ~any(strcmpi(lh, scn(:, 3)))
        try
            neuroelf_gui('openfile', lh);
        catch xfferror
            oeo = xfferror;
            try
                neuroelf_gui('scenerylist', 'set', scs);
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
            rethrow(oeo);
        end
    end
    if ~any(strcmpi(rh, scn(:, 3)))
        try
            neuroelf_gui('openfile', rh);
        catch xfferror;
            oeo = xfferror;
            try
                neuroelf_gui('scenerylist', 'set', scs);
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
            rethrow(oeo);
        end
    end

    % get updated scenery list
    scn = neuroelf_gui('scenerylist');
    for lc = 1:size(scn, 1)
        scn{lc, 3} = aft_FilenameOnDisk(scn{lc, 4});
    end
    scn(:, 3) = strrep(scn(:, 3), '\', '/');
    lh = findfirst(strcmpi(lh, scn(:, 3)));
    rh = findfirst(strcmpi(rh, scn(:, 3)));
    mh = [];

else

    % no LH/RH
    lh = [];
    rh = [];

    % iterate over surface files
    for srfc = 1:numel(opts.srf)
        if ~any(strcmpi(opts.srf{srfc}, scn(:, 3)))
            try
                neuroelf_gui('openfile', opts.srf{srfc});
            catch xfferror
                oeo = xfferror;
                try
                    neuroelf_gui('scenerylist', 'set', scs);
                catch xfferror
                    neuroelf_lasterr(xfferror);
                end
                rethrow(oeo);
            end
        end
    end

    % then re-find
    scn = neuroelf_gui('scenerylist');
    for lc = 1:size(scn, 1)
        scn{lc, 3} = aft_FilenameOnDisk(scn{lc, 4});
    end
    scn(:, 3) = strrep(scn(:, 3), '\', '/');
    for srfc = 1:numel(opts.srf)
        srfi = findfirst(strcmpi(opts.srf{srfc}, scn(:, 3)));
        if isempty(srfi)
            error('neuroelf:xff:NeuroElfGUIError', ...
                'Error loading specified surface file (%s).', opts.srf{srfc});
        end
        opts.srf{srfc} = srfi;
    end
    mh = cat(2, opts.srf{:});
end

% check axes if given
if ~isempty(opts.axes)
    opts.hFig.MLHandle = get(opts.axes, 'Parent');
    if ~strcmpi(get(opts.hFig.MLHandle, 'Type'), 'figure') || ...
        isempty(regexpi(get(opts.hFig.MLHandle, 'Tag'), 'BS[0-9a-f]+_Figure'))
        opts.axes = [];
    else
        axc = get(opts.axes, 'Children');
        if isempty(axc) || ~strcmpi(get(axc(1), 'Type'), 'patch')
            opts.axes = [];
        else
            opts.hFig = xfigure(opts.hFig.MLHandle);
        end
    end
end

% preset angles
anglez = 0;
angler = 0;

% no valid axes available
if isempty(opts.axes)

    % get current viewpoint
    viewp = neuroelf_gui('scenerytrans');

    % depending on view
    switch (opts.view)

        % frontal view
        case {'f', 'fr'}

            % viewpoint settings
            slist = [lh, rh, mh];
            anglez = 90;

        % inferior
        case {'i', 'in'}
            slist = [lh, rh, mh];
            anglez = 90;
            angler = -90;

        % left lateral view
        case {'ll'}
            slist = [lh, mh];
            anglez = 180;

        % left medial
        case {'lm'}
            slist = [lh, mh];
            anglez = 0;

        % occipetal
        case {'o', 'oc'}
            slist = [lh, rh, mh];
            anglez = -90;

        % right lateral
        case {'rl'}
            slist = [rh, mh];
            anglez = 0;

        % right medial
        case {'rm'}
            slist = [lh, mh];
            anglez = 180;

        % superior
        case {'s', 'su'}
            slist = [lh, rh, mh];
            anglez = -90;
            angler = 90;
    end

    % other surfaces altogether
    if opts.bothhemi
        slist = [slist(:); lh(:); rh(:)];
    end

    % every surface only once!
    slist = unique(slist(:));

    % set scenery viewpoint
    neuroelf_gui('scenerylist', 'set', slist);
    neuroelf_gui('scenerytrans', {[0, 0, 0], [anglez, angler], 1});

    % make sure screen is up-to-date!
    drawnow;
    pause(0.2);

    % undock
    [opts.hFig] = neuroelf_gui('undock');

    % set viewpoint back
    neuroelf_gui('scenerytrans', viewp);

    % bring figure to front
    figure(opts.hFig.MLHandle);
    drawnow;

    % get axes object
    tags = opts.hFig.TagStruct;
    tagf = fieldnames(tags);
    ax = tags.(tagf{find(~cellfun('isempty', regexpi(tagf, 'AX_.*_Slice_Zoom')))}).MLHandle;
    opts.axes = ax;
end

% perform transformation
trf = spmtrf([0, 0, 0], [0, 0, (180 / pi) * anglez]);
tpcrd = opts.pcrd * trf';

% hovering ?
if opts.hover > 0
    tpcrd(:, 1) = tpcrd(:, 1) + opts.hover;
elseif opts.hover < 0
    tpcrd(:, 1) = -opts.hover;
end

% background color ?
if ~isempty(opts.bgcolor)
    if strcmpi(get(opts.axes, 'Visible'), 'on')
        set(opts.axes, 'Color', opts.bgcolor);
    else
        set(get(opts.axes, 'Parent'), 'Color', opts.bgcolor);
    end
end

% force hold
hold(opts.axes, 'on');

% plot points
ph = scatter3(opts.axes, tpcrd(:, 1), tpcrd(:, 2), tpcrd(:, 3), ...
    opts.symsize, opts.color, opts.symbol, opts.filled{:});

% labels
th = nan(np, 1);
if opts.enum > 0
    for pc = 1:np
        th(pc) = text(tpcrd(pc, 1), tpcrd(pc, 2), tpcrd(pc, 3), ...
            sprintf('%d', p(pc, end)), 'Parent', opts.axes);
        if ~isempty(opts.labcolor)
            set(th(pc), 'Color', opts.labcolor(pc, :));
        end
    end
elseif ~isempty(opts.label)
    for pc = 1:np
        th(pc) = text(tpcrd(pc, 1), tpcrd(pc, 2), tpcrd(pc, 3), ...
            opts.label{pc}, 'Parent', opts.axes);
        if ~isempty(opts.labcolor)
            set(th(pc), 'Color', opts.labcolor(pc, :));
        end
    end
end
if ~isnan(th(1))
    if ~isempty(opts.labfont)
        set(th, 'FontName', opts.labfont);
    end
    if ~isempty(opts.labfontsz)
        set(th, 'FontSize', opts.labfontsz);
    end
end
opts.th = th;

% set information to enable updates
set(ph, 'UserData', opts);

% force redraw
axn = get(opts.axes, 'Tag');
if numel(axn) == 22
    try
        neuroelf_gui('setsurfpos', axn(4:11));
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

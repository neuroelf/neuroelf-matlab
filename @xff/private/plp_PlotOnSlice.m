function [ph, th] = plp_PlotOnSlice(xo, opts)
% PLP::PlotOnSlice  - plot points on a slice
%
% FORMAT:       [ph, tx] = plp.PlotOnSlice(opts)
%
% Input fields:
%
%       opts        options
%        .axes      axes object, if not given, create new from GUI
%        .bcolor    optional 1x3 RGB color for border of dots (if filled)
%        .bgcolor   background RGB color
%        .color     either column name/number or 1x3 RGB color
%        .cond      conditional statement, e.g.
%                   'Study >= 1 & Study <= 3 & Type == 2'
%                   whereas names are replaced by their column data
%        .conv      convention (for L/R), either or {'neuro'} or 'radio'
%        .dist      distance (fills required .x/y/zrange if not given)
%        .enum      enumerate points (can be a double as start value)
%        .filled    either logical value or column name/number
%        .frame     coordinate frame (2x3 Tal coords)
%        .imethana  anatomical interpolation method (default: linear)
%        .imethstat stat-map interpolation method (default: linear)
%        .joinstats flag whether to join or overlay different stats (false)
%        .label     either column name/number or fixed char
%        .labcolor  either column name/number or 1x3 RGB color
%        .labfont   label font name (empty for default)
%        .labfontsz label font size (empty for default)
%        .labshift  either mm value (radius from center) or X/Y screen vals
%        .oversmp   oversampling (default: 8)
%        .slice     slice coordinate (mandatory, 1x1 position or 1x3 coord)
%        .stalpha   alpha-level of stats map (unless VMP)
%        .statsmaps index of maps to display
%        .statsvar  stats object (e.g. VMP, GLM, ...)
%        .stthresh  1x2 stats threshold (unless VMP)
%        .symbol    either column name/number or 1x1 char symbol
%        .symsize   either column name/number or negative for size
%        .tpvol     only needed if .vmr object is a 4-D (VTC) object
%        .view      either of (mandatory field)
%                   'ax' - axial (transversal)
%                   'co' - coronal
%                   'sa' - sagittal
%        .vmr       alternative VMR (or img/hdr/nii/HEAD) filename
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
% Build:    16021413
% Date:     Feb-14 2016, 1:28 PM EST
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
   ~isfield(opts, 'slice') || ~ischar(opts.view) || ~isa(opts.slice, 'double') || ...
   ~any([1, 3] == numel(opts.slice)) || ...
    any(isinf(opts.slice) | isnan(opts.slice) | abs(opts.slice) > 128) || ...
   ~any(strcmpi(opts.view(:)', {'a', 'ax', 'axial', 'c', 'co', 'cor', ...
        'coronal', 's', 'sa', 'sag', 'sagittal', 'tra', 'transversal'}))
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
rtv = bc.RunTimeVars;
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
if ~isfield(opts, 'cond') || ~ischar(opts.cond) || isempty(opts.cond)
    opts.cond = '';
else
    opts.cond = ['(' lower(opts.cond(:)') ')'];
end
if ~isfield(opts, 'conv') || ~ischar(opts.conv) || isempty(opts.conv) || ...
   ~any('nr' == lower(opts.conv(1)))
    opts.conv = 'n';
else
    opts.conv = lower(opts.conv(1));
end
if ~isfield(opts, 'dist') || ~isa(opts.dist, 'double') || numel(opts.dist) ~= 1 || ...
    isinf(opts.dist) || isnan(opts.dist) || opts.dist < 0
    opts.dist = [];
end
if ~isfield(opts, 'enum') || numel(opts.enum) ~= 1 || ...
   (~islogical(opts.enum) && ~isa(opts.enum, 'double'))
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
if ~isfield(opts, 'filled') || (~isa(opts.filled, 'double') && ...
    ~ischar(opts.filled) && ~islogical(opts.filled))
    opts.filled = true;
end
if ischar(opts.filled)
    opts.filled = findfirst(strcmpi(bc.ColumnNames, opts.filled(:)'));
    if isempty(opts.filled) || opts.filled > size(bc.Points, 2)
        opts.filled = false;
    end
elseif isa(opts.filled, 'double') && (numel(opts.filled) ~= 1 || ...
    isinf(opts.filled) || isnan(opts.filled) || opts.filled < 1 || ...
    ~any((1:size(bc.Points, 2)) == opts.filled))
    opts.filled = false;
elseif islogical(opts.filled) && numel(opts.filled) ~= 1
    opts.filled = any(opts.filled(:));
end
if ~isfield(opts, 'frame') || ~isa(opts.frame, 'double') || ...
   ~isequal(size(opts.frame), [2, 3]) || ...
    any(isinf(opts.frame(:)) | isnan(opts.frame(:))) || ...
    opts.frame(2, :) <= opts.frame(1, :)
    opts.frame = [128, 128, 128; -127.99, -127.99, -127.99];
else
    opts.frame = opts.frame([2, 1], :);
end
if ~isfield(opts, 'imethana') || ~ischar(opts.imethana) || ...
   ~any(strcmpi(opts.imethana(:)', {'nearest', 'linear', 'cubic'}))
    opts.imethana = 'linear';
else
    opts.imethana = lower(opts.imethana(:)');
end
if ~isfield(opts, 'imethstat') || ~ischar(opts.imethstat) || ...
   ~any(strcmpi(opts.imethstat(:)', {'nearest', 'linear', 'cubic'}))
    opts.imethstat = 'linear';
else
    opts.imethstat = lower(opts.imethstat(:)');
end
if ~isfield(opts, 'joinstats') || ~islogical(opts.joinstats) || numel(opts.joinstats) ~= 1
    opts.joinstats = false;
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
elseif isa(opts.label, 'double') && (numel(opts.label) ~= 1 || ...
    isinf(opts.label) || isnan(opts.label) || ~any((1:size(bc.Points, 2)) == opts.label))
    opts.label = [];
elseif isa(opts.label, 'double')
    try
        opts.labeltxt = rtv.ColumnIsText.(bc.ColumnNames{labeli});
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
elseif isa(opts.labcolor, 'double') && numel(opts.labcolor) == 1 && ...
   ~any(1:size(bc.Points, 2) == opts.labcolor)
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
if ~isfield(opts, 'labshift') || ~isa(opts.labshift, 'double') || ...
   ~any([1, 2] == numel(opts.labshift)) || ...
    any(isinf(opts.labshift) | isnan(opts.labshift))
    opts.labshift = [4, -4];
elseif numel(opts.labshift) == 1
    opts.labshift = abs(opts.labshift);
end
if ~isfield(opts, 'oversmp') || ~isa(opts.oversmp, 'double') || numel(opts.oversmp) ~= 1 || ...
    isinf(opts.oversmp) || isnan(opts.oversmp) || opts.oversmp < 1 || opts.oversmp > 32
    opts.oversmp = 8;
else
    opts.oversmp = floor(opts.oversmp);
end
if numel(opts.slice) == 1
    opts.slice = opts.slice([1, 1, 1]);
end
if ~isfield(opts, 'stalpha') || ~isa(opts.stalpha, 'double') || numel(opts.stalpha) ~= 1  || ...
    isinf(opts.stalpha) || isnan(opts.stalpha) || opts.stalpha > 1
    opts.stalpha = 1;
end
if ~isfield(opts, 'statsmaps') || ~isa(opts.statsmaps, 'double') || isempty(opts.statsmaps) || ...
    any(isinf(opts.statsmaps(:)) | isnan(opts.statsmaps(:)) | ...
    opts.statsmaps(:) < 1 | opts.statsmaps(:) ~= fix(opts.statsmaps(:)))
    opts.statsmaps = [];
else
    opts.statsmaps = opts.statsmaps(:);
end
if ~isfield(opts, 'statsvar') || numel(opts.statsvar) ~= 1 || ...
   ~xffisobject(opts.statsvar, true)
    opts.statsvar = [];
end
if ~isfield(opts, 'stthresh') || ~isa(opts.stthresh, 'double') || numel(opts.stthresh) ~= 2 || ...
    any(isinf(opts.stthresh) | isnan(opts.stthresh) | opts.stthresh <= 0) || ...
    opts.stthresh(2) <= opts.stthresh(1)
    opts.stthresh = [sqrt(eps), 1];
else
    opts.stthresh = opts.stthresh(:)';
end
if ~isfield(opts, 'symbol') || (~isa(opts.symbol, 'double') && ~ischar(opts.symbol))
    opts.symbol = 'o';
end
if ischar(opts.symbol) && (numel(opts.symbol) ~= 1 || ~any('.ox+*sdv^<>ph' == opts.symbol'))
    opts.symbol = findfirst(strcmpi(bc.ColumnNames, opts.symbol(:)'));
    if isempty(opts.symbol) || opts.symbol > size(bc.Points, 2)
        opts.symbol = '.';
    end
elseif isa(opts.symbol, 'double') && (numel(opts.symbol) ~= 1 || ...
    isinf(opts.symbol) || isnan(opts.symbol) || ~any((1:size(bc.Points, 2)) ~= opts.symbol))
    opts.symbol = '.';
end
if ~isfield(opts, 'symsize') || (~isa(opts.symsize, 'double') && ~ischar(opts.symsize))
    opts.symsize = [];
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
if ~isfield(opts, 'tpvol') || ~isa(opts.tpvol, 'double') || numel(opts.tpvol) ~= 1 || ...
    isinf(opts.tpvol) || isnan(opts.tpvol) || opts.tpvol < 1
    opts.tpvol = 1;
else
    opts.tpvol = floor(opts.tpvol);
end
opts.view = lower(opts.view(:)');
if ~isempty(opts.dist)
    switch (opts.view(1))
        case {'a', 't'}
            if ~isfield(opts, 'zrange')
                opts.zrange = opts.slice(3) + [-opts.dist, opts.dist];
            end
        case 'c'
            if ~isfield(opts, 'yrange')
                opts.yrange = opts.slice(2) + [-opts.dist, opts.dist];
            end
        case 's'
            if ~isfield(opts, 'xrange')
                opts.xrange = opts.slice(1) + [-opts.dist, opts.dist];
            end
    end
end
if ~isfield(opts, 'vmr') || ~ischar(opts.vmr) || isempty(opts.vmr) || ...
    exist(opts.vmr(:)', 'file') ~= 2
    opts.vmr = '';
else
    [ia, opts.vmr] = isabsolute(opts.vmr(:)');
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
    cond = [cond sprintf(' & p(:,%d)>=%f & p(:,%d)<=%f', xc, opts.xrange(1), xc, opts.xrange(2))];
end
yc = findfirst(strcmpi(bc.ColumnNames, 'y'));
if isempty(yc)
    yc = 2;
end
if ~isempty(opts.yrange)
    cond = [cond sprintf(' & p(:,%d)>=%f & p(:,%d)<=%f', yc, opts.yrange(1), yc, opts.yrange(2))];
end
zc = findfirst(strcmpi(bc.ColumnNames, 'z'));
if isempty(zc)
    zc = 3;
end
if ~isempty(opts.zrange)
    cond = [cond sprintf(' & p(:,%d)>=%f & p(:,%d)<=%f', zc, opts.zrange(1), zc, opts.zrange(2))];
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

% try to get colors
if numel(opts.color) == 1
    try
        opts.color = bc.Colors(p(:, opts.color), :);
    catch xfferror
        neuroelf_lasterr(xfferror);
        if ~isempty(bc.Colors)
            opts.color = bc.Colors(ones(np, 1), :);
        else
            opts.color = 255 .* [ones(np, 1), zeros(np, 2)];
        end
    end
else
    opts.color = opts.color(ones(np, 1), :);
end
opts.color = (1 / 255) .* opts.color;

% get filled status
if ~islogical(opts.filled)
    try
        opts.filled = (p(:, opts.filled) > 0);
    catch xfferror
        neuroelf_lasterr(xfferror);
        opts.filled = false(np, 1);
    end
end

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

% symbols and size
if isa(opts.symbol, 'double')
    try
        opts.symbol = bc.Symbols(p(:, opts.symbol));
    catch xfferror
        neuroelf_lasterr(xfferror);
        opts.symbol = repmat({'.'}, np, 1);
    end
else
    opts.symbol = repmat({opts.symbol}, np, 1);
end
if ~isempty(opts.symsize) && opts.symsize > 0
    try
        opts.symsize = p(:, opts.symsize);
        opts.symsize(opts.symsize < 0) = bc.SymbolSizes(-opts.symsize(opts.symsize < 0));
    catch xfferror
        neuroelf_lasterr(xfferror);
        opts.symsize = [];
    end
elseif ~isempty(opts.symsize) && opts.symsize < 0
    opts.symsize = -opts.symsize(ones(np, 1));
else
    opts.symsize = [];
end

% make sure any requested VMR is loaded
if ~isempty(opts.vmr)
    tlv = neuroelf_gui('varlist');
    tlf = fieldnames(tlv);
    for tlc = 1:numel(tlf)
        if strcmpi(opts.vmr, tlv.(tlf{tlc}).F)
            try
                neuroelf_gui('openfile', tlv.(tlf{tlc}));
                opts.vmr = '';
                break;
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
        end
    end
end
if ~isempty(opts.vmr)
    try
        neuroelf_gui('openfile', opts.vmr);
        pause(0.25);
        drawnow;
    catch xfferror
        rethrow(xfferror);
    end
end
try
    opts.vmr = neuroelf_gui('slicevar');
catch xfferror
    rethrow(xfferror);
end

% check axes if given
if ~isempty(opts.axes)
    imgh = findobj('Parent', opts.axes, 'Type', 'image');
    if numel(imgh) ~= 1
        opts.axes = [];
    else
        imgs = size(get(imgh, 'CData'));
    end
end

% depending on direction
switch (opts.view(1))
    case {'a', 't'}
        drc = 3;
        angler = 90;
        anglez = 90;
    case 'c'
        drc = 2;
        angler = 180;
        anglez = 90;
    case 's'
        drc = 1;
        angler = 180;
        anglez = 180;
end
opts.frame(:, drc) = [opts.slice(drc); opts.slice(drc) - 0.01];
angler = (pi / 180) * angler;
anglez = (pi / 180) * anglez;

% no valid axes available
if isempty(opts.axes)
    try
        [idata, opts.axes] = neuroelf_gui('vismontage_create_ex', struct( ...
            'atrans', ~isempty(opts.bgcolor), 'atranscol', opts.bgcolor, 'blx', [1, 1], ...
            'brds', 0, 'drc', drc, 'drs', 256, 'filename', '', 'flp', false, ...
            'flx', (opts.conv == 'n' && opts.view(1) ~= 's'), 'fontcolor', [0, 0, 0], ...
            'fontname', 'Arial', 'fontsize', 12, 'frame', opts.frame, ...
            'imeth', opts.imethstat, 'imetha', opts.imethana, 'join', opts.joinstats, ...
            'ppv', opts.oversmp, 'showinfig', true, 'slcoord', false, 'slvar', opts.vmr, ...
            'stalp', opts.stalpha, 'stthr', opts.stthresh, 'stvar', opts.statsvar, ...
            'stvix', opts.statsmaps, 'sws', false, 'tpvol', opts.tpvol));
        imgs = size(idata);
    catch xfferror
        rethrow(xfferror);
    end
else
    opts.oversmp = imgs(1) / 256;
end

% perform transformation
atrf = spmtrf([0, 0, 0], [0, angler, 0]) * spmtrf([0, 0, 0], [0, 0, anglez]) * ...
    spmtrf([0, 0, 0], [0, 0, 0], opts.oversmp(ones(1, 3)));
tpcrd = opts.pcrd * atrf';
tpcrd(:, 1) = 160;
if opts.conv == 'n' && opts.view(1) ~= 's'
    tpcrd(:, 2) = -tpcrd(:, 2);
end

% force hold
hold(opts.axes, 'on');

% plot points
ph = scatter(opts.axes, tpcrd(:, 2) + 0.5 * imgs(1), tpcrd(:, 3) + 0.5 * imgs(2));
if ~isempty(opts.color)
    set(ph, 'CData', opts.color);
end
if numel(opts.filled) == 1 && opts.filled
    set(ph, 'MarkerFaceColor', 'flat');
end
if ~isempty(opts.symsize)
    set(ph, 'SizeData', opts.symsize(:));
end

% iterate over points
phc = get(ph, 'Children');
phc = phc(end:-1:1);
for cc = 1:numel(phc)
    if ~isempty(opts.symbol)
        set(phc(cc), 'Marker', opts.symbol{cc});
    end
    if ~isempty(opts.color) && numel(opts.filled) == numel(phc) && opts.filled(cc)
        set(phc(cc), 'MarkerFaceColor', opts.color(cc, :));
    end
    if ~isempty(opts.bcolor) && ((numel(opts.filled) == 1 && ...
         opts.filled) || (numel(opts.filled) == numel(phc) && opts.filled(phc)))
        set(phc(cc), 'MarkerEdgeColor', opts.bcolor);
    end
end

% shift for labels
if numel(opts.labshift) == 1
    tpdir = tpcrd(:, 2:3);
    tpdis = opts.labshift ./ sqrt(sum(tpdir .* tpdir, 2));
    tpcrd(:, 2) = tpcrd(:, 2) + (tpdir(:, 1) .* tpdis) + 0.5 * imgs(1);
    tpcrd(:, 3) = tpcrd(:, 3) + (tpdir(:, 2) .* tpdis) + 0.5 * imgs(2);
elseif numel(opts.labshift) == 2
    tpcrd(:, 2) = tpcrd(:, 2) + (0.5 * imgs(1) + opts.labshift(1));
    tpcrd(:, 3) = tpcrd(:, 3) + (0.5 * imgs(2) + opts.labshift(2));
else
    tpcrd(:, 2) = tpcrd(:, 2) + 0.5 * imgs(1);
    tpcrd(:, 3) = tpcrd(:, 3) + 0.5 * imgs(2);
end

% labels
th = nan(np, 1);
if opts.enum > 0
    for pc = 1:np
        th(pc) = text(tpcrd(pc, 2), tpcrd(pc, 3), sprintf('%d', p(pc, end)), 'Parent', opts.axes);
        if ~isempty(opts.labcolor)
            set(th(pc), 'Color', opts.labcolor(pc, :));
        end
    end
elseif ~isempty(opts.label)
    for pc = 1:np
        th(pc) = text(tpcrd(pc, 2), tpcrd(pc, 3), opts.label{pc}, 'Parent', opts.axes);
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

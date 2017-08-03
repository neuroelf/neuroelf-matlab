function [l, p] = tcplot(varargin)
% tcplot  - plot timecourse data
%
% FORMAT        [l, p] = tcplot([a ,] x, y [, s [, sn [, opts]]])
%
% Input fields:
%
%       a           axes handle
%       x           x ordinate values
%       y           y ordinate values
%       s           y std error estimate
%       sn          if given s is considered the positive, sn the negative
%       opts        optional settings
%        .color     [0 .. 1] -scaled RGB color for line
%        .lwidth    line width
%        .scolor    color for std patch
%        .spalpha   std patch alpha value (default: 1)
%        .spline    std patch line (default: true)
%
% Output fields:
%
%       l           optional handle of line
%       p           optional handle of std error patch

% Version:  v1.0
% Build:    14110514
% Date:     Nov-05 2014, 2:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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

% argument check
if nargin < 2 || ...
   (~isnumeric(varargin{1}) && ...
    ~ishandle(varargin{1})) || ...
    isempty(varargin{1})
    error( ...
        'neuroelf:BadArgument', ...
        'Two few arguments or bad first argument.' ...
    );
end

% first arg a valid axis object
if numel(varargin{1}) == 1 && ...
   (isa(varargin{1}, 'double') || ...
    isa(varargin{1}, 'matlab.graphics.axis.Axes')) && ...
    ishandle(varargin{1}) && ...
    strcmpi(get(varargin{1}, 'Type'), 'axes')
    ax = varargin{1};
    args = varargin(2:nargin);
else
    ax = [];
    args = varargin(1:nargin);
end
nargs = numel(args);

% get remaining arguments
if nargs < 2 || ...
   ~isnumeric(args{1}) || ...
   ~isnumeric(args{2}) || ...
    numel(args{1}) ~= numel(args{2})
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
x = double(args{1}(:));
y = double(args{2}(:));
if nargs > 2 && ...
    isnumeric(args{3}) && ...
    numel(args{3}) == numel(x)
    s = double(args{3}(:));
    if nargs > 3 && ...
        isnumeric(args{4}) && ...
        numel(args{4}) == numel(x)
        sn = double(args{4}(:));
    else
        sn = s;
    end
else
    s = [];
end
if nargs > 2 && ...
    isstruct(args{end}) && ...
    numel(args{end}) == 1
    opts = args{end};
else
    opts = struct;
end
if ~isfield(opts, 'color') || ...
   ~isa(opts.color, 'double') || ...
    numel(opts.color) ~= 3 || ...
    any(isinf(opts.color) | isnan(opts.color) | opts.color < 0 | opts.color > 1)
    opts.color = [];
else
    opts.color = opts.color(:)';
end
if ~isfield(opts, 'lwidth') || ...
   ~isa(opts.lwidth, 'double') || ...
    numel(opts.lwidth) ~= 1 || ...
    isinf(opts.lwidth) || ...
    isnan(opts.lwidth) || ...
    opts.lwidth <= 0 || ...
    opts.lwidth >= 5
    opts.lwidth = [];
end
if ~isfield(opts, 'scolor') || ...
   ~isa(opts.scolor, 'double') || ...
    numel(opts.scolor) ~= 3 || ...
    any(isinf(opts.scolor) | isnan(opts.scolor) | opts.scolor < 0 | opts.scolor > 1)
    opts.scolor = [];
else
    opts.scolor = opts.scolor(:)';
end
if ~isfield(opts, 'spalpha') || ...
   ~isa(opts.spalpha, 'double') || ...
    numel(opts.spalpha) ~= 1 || ...
    isinf(opts.spalpha) || ...
    isnan(opts.spalpha) || ...
    opts.spalpha < 0 || ...
    opts.spalpha > 1
    opts.spalpha = [];
end
if ~isfield(opts, 'spline') || ...
   ~islogical(opts.spline) || ...
    numel(opts.spline) ~= 1
    opts.spline = true;
end

% get axes to plot into
if isempty(ax)
    ax = gca;
else
    af = get(ax, 'Parent');
    while ~strcmpi(get(af, 'Type'), 'figure')
        af = get(af, 'Parent');
        if strcmpi(get(af, 'Type'), 'root')
            error( ...
                'neuroelf:UIError', ...
                'Invalid object chain detected.' ...
            );
        end
    end
    set(0, 'CurrentFigure', af);
    set(af, 'CurrentAxes', ax);
end
hold(ax, 'on');

% plot the patch first
if ~isempty(s)
    if ~isempty(opts.scolor)
        pcol = opts.scolor;
    elseif ~isempty(opts.color)
        pcol = 0.5 + 0.5 * opts.color;
    else
        pcol = [0.8, 0.8, 0.8];
    end
    p = patch([x; x(end:-1:1); x(1)], ...
        [y + s; y(end:-1:1) - sn(end:-1:1); y(1) + s(1)], pcol, 'Parent', ax);
    if ~isempty(opts.spalpha)
        set(p, 'FaceAlpha', opts.spalpha);
    end
    if ~opts.spline
        set(p, 'EdgeAlpha', 0);
    end
else
    p = [];
end

% plot the line last
l = plot(ax, x, y);
if ~isempty(opts.color)
    set(l, 'Color', opts.color);
end
if ~isempty(opts.lwidth)
    set(l, 'LineWidth', opts.lwidth);
end

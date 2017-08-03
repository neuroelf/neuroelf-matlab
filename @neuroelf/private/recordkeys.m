function rks = recordkeys(opts)
% recordkeys  - record keyboard inputs
%
% FORMAT:       rks = recordkeys([opts])
%
% Input fields:
%
%       opts        optional settings
%        .prompt    prompt string (text on figure, default: 'recording...')
%        .releases  flag, also record release events (default: false)
%        .timeout   number of seconds to timeout (default: Inf)
%        .title     figure title (default: 'NeuroElf - keyboard recorder')
%        .xchar     exit character (default: 'x')
%        .xcmod     exit character modifiers (default: {'control'})
%
% Output fields:
%
%       rks         Ex2 cell array with first column time (since start),
%                   and second column the pressed or released key

% Version:  v0.9c
% Build:    11061317
% Date:     Jun-13 2011, 4:30 PM EST
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

% global variable for recording
global ne_ui;

% check global variable
if ~isfield(ne_ui, 'rks') || ...
   ~isstruct(ne_ui.rks) || ...
    numel(ne_ui.rks) ~= 1
    ne_ui.rks = struct( ...
        'fig',   [], ...
        'nume',  0, ...
        'opts',  struct, ...
        'rks',   {cell(100, 2)}, ...
        'start', -Inf);
end

% argument check
if nargin < 1 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'prompt') || ...
   ~ischar(opts.prompt)
    opts.prompt = 'recording...';
else
    opts.prompt = opts.prompt(:)';
end
if ~isfield(opts, 'releases') || ...
   ~islogical(opts.releases) || ...
    numel(opts.releases) ~= 1
    opts.releases = false;
end
if ~isfield(opts, 'timeout') || ...
   ~isa(opts.timeout, 'double') || ...
    numel(opts.timeout) ~= 1 || ...
    isnan(opts.timeout) || ...
    opts.timeout <= 0
    opts.timeout = Inf;
end
if ~isfield(opts, 'title') || ...
   ~ischar(opts.title) || ...
    isempty(opts.title)
    opts.title = 'NeuroElf - keyboard recorder';
else
    opts.title = opts.title(:)';
end
if ~isfield(opts, 'xchar') || ...
   ~ischar(opts.xchar) || ...
    numel(opts.xchar) ~= 1 || ...
   ~any(lower(opts.xchar) == '0123456789abcdefghijklmnopqrstuvwxyz')
    opts.xchar = 'x';
    opts.xcmod = {'control'};
else
    opts.xchar = lower(opts.xchar);
end
if ~isfield(opts, 'xcmod') || ...
   ~iscell(opts.xcmod)
    opts.xcmod = {};
end
for xc = numel(opts.xcmod):-1:1
    if ~any(strcmpi(opts.xcmod{xc}(:)', {'control', 'alt', 'command', 'shift'}))
        opts.xcmod(xc) = [];
    else
        opts.xcmod{xc} = lower(opts.xcmod{xc}(:)');
    end
end
if ~isempty(opts.xcmod)
    opts.xcmod = unique(lower(opts.xcmod(:)));
end

% create figure
f = figure;
set(f, 'NumberTitle', 'off', 'Name', opts.title, 'Units', 'pixels', ...
    'Color', [0.8, 0.8, 0.8]);
fp = get(f, 'Position');

% add central prompt
u = uicontrol('Parent', f, 'Style', 'text', 'BackgroundColor', [0.8, 0.8, 0.8], ...
    'Position', [0, 0.5 * fp(4), fp(3), 20], ...
    'String', opts.prompt, 'HorizontalAlignment', 'center');

% set callbacks
set(f, 'WindowKeyPressFcn',   @rks_press);
if opts.releases
    set(f, 'WindowKeyReleaseFcn', @rks_release);
end

% recompute timeout
ne_ui.rks.fig = f;
ne_ui.rks.opts = opts;
ne_ui.rks.start = now;
timeout = ne_ui.rks.start + opts.timeout / 86400;

% wait until time out is reached or figure closed
while ishandle(f) && ...
    timeout >= now

    % update
    if ishandle(f)
        if isinf(timeout)
            set(u, 'String', sprintf('%s (%d event%s captured)', ...
                opts.prompt, ne_ui.rks.nume, plurals(ne_ui.rks.nume)));
        else
            set(u, 'String', sprintf('%s (%d event%s captured, %d seconds left)', ...
                opts.prompt, ne_ui.rks.nume, plurals(ne_ui.rks.nume), ...
                floor(86400 * (timeout - now))));
        end
    end

    % pause
    pause(0.1);
end

% delete figure?
if ishandle(f)
    delete(f);
end

% return recorded events
rks = ne_ui.rks.rks(1:ne_ui.rks.nume, :);

% and de-initialize
ne_ui.rks = [];



% UI functions



function rks_press(src, ke, varargin)
global ne_ui;
o = ne_ui.rks.opts;

% time
n = now;

% what keys
k = ke.Key;
m = ke.Modifier;
m = unique(lower(m(:)));

% quit program
if numel(m) == numel(o.xcmod) && ...
    all(strcmp(m, o.xcmod)) && ...
    strcmpi(k, o.xchar)

    % delete figure
    delete(ne_ui.rks.fig);

    % return
    return;
end

% record current key press
if isempty(m) || ...
   ~any(strcmpi(k, m))
    m{end+1} = k;
end
ne_ui.rks.nume = ne_ui.rks.nume + 1;
if size(ne_ui.rks.rks, 1) < ne_ui.rks.nume
    ne_ui.rks.rks{ne_ui.rks.nume + 99, 2} = [];
end
ne_ui.rks.rks(ne_ui.rks.nume, :) = {m(:)', 86400 * (n - ne_ui.rks.start)};


function rks_release(src, ke, varargin)
global ne_ui;

% time
n = now;

% what keys (remain)
m = ke.Modifier;
m = unique(lower(m(:)));

% record current key press
ne_ui.rks.nume = ne_ui.rks.nume + 1;
if size(ne_ui.rks.rks, 1) < ne_ui.rks.nume
    ne_ui.rks.rks{ne_ui.rks.nume + 99, 2} = [];
end
ne_ui.rks.rks(ne_ui.rks.nume, :) = {m(:)', 86400 * (n - ne_ui.rks.start)};

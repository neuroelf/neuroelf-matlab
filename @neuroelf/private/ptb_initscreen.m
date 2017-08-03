function ptbvars = ptb_initscreen(opts)
% ptb_initscreen  - initialize PsychToolbox screen
%
% FORMAT:       ptbvars = ptb_initscreen([opts])
%
% Input fields:
%
%       opts        1x1 optional struct with settings
%        .bgcolor   background color RGB (default: [0, 0, 0]);
%        .debugsize 1x4 optional size/position of "Screen" (default: [])
%        .fontname  string font name (default: don't set)
%        .fontsize  1x1 double font size (default: don't set)
%        .fontstyle 1x1 double font style (default: don't set)
%        .hcursor   if true (default) call HideCursor
%        .screen    screen number to use (default: last)
%        .skipsync  value for 'SkipSyncTests' Preference (default: 2)
%        .textures  struct with fields where values are filenames to images
%        .unifykeys if true (default) call KbName('UnifyKeyNames')
%        .verbosity value for 'Verbosity' Preference (default: 0)
%        .visdebug  value for 'VisualDebugLevel' Preference (default: 0)
%
% Output fields:
%
%       ptbvars     1x1 struct with fields
%        .screen        screen (wPtr) pointer
%        .screenrect    screen rectangle
%        .texturenames  texture names
%        .textures      textures (struct)
%        .texturesizes  texture sizes (struct with 1x2 values)
%
% Note: for fontstyles, check Screen('TextStyle?')

% Version:  v0.9c
% Build:    12051712
% Date:     May-01 2012, 4:43 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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

% available screens
try
    as = Screen('Screens');
catch ne_eo;
    rethrow(ne_eo);
end

% arguments
if nargin < 1 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bgcolor') || ...
   ~isa(opts.bgcolor, 'double') || ...
    numel(opts.bgcolor) ~= 3 || ...
    any(isinf(opts.bgcolor) | isnan(opts.bgcolor) | opts.bgcolor < 0 | opts.bgcolor > 255)
    opts.bgcolor = [0, 0, 0];
else
    opts.bgcolor = round(opts.bgcolor(:)');
end
if ~isfield(opts, 'debugsize') || ...
   ~isa(opts.debugsize, 'double') || ...
   ~isequal(size(opts.debugsize), [1, 4]) || ...
    any(isinf(opts.debugsize) | isnan(opts.debugsize) | opts.debugsize < 0 | opts.debugsize ~= fix(opts.debugsize)) || ...
    any(opts.debugsize(3:4) - opts.debugsize(1:2) < 1)
    opts.debugsize = [];
end
if ~isfield(opts, 'fontname') || ...
   ~ischar(opts.fontname) || ...
    isempty(opts.fontname)
    opts.fontname = '';
end
if ~isfield(opts, 'fontsize') || ...
   ~isa(opts.fontsize, 'double') || ...
    numel(opts.fontsize) ~= 1 || ...
    isinf(opts.fontsize) || ...
    isnan(opts.fontsize) || ...
    opts.fontsize < 1 || ...
    opts.fontsize > 400 || ...
    opts.fontsize ~= fix(opts.fontsize)
    opts.fontsize = [];
end
if ~isfield(opts, 'fontstyle') || ...
   ~isa(opts.fontstyle, 'double') || ...
    numel(opts.fontstyle) ~= 1 || ...
    isinf(opts.fontstyle) || ...
    isnan(opts.fontstyle) || ...
    opts.fontstyle ~= fix(opts.fontstyle) || ...
    opts.fontstyle < 0 || ...
    opts.fontstyle > 255
    opts.fontstyle = [];
end
if ~isfield(opts, 'hcursor') || ...
   ~islogical(opts.hcursor) || ...
    numel(opts.hcursor) ~= 1
    opts.hcursor = true;
end
if ~isfield(opts, 'screen') || ...
   ~isa(opts.screen, 'double') || ...
    numel(opts.screen) ~= 1 || ...
    isinf(opts.screen) || ...
    isnan(opts.screen) || ...
   ~any(as == opts.screen)
    opts.screen = max(as(:));
end
if ~isfield(opts, 'skipsync') || ...
   ~isa(opts.skipsync, 'double') || ...
    numel(opts.skipsync) ~= 1 || ...
    isinf(opts.skipsync) || ...
    isnan(opts.skipsync) || ...
   ~any(0:2 == opts.skipsync)
    opts.skipsync = 2;
end
if isfield(opts, 'textures') && ...
    iscell(opts.textures) && ...
   ~isempty(opts.textures)
    t = struct;
    d = ceil(log10(numel(opts.textures) + 1));
    d = sprintf('t%%0%dd', d);
    for tc = 1:numel(opts.textures)
        if ischar(opts.textures{tc}) && ...
           ~isempty(opts.textures{tc}) && ...
            exist(opts.textures{tc}, 'file') == 2
            try
                t.(sprintf(d, tc)) = opts.textures{tc};
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
    opts.textures = t;
end
if ~isfield(opts, 'textures') || ...
   ~isstruct(opts.textures) || ...
    numel(opts.textures) ~= 1
    opts.textures = struct;
end
if ~isfield(opts, 'unifykeys') || ...
   ~islogical(opts.unifykeys) || ...
    numel(opts.unifykeys) ~= 1
    opts.unifykeys = true;
end
if ~isfield(opts, 'verbosity') || ...
   ~isa(opts.verbosity, 'double') || ...
    numel(opts.verbosity) ~= 1 || ...
    isinf(opts.verbosity) || ...
    isnan(opts.verbosity) || ...
   ~any(0:2 == opts.verbosity)
    opts.verbosity = 0;
end
if ~isfield(opts, 'visdebug') || ...
   ~isa(opts.visdebug, 'double') || ...
    numel(opts.visdebug) ~= 1 || ...
    isinf(opts.visdebug) || ...
    isnan(opts.visdebug) || ...
   ~any(0:2 == opts.visdebug)
    opts.visdebug = 0;
end

% make calls
try
    Screen('Preference', 'SkipSyncTests',    opts.skipsync);
    Screen('Preference', 'Verbosity',        opts.verbosity);
    Screen('Preference', 'VisualDebugLevel', opts.visdebug);
catch ne_eo;
    rethrow(ne_eo);
end

% hide cursor
if opts.hcursor
    HideCursor;
end

% also unify keyboard names
if opts.unifykeys
    KbName('UnifyKeyNames');
end

% start by closing prior sessions
Screen('CloseAll');

% open screen
if isempty(opts.debugsize)
    [s, r] = Screen('OpenWindow', opts.screen, opts.bgcolor);
else
    [s, r] = Screen('OpenWindow', opts.screen, opts.bgcolor, opts.debugsize);
end

% font settings
if ~isempty(opts.fontname)
    try
        Screen('TextFont', s, opts.fontname);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end
if ~isempty(opts.fontsize)
    try
        Screen('TextSize', s, opts.fontsize);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end
if ~isempty(opts.fontstyle)
    try
        Screen('TextStyle', s, opts.fontstyle);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% return textures
tn = opts.textures;
t = tn;
ts = t;

% load textures
tf = fieldnames(t);
for tc = 1:numel(tf)
    tl = [];
    if ischar(t.(tf{tc})) && ...
       ~isempty(t.(tf{tc})) && ...
        exist(t.(tf{tc}), 'file') == 2
        try
            ti = imread(t.(tf{tc}));
            tl = Screen('MakeTexture', s, ti);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
    if ~isempty(tl)
        t.(tf{tc}) = tl;
        ts.(tf{tc}) = [size(ti, 2), size(ti, 1)];
    else
        tn = rmfield(tn, tf{tc});
        t = rmfield(t, tf{tc});
        ts = rmfield(ts, tf{tc});
    end
end

% get names resolution array
tn = [struct2cell(tn), fieldnames(tn)];

% create output
ptbvars = struct( ...
    'screen', s, ...
    'screenrect', r, ...
    'texturenames', {tn}, ...
    'textures', t, ...
    'texturesizes', ts);

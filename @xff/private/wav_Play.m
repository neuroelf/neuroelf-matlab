function xo = wav_Play(xo, opts)
% WAV::Play  - play wave file
%
% FORMAT:       wav.Play([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .freq      frequency
%        .range     1x2 range (if fractional, in seconds)
%        .tfunc     timer function (handle or cell array)
%
% No output fields.
%
% TYPES: WAV

% Version:  v1.1
% Build:    16012915
% Date:     Jan-29 2016, 3:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, 2016, Jochen Weber
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

% using NeuroElf functions
using(neuroelf, {'limitrangec'});

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'wav')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% options
if nargin > 1 && isa(opts, 'double') && numel(opts) == 1
    opts = struct('freq', opts);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'freq') || ~isa(opts.freq, 'double') || numel(opts.freq) ~= 1 || ...
    isinf(opts.freq) || isnan(opts.freq) || opts.freq < 1
    opts.freq = xo.C.FormatChunk.SampleRate;
end
if ~isfield(opts, 'range') || ~isa(opts.range, 'double') || numel(opts.range) < 2 || ...
    any(isnan(opts.range(:)) | opts.range(:) < 0)
    opts.range = [1, size(xo.C.Data, 2)];
elseif numel(opts.range) > 2
    opts.range = limitrangec(opts.range(:), 0, size(xo.C.Data, 2), 1);
end
if any(opts.range < 1 | opts.range ~= fix(opts.range))
    opts.range = fix(limitrangec(1 + xo.C.FormatChunk.SampleRate .* opts.range, ...
        1, size(xo.C.Data, 2), 1));
end

% simple replay?
if isfield(xo.H, 'Player') && numel(xo.H.Player) == 1 && isa(xo.H.Player, 'audioplayer')
    try
        apd = get(xo.H.Player, 'UserData');
        if isstruct(apd) && numel(apd) == 1 && isfield(apd, 'freq') && ...
            isequal(apd.freq, opts.freq) && isfield(apd, 'range') && ...
            isequal(apd.range, opts.range)
            stop(xo.H.Player);
            play(xo.H.Player);
            return;
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

% handle set?
if isfield(xo.H, 'Player') && ~isempty(xo.H.CleanUp)
   
    % delete and re-set
    for hc = numel(xo.H.CleanUp):-1:1
        if isa(xo.H.CleanUp{hc}, 'audioplayer')
            try
                delete(xo.H.CleanUp{hc});
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
            xo.H.CleanUp(hc) = [];
        end
    end
    xo.H.Player = [];
end

% create player
if numel(opts.range) > 2
    data = xo.C.Data(:, :);
    p = audioplayer(data(:, opts.range), opts.freq);
elseif opts.range(2) >= opts.range(1)
    p = audioplayer(xo.C.Data(:, opts.range(1):opts.range(2))', opts.freq);
else
    p = audioplayer(xo.C.Data(:, opts.range(1):-1:opts.range(2))', opts.freq);
end

% set fields with options
set(p, 'UserData', struct('freq', opts.freq, 'range', opts.range));

% add to cleanup
xo.H.Player = p;
xo.H.CleanUp{end+1} = p;

% play
play(xo.H.Player);

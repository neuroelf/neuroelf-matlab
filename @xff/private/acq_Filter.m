function xo = acq_Filter(xo, channels, opts)
% ACQ::Filter  - filter channel data using prefilter
%
% FORMAT:       [acq = ] acq.Filter([channels [, opts]])
%
% Input fields:
%
%       channels    list of channels (default: all)
%       opts        optional settings
%        .dest      destination channels (must match number of inputs)
%        .freq      signal frequency (default: from file)
%        .fwin      filtering window seconds (default: 2)
%        .fwincut   filtering frequency cutoff in seconds (default: 0.5)
%        .fwslide   filtering window slide in seconds (default: fwin/16)
%        .kern      smoothing kernel in seconds (default: 0.05)
%        .post      post-diff smoothing kernel in seconds (default: 0.1)
%        .stdt      std(abs(diff)) change threshold (default: 0.01)
%
% Output fields:
%
%       acq         altered object
%
% Using: prefilter.

% Version:  v1.1
% Build:    16020214
% Date:     Feb-02 2016, 2:18 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'acq') || ...
   ~isa(channels, 'double') || isempty(channels) || numel(channels) ~= max(size(channels)) || ...
    any(isinf(channels) | isnan(channels) | channels <= 0 | channels ~= fix(channels))
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
channels = unique(min(numel(bc.Channel), channels));
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dest') || ~isa(opts.dest, 'double') || ...
    numel(opts.dest) ~= numel(channels) || numel(opts.dest) ~= max(size(opts.dest)) || ...
    any(isinf(opts.dest) | isnan(opts.dest) | opts.dest <= 0 | opts.dest ~= fix(opts.dest))
    opts.dest = channels;
else
    opts.dest = unique(opts.dest(:)');
    odav = false(1, max(opts.dest));
    odav(1:numel(bc.Channel)) = true;
    odav(opts.dest) = true;
    if ~all(odav)
        error('neuroelf:xff:badArgument', ...
            'The destination channels must not skip and create an empty channel!');
    end
    if numel(odav) > numel(bc.Channel)
        bc.Channel(end+1:numel(odav)) = bc.Channel(1);
        bc.NrOfChannels = numel(bc.Channel);
    end
    olc = bc.Channel;
    for dcc = 1:numel(opts.dest)
        dc = opts.dest(dcc);
        bc.Channel(dc).ChannelNumber = dc;
        bc.Channel(dc).Comment = sprintf('Filtered channel %d', channels(dcc));
        bc.Channel(dc).ColorRGB = olc(channels(dcc)).ColorRGB;
        bc.Channel(dc).DisplayChannel = olc(channels(dcc)).DisplayChannel;
        bc.Channel(dc).VoltageOffset = olc(channels(dcc)).VoltageOffset;
        bc.Channel(dc).VoltageFactor = olc(channels(dcc)).VoltageFactor;
        bc.Channel(dc).Unit = olc(channels(dcc)).Unit;
        bc.Channel(dc).AmplitudeScale = olc(channels(dcc)).AmplitudeScale;
        bc.Channel(dc).AmplitudeOffset = olc(channels(dcc)).AmplitudeOffset;
    end
end
if ~isfield(opts, 'freq') || ~isa(opts.freq, 'double') || numel(opts.freq) ~= 1 || ...
    isinf(opts.freq) || isnan(opts.freq) || opts.freq < 0.1 || opts.freq > 1e6
    opts.freq = 1000 / bc.MillisecsPerSample;
end
if ~isfield(opts, 'kern') || ~isa(opts.kern, 'double') || numel(opts.kern) ~= 1 || ...
    isinf(opts.kern) || isnan(opts.kern) || opts.kern < 0
    opts.kern = 0.05;
end
opts.kern = opts.kern * opts.freq;
if ~isfield(opts, 'post') || ~isa(opts.post, 'double') || numel(opts.post) ~= 1 || ...
    isinf(opts.post) || isnan(opts.post) || opts.post < 0
    opts.post = 0.1;
end
opts.post = opts.post * opts.freq;
if ~isfield(opts, 'stdt') || ~isa(opts.stdt, 'double') || numel(opts.stdt) ~= 1 || ...
    isinf(opts.stdt) || isnan(opts.stdt) || opts.stdt <= 0
    opts.stdt = 0.02;
else
    opts.stdt = max(1e-6, opts.stdt - floor(opts.stdt));
end

% iterate over channels
for cc = 1:numel(channels)

    % set with filtering result
    bc.Channel(opts.dest(cc)).Data = ...
        ne_methods.prefilter(acq_ChannelData(xo, channels(cc)), opts);
end

% set back in memory
xo.C = bc;

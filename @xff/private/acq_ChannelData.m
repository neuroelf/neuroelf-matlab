function cdata = acq_ChannelData(xo, channel, idx)
% ACQ::ChannelData  - get channel data from RawData/Channel(C).Data
%
% FORMAT:       cdata = acq.ChannelData(channel [, idx])
%
% Input fields:
%
%       channel     channel number
%       idx         indices (passed on!)
%
% Output fields:
%
%       cdata       channel data

% Version:  v1.1
% Build:    16020214
% Date:     Feb-02 2016, 2:19 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'acq') || ...
   ~isa(channel, 'double') || numel(channel) ~= 1 || isinf(channel) || isnan(channel) || ...
    channel < 1 || channel ~= fix(channel)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if channel > bc.NrOfChannels
    error('neuroelf:xff:badArgument', 'Channel number out of bounds.');
end

% if channel data is available
if ~isempty(bc.Channel(channel).Data)
    cdata = bc.Channel(channel).Data(:);

% otherwise
else
    f = double(bc.Channel(channel).AmplitudeScale);
    o = double(bc.Channel(channel).AmplitudeOffset);
    if isinf(f) || isnan(f) || f == 0
        f = 1;
    end
    if isinf(o) || isnan(o)
        o = 0;
    end
    cdata = o + f .* double(bc.RawData(channel, :)');
    cdn = isnan(cdata);
    if any(cdn)
        cdata(cdn) = sum(cdata(~cdn)) / sum(~cdn);
    end

    % and store
    bc.Channel(channel).Data = cdata;
    xo.C = bc;
end

% indices
if nargin > 2 && (islogical(idx) || isa(idx, 'double'))
    try
        cdata = cdata(idx);
    catch xfferror
        rethrow(xfferror);
    end
end

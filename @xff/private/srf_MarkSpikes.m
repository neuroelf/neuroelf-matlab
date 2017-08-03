function [xo, mi, mn1, mn2, mn3] = srf_MarkSpikes(xo, xo2, cutoff, varargin)
% SRF::MarkSpikes  - find spiked vertices (1-based indices)
%
% FORMAT:       [srf, sps, npmin, vdmin, vdmean] = srf.MarkSpikes([spikesrf, co]);
%
% Input fields:
%
%       spikesrf    srf object to determine spikes from (smoothed SRF)
%       co          simply passed as cutoff to SRF::FindSpikes
%
% Output fields:
%
%       srf         colored SRF
%       sps         Vx1 vertex list of found spikes
%       npmin       Vx1 minimal product list per vertex
%       vdmin       Vx1 min vertex-to-neighbors distance
%       vdmean      Vx1 mean vertex-to-neighbors distance
%
% See also @xff/private/srf_FindSpikes

% Version:  v1.1
% Build:    16021119
% Date:     Feb-11 2016, 7:46 PM EST
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% special case: information given!
numvert = bc.NrOfVertices;
if nargin > 4 && isa(xo2, 'double') && isa(cutoff, 'double') && ...
    isa(varargin{1}, 'double') && isa(varargin{2}, 'double') && ...
    numel(xo2) == numvert && numel(cutoff) == numvert && ...
    numel(varargin{1}) == numvert && numel(varargin{2}) == 1 && ...
   ~isinf(varargin{2}) && ~isnan(varargin{2}) && varargin{2} > 0 && varargin{2} < 1
    mn1 = xo2;
    mn2 = cutoff;
    mn3 = varargin{1};
    cutoff = varargin{2};
    mi = find((mn1 .* (mn2 / mean(mn2)) .* (mn3 / mean(mn3))) <= cutoff);

% continue normally
else

    % second SRF given
    if nargin < 2 || numel(xo2) ~= 1 || ~xffisobject(xo2, true, 'srf')
        xo2 = xo;
    else
        bc2 = xo2.C;
        if bc2.NrOfVertices ~= bc.NrOfVertices
            error('neuroelf:xff:badArgument', 'For multiple SRF objects, NrOfVertices must match.');
        end
    end
    if nargin < 3
        cutoff = Inf;
    end

    % get information on secondary SRF
    [mi, mn1, mn2, mn3] = srf_FindSpikes(xo2, cutoff);
end

% normalize mn1, mn2, mn3 between [0..1]
r1 = [min(mn1), max(mn1) + eps];
r2 = [min(mn2), max(mn2) + eps];
r3 = [min(mn3), max(mn3) + eps];

% color srf
vc = zeros(numel(mn1), 4);
vc(mi, 1) = NaN;
vc(:, 2) = fix(32 + 200 * (mn1 - r1(1)) / (r1(2) - r1(1)));
vc(:, 3) = fix(32 + 200 * (mn2 - r2(1)) / (r2(2) - r2(1)));
vc(:, 4) = fix(32 + 200 * (mn3 - r3(1)) / (r3(2) - r3(1)));
bc.VertexColor = vc;
bc.ConvexRGBA = [1/8, 1/8, 1/8, 1];
xo.C = bc;

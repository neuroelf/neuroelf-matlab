function xo = srf_ShrinkSpikes(xo, factor)
% SRF::ShrinkSpikes  - shrink marked spikes
%
% FORMAT:       [srf] = srf.ShrinkSpikes([factor]);
%
% Input fields:
%
%       factor      shrinking factor (default: 1)
%
% Output fields:
%
%       srf         altered SRF

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:52 PM EST
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
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isa(factor, 'double') || numel(factor) ~= 1 || ...
    isinf(factor) || isnan(factor) || factor <= 0
    factor = 1;
else
    factor = real(factor);
end
factpp = factor + 1;

% get neighbors and colors
nb = bc.Neighbors;
co = bc.VertexColor;
cm = co(:, 1);

% iterate over spike vertices
for vc = find(isnan(cm(:)'))

    % first check that at least one neighbor is also marked
    nl = nb{vc, 2};
    if ~any(isnan(cm(nl)))
        co(vc, 1) = 0;

    % if good then combine color
    else
        co(vc, 2:4) = (factor * co(vc, 2:4) + sum(co(nl, 2:4))) / factpp;
    end
end

% store back
bc.VertexColor = co;
xo.C = bc;

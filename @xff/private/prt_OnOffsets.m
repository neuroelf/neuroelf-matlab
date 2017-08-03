function oo = prt_OnOffsets(xo)
% PRT::OnOffsets  - retrieve all on- and offsets with a condition marker
%
% FORMAT:       oo = prt.OnOffsets
%
% No input fields.
%
% Output fields:
%
%       oo          Ox3 (or wider, with parameters) on- and offsets
%
% Note: the first column in oo contains the condition number.

% Version:  v1.1
% Build:    16021017
% Date:     Feb-10 2016, 5:26 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% retrieve conditions
bc = xo.C;
c = bc.Cond;

% number of onsets and parameters
noo = sum(cellfun('size', {c.OnOffsets}, 1));
npar = max(cellfun('size', {c.Weights}, 2));

% create outputs and preset parameters to NaN
oo = zeros(noo, 3 + npar);
if npar > 0
    oo(:, 4:end) = NaN;
end

% iterate over conditions
oc = 1;
for cc = 1:numel(c)
    ot = oc + size(c(cc).OnOffsets, 1) - 1;
    oo(oc:ot, 1) = cc;
    oo(oc:ot, 2:3) = c(cc).OnOffsets;
    if npar > 0 && ~isempty(c(cc).Weights)
        oo(oc:ot, 4:(3+size(c(cc).Weights, 2))) = c(cc).Weights;
    end
    oc = ot + 1;
end

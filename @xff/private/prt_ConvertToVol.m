function xo = prt_ConvertToVol(xo, tr)
% PRT::ConvertToVol  - convert a ms-based PRT to volume-based
%
% FORMAT:       [prt] = prt.ConvertToVol(tr)
%
% Input fields:
%
%       tr          TR
%
% Output fields:
%
%       prt         altered protocol

% Version:  v1.1
% Build:    16021017
% Date:     Feb-10 2016, 5:36 PM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt') || ...
   ~isa(tr, 'double') || numel(tr) ~= 1 || isinf(tr) || isnan(tr) || tr < 1 || tr ~= fix(tr)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% return if not in volumes
if isempty(bc.ResolutionOfTime) || lower(bc.ResolutionOfTime(1)) ~= 'm'
    return;
end

% iterate over conditions
c = bc.Cond;
for cc = 1:numel(c)
    c(cc).OnOffsets = [round(c(cc).OnOffsets(:, 1) / tr) + 1, round(c(cc).OnOffsets(:, 2) / tr)];
    bv = find(c(cc).OnOffsets(:, 2) < c(cc).OnOffsets(:, 1));
    if ~isempty(bv)
        c(cc).OnOffsets(bv, 2) = c(cc).OnOffsets(bv, 1);
    end
end
bc.Cond = c;

% set to Volumes
bc.ResolutionOfTime = 'Volumes';
xo.C = bc;

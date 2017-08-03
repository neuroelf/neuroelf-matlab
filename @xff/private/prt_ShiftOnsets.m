function xo = prt_ShiftOnsets(xo, shiftval, conds)
% PRT::ShiftOnsets  - shift onsets in conditions of a PRT
%
% FORMAT:       [prt =] prt.ShiftOnsets(shiftval [, conds])
%
% Input fields:
%
%       shiftval    1x1 double value that will be added to OnOffsets
%       conds       optional condition selection (pattern or list, def: all)
%
% Output fields:
%
%       prt         altered PRT
%
% Examples:
%
%   prt.ShiftOnsets(4000);
%   prt.ShiftOnsets(2000, '.*response');
%

% Version:  v1.1
% Build:    16021016
% Date:     Feb-10 2016, 4:42 PM EST
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
   ~isa(shiftval, 'double') || numel(shiftval) ~= 1 || isinf(shiftval) || isnan(shiftval)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if shiftval == 0
    return;
end
bc = xo.C;

% get number of conditions
ncon = numel(bc.Cond);

% calc new color ?
if nargin < 3 || (~ischar(conds) && ~isa(conds, 'double'))
    conds = 1:ncon;
elseif ischar(conds)
    tconds = false(1, ncon);
    conds = conds(:)';
    for cc = 1:ncon
        if ~isempty(regexpi(bc.Cond(cc).ConditionName{1}, conds))
            tconds(cc) = true;
        end
    end
    conds = find(tconds);
else
    conds = conds(:)';
    conds(isinf(conds) | isnan(conds) | conds < 1 | conds > ncon) = [];
    conds = unique(round(conds));
end

for cc = conds
    bc.Cond(cc).OnOffsets = bc.Cond(cc).OnOffsets + shiftval;
end

% set new content
xo.C = bc;

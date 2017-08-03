function xo = glm_ReorderPredictors(xo, order)
% GLM::ReorderPredictors  - reorder predictors (RFX-GLM only!)
%
% FORMAT:       [glm = ] glm.ReorderPredictors(order);
%
% Input fields:
%
%       order       either 1xP double or cell array with new order
%
% No output fields.
%
% Using: findfirst, multimatch.

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:54 AM EST
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

% neuroelf library
global ne_methods;
findfirst  = ne_methods.findfirst;
multimatch = ne_methods.multimatch;

% check arguments
if nargin ~= 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
   (~isa(order, 'double') && ~iscell(order)) || isempty(order)
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end
bc = xo.C;
if bc.ProjectTypeRFX ~= 1
    error('neuroelf:xff:badObject', 'Only valid for RFX-GLMs.');
end
order = order(:);
if isa(order, 'double')
    if any(isinf(order) | isnan(order) | order < 1 | order > numel(order) | order ~= round(order)) || ...
        numel(order) ~= numel(unique(order))
        error('neuroelf:xff:badArgument', ...
            'Invalid numeric order, must be unique integer numbers >= 1.');
    end
else
    if ~all(cellfun(@ischar, order)) || any(cellfun('isempty', order)) || ...
        numel(order) ~= numel(unique(order))
        error('neuroelf:xff:badArgument', ...
            'Invalid string order, must be unique strings.');
    end
end

% number of subject predictors (without constant!)
nsp = bc.NrOfSubjectPredictors - 1;
if ~any(numel(order) == (nsp + [0, 1]))
    error('neuroelf:xff:badArgument', 'Insufficient order argument.');
end

% for strings
if iscell(order)

    % linearize (to be sure)
    for pc = 1:numel(order)
        order{pc} = order{pc}(:)';
    end

    % get first subject predictor names (see GLM::SubjectPredictors)
    p = {bc.Predictor(:).Name2};
    pl = p{end};
    fs = p{1}(1:findfirst(p{1} == ':'));
    p = strrep(p, [fs ' '], '');
    keep = cellfun('isempty', regexp(p, '^Subject\s+.*:\s+')) & ~strcmpi(p, pl);
    p = p(keep);
    p = p(:);

    % test
    order = multimatch(order, p);
    if any(order == 0) || (numel(order) < numel(p) && any(order == numel(p))) || ...
        numel(order) ~= numel(unique(order))
        error('neuroelf:xff:badArgument', 'Invalid strings for reordering.');
    end
end

% remove constant (must always be the last predictor)
if numel(order) > nsp
    order(order == numel(order)) = [];
end
order = order(:)';

% number of BetaMap spatial dims
nsd = ndims(bc.GLMData.RFXGlobalMap);
spi = repmat({':'}, 1, nsd);
spi{1, end+1} = [order, nsp + 1];

% now do the work, start with Predictors
p = 1:numel(bc.Predictor);
for sc = 1:numel(bc.GLMData.Subject);

    % change that range
    p((sc-1)*nsp+1:sc*nsp) = order + (sc - 1) * nsp;

    % then alter BetaMaps
    bm = bc.GLMData.Subject(sc).BetaMaps;
    if istransio(bm)
        bm = resolve(bm);
    end
    bc.GLMData.Subject(sc).BetaMaps = bm(spi{:});
end
bc.Predictor = bc.Predictor(:)';
bc.Predictor = bc.Predictor(1, p);

% finally, change contrasts
rtv = bc.RunTimeVars;
if isfield(rtv, 'Contrasts') && iscell(rtv.Contrasts) && ~isempty(rtv.Contrasts)
    for cc = 1:size(rtv.Contrasts, 1)
        rtv.Contrasts{cc, 2} = rtv.Contrasts{cc, 2}(:);
        if numel(rtv.Contrasts{cc, 2}) == nsp
            rtv.Contrasts{cc, 2} = rtv.Contrasts{cc, 2}(order, 1);
        elseif numel(rtv.Contrasts{cc, 2}) == (nsp + 1)
            rtv.Contrasts{cc, 2} = rtv.Contrasts{cc, 2}([order, nsp + 1], 1);
        end
    end
    bc.RunTimeVars.Contrasts = rtv.Contrasts;
end

% set back in GLM
xo.C = bc;

function xo2 = prt_ShuffleLabels(xo, labels, minons)
% PRT::ShuffleLabels  - shuffle labels within specified conditions
%
% FORMAT:       sprt = prt.ShuffleLabels(labels [, minons]);
%
% Input fields:
%
%       labels      1xL cell array with labels (condition names) to shuffle
%       minons      1x1/1xL minimum number of onsets per label (default: 1)
%
% Output fields:
%
%       sprt        altered PRT with shuffled labels
%
% Using: findfirst.

% Version:  v1.1
% Build:    16021016
% Date:     Feb-10 2016, 4:41 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt') || ...
   ~iscell(labels) || isempty(labels) || ~ischar(labels{1}) || isempty(labels{1})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
ncon = numel(bc.Cond);
if ncon == 0
    return;
end
connames = {bc.Cond.ConditionName};
connames = cat(1, connames{:});
labels = labels(:)';
nlab = numel(labels);
labcond = zeros(1, nlab);
labcondn = labcond;
for lc = 1:nlab
    if ~ischar(labels{lc}) || isempty(labels{lc}) || ~any(strcmpi(labels{lc}, connames))
        error('neuroelf:xff:badArgument', 'Invalid list of labels in call.');
    end
    labcond(lc) = ne_methods.findfirst(strcmpi(labels{lc}, connames));
    labcondn(lc) = size(bc.Cond(labcond(lc)).OnOffsets, 1);
end
if nlab ~= numel(unique(lower(labels)))
    error('neuroelf:xff:badArgument', 'Labels must be unique.');
end
if nargin < 3 || ~isa(minons, 'double') || ~any([1, nlab] == numel(minons)) || ...
    any(isinf(minons(:)) | isnan(minons(:)))
    minons = 1;
end
if numel(minons) == 1 && nlab > 1
    minons = minons .* ones(1, nlab);
else
    minons = minons(:)';
end
minons = min(minons, labcondn);

% get list of labels to use
numlab = sum(labcondn);
labnum = zeros(1, numlab);
lt = 1;
for lc = 1:nlab
    labnum(lt:lt+minons(lc)-1) = labcond(lc);
    lt = lt + minons(lc);
end

% add random entries if desired
if lt <= numlab
    labnum(lt:end) = labcond(ceil(nlab .* rand(1, numlab + 1 - lt)));
end

% shuffle labels
[trash, labshuffle] = sort(rand(numlab, 1));
labnum = labnum(labshuffle);

% get on-offsets of selected conditions
try
    labcondo = cat(1, bc.Cond(labcond).OnOffsets);
    labcondw = cat(1, bc.Cond(labcond).Weights);
    if ~isempty(labcondw) && size(labcondw, 1) ~= size(labcondo, 1)
        error('BAD_PARAMS');
    end
catch xfferror
    neuroelf_lasterr(xfferror);
    error('neuroelf:xff:badArgument', ...
        'Conditions with swapped labels must have the same number of parameters.');
end

% create copy
xo2 = aft_CopyObject(xo);

% assign new on/offsets
for lc = 1:nlab
    bc.Cond(labcond(lc)).OnOffsets = labcondo(labnum == labcond(lc), :);
    if ~isempty(labcondw)
        bc.Cond(labcond(lc)).Weights = labcondw(labnum == labcond(lc), :);
    end
    bc.Cond(labcond(lc)).NrOfOnOffsets = sum(labnum == labcond(lc));
end

% set new content
xo2.C = bc;

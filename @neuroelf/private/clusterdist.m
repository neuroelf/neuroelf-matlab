function pd = clusterdist(p, s)
% clusterdist  - compute clustering distance between groups of points
%
% FORMAT:       pcd = clusterdist(p, s)
%
% Input fields:
%
%       p           points (PxD double)
%       s           labels (Px1 cell array of strings with subset IDs)
%
% Output fields:
%
%       pcd          point-clustering-distance measure
%

% Version:  v0.9b
% Build:    11051114
% Date:     Apr-08 2011, 10:18 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
if nargin < 2 || ...
   ~isa(p, 'double') || ...
    ndims(p) > 2 || ...
    isempty(p) || ...
    any(isinf(p(:)) | isnan(p(:))) || ...
   ~iscell(s) || ...
    numel(s) ~= size(p, 1)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
s = s(:);

% number of points
np = size(p, 1);

% get unique IDs (titles, ...)
us = uunion(lower(s), {});
ns = numel(us);
if ns == 1
    error( ...
        'neuroelf:BadArgument', ...
        'Not enough subsets in set given.' ...
    );
end

% find rows that match id
usrows = cell(ns, 1);
for sc = 1:ns
    usrows{sc} = find(strcmpi(s, us{sc}));
end

% initialize output
pd = zeros(np, 1);

% for each point
for pc = 1:size(p, 1)

    % initialize the metric
    pdc = 0;

    % reverse set
    sids = 1:ns;
    sids(strcmpi(us, s{pc})) = [];

    % iterate over sets (other than own)
    for sc = sids

        % compute distances between points
        pds = psetdists(p(usrows{sc}, :), p(pc, :));

        % add minimum to metric
        pdc = pdc - log(1 ./ (1 + min(pds(:))));
    end

    % store in output
    pd(pc) = pdc;
end

% divide by number of studies - 1
pd = pd ./ (ns - 1);

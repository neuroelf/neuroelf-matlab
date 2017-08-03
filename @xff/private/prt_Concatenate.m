function xo = prt_Concatenate(xo, xo2, nvol, tr)
% PRT::Concatenate  - concatenate two PRTs in time
%
% FORMAT:       [prt] = prt.Concatenate(prt2, tvol, tr)
%
% Input fields:
%
%       prt2        PRT object used for concatenation
%       nvol        number of volumes to shift
%       tr          tr
%
% Output fields:
%
%       prt         altered PRT object
%
% Note: conditions that are only present in 2nd PRT will be added

% Version:  v1.1
% Build:    16021017
% Date:     Feb-10 2016, 5:41 PM EST
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
if nargin < 4 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt') || ...
    numel(xo2) ~= 1 || ~xffisobject(xo2, true, 'prt') || ...
    numel(nvol) ~= 1 || ~isa(nvol, 'double') || isinf(nvol) || isnan(nvol) || nvol <= 0 || ...
    numel(tr) ~= 1 || ~isa(tr, 'double') || isinf(tr) || isnan(tr) || tr <= 0
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc1 = xo.C;
bc2 = xo2.C;
if ~strcmpi(bc1.ResolutionOfTime, bc2.ResolutionOfTime)
    error('neuroelf:xff:badArgument', 'Protocols must have same ResolutionOfTime setting.');
end
if lower(bc1.ResolutionOfTime(1)) == 'm'
    tdiff = nvol * tr;
else
    tdiff = nvol;
end
ncon1 = numel(bc1.Cond);
ncon2 = numel(bc2.Cond);

% get names
names1 = cell(1, ncon1);
names2 = cell(1, ncon2);
for nc = 1:ncon1
    names1{nc} = bc1.Cond(nc).ConditionName{1};
end
for nc = 1:ncon2
    names2{nc} = bc2.Cond(nc).ConditionName{1};
end

% iterate over conditions of 2nd PRT
for nc = 1:ncon2

    % find match
    match = strcmpi(names2{nc}, names1);

    % if match
    if any(match)
        match = find(match);
        match = match(1);

        % take care of weights
        if size(bc1.Cond(match).Weights, 1) ~= size(bc1.Cond(match).OnOffsets, 1)
            bc1.Cond(match).Weights = zeros(size(bc1.Cond(match).OnOffsets, 1), 0);
        end
        if size(bc2.Cond(nc).Weights, 1) ~= size(bc2.Cond(nc).OnOffsets, 1)
            bc2.Cond(nc).Weights = zeros(size(bc2.Cond(nc).OnOffsets, 1), 0);
        end
        if size(bc1.Cond(match).Weights, 2) < size(bc2.Cond(nc).Weights, 2)
            bc1.Cond(match).Weights(:, end+1:size(bc2.Cond(nc).Weights, 2)) = 0;
        elseif size(bc1.Cond(match).Weights, 2) > size(bc2.Cond(nc).Weights, 2)
            bc2.Cond(nc).Weights(:, end+1:size(bc1.Cond(match).Weights, 2)) = 0;
        end

        % add to equally named condition with offset
        bc1.Cond(match).OnOffsets = [bc1.Cond(match).OnOffsets; bc2.Cond(nc).OnOffsets + tdiff];
        bc1.Cond(match).NrOfOnOffsets = size(bc1.Cond(match).OnOffsets, 1);
        bc1.Cond(match).Weights = [bc1.Cond(match).Weights; bc2.Cond(nc).Weights];

    % otherwise
    else

        % add new condition and shift onsets
        bc1.Cond(end+1) = bc2.Cond(nc);
        bc1.Cond(end).OnOffsets = bc1.Cond(end).OnOffsets + tdiff;
    end
end

% store into old structure
bc1.NrOfConditions = numel(bc1.Cond);
xo.C = bc1;

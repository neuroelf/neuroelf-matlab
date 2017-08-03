function [cs, tm, cc, crc, cst] = clustermeshmapbin(m, nei, crd, tra, trb, st, cm)
% clustermeshmapbin  - cluster binary mesh-based map
%
% FORMAT:       [cs, tm, vl, cl, c] = clustermeshmapbin(m, nei, c, ta, tb, st, cm)
%
% Input fields:
%
%       m           binary map (Vx1)
%       nei         neighborhood information
%       crd         coordinates
%       tra         triangle areas
%       trb         vertex->triangle reference
%       st          size threshold (in coordinate units^2)
%       cm          clustering method: default: 1
%                   1 vertex-connectivity
%                   2 edge-connectivity
%
% Output fields:
%
%       cs          cluster sizes
%       tm          clustered map (uint32 with 0's for unclustered vertices)
%       vl          vertex list [index, clusterindex]
%       cl          Vx4 list of cluster coordinates [xyz, clusterindex]
%       c           Cx1 struct array with
%        .coords    Vx3 coordinates
%        .vertices  Vx1 vertex indices

% Version:  v0.9a
% Build:    11050223
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
if nargin < 5 || ...
   ~islogical(m) || ...
    numel(m) ~= size(m, 1) || ...
   ~iscell(nei) || ...
    numel(nei) ~= numel(m) || ...
   ~isa(crd, 'double') || ...
   ~isequal(size(crd), [numel(m), 3]) || ...
    any(isinf(crd(:)) | isnan(crd(:))) || ...
   ~isa(tra, 'double') || ...
    numel(tra) ~= size(tra, 1) || ...
    any(isinf(tra) | isnan(tra) | tra < 0) || ...
   ~iscell(trb) || ...
    numel(trb) ~= numel(m)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 6 || ...
   ~isa(st, 'double') || ...
    numel(st) ~= 1 || ...
    isinf(st) || ...
    isnan(st) || ...
    st < 0
    st = 0;
end
if nargin < 7 || ...
   ~isnumeric(cm) || ...
    numel(cm) ~= 1 || ...
   ~any((1:2) == cm)
    cm = 1;
else
    cm = double(cm);
end
if cm == 2
    warning( ...
        'neuroelf:NotYetImplemented', ...
        'Not yet implemented: Edge-connectivity; falling back...' ...
    );
    cm = 1;
end

% find indices of "above threshold" vertices
cc = find(m);
total = numel(cc);

% empty map
if total == 0
    cs = zeros(0, 1);
    tm = uint32(0);
    tm(numel(m), 1) = 0;
    cc = zeros(0, 2);
    crc = zeros(0, 4);
    cst = emptystruct({'coords', 'vertices'}, [0, 1]);
    return;
end

% otherwise
cc = [cc(:), zeros(total, 1)];

% copy m into numeric array
m = int32(m) - 1;

% keep track of "clustered vertices"
done = 0;
nextc = 1;
nextv = 1;

% start with first index
if cm == 1
    while done < total

        % continue nextv til next to cluster
        while (m(cc(nextv)) ~= 0)
            nextv = nextv + 1;
        end

        % mark first index and then traverse neighbors
        thisc = 0;
        nextcv = cc(nextv);

        % repeat until no more vertices found
        while ~isempty(nextcv)
            m(nextcv) = nextc;
            thisc = thisc + numel(nextcv);
            lastcv = nextcv;
            lastnei = nei(lastcv);
            nextcv = cat(2, lastnei{:});
            nextcv = unique(nextcv(m(nextcv) == 0));
        end

        % increase nextc and done
        nextc = nextc + 1;
        done = done + thisc;
    end
else
    error( ...
        'neuroelf:NotYetImplemented', ...
        'Not yet implemented feature: EDGE-CONN.' ...
    );
end

% make good m,tm (first step)
m(m < 0) = 0;
tm = uint32(m);

% create list
nextc = nextc - 1;
cl = cell(1, nextc);

% sort list of marked vertices
[sm, si] = sort(m(cc(:, 1)));

% get "onsets"
io = [1; 1 + find(diff(sm(:))); total + 1];

% initialize cluster sizes
cs = zeros(nextc, 1);

% put numbers into cc and cl
athird = 1 / 3;
cst = emptystruct({'coords', 'vertices'}, [nextc, 1]);
csr = false(nextc, 1);
for clc = 1:nextc
    cc(si(io(clc):io(clc+1)-1), 2) = clc;
    cl{clc} = cc(si(io(clc):io(clc+1)-1));
    cs(clc) = athird * sum(tra(cat(2, trb{cl{clc}})));
    cst(clc).vertices = cl{clc};
    cst(clc).coords = crd(cl{clc}, :);

    % apply size treshold
    if cs(clc) < st
        cc(si(io(clc):io(clc+1)-1), 2) = -1;
        csr(clc) = true;
    end
end

% remove from thresholded map
if any(csr)
    tm(cc(cc(:, 2) < 0)) = 0;

    % and all other outputs
    cs(csr) = [];
    cc(cc(:, 2) < 0, :) = [];
    cst(csr) = [];
end

% create coordinate list
crc = [crd(cc(:, 1), :), cc(:, 2)];

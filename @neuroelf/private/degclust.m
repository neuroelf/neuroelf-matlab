function d = degclust(vol, k, clconn, t)
% degclust  - degree of clustering of data
%
% FORMAT:       d = degclust(vol, k [, clconn, [ t]])
%
% Input fields:
%
%       vol         either binary (logical) or double 3D data
%       k           size threshold (in voxels)
%       clconn      clustering connectivity (see clustercoordsc, default: 2)
%       t           threshold (required for double data)
%
% Output fields:
%
%       d           degree of clustering ([0 .. 1])

% Version:  v0.9a
% Build:    10051716
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
if nargin < 2 || ...
   (~isnumeric(vol) && ...
    ~islogical(vol)) || ...
    ndims(vol) ~= 3 || ...
   ~isa(k, 'double') || ...
    numel(k) ~= 1 || ...
    isinf(k) || ...
    isnan(k) || ...
    k < 1 || ...
    k ~= fix(k)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing vol/k argument to cluster.' ...
    );
end
if nargin < 3 || ...
   ~isa(clconn, 'double') || ...
    numel(clconn) ~= 1 || ...
    isinf(clconn) || ...
    isnan(clconn) || ...
   ~any(1:5 == clconn)
    clconn = 2;
end
if ~islogical(vol)
    if nargin < 4 || ...
       ~isa(t, 'double') || ...
        numel(t) ~= 1 || ...
        t <= 0
        error( ...
            'neuroelf:BadArgument', ...
            'Bad or missing threshold for clustering.' ...
        );
    end
    if ~isa(vol, 'double')
        vol = double(vol(:, :, :));
    end
end

% for binary only one way
if ~isa(vol, 'double')

    % cluster
    cs = clustercoordsc(vol, clconn, k);

    % compute degree of clustering
    d = sum(cs) ./ sum(vol(:));

% for double values, allow both tails to be clustered separately
else
    pvol = (vol >= t);
    nvol = ((-vol) >= t);
    pcs = clustercoordsc(pvol, clconn, k);
    ncs = clustercoordsc(nvol, clconn, k);

    % compute as sums
    d = (sum(pcs) + sum(ncs)) ./ (sum(pvol(:)) + sum(nvol(:)));
end

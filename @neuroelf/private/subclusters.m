function sc = subclusters(vox, val, sep)
% subclusters  - separate mega cluster into sub-clusters
%
% FORMAT:       sc = subclusters(vox [, val], sep)
%
% Input fields:
%
%       vox         either a 3D volume with values or a Vx3 list
%       val         required if vox is a list of coordinates
%       sep         cluster separation (in voxel resolution)
%
% Output fields:
%
%       sc          1xC sub clusters structure
%        .peak      peak coordinate
%        .val       values
%        .vox       voxel coordinates

% Version:  v0.9a
% Build:    10062206
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
if nargin < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

if numel(val) == 1
    sep = val;
    msk = double(vox > 0);
    [vox, voy, voz] = ind2sub(size(msk), find(msk));
    vox = [vox, voy, voz];
else
    vox = round(vox);
    vox = 1 + vox - ones(size(vox, 1), 1) * min(vox);
    msk = zeros(max(vox));
    msk(sub2ind(max(vox), vox(:, 1), vox(:, 2), vox(:, 3))) = 1;
end

% sort values
[sv, svi] = sort(val);

% start with marking highest value as first cluster
cc = 2;
msk(vox(svi(1), 1), vox(svi(1), 2), vox(svi(1), 3)) = cc;

% iterate until all values are used

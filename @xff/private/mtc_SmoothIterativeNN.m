function xo = mtc_SmoothIterativeNN(xo, srf, nweight, numiter)
% MTC::SmoothIterativeNN  - iteratively smooth MTC (nearest neighbor)
%
% FORMAT:       [mtc] = mtc.SmoothIterativeNN(srf, nweight, numiter);
%
% Input fields:
%
%       srf         matching surface file (number of vertices)
%       nweight     neighbor weighting (divided by num. of neighbors)
%       numiter     number of iterations

% Version:  v1.1
% Build:    16020917
% Date:     Feb-09 2016, 5:06 PM EST
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
if nargin < 4 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mtc') || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, 'srf') || ...
    numel(nweight) ~= 1 || ~isa(nweight, 'double') || ...
    isinf(nweight) || isnan(nweight) || nweight <= 0 || ...
    numel(numiter) ~= 1 || ~isa(numiter, 'double') || ...
    isinf(numiter) || isnan(numiter) || numiter ~= fix(numiter)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
srfc = srf.C;
if bc.NrOfVertices ~= srfc.NrOfVertices
    error('neuroelf:xff:invalidObject', 'NrOfVertices mismatch between MTC and SRF.');
end

% get neighbors
nei = srfc.Neighbors;
nv = size(nei, 1);
normw = nweight + 1;

% get initial MTC data (resolve transio if needed!)
mtcd = double(bc.MTCData(:, :));

% iterate
for ic = 1:numiter

    % copy data
    mtcdo = mtcd;

    % iterate over vertices
    for vc = 1:nv

        % calculus
        mtcd(:, vc) = (mtcdo(:, vc) + (nweight / nei{vc, 1}) * sum(mtcdo(:, nei{vc, 2}), 2)) / normw;
    end
end

% store result in MTC
bc.MTCData = single(mtcd);
xo.C = bc;

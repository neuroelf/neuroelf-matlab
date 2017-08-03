function ovecs = orthvecs(ivecs, vorder, detr)
% orthvecs  - orthogonalize vectors against each other
%
% FORMAT:       ovecs = orthvecs(ivecs [, vorder, detr])
%
% Input fields:
%
%       ivecs       NxV matrix with V number of vectors
%       vorder      must contain all numbers from 1..V, default 1:V
%       detr        if given and evaluates to true, detrend vectors
%
% Output fields:
%
%       ovecs       orthogonalized vectors
%
% Note: for two vectors, the "cheaper" orthvec function can be used

% Version:  v0.9a
% Build:    11111518
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
if nargin < 1 || ...
   ~isa(ivecs, 'double') || ...
    isempty(ivecs) || ...
    ndims(ivecs) > 2 || ...
    size(ivecs, 1) < size(ivecs, 2) || ...
    any(isinf(ivecs(:)) | isnan(ivecs(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Valid matrix of vectors without Inf/Nans are required.' ...
    );
end
numvals = size(ivecs, 1);
numvecs = size(ivecs, 2);
if nargin > 1 && ...
    isa(vorder, 'double') && ...
    numel(vorder) == numvecs && ...
    numel(vorder) == max(size(vorder)) && ...
    all(vorder > 0) && ...
    numel(unique(fix(vorder))) == numvecs
    vorder = fix(vorder(:)');
else
    vorder = 1:numvecs;
end
if nargin > 2 && ...
   ~isempty(detr) && ...
   (isnumeric(detr) || islogical(detr)) && ...
    detr(1)
    detr = true;
else
    detr = false;
end

% generate output
ovecs = zsz(ivecs);
tcol = vorder(1);
if detr
    ovecs(:, tcol) = ivecs(:, tcol) - mean(ivecs(:, tcol));
else
    ovecs(:, tcol) = ivecs(:, tcol);
end

% iterate in order of importance
for vc = vorder(2:numvecs)

    % get vector which is to be tested
    testvec = ivecs(:, vc);
    if detr
        testvec = detrend(testvec);
    end

    % regress on partial matrix
    rmat = ovecs(:, tcol);
    rmat = rmat - repmat(mean(rmat), [numvals, 1]);
    tmat = rmat';
    b = (tmat * rmat) \ (tmat * testvec);

    % fitted -> residual
    r = testvec - rmat * b;

    % residual goes into output !
    ovecs(:, vc) = r;
    tcol = [tcol, vc];
end

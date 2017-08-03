function tstat = glmtstat(c, b, iXX, se, rw)
% tstat  - compute GLM t statistic
%
% FORMAT:       tstat = glmtstat(c, b, iXX, se [, rw])
%
% Input fields:
%
%       c       1xP contrast vector (right padded with zeros)
%       b       beta estimates (double or single)
%       iXX     pinv of X'*X
%       se      standard error of estimation (double or single)
%       rw      robust weights (to adjust d.f.)
%
% Output fields:
%
%       tstat   t-statistic
%
% See also @xff/private/sdm_CalcBetas

% Version:  v1.0
% Build:    15122921
% Date:     Dec-29 2015, 9:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2013, 2015, Jochen Weber
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
if nargin < 4 || ...
   ~isa(c, 'double') || ...
   (~isa(b, 'double') && ...
    ~isa(b, 'single')) || ...
   ~isa(iXX, 'double') || ...
   (~isa(se, 'double') && ...
    ~isa(se, 'single')) || ...
    isempty(c) || ...
    isempty(b) || ...
    ndims(iXX) > 3 || ...
    size(iXX, 1) ~= size(iXX, 2) || ...
   (size(iXX, 1) > 1 && ...
    ~any(size(b) == size(iXX, 1))) || ...
    numel(se) ~= (numel(b) / size(iXX, 1)) || ...
   (size(iXX, 3) ~= 1 && ...
    size(iXX, 3) ~= numel(se))
    error( ...
        'neuroelf:BadArgument', ...
        'Missing, invalid or bad sized argument given.' ...
    );
end

% check sizes
csize = size(iXX, 1);
c = c(:)';
if numel(c) < csize
    c = [c, zeros(1, csize - numel(c))];
else
    c = c(1:csize);
end

% determine beta dimension
bd = find(size(b) == csize);
if isempty(bd)
    bd = ndims(b) + 1;
else
    bd = bd(end);
end

% if needed, reorder beta
if bd(1) ~= 1
    b = reshape(permute(b, [bd(1), setdiff(1:numel(size(b)), bd(1))]), ...
        [csize, numel(se)]);
end

% calc stats
if size(iXX, 3) == 1
    
    % ensure iXX's bad entries don't invalidate stats
    ciXX = iXX(c ~= 0, c ~= 0);
    cc = c(1, c ~= 0);
    ciXXc = cc * ciXX * cc';
    
    % compute t-stat
    tstat = reshape((c * b) ./ (se(:)' * sqrt(ciXXc)), size(se));
else
    tstat = reshape((c * b) ./ (se(:)' .* lsqueeze(sqrt( ...
        mtimesnd(mtimesnd(repmat(c, [1, 1, numel(se)]), iXX), repmat(c', [1, 1, numel(se)]))))'), size(se));
end
tstat(isinf(tstat) | isnan(tstat)) = 0;

% adjust d.f.
if nargin > 4 && ...
    isa(rw, 'double')

    % set size
    srw = size(rw);
    drw = numel(srw);
    sse = size(se);
    while(sse(end) == 1)
        sse(end) = [];
    end

    % single set of weights
    if numel(rw) == max(size(rw))
        tstat = sdist('tinv', sdist('tcdf', tstat, sum(rw(:)) - csize), numel(rw) - csize);

    % otherwise sizes must match
    elseif isequal(sse, srw(1:end-1))

        tstat = sdist('tinv', sdist('tcdf', tstat, sum(rw, drw) - csize), ...
            size(rw, drw) - csize);
    end
end

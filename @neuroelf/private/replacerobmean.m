function data = replacerobmean(data, dim, kin)
% noinfnan  - replace valid values with more robust mean estimate
%
% FORMAT:       data = replacerobmean(data [, dim [, kin]])
%
% Input fields:
%
%       data        N-d single/double input
%       dim         optional dimension, default: first non-singleton
%       kin         keep Inf/NaN values (default: false)
%
% Output fields:
%
%       data        data with amended content

% Version:  v0.9b
% Build:    11051718
% Date:     Mar-17 2011, 8:50 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
    isempty(data) || ...
   (~isa(data, 'double') && ...
    ~isa(data, 'single'))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
    isinf(dim) || ...
    isnan(dim) || ...
    dim < 1 || ...
    dim > ndims(data) || ...
    dim ~= fix(dim)
    dim = findfirst(size(data) > 1);
    if isempty(dim)
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid data (size).' ...
        );
    end
end
if nargin < 3 || ...
   ~islogical(kin) || ...
    numel(kin) ~= 1
    kin = false;
end

% first, compute mean and replace bad elements
[md, ge] = meannoinfnan(data, dim);
rma = ones(1, ndims(data));
rma(dim) = size(data, dim);
mdd = repmat(md, rma);
badd = false;
if ~all(ge(:))
    badd = true;
    if kin
        kinv = data(~ge);
    end
    data(~ge) = mdd(~ge);
end

% compute difference to mean, relative to STD
wmax = 1 ./ sdist('normpdf', 0, 0, 2);
w = limitrangec(wmax .* sdist('normpdf', ...
    repmat(1 ./ limitrangec(sqrt(varc(data, dim, true)), 1e-10, 1e10, 1e10), rma) .* ...
    abs(data - mdd), 0, 2), 0, 1, 0);

% robust needed?
maxiter = 30;
while any(w(:) < 0.5) && ...
    maxiter > 0

    % compute weights (as a factor)
    data = w .* data + (1 - w) .* mdd;
    md = sum(data, dim) ./ size(data, dim);
    mdd = repmat(md, rma);
    if badd
        data(~ge) = mdd(~ge);
    end
    w = limitrangec(wmax .* sdist('normpdf', ...
        repmat(1 ./ limitrangec(sqrt(varc(data, dim, true)), 1e-10, 1e10, 1e10), rma) .* ...
        abs(data - mdd), 0, 2), 0, 1, 0);
    maxiter = maxiter - 1;
end
if kin && ...
    badd
    data(~ge) = kinv;
end

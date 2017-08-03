function s = spatent(m, n, rz)
% spatent  - spatial entropy of map values
%
% FORMAT:       s = spatent(m [, n [, rz]])
%
% Input fields:
%
%       m           3D map data (0-values are removed)
%       n           optional number of bins (default: number of values/100)
%       rz          flag, remove zeros (default: false)
%
% Output fields:
%
%       s           spatial entropy

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
if nargin < 1 || ...
   ~isnumeric(m) || ...
    numel(m) < 4 || ...
    any(isinf(m(:)) | isnan(m(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% get values of interest
if ~isa(m, 'double')
    m = double(m);
end
m = m(:);
if nargin > 2 && ...
    islogical(rz) && ...
    numel(rz) == 1 && ...
    rz
    m = m(m ~= 0);
end
nm = numel(m);

% number of bins given?
if nargin < 2 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n < 3 || ...
    n > (nm / 2)
    n = max(3, round(nm / 100));
else
    n = round(n);
end

% compute histogram
mm = minmaxmean(m);
md = mm(2) - mm(1);
hm = histc(m, mm(1):md/n:mm(2));

% remove empty bins
hm(hm == 0) = [];

% then compute "relative" histogram (where sum = 1)
hm = (1 / numel(m)) .* hm;

% and finally compute spatial entropy
s = sum(hm .* log2(hm));

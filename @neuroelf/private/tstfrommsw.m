function [t, ut, sd] = tstfrommsw(m1, m2, s1, s2, w1, w2)
% tstfrommsw  - two-sample t-test from means, stds, and weights
%
% FORMAT:       t = tstfrommsw(m1, m2, s1, s2, w1, w2)
%
% Input fields:
%
%       m1, m2      means for samples 1 and 2 (all must match in size)
%       s1, s2      STDs for samples 1 and 2
%       w1, w2      weights (weighted sample sizes) for samples 1 and 2
%       opts        optional settings
%        .olsalso   flag, do OLS regression as well (default: false)
%        .olsonly   only OLS (for other functions, default: false)
%        .olstpat   OLS t-vol file pattern (def.: [pwd/G%d-G%d_ot.img])
%        .tvolpat   t-volume filename pattern (def.: [pwd/G%d-G%d_t.img])
%        .usenulls  flag, accept 0 values as real (default: false);
%        .wvols     create weight volumes
%        .wvolsinv  invert weight volumes (so that data is 1 - w)
%        .wvolspat  w-volume filename pattern (def.: [pwd/sub%03d_w.img])
%
% Output fields:
%
%       tv          either set of pairwise t-test volumes or XxYxZxP data
%       wv          either set of weighting volumes or XxYxZxV data
%       ot          either set of pairwise OLS t-test vols or data
%
% Note: t-stats will be re-computed to nominal d.f. !

% Version:  v0.9b
% Build:    13041418
% Date:     Apr-08 2011, 9:16 PM EST
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
if nargin < 6 || ...
   (~isa(m1, 'single') && ...
    ~isa(m1, 'double')) || ...
   (~isa(m2, 'single') && ...
    ~isa(m2, 'double')) || ...
   (~isa(s1, 'single') && ...
    ~isa(s1, 'double')) || ...
   (~isa(s2, 'single') && ...
    ~isa(s2, 'double')) || ...
   (~isa(w1, 'single') && ...
    ~isa(w1, 'double')) || ...
   (~isa(w2, 'single') && ...
    ~isa(w2, 'double')) || ...
   ~isequal(size(m1), size(m2)) || ...
   ~isequal(size(m1), size(s1)) || ...
   ~isequal(size(m1), size(s2)) || ...
   ~isequal(size(m1), size(w1)) || ...
   ~isequal(size(m1), size(w2))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument.' ...
    );
end

% ensure double precision
if ~isa(m1, 'double')
    m1 = double(m1);
end
if ~isa(m2, 'double')
    m2 = double(m2);
end
if ~isa(s1, 'double')
    s1 = double(s1);
end
if ~isa(s2, 'double')
    s2 = double(s2);
end
if ~isa(w1, 'double')
    w1 = double(w1);
end
if ~isa(w2, 'double')
    w2 = double(w2);
end

% convert values according to required formula
s1 = s1 .* s1 ./ w1;
s2 = s2 .* s2 ./ w2;

% compute raw t-score
md = m1 - m2;
ut = md ./ sqrt(s1 + s2);

% compute DF correction factors
df = w1 + w2 - 2;
df12 = ((s1 + s2) .^ 2) ./ (s1 .* s1 ./ (w1 - 1) + s2 .* s2 ./ (w2 - 1));

% find and remove bad values
badv = find(isinf(ut) | isnan(ut) | isnan(df12) | df12 < 1);
ut(badv) = 0;
df12(badv) = 1;

% then re-compute t-stats with correction factor
t = -sign(ut) .* sdist('tinv', sdist('tcdf', -abs(ut), df12), df);

% compute SD
if nargout > 2
    sd = (md ./ t) .* sqrt(df);
    sd(isnan(sd)) = Inf;
end

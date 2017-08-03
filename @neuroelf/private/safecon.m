function bc = safecon(b, w)
% safecon  - safe contrast computation taking missing values into account
%
% FORMAT:       bc = safecon(b, w)
%
% Input fields:
%
%       b           beta estimates
%       w           contrast weighting vector
%
% Output fields:
%
%       bc          beta contrast

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
if nargin < 2 || ...
   ~isnumeric(b) || ...
   ~isa(w, 'double') || ...
    size(b, ndims(b)) ~= numel(w) || ...
    any(isinf(w(:)) || isnan(w(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% find weights that are ~= 0
w = w(:)';
wn = find(w < 0);
wp = find(w > 0);

% check contrast
if ~isempty(wn) && ...
   ~isempty(wp)
    if sum(w([wn, wp])) ~= 0
        warning( ...
            'neuroelf:Warning', ...
            'Differential contrast sum ~= 0!' ...
        );
    end
end

% re-weigh both parts to one
w(wn) = (-1 / sum(w(wn))) .* w(wn);
w(wp) = (1 / sum(w(wp))) .* w(wp);

% build contrast output
sb = size(b);
sb(end) = 1;
bn = zeros(sb);
bp = zeros(sb);
bs = zeros(sb);

% build subsref argument
sr = struct('type', '()', 'subs', {repmat({':'}, 1, numel(sb))});

% start with negative weights
for wc = wn
    sr.subs{end} = wc;
    bt = double(subsref(b, sr));
    bt(isinf(bt) | isnan(bt)) = 0;
    bn = bn + w(wc) .* bt;
    bs = bs + w(wc) .* (bt ~= 0);
end
bn = bn ./ abs(bs);
bn(bs == 0) = 0;

% continue with positive weights
bs = zeros(sb);
for wc = wp
    sr.subs{end} = wc;
    bt = double(subsref(b, sr));
    bt(isinf(bt) | isnan(bt)) = 0;
    bp = bp + w(wc) .* bt;
    bs = bs + w(wc) .* (bt ~= 0);
end
bp = bp ./ abs(bs);
bp(bs == 0) = 0;

% combine
bc = bn + bp;

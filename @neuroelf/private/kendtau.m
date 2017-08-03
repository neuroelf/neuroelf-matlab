function [t, np, nn] = kendtau(v, m, d)
% kendtau  - returns the Kendall tau statistic for two vectors
%
% FORMAT:       t = kendtau(v, m [, d])
%
% Input fields:
%
%       v           vector (Vx1 or 1xV)
%       m           matrix (either VxM or MxV)
%       d           dim of m to use to build pairs, default: first match
%
% Output fields
%
%       t           kendall tau statistic (for each "variable pairing")
%       np          number of pairs
%       nn          nominal number of pairs
%
% Note: see http://www.statsdirect.com/help/nonparametric_methods/kend.htm

% Version:  v1.1
% Build:    16041315
% Date:     Apr-13 2016, 3:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2013, 2016, Jochen Weber
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
if nargin < 2 || ~isnumeric(v) || ~isnumeric(m) || isempty(v) || ...
    max(size(v)) ~= numel(v) || ~any(size(m) == numel(v))
    error('neuroelf:general:badArgument', 'Invalid or missing argument.');
end
n = numel(v);

% check for ties in v
sp = kendtaupairsign(v(:), 1);
tv = sum(sp == 0);

% for a simple case
nn = 0.5 * n * (n - 1);
if n == numel(m)

    % get sign of pairs
    sm = kendtaupairsign(m(:), 1);
    tm = sum(sm == 0);

    % compute simple statistic
    if tv == 0 && tm == 0
        np = nn;
        t = sum(2 * (sp == sm) - 1) / np;

    % or yet a bit more complicated
    else
        ui = (sp ~= 0 & sm ~= 0);
        np = sum(ui);
        t = (sum(sp(ui) == sm(ui)) - sum(sp(ui) ~= sm(ui))) / np;
        if isinf(t) || isnan(t)
            t = 0;
        end
        if abs(t) > 1
            t = sign(t);
        end
    end

% otherwise
else
    if nargin < 3 || ~isa(d, 'double') || numel(d) ~= 1 || ...
        isinf(d) || isnan(d) || d < 1 || d > ndims(m)
        d = findfirst(size(m) == n);
    else
        d = round(real(d));
    end
    if size(m, d) ~= n
        error('neuroelf:general:badArgument', 'Invalid dim argument.');
    end

    % resolve transio if necessary
    if istransio(m)
        m = resolve(m);
    end

    % get target size
    s = size(m);
    ud = setdiff(1:numel(s), d);
    ts = s(ud);
    if numel(ts) < 2
        ts(2) = 1;
    end

    % make sure m is in "good shape"
    m = reshape(permute(double(m), [d, ud]), n, prod(ts));

    % build index of good pairs in vector
    uip = sp ~= 0;

    % make useful estimate of how many items in m to take at the same time
    ml = max(1, floor(1e7 / numel(sp)));

    % prepare output
    t = zeros(ts);
    np = zeros(ts);
    tt = numel(t);

    % iterate over m
    mf = 1;
    to = ones(1, ml);
    while mf <= tt

        % get upper limit of next block
        mt = min(mf + ml - 1, tt);

        % compute signs of equal pairs
        sm = kendtaupairsign(m(:, mf:mt), 1);

        % get ties
        if (mt + 1 - mf) ~= ml
            to = ones(1, mt + 1 - mf);
        end
        ui = (uip(:, to) & sm ~= 0);

        % compute statistic
        npp = sum(ui);
        tp = ((sum((sp(:, to) == sm) & ui) - sum((sp(:, to) ~= sm) & ui))) ./ npp;

        % reject errors
        tp(isinf(tp) | isnan(tp)) = 0;
        tw = (abs(tp) > 1);
        if any(tw)
            tp(tw) = sign(tp(tw));
        end

        % store in output
        t(mf:mt) = tp;
        np(mf:mt) = npp;

        % increase lower bound
        mf = mt + 1;
    end
end

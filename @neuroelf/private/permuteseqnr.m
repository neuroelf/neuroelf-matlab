function s = permuteseqnr(s, n)
%PERMUTESEQNR  Permute a sequence without repeats.
%   S = PERMUTESEQNR(S) permutes the sequence SINPUT into a random order,
%   not allowing repeats. If one of the elements with a frequency of at
%   least 0.5, an error will be thrown.
%
%   S = PERMUTESEQNR(S, N) will generate N separate sequences, such that
%   the output S is a cell array with separate sequences

% Version:  v1.1
% Build:    16053121
% Date:     May-31 2016, 9:41 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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
if nargin < 1 || isempty(s)
    error('neuroelf:general:missingArgument', 'Missing argument: s.');
end
if ~isnumeric(s) && ~ischar(s)
    error('neuroelf:general:badArgument', 'Sequence must be numeric or char.');
end

% keep track of input class
cs = class(s);
if strcmpi(cs, 'double')
    cs = [];
else
    eval(['cs=@' cs ';']);
end

% then make sure no invalid values in input
s = double(s(:));
if any(s ~= fix(s) | isinf(s))
    error('neuroelf:general:badArgument', 'Data must be integer.');
end
ns = numel(s);

% how many sequences
if nargin < 2 || ~isa(n, 'double') || numel(n) ~= 1 || isinf(n) || isnan(n) || n < 1
    n = 1;
end

% all elements within range
e = unique(s);
mins = e(1);
maxs = e(end);
ne = numel(e);
if mins ~= 1 || ne ~= maxs
    so = zeros(1, maxs + max(0, 1 - mins));
    if mins > 0
        so(e) = 1:ne;
        s = so(s);
    else
        so(e - (mins - 1)) = 1:ne;
        s = so(s - (mins - 1));
    end
else
    so = [];
end

% compute frequencies
ft = histc(s, 1:ne);
f = ft ./ ns;

% invalid?
if any(f >= 0.5)
    error('neuroelf:general:badArgument', 'Sequence has elements occurring with more than .5 frequency.');
end

% generate corrected sampling frequencies (for non-repeats)
cf = f ./ (1 - f);
cf = cf ./ sum(cf);
cfo = f;
while any(abs(cfo - cf) > sqrt(eps))
    cfo = cf;
    cf = f ./ (1 - cf);
    cf = cf ./ sum(cf);
end
ff = cf ./ f;

% sequences
s = cell(n, 1);
for sc = 1:n

    % copy tally
    t = ft;

    % generate sequence
    ss = zeros(ns, 1);

    % begin with initial draw
    r = rand(ns, 1);
    ss(1) = 1 + sum(r(1) > cumsum(cf));
    t(ss(1)) = t(ss(1)) - 1;
    nd = ns - 1;

    % keep drawing
    td = 1;
    while td < ns
        f = ff .* (t ./ nd);
        f(ss(td)) = 0;
        f = f ./ sum(f);
        td = td + 1;
        ss(td) = 1 + sum(r(td) > cumsum(f));
        t(ss(td)) = t(ss(td)) - 1;
    end
    
    % re-value
    if ~isempty(so)
        ss = e(ss);
    end

    % change class?
    if ~isempty(cs)
        ss = cs(ss);
    end

    % store
    s{sc} = ss;
end

% single sequence
if n == 1
    s = s{1};
end

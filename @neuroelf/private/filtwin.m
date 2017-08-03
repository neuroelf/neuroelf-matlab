function w = filtwin(n, t, p)
% filtwin  - compute filtering window function
%
% FORMAT:       w = filtwin(n [, t [, p]])
%
% Input fields:
%
%       n           window length (N > 2)
%       t           type, either one of
%                   'gauss', {'hamming'}, 'hann', 'rect', 'tri', 'tukey'
%                   or a set of parameters (1xP double) for equation
%                   w(x[0:n-1]) = t(1) - t(2) * cos(2*pi*x/(n-1)) + ...
%                   being a higher-order generalized cosine window
%                   for the gauss window, an alpha parameter of 0.4 is used
%                   for the hamming window, a default of [0.53836] is used
%                   for the tukey window, a default of 0.5 is used
%       p           additional parameter for named function (see above)
%
% Output fields:
%
%       w           window function 
%
% Note: see http://en.wikipedia.org/wiki/Window_function for details

% Version:  v1.0
% Build:    15112720
% Date:     Nov-27 2015, 8:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2015, Jochen Weber
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
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n < 3
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument n.' ...
    );
end
n = ceil(n);
if nargin < 3 || ...
   ~isa(p, 'double') || ...
    isempty(p) || ...
    any(isinf(p(:)) | isnan(p(:)))
    p = [];
end
if nargin < 2 || ...
    isempty(t) || ...
   ((~ischar(t) || ...
     ~any(strcmpi(t(:)', {'blackman', 'gauss', 'hamm', 'hamming', 'hann', ...
         'nutall', 'rect', 'rectangle', 'tri', 'triangle', 'tukey', 'welch'}))) && ...
    (~isa(t, 'double') || ...
      numel(t) < 2 || ...
      length(t) ~= numel(t) || ...
      any(isinf(t) | isnan(t) | abs(t) > 4)))
    t = [0.53836, 0.46164];
elseif ischar(t)
    switch (lower(t(:)'))
        case {'blackman'}
            t = [7938, 9240, 1430] ./ 18608;
        case {'gauss'}
            if isempty(p)
                p = 0.4;
            end
            n = n - 1;
            nc = (0:n)';
            w = exp(-0.5 .* (((2 / (p(1) * n)) .* (nc - (0.5 * n))) .^ 2));
            return;
        case {'hamm', 'hamming'}
            t = [0.53836, 0.46164];
        case {'hann'}
            t = [0.5, 0.5];
        case {'nutall'}
            t = [0.3635819, 0.4891775, 0.1365995, 0.0106411];
        case {'rect', 'rectangle'}
            w = ones(n, 1);
            return;
        case {'tri', 'triangle'}
            n = 0.5 * n;
            if n == fix(n)
                w = (1/n):(1/n):1;
                w = [w'; w(1, end:-1:1)'];
            else
                w = (0.5/n):(1/n):((n-0.25)/n);
                w = [w'; 1; w(1, end:-1:1)'];
            end
            return;
        case {'tukey'}
            if isempty(p)
                p = 0.5;
            end
            rectlen = round(n * (1 - p / 2));
            sinclen = round(n * p / 2);
            sfunc = sin((pi / sinclen) * (0:sinclen)');
            w = convones(sfunc, rectlen);
            w = (1 / max(w)) .* w;
            return;
        case {'welch'}
            w = 1 - (((2 / (n - 1)) .* ((0:(n-1))' - (n - 1) / 2)) .^ 2);
            return;
    end
end

% numeric window
w = zeros(n, 1);
nc = (2 * pi / (n - 1)) .* (0:(n-1))';
for c = 1:numel(t)
    w = w + ((2 * mod(c, 2) - 1) .* t(c)) .* cos((c - 1) .* nc);
end

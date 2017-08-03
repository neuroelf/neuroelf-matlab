function z = fisherr2z(r, d, n)
% fisherr2z  - convert Pearson correlation into z values
%
% FORMAT:       z = fisherr2z(r [, inverse [, n]])
%
% Input fields:
%
%       r           correlation value(s) (or z for inverse)
%       inverse     if given and 1x1 logical true, inverse operation
%       n           number of observations in correlation
%
% Output fields:
%
%       z           z values (or r for inverse)
%
% Note: if n given and > 3, z-scores will be (approximately) normally
%       distributed with SD:=1

% Version:  v0.9a
% Build:    11052014
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
   (~isa(r, 'double') && ...
    ~isa(r, 'single')) || ...
    any(isinf(r(:)) | isnan(r(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing r argument.' ...
    );
end
if ~isa(r, 'double')
    r = double(r);
end

% correct for SD
if nargin > 2 && ...
    isa(n, 'double') && ...
    numel(n) == 1 && ...
   ~isinf(n) && ...
   ~isnan(n) && ...
    n > 3
    sd = (1 ./ (n - 3)) .^ .5;
else
    sd = 1;
end

% compute the desired direction
if nargin < 2 || ...
   ~islogical(d) || ...
   ~d(1)
    z = 0.5 * log((1 + r) ./ (1 - r));
    if ~isreal(z)
        znr = (imag(z) ~= 0);
        z(znr) = Inf * sign(real(z(znr)));
    end
    if sd ~= 1
        z = (1 ./ sd) * z;
    end
else
    if sd ~= 1
        r = sd .* r;
    end
    z = exp(1) .^ (2 .* r);
    z = (z - 1) ./ (z + 1);
end

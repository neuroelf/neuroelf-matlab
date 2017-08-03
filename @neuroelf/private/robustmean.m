function [m, w] = robustmean(x, d)
% robustmean  - robust mean estimation
%
% FORMAT:       [m, w] = robustmean(x [, dim])
%
% Input fields:
%
%       x           data
%       dim         dim (default: 1st non singleton)
%
% Output fields:
%
%       m           robust mean estimate along dim
%       w           weights

% Version:  v0.9b
% Build:    13110616
% Date:     Jun-15 2010, 12:37 PM EST
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
   ~isa(d, 'double') || ...
    numel(d) ~= 1 || ...
    isinf(d) || ...
    isnan(d) || ...
    d < 1 || ...
    d > ndims(x) || ...
    d ~= fix(d)
    d = findfirst(size(x) > 1);
end
sd = size(x, d);
if nargout == 1
    if numel(x) == sd
        m = fitrobustbisquare(ones(sd, 1), x(:));
    else
        m = fitrobustbisquare_img(ones(sd, 1), x);
    end
else
    if numel(x) == sd
        [m, r, w] = fitrobustbisquare(ones(sd, 1), x(:));
    else
        [m, w] = fitrobustbisquare_img(ones(sd, 1), x);
    end
end

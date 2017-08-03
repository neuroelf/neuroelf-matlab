function lta = ltriasc(m, p)
% ltriasc  - ASCII form of a lower triangle matrix (for Mx)
%
% FORMAT:       lta = ltriasc(m [, p])
%
% Input fields:
%
%       m           double matrix
%       p           optional precision (default: 8)
%
% Output fields:
%
%       lta         ASCII form of lower triangle matrix

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
   ~isa(m, 'double') || ...
   ~isreal(m) || ...
    isempty(m) || ...
    ndims(m) ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid matrix argument supplied.' ...
    );
end
if nargin < 2 || ...
   ~isa(p, 'double') || ...
   ~isreal(p) || ...
    numel(p) ~= 1 || ...
    isinf(p) || ...
    isnan(p) || ...
    p < 1 || ...
    p > 15
    pstr = ' %11.6f';
else
    pstr = sprintf(' %%%d.%df', ceil(p + 5), floor(p));
end

% initialize cell array
rc = size(m, 1);
cc = size(m, 2);
lta = cell(1, rc);

% iterate over rows
for rc = 1:rc

    % build row
    lta{rc} = sprintf(pstr, m(rc, 1:min(cc, rc)));
end

% build output
lta = gluetostring(lta, char(10));

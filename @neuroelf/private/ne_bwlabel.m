function [l, num] = ne_bwlabel(bw, n)
% ne_bwlabel  - NeuroElf implementation of bwlabel
%
% FORMAT:       [l, num] = ne_bwlabel(bw [, n])
%
% Input fields:
%
%       bw          binary 2D matrix (image) to label (into clusters)
%       n           connectivity scheme, either 4 or {8}
%
% Output fields:
%
%       l           labels (double matrix of size of input matrix)
%       num         number of labels
%
% Note: this functions passes the input on to clustercoordsc

% Version:  v0.9c
% Build:    11100513
% Date:     Oct-05 2011, 1:16 PM EST
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
   ~isnumeric(bw) || ...
    ndims(bw) > 2 || ...
    isempty(bw) || ...
    any(size(bw) < 3)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if ~islogical(bw)
    bw = (bw >= 0);
end
if nargin < 2 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
   ~any([4, 8] == n)
    n = 5;
elseif n == 8
    n = 5;
end

% pass on to clustercoordsc
[cs, l] = clustercoordsc(bw, n);

% second output
num = numel(cs);

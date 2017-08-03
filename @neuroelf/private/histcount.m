function hc = histcount(v, from, to, step, v2, from2, to2, step2)
% histcount  - histogram across equidistant bins (1D/2D)
%
% FORMAT:       hc = histcount(v, from, to, [step, [v2, from2, to2, step2]]);
%
% Input fields:
%
%       v, v2       value array
%       from, from2 range definition (begin)
%       to, to2     range definition (end)
%       step, step2 range definition (stepsize, default: 1)
%
% Output fields:
%
%       hc          histogram count
%
% Note: if v2 is set to a 1x1 double, it gives the dimension along which
%       v is histogramed
%
% Note: if v2 is a double array of the size of v, and from2 is empty
%       then the histogram is weighed by the information in v2; to combine
%       dimension and weighting, the order must be dimension first
%
% Example:
%
%    % create random numbers (5 columns with 40 values each)
%    r = randn(40, 5);
%
%    % create random weights
%    w = rand(size(r));
%
%    % compute weighted histogram along column dimension
%    wh = histcount(r, -3, 3, 0.1, 1, w);
%
%    % smooth histogram
%    swh = flexinterpn(wh, [Inf, Inf; 1, 1; 1, 1; size(wh)], ...
%          {smoothkern(10), [0; 1; 0]}, {1, 1});

% Version:  v0.9d
% Build:    14052616
% Date:     May-26 2014, 4:49 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010 - 2014, Jochen Weber
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

% simple error
error( ...
    'neuroelf:MissingMEXFile', ...
    'This is a compiled function, but the MEX file is missing.' ...
);

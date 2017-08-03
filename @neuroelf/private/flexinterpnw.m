function [idata, iw] = flexinterpnw(d, dw, coords, varargin)
% flexinterpnw  - flexible data interpolation (up to 4D) with weights
%
% FORMAT:       [id, iw] = flexinterpnw(d, dw, coords [, k, ks [, hold [, qt]]])
%
% Input fields:
%
%       d           up to 4D data
%       dw          data weights
%       coords      CxD (number of coords -by- number of dims) coordinates
%               or  4xD range specification with
%                   1st row: Inf
%                   2nd row: from
%                   3rd row: step
%                   4th row: to, such that
%                   a given dim's range is sample at from:step:to
%       kernel      Kx1 interpolation kernel (e.g. from sinc(x)/sinc(x/a))
%       ksampling   kernel sampling rate (in samples)
%       hold        default sample at invalid coordinates
%       qtrf        quaternion transformation applied to a given range
%
% Output fields:
%
%       id          interpolated data
%       iw          interpolated weights
%
% Note: This function uses flexinterpn, which is a compiled function!

% Version:  v0.9c
% Build:    11120715
% Date:     Dec-07 2011, 3:25 PM EST
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

% basic argument check
if nargin < 3 || ...
    nargout ~= 2 || ...
   ~isnumeric(d) || ...
   ~isa(dw, 'double') || ...
   ~isequal(size(d), size(dw)) || ...
   ~isa(coords, 'double')
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% try to use flexinterpn
try
    idata = double(d);
    idata = idata .* dw;
    idata = flexinterpn(idata, coords, varargin{:});
    iw = flexinterpn(dw, coords, varargin{:});
    idata = idata ./ iw;
    idata(iw == 0) = 0;
catch ne_eo;
    rethrow(ne_eo);
end

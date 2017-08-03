function idata = flexinterpn(data, coords, kernel, ksampling, hold, qtrf)
% flexinterpn  - flexible data interpolation (up to 4D)
%
% FORMAT:       idata = flexinterpn(data, coords [, k, ks [, hold [, qt]]])
%
% Input fields:
%
%       data        up to 4D data
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
%       idata       interpolated data
%
% Example:
%
%       data = randn(16, 16, 16);
%       kernel = sinc(-3:0.001:3) .* sinc(-1:1/3000:1);
%       kernel([1,end]) = 0;
%       meshing = 0.625:0.25:16.375;
%       [x, y, z] = meshgrid(meshing, meshing, meshing);
%       idata = reshape( ...
%           flexinterpn(data, [x(:), y(:), z(:)], kernel(:), 1000), ...
%           [64, 64, 64]);
%
% Observe that the sample result is obtained by using...
%
%       idata = flexinterpn(data, repmat([Inf; 0.625; 0.25; 16.375], 1, 3), ...
%           kernel(:), 1000);
%
% Note: This is a compiled function.
%       Also, make sure the kernel has a 0 value at the beginning and end
%       e.g. by issuing k([1, end]) = 0;
%       Lastly, if different kernels are to be used (only valid for
%       coordinate ranges!), both the kernel and ksampling arguments can
%       be given as a 1-by-ndims(data) cell array

% Version:  v0.9b
% Build:    10062206
% Date:     Jun-21 2010, 12:11 AM EST
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

% bail out
error( ...
    'neuroelf:MEXFileMissing', ...
    'This is a MEX compiled function, but the MEX file is missing.' ...
);

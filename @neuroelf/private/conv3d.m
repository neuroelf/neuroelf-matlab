function v = conv3d(v, k, o, t)
% conv3d  - perform 3D convolution
%
% FORMAT:       v = conv3d(v, k [, o, t])
%
% Input fields:
%
%       v           3D volume (logical or numeric)
%       k           3D convolution kernel (must be double)
%                   -  special 1x1 kernel values:
%                      0 - only apply == threshold (values ~= thresh := 0)
%                      1 - binary erode (full erosion!)
%                      2 - simple binary dilate (only main directions, faces)
%                      3 - complex binary dilate (edge neighbors)
%                      4 - full binary dilate (vertex neighbors)
%                          if input array is not logical, threshold must be
%                          given, by default uses == threshold, add 64 to use
%                          values >= thresh and -64 to use values <= thresh
%                     32 - edge enhancer
%                     48 - gradient magnitude
%                     49 - mean gradient cosine/colinearity in 3x3x3 cube
%                   - special 1x2 kernel values:
%                     [fwhm, func] - build kernel with FWHM and func
%                     func 0: gaussian
%                     func 1: linear
%                   - special 1x3 kernel values:
%                     [ g1 ,  g2 ,  g3 ] - free spatial gradient
%                     [+/-1,  0  ,  0  ] - first dim gradient
%                     [ 0  , +/-1,  0  ] - second dim gradient
%                     [ 0  ,  0  , +/-1] - third dim gradient
%                   - special 1x4 kernel values:
%                     [k1, k2, k3, func] - build kernel with FWHM
%       o           output type (1x1 double, if not given or 0, as input)
%                      1 - logical volume
%       t           optional threshold value (see above)
%
% Output fields:
%
%       v           convolved and/or thresholded/computed volume
%
% Note: for 1x1 kernel arguments 0 to 4 , a logical array with the
%       requested threshold is produced
%
% Note: this is a MEX compiled function !

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

% bogus code
error( ...
    'neuroelf:MEXMissing', ...
    'This is a compiled function, but the MEX file is missing.' ...
);

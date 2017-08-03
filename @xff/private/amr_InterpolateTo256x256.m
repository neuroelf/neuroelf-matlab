function xo = amr_InterpolateTo256x256(xo, method)
% AMR::InterpolateTo256x256  - interpolate the AMR to a 256x256 matrix
%
% FORMAT:       amr.InterpolateTo256x256([method])
%
% Input fields:
%
%       method      string: 'neareast', {'linear'}, 'cubic', 'lanczos3'
%
% No output fields.

% Version:  v1.1
% Build:    16012715
% Date:     Jan-27 2016, 3:32 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'amr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~ischar(method) || isempty(method) || ...
    ~any(strcmpi(method(:)', {'nearest', 'linear', 'cubic', 'lanczos3'}))
    method = 'linear';
end
method = lower(method(:)');

% check AMR size
bc = xo.C;
asz = size(bc.Slice(1).AMRData);
nsl = numel(bc.Slice);
for sc = 1:nsl
    if numel(size(bc.Slice(sc).AMRData)) ~= 2 || any(asz ~= size(bc.Slice(sc).AMRData))
        error('neuroelf:xff:badObject', 'AMR slices are not equally sized.');
    end
end

% interpolation needed?
if all(asz == 256)
    return;
end

% build grid coordinates
[nx, ny] = meshgrid(1:((asz(1) - 1) / 255):asz(1), 1:((asz(2) - 1) / 255):asz(2));

% interpolate data and new settings
for sc = 1:nsl
    slc = bc.Slice(sc);
    slcd= uint8(interp2(double(slc.AMRData), nx, ny, method));
    slcd(slcd > 225) = 225;
    slc.AMRData = slcd;
    slc.BITMAPFILEHEADER.bfSize = 65536 + slc.BITMAPFILEHEADER.bfOffBits;
    slc.BITMAPINFOHEADER.biWidth = 256;
    slc.BITMAPINFOHEADER.biHeight = 256;
    slc.BITMAPINFOHEADER.biImageSize = 65536;
    bc.Slice(sc) = slc;
end

% set back to object and clear filename
xo.C = bc;
xo.F = '';

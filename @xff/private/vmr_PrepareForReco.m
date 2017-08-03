function xo = vmr_PrepareForReco(xo, sval, bval)
% VMR::PrepareForReco  - prepare border for reconstruction
%
% FORMAT:       [vmr] = vmr.PrepareForReco([sval, bval])
%
% Input fields:
%
%       sval        an optional segmentation value (default: 240)
%       bval        an optional border value (default: 235)
%
% Output fields:
%
%       vmr         prepared VMR
%
% NOTE: this function does NOT yet work as expected !!
%
% Using: conv3d.

% Version:  v1.1
% Build:    16021316
% Date:     Feb-13 2016, 4:40 PM EST
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

% neuroelf library
global ne_methods;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isa(sval, 'double') || numel(sval) ~= 1 || isnan(sval) || sval < 1 || sval > 255
    sval = 240;
else
    sval = fix(real(sval));
end
if nargin < 3 || ~isa(bval, 'double') || numel(bval) ~= 1 || isnan(bval) || bval < 1 || bval > 255
    bval = 235;
else
    bval = fix(real(bval));
end

% get VMRdata
vmrd = bc.VMRData(:, :, :);

% replace bval with sval
vmrd(vmrd == bval) = sval;

% create boolean matrix with non matching color
vmrb = (vmrd ~= sval);

% find border by notting array with dilation operator
vmrd(~vmrb & ne_methods.conv3d(vmrb, 2)) = bval;

% back in VMR
bc.VMRData = vmrd;
xo.C = bc;

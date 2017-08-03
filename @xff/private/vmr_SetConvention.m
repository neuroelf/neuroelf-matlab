function xo = vmr_SetConvention(xo, cnv)
% VMR::SetConvention  - set radiological/neurological convention
%
% FORMAT:       [vmr = ] vmr.SetConvention(cnv)
%
% Input fields:
%
%       cnv         either 0 / 'n' for neurological
%                   or 1 / 'r' for radiological convention
%
% Output fields:
%
%       vmr         (possibly) altered object with desired convention
%
% Note: the VMR *must* be in sagittal orientation

% Version:  v1.1
% Build:    16021316
% Date:     Feb-13 2016, 4:51 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ...
    numel(cnv) ~= 1 || ((~isa(cnv, 'double') || ~any([0, 1] == cnv)) && ...
    (~ischar(cnv) || ~any('nr' == cnv)))
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;

% what target convention
if ischar(cnv)
    cnv = find('nr' == cnv) - 1;
end

% same as current convention
if cnv == bc.Convention

    % do nothing
    return;
end

% swap data
bc.VMRData = bc.VMRData(:, :, end:-1:1);
if ~isempty(bc.VMRData16) && numel(size(bc.VMRData16)) == numel(size(bc.VMRData)) && ...
    all(size(bc.VMRData16) == size(bc.VMRData))
    bc.VMRData16 = bc.VMRData16(:, :, end:-1:1);
end

% patch offset
bc.OffsetX = bc.FramingCube - (bc.OffsetX + size(bc.VMRData, 1));

% set convention
bc.Convention = cnv;

% add to transformations
t = bc.Trf;
if cnv == 0
    t(end+1).NameOfSpatialTransformation = 'Changing from Radiological to Neurological convention';
else
    t(end+1).NameOfSpatialTransformation = 'Changing from Neurological to Radiological convention';
end
t(end).TypeOfSpatialTransformation = 2;
t(end).SourceFileOfSpatialTransformation = xo.F;
t(end).NrOfSpatialTransformationValues = 16;
t(end).TransformationValues = diag([1, 1, -1, 1]);
bc.Trf = t;
bc.NrOfPastSpatialTransformations = numel(t);
xo.C = bc;

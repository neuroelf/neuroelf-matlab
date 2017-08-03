function xo = voi_ToISO(xo, raw, iso, iflag)
% VOI::ToISO  - apply transformation to (and form) ISO-space
%
% FORMAT:       [voi = ] voi.ToISO(raw, iso [, iflag]))
%
% Input fields:
%
%       raw         VMR object of original (scanner-space) VMR
%       iso         ISO-(1mm-)VMR object
%       iflag       if given and true, inverse transformation
%
% Output fields:
%
%       voi         VOI with transformed coordinates

% Version:  v1.1
% Build:    16021016
% Date:     Feb-10 2016, 4:22 PM EST
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
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'voi') || ...
    numel(raw) ~= 1 || ~xffisobject(raw, true, 'vmr') || ...
    numel(iso) ~= 1 || ~xffisobject(iso, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 4 || ~islogical(iflag) || numel(iflag) ~= 1
    iflag = false;
end

% get contents of objects
vc = xo.C;
rc = raw.C;
ic = iso.C;

% check whether ic is ISO VMR content
if isempty(ic.Trf) || numel(ic.Trf) > 2 || ...
    isempty(regexpi(ic.Trf(1).NameOfSpatialTransformation, '^isovoxel')) || ...
   (numel(ic.Trf) > 1 && isempty(regexpi(ic.Trf(2).NameOfSpatialTransformation, 'sagittal')))
    error('neuroelf:xff:badArgument', 'Invalid ISO-Voxel VMR object.');
end

% compute offset
of = 0.5 .* ([rc.DimX * rc.VoxResX, rc.DimY * rc.VoxResY, rc.DimZ * rc.VoxResZ] - 256);

% and get transformation
trf = {ic.Trf.TransformationValues};
if numel(trf) > 1
    trf = trf{1} * trf{2};
else
    trf = trf{1};
end

% inverse
if iflag

% forward
else

    % inverse sampling transformation
    trf = inv(trf)';

    % compute difference for two coordinates
    cdiff = abs((trf * [1;1;1;1]) - (trf * [2;2;2;1]));

    % create oversampling steps
    ost = ceil(2 * cdiff(1:3));
    osx = (-1 / 3):(2 / (ost(1) * 3)):(1 / 3);
    osy = (-1 / 3):(2 / (ost(2) * 3)):(1 / 3);
    osz = (-1 / 3):(2 / (ost(3) * 3)):(1 / 3);
    osn = numel(osx) * numel(osy) * numel(osz);

    % iterate over VOIs
    for voic = 1:numel(vc.VOI)

        % get coordinates
        crd = vc.VOI(voic).Voxels;
        ncrd = size(crd, 1);
        ocrd = ones(ncrd, 1);

        % create target coordinates
        tcrd = zeros(osn * ncrd, 4);

        % for oversampling steps
        tc = 1;
        for xc = osx
            for yc = osy
                for zc = osz
                    tcrd(tc:tc+ncrd-1, :) = ...
                        [crd(:, 1) + xc, crd(:, 2) + yc, crd(:, 3) + zc, ocrd] * trf;
                    tc = tc + ncrd;
                end
            end
        end

        % extract coordinate and add offset
        tcrd = unique(round(tcrd(:, 1:3)), 'rows');
        tcrd = [tcrd(:, 1) - of(1), tcrd(:, 2) - of(2), tcrd(:, 3) - of(3)];

        % store in VOI
        vc.VOI(voic).Voxels = tcrd;
        vc.VOI(voic).NrOfVoxels = size(tcrd, 1);
    end

    % patch header
    vc.OriginalVMRResolutionX = ic.VoxResX;
    vc.OriginalVMRResolutionY = ic.VoxResY;
    vc.OriginalVMRResolutionZ = ic.VoxResZ;
    vc.OriginalVMROffsetX = ic.OffsetX;
    vc.OriginalVMROffsetY = ic.OffsetY;
    vc.OriginalVMROffsetZ = ic.OffsetZ;
    vc.OriginalVMRFramingCubeDim = ic.FramingCube;
    vc.Convention = ic.Convention;
end

% set back
xo.C = vc;

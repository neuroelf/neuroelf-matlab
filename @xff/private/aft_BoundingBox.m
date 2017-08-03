function bbox = aft_BoundingBox(xo)
% AFT::BoundingBox  - get bounding box
%
% FORMAT:       bbox = obj.BoundingBox;
%
% No input fields
%
% Output fields:
%
%       bbox        struct with fields
%        .BBox      2x3 offset and offset + size - 1
%        .FCube     framing cube
%        .DimXYZ    data dimensions
%        .ResXYZ    data resolution
%        .QuatB2T   BV2Tal quaternion
%        .QuatT2B   Tal2BV quaternion
%
% TYPES: AVA, CMP, DDT, GLM, HDR, HEAD, MGH, MSK, NLF, SRF, TVL, VDW, VMP, VMR, VTC
%
% Note: output is in BV's *internal* notation (axes not in TAL order!)

% Version:  v1.1
% Build:    17071010
% Date:     Jul-10 2017, 10:54 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, 2017, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% switch on filetype
bc = xo.C;
ft = lower(xo.S.Extensions{1});
try
    switch (ft)
        case {'ava', 'cmp', 'ddt', 'glm', 'msk', 'vdw', 'vmp', 'vtc'}
            fc = [256, 256, 256];
            ons = [bc.XStart, bc.YStart, bc.ZStart];
            offs = [bc.XEnd, bc.YEnd, bc.ZEnd];
            resxyz = bc.Resolution * ones(1, 3);
            if bc.Resolution > 1
                bbox = [floor(ons); ceil(offs - 1)];
            else
                if ~strcmp(ft, 'vmp') || bc.NativeResolutionFile == 1
                    bbox = [floor(ons); offs];
                else
                    bbox = [floor(ons); offs + 1];
                end
            end
            dimxyz = ceil(diff(bbox) ./ resxyz);
            rcnv = true;
        case {'hdr', 'head', 'mgh', 'nlf'}
            fc = [256, 256, 256];
            switch (ft)
                case 'hdr'
                    cfr = hdr_CoordinateFrame(xo);
                case 'head'
                    cfr = head_CoordinateFrame(xo);
                case 'mgh'
                    cfr = mgh_CoordinateFrame(xo);
                case 'nlf'
                    cfr = nlf_CoordinateFrame(xo);
            end
            ctr = cfr.Trf;
            dimxyz = cfr.Dimensions(1:3);
            resxyz = cfr.Resolution;
            % regular cube ?
            if abs(max(abs(ctr(1:3, 1))) - resxyz(1)) < 1e-4 && ...
                abs(max(abs(ctr(1:3, 2))) - resxyz(2)) < 1e-4 && ...
                abs(max(abs(ctr(1:3, 3))) - resxyz(3)) < 1e-4
                ofmc = 128 - [ ctr * [0.5 + zeros(3, 1); 1], ctr * [0.5 + dimxyz(:); 1]]';
            else
                ofmc = 128 - [ ...
                    ctr * [    1    ;    1    ;    1    ; 1], ...
                    ctr * [dimxyz(1);    1    ;    1    ; 1], ...
                    ctr * [    1    ;dimxyz(2);    1    ; 1], ...
                    ctr * [dimxyz(1);dimxyz(2);    1    ; 1], ...
                    ctr * [    1    ;    1    ;dimxyz(3); 1], ...
                    ctr * [dimxyz(1);    1    ;dimxyz(3); 1], ...
                    ctr * [    1    ;dimxyz(2);dimxyz(3); 1], ...
                    ctr * [dimxyz(1);dimxyz(2);dimxyz(3); 1]]';
            end
            [xd{1:2}] = max(abs(ctr(1:3, 2)));
            [yd{1:2}] = max(abs(ctr(1:3, 3)));
            [zd{1:2}] = max(abs(ctr(1:3, 1)));
            ao = [xd{2}, yd{2}, zd{2}];
            ons = floor(min(ofmc(:, ao)) + 1e-2);
            offs = ceil(max(ofmc(:, ao)) - 1e-2);
            bbox = [ons; offs - 1];
            rcnv = cfr.IsRadiological;
        case 'srf'
            fc = 2 * bc.MeshCenter;
            bbox = [[0, 0, 0]; fc - 1];
            dimxyz = fc;
            rcnv = true;
            resxyz = [1, 1, 1];
        case 'vmr'
            dimxyz = size(bc.VMRData);
            rcnv = bc.Convention;
            resxyz = [bc.VoxResX, bc.VoxResY, bc.VoxResZ];
            if all(resxyz == resxyz(1))
                fc = resxyz .* bc.FramingCube(ones(1, 3));
            else
                fc = [256, 256, 256];
            end
            offs = resxyz.* [bc.OffsetX, bc.OffsetY, bc.OffsetZ];
            bbox = [offs; offs + resxyz .* (dimxyz - 1)];
        otherwise
            error('BAD_FILETYPE');
    end
catch xfferror
    error('neuroelf:xff:internalError', ...
        'BoundingBox is not a valid method of type %s: %s.', upper(ft), xfferror.message ...
    );
end

% generate quaternions for TAL2BVC and inverse
fch = 0.5 * fc;
denxyz = 1 ./ resxyz;
b2t = [0, 0, -resxyz(1), fch(3) + resxyz(3) - bbox(1, 3); ...
    -resxyz(2), 0, 0, fch(1) + resxyz(1) - bbox(1, 1); ...
    0, -resxyz(3), 0, fch(2) + resxyz(2) - bbox(1, 2); 0, 0, 0, 1];
t2b = [0, -denxyz(2), 0, denxyz(2) * b2t(2, 4); 0, 0, -denxyz(3), denxyz(3) * b2t(3, 4); ...
    -denxyz(1), 0, 0, denxyz(1) * b2t(1, 4); 0, 0, 0, 1];

% generate output struct
bbox = struct('BBox', bbox, 'DimXYZ', dimxyz, 'FCube', fc, 'RadCnv', rcnv, ...
    'ResXYZ', resxyz, 'QuatB2T', b2t, 'QuatT2B', t2b);

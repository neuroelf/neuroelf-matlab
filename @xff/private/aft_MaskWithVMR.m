function xo = aft_MaskWithVMR(xo, maskvmr, thold)
% AFT::MaskWithVMR  - zero values where VMR value beyond threshold
%
% FORMAT:       [obj = ] obj.MaskWithVMR(maskvmr [, threshold])
%
% Input fields:
%
%       maskvmr     masking VMR object
%       threshold   threshold value(s), default: 11 (no upper boundary)
%
% Output fields:
%
%       obj         masked object
%
% TYPES: CMP, DDT, GLM, VDW, VMP, VMR, VTC
%
% Using: bvcoordconv.

% Version:  v1.1
% Build:    16020214
% Date:     Feb-02 2016, 2:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true) || numel(maskvmr) ~= 1 || ...
   ~xffisobject(maskvmr, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
ext = lower(xo.S.Extensions{1});
if ~any(strcmp(ext, {'amr', 'ava', 'cmp', 'ddt', 'dmr', 'fmr', 'glm', 'head', ...
    'map', 'msk', 'vdw', 'vmp', 'vmr', 'vtc'}))
    error('neuroelf:xff:badArgument', 'Invalid object for MaskWithVMR method.');
end
bc1 = xo.C;
if nargin < 3 || ~isnumeric(thold) || ~any(1:2 == numel(thold)) || ...
    any(isinf(thold) | isnan(thold) | thold < 0 | thold > 32767)
    thold = 11;
    if isfield(maskvmr.C.RunTimeVars, 'ScalingWindow')
        thold = max(thold, maskvmr.C.RunTimeVars.ScalingWindow(1));
    end
end
thold = fix(thold);

% pass on errors
try
    % get bounding box
    bbox = aft_BoundingBox(xo);
    vbox = aft_BoundingBox(maskvmr);

    % get final transformation matrix
    trf = ne_methods.bvcoordconv(zeros(0, 3), 'tal2bvc', vbox) * ...
          ne_methods.bvcoordconv(zeros(0, 3), 'bvc2tal', bbox);

    % create argument for sampling
    range = [Inf, Inf, Inf; ones(2, 3); size(aft_GetVolume(xo, 1))];

    % get masking data
    mskd = aft_SampleData3D(maskvmr, range, ...
        struct('method', 'nearest', 'space', 'bvc', 'trans', trf));

    % threshold data
    if numel(thold) == 1
        mskd = (mskd < thold);
    else
        mskd = ((mskd < thold(1)) | (mskd > thold(2)));
    end
catch xfferror
    rethrow(xfferror);
end

% depending on filetype
rsetc = true;
switch (ext)
    case 'cmp'
        for mc = 1:numel(bc1.Map)
            bc1.Map(mc).CMPData(mskd) = 0;
            if ~isempty(bc1.Map(mc).CMPDataCT)
                bc1.Map(mc).CMPDataCT(mskd) = 0;
            end
        end
    case 'ddt'
        bc1.TensorEigenVs(repmat(shiftdim(mskd, -1), 12, 1)) = 0;
    case 'glm'
        if bc1.ProjectType ~= 1
            warning('neuroelf:xff:badArgument', 'Only valid for VTC-based GLMs.');
            return;
        end
        if bc1.ProjectTypeRFX <= 0
            np = size(bc1.GLMData.BetaMaps, 4);
            bc1.GLMData.MultipleRegressionR(mskd) = 0;
            bc1.GLMData.MCorrSS(mskd) = 0;
            bc1.GLMData.BetaMaps(repmat(mskd, [1, 1, 1, np])) = 0;
            bc1.GLMData.XY(repmat(mskd, [1, 1, 1, np])) = 0;
            bc1.GLMData.TimeCourseMean(mskd) = 0;
        else
            np = size(bc1.GLMData.Subject(1).BetaMaps, 4);
            for sc = 1:numel(bc1.GLMData.Subject)
                bc1.GLMData.Subject(sc).BetaMaps(repmat(mskd, [1, 1, 1, np])) = 0;
            end
        end
    case 'vdw'
        nvol = size(bc1.VDWData, 1);
        bc1.VDWData(repmat(shiftdim(mskd, -1), nvol, 1)) = 0;
    case 'vmp'
        for mc = 1:numel(bc1.Map)
            bc1.Map(mc).VMPData(mskd) = 0;
            if ~isempty(bc1.Map(mc).VMPDataCT)
                bc1.Map(mc).VMPDataCT(mskd) = 0;
            end
        end
    case 'vmr'
        if ~isempty(bc1.VMRData)
            bc1.VMRData(mskd) = 0;
        else
            bc1.VMRData16(mskd) = 0;
        end
    case 'vtc'
        nvol = size(bc1.VTCData, 1);
        bc1.VTCData(repmat(shiftdim(mskd, -1), nvol, 1)) = 0;
    otherwise
        rsetc = false;
        warning('neuroelf:xff:notYetImplemented', ...
            'Function not yet implemented for this type of object.');
end

% reset content
if rsetc
    xo.C = bc1;
end

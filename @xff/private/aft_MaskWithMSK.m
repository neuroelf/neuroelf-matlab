function xo = aft_MaskWithMSK(xo, mask)
% AFT::MaskWithMSK  - zero values where MSK value is 0
%
% FORMAT:       [obj = ] obj.MaskWithMSK(mask)
%
% Input fields:
%
%       mask        masking MSK object
%
% Output fields:
%
%       obj         masked object
%
% TYPES: AVA, CMP, DDT, GLM, MSK, VDW, VMP, VTC

% Version:  v1.1
% Build:    16020214
% Date:     Feb-02 2016, 2:30 PM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true) || numel(mask) ~= 1 || ...
   ~xffisobject(mask, true, 'msk')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
ext = lower(xo.S.Extensions{1});
if ~any(strcmp(ext, {'ava', 'cmp', 'ddt', 'glm', 'msk', 'vdw', 'vmp', 'vtc'}))
    error('neuroelf:xff:badArgument', 'Invalid object for MaskWithMSK method.');
end
bc = xo.C;

% pass on errors
try
    % get bounding box
    bbox = aft_BoundingBox(xo);
    vbox = aft_BoundingBox(mask);

    % check settings
    if ~isequal(bbox.BBox, vbox.BBox) || ~isequal(bbox.ResXYZ, vbox.ResXYZ)
        error('neuroelf:xff:badArgument', 'Bounding boxes mismatch.');
    end
catch xfferror
    rethrow(xfferror);
end

% get mask as logical array
mask = mask.C;
mask = (double(mask.Mask) == 0);

% depending on filetype
switch (ext)
    case 'ava'
        f = fieldnames(bc.Maps);
        for fc = 1:numel(f)
            bc.Maps.(f{fc})(repmat(mask, [1, 1, 1, size(bc.Maps.(f{fc}), 4)])) = 0;
        end
    case 'cmp'
        for mc = 1:numel(bc.Map)
            bc.Map(mc).CMPData(mask) = 0;
            bc.Map(mc).CMPDataCT = [];
        end
    case 'ddt'
        bc.TensorEigenVs(repmat(shiftdim(mask, -1), 12, 1)) = 0;
    case 'glm'
        if bc.ProjectType ~= 1
            warning('neuroelf:xff:badArgument', 'Only valid for VTC-based GLMs.');
            return;
        end
        if bc.ProjectTypeRFX <= 0
            np = size(bc.GLMData.BetaMaps, 4);
            bc.GLMData.MultipleRegressionR(mask) = 0;
            bc.GLMData.MCorrSS(mask) = 0;
            bc.GLMData.BetaMaps(repmat(mask, [1, 1, 1, np])) = 0;
            bc.GLMData.XY(repmat(mask, [1, 1, 1, np])) = 0;
            bc.GLMData.TimeCourseMean(mask) = 0;
        else
            np = size(bc.GLMData.Subject(1).BetaMaps, 4);
            for sc = 1:numel(bc.GLMData.Subject)
                bc.GLMData.Subject(sc).BetaMaps(repmat(mask, [1, 1, 1, np])) = 0;
            end
        end
    case 'msk'
        bc.Mask(mask) = 0;
    case 'vdw'
        nvol = size(bc.VDWData, 1);
        bc.VDWData(repmat(shiftdim(mask, -1), nvol, 1)) = 0;
    case 'vmp'
        for mc = 1:numel(bc.Map)
            bc.Map(mc).VMPData(mask) = 0;
            bc.Map(mc).VMPDataCT = [];
            if isfield(bc.Map(mc), 'RunTimeVars') && isfield(bc.Map(mc).RunTimeVars) && ...
                isfield(bc.Map(mc).RunTimeVars, 'FWHMResImg') && ...
                isequal(size(bc.Map(mc).RunTimeVars.FWHMResImg), size(bc.Map(mc).VMPData))
                bc.Map(mc).RunTimeVars.FWHMResImg(mask) = 0;
            end
        end
    case 'vtc'
        nvol = size(bc.VTCData, 1);
        bc.VTCData(repmat(shiftdim(mask, -1), nvol, 1)) = 0;
    otherwise
        warning('neuroelf:xff:notYetImplemented', ...
            'Function not yet implemented for this type of object.');
        return;
end

% reset content
xo.C = bc;

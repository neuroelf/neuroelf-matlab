function [xo, c] = dmr_MaskWithVMR(xo, vmr, trf, cc)
% DMR::MaskWithVMR  - mask DMR volumes with a segmentation result
%
% FORMAT:       [dmr, c] = dmr.MaskWithVMR(vmr [, trf, cc])
%
% Input fields:
%
%       vmr         VMR dataset
%       trf         optional 1xT cell array with IA, FA alignment files
%       cc          color codes used for mask (default: [235, 240])
%
% Output fields:
%
%       dmr         altered DMR/DWI (Prefix adapted!)
%       c           optional output: DMR volume coordinate list in mask
%
% Using: samplefmrspace.

% Version:  v1.1
% Build:    16020110
% Date:     Feb-01 2016, 10:56 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'dmr') || ...
    numel(vmr) ~= 1 || ~xffisobject(vmr, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~iscell(trf)
    trf = {};
else
    trf = trf(:)';
end
if nargin < 4 || ~isnumeric(cc) || isempty(cc) || any(isinf(cc(:)) | isnan(cc(:)))
    cc = [235, 240];
else
    cc = min(255, max(0, cc(:)'));
end
cc = uint8(cc);

% get dmr data
bc = xo.C;
if istransio(bc.DWIData)
    bc.DWIData = resolve(bc.DWIData);
end
vc = vmr.C;
vd = vc.VMRData(:, :, :);

% get coordinates
vc = false(size(vd));
for ccc = 1:numel(cc)
    vc(vd == cc(ccc)) = true;
end
[vcc{1:3}] = ind2sub(size(vc), find(vc));
vc = [vcc{3}(:), vcc{1}(:), vcc{2}(:)];

% get coordinates *in mask*
[vcc{1:2}] = ne_methods.samplefmrspace(bc.DWIData, vc, xo, trf, 'linear', true);
vc = round(vcc{2});
vc( vc(:, 1) < 1 | vc(:, 1) > bc.ResolutionX | ...
    vc(:, 2) < 1 | vc(:, 2) > bc.ResolutionY | ...
    vc(:, 3) < 1 | vc(:, 3) > bc.NrOfSlices, :) = [];
vc = unique(vc, 'rows');
c = false([bc.ResolutionX, bc.ResolutionY, bc.NrOfSlices]);
nd = numel(c);
c(sub2ind(size(c), vc(:,1), vc(:,2), vc(:,3))) = 1;
vcc = find(~c);

% mask false values
if bc.DataStorageFormat == 3
    for vlc = 0:(bc.NrOfVolumes - 1)
        bc.DWIData(vcc + nd * vlc) = 0;
    end
else
    vcc = bc.NrOfVolumes * (c - 1);
    for vlc = 1:bc.NrOfVolumes
        bc.DWIData(vcc + vlc) = 0;
    end
end
bc.Prefix = [bc.Prefix '_masked'];

% AMR loaded
amrf = bc.LoadAMRFile;
if ischar(amrf) && ~isempty(amrf) && exist([fileparts(xo.F) '/' amrf], 'file') == 2
    amrf = [fileparts(xo.F) '/' amrf];
    try
        amr = [];
        amr = xff(amrf);
        if ~xffisobject(amr, true, 'amr')
            error('BAD_REFERENCED_AMR');
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
        if xffisobject(amr, true)
            delete(amr);
        end
        warning('neuroelf:xff:internalError', 'Error applying mask to AMR.');
    end
    amrc = amr.C;
    xf = size(amrc.Slice(1).AMRData, 1) / bc.ResolutionX;
    yf = size(amrc.Slice(1).AMRData, 2) / bc.ResolutionY;
    if xf ~= round(xf) || yf ~= round(yf)
        warning('neuroelf:xff:internalError', ...
            'Invalid resolution ratio between DMR/AMR.');
    else
        for slc = 1:bc.NrOfSlices
            amrc.Slice(slc).AMRData(~c(floor(1:(1/xf):(bc.ResolutionX + 0.99)), ...
                ceil(bc.ResolutionY:(-1/yf):0.01), slc)) = 0;
        end
        amr.C = amrc;
        amrf = [bc.LoadAMRFile(1:end-4) '_masked.amr'];
        try
            aft_SaveAs(amr, [fileparts(xo.F) '/' amrf]);
            bc.LoadAMRFile = amrf;
        catch xfferror
            neuroelf_lasterr(xfferror);
            warning('neuroelf:xff:internalError', 'Error saving masked AMR file.');
        end
    end
    delete(amr);
end

% set back
xo.C = bc;

% copy coordinate list
c = vc;

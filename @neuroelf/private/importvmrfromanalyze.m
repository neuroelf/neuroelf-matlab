function vmr = importvmrfromanalyze(img, imeth, hl, res, trf, bbox)
% importvmrfromanalyze  - import a VMR from an Analzye file
%
% FORMAT:       vmr = importvmrfromanalyze(img [, imeth [, hlim, res, trf]])
%
% Input fields:
%
%       img         HDR/IMG/NII filename
%       imeth       interpolation method (optional, see flexinterpn_method)
%       hlim        limit values on histogram (1x2 double, [0.001 0.999])
%       res         resolution (default: 1, optional 0.5)
%       trf         additional TRF to apply (e.g. from coregistration)
%       bbox        bounding box
%
% Output fields:
%
%       vmr         created VMR object
%
% Note: this function requires the MEX file flexinterpn.

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:22 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% argument check
if nargin < 1 || ...
   ~ischar(img) || ...
    isempty(img) || ...
    exist(img(:)', 'file') ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad input argument.' ...
    );
end
if nargin < 2 || ...
   ~ischar(imeth) || ...
    isempty(imeth)
    imeth = 'cubic';
elseif ~any(strcmpi(imeth(:)', {'cubic', 'lanczos3', 'lanczos8', 'linear', 'nearest'}))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid interpolation method.' ...
    );
else
    imeth = lower(imeth(:)');
end
if nargin < 3 || ...
   ~isa(hl, 'double')
    hl = [0.001, 0.999];
elseif numel(hl) ~= 2 || ...
    any(isinf(hl) | isnan(hl) | hl < 0 | hl > 4095)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid histogram limits.' ...
    );
else
    hl = hl(:)';
end
if nargin < 4 || ...
   ~isa(res, 'double') || ...
    numel(res) ~= 1
    res = 1;
elseif isinf(res) || ...
    isnan(res) || ...
   ~any(res == [0.5, 1])
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid resolution specified.' ...
    );
end
if nargin < 5 || ...
   ~isa(trf, 'double') || ...
   ~isequal(size(trf), [4, 4])
    trf = {};
elseif any(isinf(trf(:)) | isnan(trf(:))) || ...
    any(trf(4, :) ~= [0, 0, 0, 1])
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid transformation matrix.' ...
    );
else
    trf = {trf};
end
if nargin < 6 || ...
   ~isa(bbox, 'double') || ...
   ~isequal(size(bbox), [2, 3]) || ...
    any(isinf(bbox(:)) | isnan(bbox(:)) | bbox(:) < 0 | bbox(:) > 256) || ...
    any(diff(bbox) <= 0)
    bbox = [0, 0, 0; 256, 256, 256];
else
    bbox = [floor(bbox(1, :)); ceil(bbox(2, :))];
end
try
    himg = cell(1, 1);
    img = strrep(strrep(img, '.img', '.hdr'), '.IMG', '.HDR');
    himg{1} = xff(img);
    if numel(himg{1}) ~= 1 ||...
       ~isxff(himg{1}, {'hdr', 'head'})
        error('BAD_IMAGECONT');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects(himg);
    error( ...
        'neuroelf:BadArgument', ...
        'Error loading image.' ...
    );
end

% try to sample volume
off = bbox(1, :);
tsc = 0;
try
    if res == 1
        rbbox = bbox;
        vmr16 = single(0);
        vmr16(diff(bbox(:, 1)), diff(bbox(:, 2)), diff(bbox(:, 3))) = 0;
        himg{1}.LoadTransIOData;
        for sc = bbox(1, 3):(bbox(2, 3) - 1)
            rbbox(:, 3) = [sc; sc + res];
            tsc = tsc + 1;
            vmr16(:, :, tsc) = ...
                himg{1}.SampleBVBox(struct('BBox', rbbox, ...
                 'ResXYZ', res), 1, imeth, trf{:});
            drawnow;
        end
        mmm = minmaxmean(vmr16, 4);
        if mmm(1) ~= 0
            vmr16 = vmr16 - mmm(1);
            mmm = mmm - mmm(1);
        end
        vmr16 = floor((4095 / mmm(2)) .* vmr16);
        vmr16(isinf(vmr16) | isnan(vmr16)) = 0;
        vmr16 = uint16(vmr16);
    else
        rbbox = bbox;
        rbsz = 2 .* diff(rbbox);
        off = 2 .* off;
        himg{1}.LoadTransIOData;
        mmm = minmaxmean(himg{1}.VoxelData);
        if mmm(3) < 64 || ...
            mmm(3) > 16384
            mmf = 4096 / (mmm(3) - mmm(1));
        else
            mmf = 1;
        end
        vmr16 = uint16(0);
        vmr16(rbsz(1), rbsz(2), rbsz(3)) = 0;
        for sc = bbox(1, 3):res:(bbox(2, 3) - 0.5 * res)
            rbbox(:, 3) = [sc; sc + res];
            tsc = tsc + 1;
            vmr16(:, :, tsc) = uint16(round(min(65535, mmf .* ...
                himg{1}.SampleBVBox(struct('BBox', rbbox, ...
                 'ResXYZ', res), 1, imeth, trf{:}))));
            drawnow;
        end
    end
catch ne_eo;
    clearxffobjects(himg);
    rethrow(ne_eo);
end
vsz = size(vmr16);

% create object
vmr = xff('new:vmr');

% set object properties
xffroot = xff();
updstate = xffroot.UpdateState('vmr', false);
vmr.VMRData = uint8([]);
vmr.DimX = vsz(1);
vmr.DimY = vsz(2);
vmr.DimZ = vsz(3);
vmr.VMRData(vsz(1), vsz(2), vsz(3)) = 0;
if res == 1
    vmr.FramingCube = 256;
    vmr.NRows = 256;
    vmr.NCols = 256;
else
    vmr.FramingCube = 512;
    vmr.NRows = 512;
    vmr.NCols = 512;
end
vmr.VMRData16 = vmr16;
vmr.PosInfoVerified = 1;
vmr.OffsetX = off(1);
vmr.OffsetY = off(2);
vmr.OffsetZ = off(3);
vmr.Slice1CenterX = -128 + 0.5 * res;
vmr.Slice1CenterY = 0;
vmr.Slice1CenterZ = 0;
vmr.SliceNCenterX = 128 - 0.5 * res;
vmr.SliceNCenterY = 0;
vmr.SliceNCenterZ = 0;
vmr.RowDirX = 0;
vmr.RowDirY = res;
vmr.RowDirZ = 0;
vmr.ColDirX = 0;
vmr.ColDirY = 0;
vmr.ColDirZ = -res;
vmr.FoVRows = 256;
vmr.FoVCols = 256;
vmr.SliceThickness = res;
vmr.GapThickness = 0;
vmr.Convention = 1;
vmr.VoxResX = res;
vmr.VoxResY = res;
vmr.VoxResZ = res;
vmr.VoxResInTalairach = 1;
vmr.VoxResVerified = 1;

% reset updatestate flag
xffroot.UpdateState('vmr', updstate);

% try to perform good limitting
if ~isequal(hl, [0, 1])
    vmr.LimitVMR(struct('recalc8b', true, 'range', hl));
else
    vmr.VMRData = uint8(floor((225.999 / double(max(vmr16(:)))) .* double(vmr16)));
end
vmr.MeanOriginalValue = floor(mean(vmr.VMRData16(:)));
vmr.MaxOriginalValue = double(max(vmr.VMRData16(:)));

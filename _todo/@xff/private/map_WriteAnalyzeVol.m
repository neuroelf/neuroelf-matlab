function rvalue = map_WriteAnalyzeVol(xo, fmr, filename, flip)
% MAP::WriteAnalyzeVol  - write an Analyze image from one volume
%
% FORMAT:       [success] = fmr.WriteAnalyzeVol(fmr, filename [, flip]);
%
% Input fields:
%
%       fmr         FMR object (required for resolution, etc.)
%       filename    analyze filename
%       flip        char string for flipping (e.g. 'xy', default: '')
%
% Output fields:
%
%       success     true if write was successful

% Version:  v1.1
% Build:    16020516
% Date:     Feb-05 2016, 4:02 PM EST
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

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'map') || ...
    numel(fmr) ~= 1 || ~xffisobject(fmr, true, 'fmr') || ...
   ~ischar(filename) || isempty(filename)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% flipping
if nargin < 4 || ~ischar(flip) || numel(flip) > 3
    flip = '';
else
    flip = lower(flip(:)');
end

% default rvalue
rvalue = false;

% get reference
bc = xo.C;
fc = fmr.C;
ofile = fmr.F;
[opath{1:2}] = fileparts(ofile);
ofile = opath{2};

% check filename
if nargin < 3 || ~ischar(filename) || isempty(filename) || ...
    numel(filename) ~= length(filename) || length(filename) < 5 || ...
  (~strcmpi(filename(end-3:end), '.hdr') && ~strcmpi(filename(end-3:end), '.img'))
    error('neuroelf:xff:badArgument', ...
        'Invalid or missing filename argument for call to %s.', mfilename);
end

% filename
filename = filename(:)';
filename = [filename(1:end-3) 'img'];
hdrfname = [filename(1:end-3) 'hdr'];
matfname = [filename(1:end-3) 'mat'];

% dimension and datatype
xpix = fc.ResolutionX;
ypix = fc.ResolutionY;
zpix = fc.NrOfSlices;
if xpix ~= bc.DimX || ypix ~= bc.DimY || zpix ~= numel(bc.Map)
    error('neuroelf:xff:badArgument', 'FMR and MAP objects mismatch.');
end
isiz = [xpix, ypix, zpix, 1];

% X/Y resolution
xres = fc.InplaneResolutionX;
yres = fc.InplaneResolutionY;
zres = fc.SliceThickness + fc.SliceGap;
tres = fc.TR / 1000;

% positional information
xvec = xres .* [fc.RowDirX; fc.RowDirY; fc.RowDirZ];
yvec = yres .* [fc.ColDirX; fc.ColDirY; fc.ColDirZ];
ovec = [ fc.Slice1CenterX; fc.Slice1CenterY; fc.Slice1CenterZ];
lvec = [ fc.SliceNCenterX; fc.SliceNCenterY; fc.SliceNCenterZ];
zspan = fc.NrOfSlices - 1;
zvec = (1 / zspan) .* (lvec - ovec);
svec = ovec - ((xpix + 1) / 2) * xvec - ((ypix + 1) / 2) * yvec;

% build mat
tmat = [[xvec, yvec, zvec], svec;  0, 0, 0, 1];

% radiological convention ?
if isfield(fc, 'Convention') && lower(fc.Convention(1)) == 'r'
    tfmat = eye(4);
    tfmat(2, 2) = -1;
    tmat = tfmat * tmat;
end

% try volume creation
try
    hdr = xff('new:hdr');
    hdrc = hdr.C;
catch xfferror
    error('neuroelf:xff:internalError', ...
        'Error creating Analyze header object: %s.', xfferror.message);
end

% set dims and data
hdrc.ImgDim.Dim(1:5) = [4, isiz];
hdrc.ImgDim.DataType = 16;
hdrc.ImgDim.BitsPerPixel = 32;
hdrc.VoxelData = single(zeros(isiz));
hdrc.ImgDim.PixSpacing(2:5) = [xres, yres, zres, tres];
hdrc.ImgDim.CalMaxDisplay = 32767;
hdrc.ImgDim.CalMinDisplay = 0;
hdrc.ImgDim.GLMax = 32767;
hdrc.ImgDim.GLMin = 0;
hdrc.DataHist.Description = sprintf('Functional map of FMR %s', ofile);
hdrc.DataHist.ScanNumber = '1';

% get map data
maps = bc.Map;
mapdata = cat(3, maps(:).Data);

% which format
switch bc.Type

    % corr-Map
    case 1
        mapdata = sign(mapdata) .* (1 - abs(mapdata));
        mapdata(abs(mapdata) == 1) = 0;

    % cross-corr-Map
    case 2
        lag = floor(mapdata);
        mapdata = 1 - (mapdata - lag);
        mapdata(mapdata == 1) = 0;

end
hdrc.VoxelData = single(mapdata);
hdr.C = hdrc;

% save hdr/img
try
    rvalue = false;
    aft_SaveAs(hdr, hdrfname);
    rvalue = hdr_SaveVoxelData(hdr, flip);

    % for CC-maps
    if bc.Type == 2

        % also save lag information
        hdrc.VoxelData = single(lag);
        hdr.C = hdrc;
        aft_SaveAs(hdr, [hdrfname(1:end-4) '_lag.hdr']);
        rvalue = hdr_SaveVoxelData(hdr, flip);
    end
catch xfferror
    warning('neuroelf:xff:internalError', ...
        'Error writing Analyze header/image %s: %s.', hdrfname, xfferror.message);
end
delete(hdr);

% leave early ?
if ~rvalue
    return;
end

% save mat
try
    eval('M=tmat;mat=tmat;save(matfname,''M'',''mat'',''-v6'');');
catch xfferror
    neuroelf_lasterr(xfferror);
end
if exist(matfname, 'file') ~= 2
    rvalue = false;
end
if bc.Type == 2
    matfname = [matfname(1:end-4) '_lag.mat'];
    try
        eval('save(matfname,''M'',''mat'',''-v6'');');
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
    if exist(matfname, 'file') ~= 2
        rvalue = false;
    end
end

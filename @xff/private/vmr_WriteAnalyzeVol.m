function rvalue = vmr_WriteAnalyzeVol(xo, filename)
% VMR::WriteAnalyzeVol  - write an Analyze image from the VMR volume
%
% FORMAT:       [success = ] vmr.WriteAnalyzeVol(filename)
%
% Input fields:
%
%       filename    volume filename
%
% Output fields:
%
%       success     true if successful

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:04 PM EST
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
if nargin ~= 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ...
   ~ischar(filename) || isempty(filename) || numel(filename) ~= length(filename) || length(filename) < 5 || ...
  (~strcmpi(filename(end-3:end), '.hdr') && ~strcmpi(filename(end-3:end), '.img'))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% default rvalue
rvalue = false;

% create struct version
[opath{1:2}] = fileparts(xo.F);
ofile = opath{2};

% filename
filename = filename(:)';
filename = [filename(1:end-3) 'img'];
hdrfname = [filename(1:end-3) 'hdr'];

% create HDR
hdr = xff('new:hdr');
hdrbc = hdr.C;

% set dims and data
hdrbc.FileMagic = 'ni1';
hdrbc.NIIFileType = 1;
hdrbc.ImgDim.Dim(1:4) = [3, 256, 256, 256];
hdrbc.ImgDim.PixSpacing(2:4) = [1, 1, 1];
vxd = aft_SampleTalBox(xo, struct('BBox', [-128, -128, -128; 127.5, 127.5, 127.5], 'ResXYZ', 1), 2);
mxv = max(vxd(:));
if mxv <= 32767
    vxd = int16(vxd);
    hdrbc.ImgDim.DataType = 4;
    hdrbc.ImgDim.BitsPerPixel = 16;
else
    vxd = single(vxd);
    hdrbc.ImgDim.DataType = 16;
    hdrbc.ImgDim.BitsPerPixel = 32;
end
hdrbc.ImgDim.CalMaxDisplay = double(mxv);
hdrbc.ImgDim.CalMinDisplay = 0;
hdrbc.ImgDim.GLMax = hdrbc.ImgDim.CalMaxDisplay;
hdrbc.ImgDim.GLMin = 0;
hdrbc.DataHist.Description = sprintf('VMR Volume of %s', ofile);
hdrbc.DataHist.Originator(1:3) = 128;
hdrbc.DataHist.NIftI1.QFormCode = 2;
hdrbc.DataHist.NIftI1.SFormCode = 2;
hdrbc.DataHist.NIftI1.QuatOffsetX = -128;
hdrbc.DataHist.NIftI1.QuatOffsetY = -128;
hdrbc.DataHist.NIftI1.QuatOffsetZ = -128;
hdrbc.DataHist.NIftI1.AffineTransX = [1, 0, 0, -128];
hdrbc.DataHist.NIftI1.AffineTransY = [0, 1, 0, -128];
hdrbc.DataHist.NIftI1.AffineTransZ = [0, 0, 1, -128];
hdrbc.VoxelData = vxd;
hdr.C = hdrbc;

% save hdr/img
try
    aft_SaveAs(hdr, hdrfname);
    rvalue = hdr_SaveVoxelData(hdr);
catch xfferror
    warning('neuroelf:xff:internalError', ...
        'Error writing Analyze header/image %s (%s).', hdrfname, xfferror.message);
    return;
end
if ~rvalue
    return;
end

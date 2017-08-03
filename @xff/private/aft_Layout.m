function layout = aft_Layout(xo)
% AFT::Layout  - get layout signature
%
% FORMAT:       layout = obj.Layout;
%
% No input fields
%
% Output fields:
%
%       layout      1x14 vector representing the layout of the dataset
%        (1:3)      XYZ-dim ([NrOfVertices, 1, 1] for MTC/SRF)
%        (4)        T-dim (1 for 3D datasets)
%        (5:7)      XYZ-start values ([1, 1, 1] for FMR, [0, 0, 0] for SRF)
%        (8:10)     XYZ-End values ((1:3) for FMR, [256, 256, 256] for SRF)
%        (11:13)    spatial resolution ([1, 1, 1] for SRF/MTC)
%        (14)       temporal resolution (1 for 3D datasets)
%        (15)       type signature
%
% TYPES: AMR, AVA, CMP, DDT, DMR, FMR, GLM, HDR, HEAD, MAP, MGH, MTC, MSK, SMP, SRF, TVL, VDW, VMP, VMR, VTC
%
% Note: output is in BV's *internal* notation (axes not in TAL order!)
%       also, this function is only valid for native BrainVoyager formats
%       with at least either temporal or spatial information!

% Version:  v1.1
% Build:    16031615
% Date:     Mar-16 2016, 3:40 PM EST
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
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% create default layout
layout = ones(1, 15);
layout(8:10) = 256;

% switch on filetype
bc = xo.C;
ft = lower(xo.S.Extensions{1});
layout(end) = sum(double(ft) .* (2 .^ ((8 * numel(ft) - 1):-8:0)));
try
    switch (ft)
        case 'amr'
            try
                layout(1:2) = size(bc.Slice(1).AMRData);
            catch xfferror
                neuroelf_lasterr(xfferror);
                layout(1:2) = 256;
            end
            layout(3) = numel(bc.Slice);
            layout(8:10) = layout(5:7) + layout(1:3) - 1;
        case {'ava', 'cmp', 'ddt', 'glm', 'msk', 'vdw', 'vmp', 'vtc'}
            if strcmp(ft, 'glm') && bc.ProjectType == 2
                layout(1) = bc.NrOfVertices;
                layout(5:7) = 0;
                return;
            end
            if strcmp(ft, 'vtc')
                layout(4) = bc.NrOfVolumes;
                layout(14) = 0.001 * bc.TR;
            end
            layout(5:10) = [bc.XStart, bc.YStart, bc.ZStart, bc.XEnd, bc.YEnd, bc.ZEnd];
            layout(11:13) = bc.Resolution;
            try
                layout(1:3) = aft_GetVolumeSize(xo);
            catch xfferror
                neuroelf_lasterr(xfferror);
                layout(1:3) = floor((layout(8:10) - layout(5:7)) / bc.Resolution);
            end
        case 'ctc'
            layout(1) = bc.NrOfChannels;
            layout(4) = bc.NrOfSamples;
            layout(14) = 1 / bc.SamplingFrequency;
        case {'dmr', 'fmr'}
            layout(1:3) = [bc.ResolutionX, bc.ResolutionY, bc.NrOfSlices];
            layout(4) = bc.NrOfVolumes;
            layout(8:10) = layout(5:7) + layout(1:3) - 1;
            layout(11:12) = [bc.InplaneResolutionX, bc.InplaneResolutionY];
            layout(13) = bc.SliceThickness + bc.GapThickness;
            layout(14) = 0.001 * bc.TR;
        case 'eig'
            layout(4) = numel(bc.EigenValues);
        case 'hdr'
            cfr = hdr_CoordinateFrame(xo);
            c1 = cfr.Trf * ones(4, 1);
            bcidd24 = bc.ImgDim.Dim(2:4);
            cl = cfr.Trf * [bcidd24(:) + 1; 1];
            layout(1:4) = cfr.Dimensions;
            layout(5:7) = 128 - round(cl([2, 3, 1]));
            layout(8:10) = 128 - round(c1([2, 3, 1]));
            layout(11:13) = cfr.Resolution;
            layout(14) = bc.ImgDim.PixSpacing(5);
        case 'head'
            cfr = head_CoordinateFrame(xo);
            c1 = cfr.Trf * ones(4, 1);
            bcidd24 = aft_GetVolumeSize(xo);
            cl = cfr.Trf * [bcidd24(:) + 1; 1];
            layout(1:4) = cfr.Dimensions;
            layout(5:7) = 128 - round(cl([2, 3, 1]));
            layout(8:10) = 128 - round(c1([2, 3, 1]));
            layout(11:13) = cfr.Resolution;
            layout(14) = 0;
        case 'map'
            try
                layout(1:2) = size(bc.Map(1).Data);
            catch xfferror
                neuroelf_lasterr(xfferror);
                layout(1:2) = [bc.DimX, bc.DimY];
            end
            layout(3) = numel(bc.Map);
            layout(8:10) = layout(5:7) + layout(1:3) - 1;
        case 'mgh'
            cfr = mgh_CoordinateFrame(xo);
            c1 = cfr.Trf * ones(4, 1);
            bcidd24 = aft_GetVolumeSize(xo);
            cl = cfr.Trf * [bcidd24(:) + 1; 1];
            layout(1:4) = cfr.Dimensions;
            layout(5:7) = 128 - round(cl([2, 3, 1]));
            layout(8:10) = 128 - round(c1([2, 3, 1]));
            layout(11:13) = cfr.Resolution;
            layout(14) = 0;
        case 'mtc'
            layout(1) = size(bc.MTCData, 2);
            layout(4) = size(bc.MTCData, 1);
            layout(5:7) = 0;
            layout(14) = 0.001 * bc.TR;
        case 'smp'
            layout(1) = bc.NrOfVertices;
            layout(5:7) = 0;
        case 'srf'
            layout(1) = size(bc.VertexCoordinate, 1);
            layout(5:7) = 0;
        case 'vmr'
            try
                layout(1:3) = size(bc.VMRData);
            catch xfferror
                neuroelf_lasterr(xfferror);
                layout(1:3) = 0;
            end
            layout(5:7) = [bc.OffsetX, bc.OffsetY, bc.OffsetZ];
            layout(8:10) = layout(5:7) + layout(1:3) - 1;
            layout(11:13) = [bc.VoxResX, bc.VoxResY, bc.VoxResZ];
        otherwise
            error('BAD_FILETYPE');
    end
catch xfferror
    error('neuroelf:xff:internalError', ...
        'Layout is not a valid method of type %s: %s.', upper(ft), xfferror.message);
end

function mtc = vtc_CreateMTC(xo, srf, opts)
% VTC::CreateMTC  - sample an MTC from a VTC
%
% FORMAT:       mtc = vtc.CreateMTC(srf, opts)
%
% Input fields:
%
%       srf         required surface file
%       opts        1x1 struct with optional fields
%        .interp    method ('nearest', {'linear'}, 'cubic')
%        .ipfrom    interpolate from P + n * normal vector, default: -3
%        .ipstep    interpolation stepsize, default: 1 (in normal vectors)
%        .ipto      interpolate to P + n * normal vector, default: 1
%        .method    method to get value ({'mean'}, 'median')
%        .recalcn   boolean, recalc normals before sampling, default: false
%
% Output fields:
%
%       mtc         MTC object
%
% Note: results slightly differ from BV's results (sampling method)
%
% Using: flexinterpn_method.

% Version:  v1.1
% Build:    16051916
% Date:     May-19 2016, 4:14 PM EST
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
flexinterpn_method = ne_methods.flexinterpn_method;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ...
    numel(srf) ~= 1 || ~xffisobject(srf, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.%s', mfilename);
end
bc = xo.C;
srfc = srf.C;
srff = srf.F;
srfh = srf.H;
uselink = ~isempty(srfc.AutoLinkedSRF) && exist(srfc.AutoLinkedSRF, 'file') == 2 && ...
   (isempty(srff) || ~strcmpi(srff, srfc.AutoLinkedSRF));
useorig = (~uselink && isfield(srfh, 'VertexCoordinateOrig') && ...
    isequal(size(srfc.VertexCoordinate), size(srfh.VertexCoordinateOrig)));
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'interp') || ~ischar(opts.interp) || ~any(strcmpi(opts.interp(:)', ...
    {'cubic', 'lanczos3', 'lanczos4', 'lanczos5', 'lanczos8', 'linear', 'nearest'}))
    opts.interp = 'linear';
else
    opts.interp = lower(opts.interp(:))';
end
if ~isfield(opts, 'ipfrom') || ~isa(opts.ipfrom, 'double') || numel(opts.ipfrom) ~= 1 || ...
    isnan(opts.ipfrom) || opts.ipfrom < -10 || opts.ipfrom > 10
    opts.ipfrom = -3;
end
if ~isfield(opts, 'ipstep') || ~isa(opts.ipstep, 'double') || numel(opts.ipstep) ~= 1 || ...
    isnan(opts.ipstep) || opts.ipstep <= 0 || opts.ipstep > 10
    opts.ipstep = 1;
end
if ~isfield(opts, 'ipto') || ~isa(opts.ipto, 'double') || numel(opts.ipto) ~= 1 || ...
    isnan(opts.ipto) || opts.ipto < opts.ipfrom || opts.ipto > 10
    opts.ipto = 1;
end
if ~isfield(opts, 'method') || ~ischar(opts.method) || ...
   ~any(strcmpi(opts.method(:)', {'mean', 'median'}))
    opts.method = 'mean';
else
    opts.method = lower(opts.method(:))';
end
if ~isfield(opts, 'recalcn') || ~islogical(opts.recalcn) || numel(opts.recalcn) ~= 1
    opts.recalcn = false;
end
vtcres  = bc.Resolution;
if uselink
    try
        tsrf = xff(srfc.AutoLinkedSRF);
        if ~xffisobject(tsrf, true, 'srf')
            if xffisobject(tsrf, true)
                delete(tsrf);
            end
            error('xff:objectLoadError', 'Error loading linked SRF.');
        end
        if opts.recalcn
            srf_RecalcNormals(tsrf);
        end
        crd = tsrf.C.VertexCoordinate;
        nrm = (1 / vtcres) .* tsrf.C.VertexNormal;
        delete(tsrf);
    catch xfferror
        neuroelf_lasterr(xfferror);
        fprintf('Error loading linked SRF for sampling: %s\n', srfc.AutoLinkedSRF);
        uselink = false;
    end
end
if ~uselink
    if opts.recalcn
        if useorig
            srfcc = srfc.VertexCoordinate;
            srfcn = srfc.VertexNormal;
            srf.C.VertexCoordinate = srfh.VertexCoordinateOrig;
            srf_RecalcNormals(srf);
            crd = srf.C.VertexCoordinate;
            nrm = (1 / vtcres) .* srf.C.VertexNormal;
            srf.C.VertexCoordinate = srfcc;
            srf.C.VertexNormal = srfcn;
        else
            srf_RecalcNormals(srf);
            crd = srf.C.VertexCoordinate;
            nrm = (1 / vtcres) .* srf.C.VertexNormal;
        end
    else
        if useorig
            crd = srfh.VertexCoordinateOrig;
            nrm = (1 / vtcres) .* srfh.VertexNormalOrig;
        else
            crd = srfc.VertexCoordinate;
            nrm = (1 / vtcres) .* srfc.VertexNormal;
        end
    end
end
ipsamp  = opts.ipfrom:opts.ipstep:opts.ipto;
ipsnum  = numel(ipsamp);
numvols = size(bc.VTCData, 1);

% number of coordinates
numv = size(crd, 1);

% create output
mtc = xff('new:mtc');
mtcc = mtc.C;
mtcc.NrOfVertices = numv;
mtcc.NrOfTimePoints = numvols;
mtcc.SourceVTCFile = xo.F;
mtcc.LinkedPRTFile = bc.NameOfLinkedPRT;
mtcc.HemodynamicDelay = bc.HemodynamicDelay;
mtcc.TR = bc.TR;
mtcc.HRFDelta = bc.HrfDelta;
mtcc.HRFTau = bc.HrfTau;
mtcc.ProtocolSegmentSize = bc.SegmentSize;
mtcc.ProtocolSegmentOffset = bc.SegmentOffset;
mtcc.MTCDataType = 1;
mtcc.MTCData = single(0);
mtcc.MTCData(numvols, numv) = 0;
mtcc.RunTimeVars.Map.DF1 = numvols - 2;
mtcc.RunTimeVars.Map.BonferroniValue = numv;

% subtract XStart, YStart, ZStart
crd = 1 + (1 / vtcres) .* [crd(:, 1) - bc.XStart, crd(:, 2) - bc.YStart, crd(:, 3) - bc.ZStart];

% prepare interpolation array
ipvals = zeros(numv, ipsnum);

% prepare interpolation coordinates
ipcrds = crd(:, :, ones(1, ipsnum));
for ipc = 1:ipsnum
    ipcrds(:, :, ipc) = crd + ipsamp(ipc) * nrm;
end

% iterate over maps
for vc = 1:numvols

    % get source volume
    mapval = permute(bc.VTCData(vc, :, :, :), [2, 3, 4, 1]);

    % interpolate according to method
    for ipc = 1:ipsnum
        ipvals(:, ipc) = flexinterpn_method(mapval, ipcrds(:, :, ipc), 0, opts.interp);
    end

    % what summary method
    switch opts.method
        case 'mean'
            mtcc.MTCData(vc, :) = single(mean(ipvals, 2));
        case 'median'
            mtcc.MTCData(vc, :) = single(median(ipvals, 2));
    end
end

% take over general fields from RunTimeVars
rtv = bc.RunTimeVars;
rtf = fieldnames(rtv);
for fc = 1:numel(rtf)
    if ~any(strcmp(rtf{fc}, {'xffID', 'TrfPlus', 'AvgVTC'}))
        mtcc.RunTimeVars.(rtf{fc}) = rtv.(rtf{fc});
    end
end

% average VTC
if isfield(rtv, 'AvgVTC')
    mtcc.RunTimeVars.AvgMTC = rtv.AvgVTC;
end

% put content in new object
mtc.C = mtcc;

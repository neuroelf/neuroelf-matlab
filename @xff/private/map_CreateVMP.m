function vmp = map_CreateVMP(xo, fmr, afs, vmpfile, meth, res, bbox, o)
% MAP::CreateVMP  - create VMP from FMR based MAP object
%
% FORMAT:       vmp = map.CreateVMP(fmr, afs [, vmpfile, meth, res, bb, o])
%
% Input fields:
%
%       fmr         FMR object matching MAP
%       afs         cell array with alignment files, {ia, fa, acpc, tal}
%       vmpfile     target VMP filename, if not given: unsaved
%       meth        interpolation, 'nearest', 'lanczos3', {'linear'}, 'cubic'
%       res         VMP resolution, default: 1
%       bb          bounding box, default: small TAL box
%       o           optional struct with settings
%        .df        degrees of freedom, 1x1 or 1x2 DF setting
%
% Output fields:
%
%       vmp         resulting VMP object
%
% Using: applyfdr, samplefmrspace.

% Version:  v1.1
% Build:    16020516
% Date:     Feb-05 2016, 4:08 PM EST
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
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'map') || ...
    numel(fmr) ~= 1 || ~xffisobject(fmr, true, 'fmr') || ...
   ~iscell(afs) || isempty(afs) || numel(afs{1}) ~= 1 || ...
   ~xffisobject(afs{1}, true, 'trf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
mapc = xo.C;
if ~any([0, 1, 2, 3] == mapc.Type)
    error('neuroelf:xff:invalidObject', 'Unsupported map type: %d.', mapc.Type);
end
if nargin < 4 || ~ischar(vmpfile) || isempty(vmpfile)
    vmpfile = '';
else
    vmpfile = vmpfile(:)';
end
fmrc = fmr.C;
if numel(mapc.Map) ~= fmrc.NrOfSlices || numel(size(mapc.Map(1).Data)) ~= 2 || ...
    any(size(mapc.Map(1).Data) ~= [fmrc.ResolutionX, fmrc.ResolutionY])
    error('neuroelf:xff:invalidObject', 'FMR and MAP dimension(s) mismatch.');
end
if nargin < 5 || ~ischar(meth) || isempty(meth)
    meth = 'linear';
else
    meth = lower(meth(:)');
end
if nargin < 6 || ~isa(res, 'double') || numel(res) ~= 1 || ~any([1, 2, 3] == res)
    res = 1;
end
if nargin < 7 || ~isa(bbox, 'double') || numel(size(bbox)) ~= 2 || ...
    any(size(bbox) ~= [2, 3]) || ...
    any(isinf(bbox(:)) | isnan(bbox(:)) | bbox(:) < 1 | bbox(:) > 254 | ...
        bbox(:) ~= fix(bbox(:))) || any((bbox(2, :) - res) <= bbox(1, :))
    bbox = [59, 57, 52; 197, 231, 172];
end
if nargin < 8 || ~isstruct(o) || numel(o) ~= 1
    o = struct;
end
if ~isfield(o, 'df') || ~isa(o.df, 'double') || numel(o.df) > 2 || ...
    any(isnan(o.df) | isinf(o.df) | o.df < 1 | o.df ~= fix(o.df))
    o.df = [];
end
if res ~= 1
    bbox(2, :) = bbox(1, :) + res .* round(diff(bbox) ./ res);
end

% generate target coordinate list and data field
if res == 1
    c = {bbox(1, 1):bbox(2, 1), bbox(1, 2):bbox(2, 2), bbox(1, 3):bbox(2, 3)};
else
    c = {bbox(1, 1):res:(bbox(2, 1) - 1), bbox(1, 2):res:(bbox(2, 2) - 1), ...
         bbox(1, 3):res:(bbox(2, 3) - 1)};
end

% check alignment files
hasia = 0;
clearia = cell(1, 1);
for ac = numel(afs):-1:1
    if numel(afs{ac}) == 1 && xffisobject(afs{ac}, true, 'trf')
        if hasia > 0
            afs(hasia) = [];
        end
        trfc = afs{ac}.C;
        if trfc.TransformationType == 1 && trfc.AlignmentStep == 1
            hasia = ac;
        end
    elseif numel(afs{ac}) ~= 1 || ~xffisobject(afs{ac}, true, 'tal')
        error('neuroelf:xff:badArgument', ...
            'Invalid alignment file object list.');
    end
end
if hasia == 0
    warning('neuroelf:xff:badArgument', ...
        'Initial alignment missing; using identity matrix instead.');
    afs = [{xff('new:trf')}, afs(:)'];
    clearia = afs{1};
end

% generate sample output with target coordinates
nslc = fmrc.NrOfSlices;
smpo = zeros(fmrc.ResolutionX, fmrc.ResolutionY, nslc);

% fill with map content
for mc = 1:nslc
    smpo(:, :, mc) = mapc.Map(mc).Data(:, :);
end
if mapc.Type == 1
    smpo = sign(smpo) .* (1 - abs(smpo));
    smpo(abs(smpo) == 1) = 0;
elseif mapc.Type == 2
    lag = floor(smpo);
    smpo = 1 - (smpo - lag);
    smpo(smpo == 1) = 0;
end
bfv = sum(smpo(:) ~= 0);

% sample at VMP coordinates
try
    smpd = permute(single(ne_methods.samplefmrspace(smpo, c, fmr, afs, meth)), [2, 3, 1]);
    if mapc.Type == 2
        smpd = smpd + round( ...
            permute(single(ne_methods.samplefmrspace(lag, c, fmr, afs, meth)), [2, 3, 1]));
    end
    smpd(isnan(smpd) | isinf(smpd)) = 0;
catch xfferror
    clearxffobjects(clearia);
    error('neuroelf:xff:internalError', 'Error sampling MAP data (%s).', xfferror.message);
end

% create and fill VMP
vmp = xff('new:vmp');
vmpc = vmp.C;
ofv = vmpc.FileVersion;
if res == 1
    vmpc.NativeResolutionFile = 0;
    vmpc.FileVersion = 4;
else
    vmpc.NativeResolutionFile = 1;
    vmpc.DocumentType = 1;
    vmpc.FileVersion = 5;
end
vmpc.Resolution = res;
vmp.C = vmpc;
vmp_Update(vmp, 'FileVersion', struct('type', '.', 'subs', 'FileVersion'), ofv);
vmpc = vmp.C;
vmpc.NrOfMaps = 1;
vmpc.XStart = bbox(1, 2);
vmpc.XEnd = bbox(2, 2);
vmpc.YStart = bbox(1, 3);
vmpc.YEnd = bbox(2, 3);
vmpc.ZStart = bbox(1, 1);
vmpc.ZEnd = bbox(2, 1);
if res ~= 1
    vmpc.LinkedPRT = fmrc.ProtocolFile;
end
vmpc.Map.Type = mapc.Type + 1;
[mapf{1:2}] = fileparts(xo.F);
vmpc.Map.Name = sprintf('VMP from MAP: %s', mapf{2});
vmpc.Map.LowerThreshold = mapc.LowerThreshold;
vmpc.Map.UpperThreshold = mapc.UpperThreshold;
if ~isempty(mapc.NrOfLags)
    vmpc.Map.NrOfLags = mapc.NrOfLags;
    vmpc.Map.MinLag = 0;
    vmpc.Map.MaxLag = mapc.NrOfLags - 1;
    vmpc.Map.CCOverLay = 1;
end
vmpc.Map.ClusterSize = round((fmrc.InplaneResolutionX * fmrc.InplaneResolutionY * ...
    (fmrc.SliceThickness + fmrc.GapThickness) * mapc.ClusterSize) / (res ^ 3));
vmpc.Map.EnableClusterCheck = 1;

% guessing here, true for r !
vmpc.Map.DF2 = 0;
if isempty(o.df)
    if map.FileVersion < 3
        if vmpc.Map.Type ~= 2
            vmpc.Map.DF1 = fmrc.NrOfVolumes - 8;
        else
            vmpc.Map.DF1 = fmrc.NrOfVolumes - 1;
        end
        if vmpc.Map.Type == 4
            vmpc.Map.DF2 = 1;
        end
    else
        vmpc.Map.DF1 = mapc.DF1;
        vmpc.Map.DF2 = mapc.DF2;
    end

% df are provided
else
    vmpc.Map.DF1 = o.df(1);
    if numel(o.df) > 1
        vmpc.Map.DF2 = o.df(2);
    end
end

% adapt Bonferroni value
bfv2 = sum(smpd(:) ~= 0);
if bfv2 < bfv
    bfv = bfv2;
end
vmpc.Map.BonferroniValue = bfv;
if res ~= 1
    vmpc.Map.FDRThresholds = [0, 1e5, 1e5];
    if any([1, 2, 4] == vmpc.Map.Type)
        mapv = smpo(~isnan(smpo(:)) & ~isinf(smpo(:)) & (smpo(:) ~= 0));
        fdrt = [0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001];
        switch (vmpc.Map.Type)
            case 1 % t-score
                vmpc.Map.FDRThresholds = [fdrt(:), ...
                    ne_methods.applyfdr(double(mapv), 't', ...
                    fdrt(:), vmpc.Map.DF1, [], true)];
            case 2 % r-value
                vmpc.Map.FDRThresholds = [fdrt(:), ...
                    ne_methods.applyfdr(double(mapv), 'r', ...
                    fdrt(:), vmpc.Map.DF1, [], true)];
            case 4 % F-score
                vmpc.Map.FDRThresholds = [fdrt(:), ...
                    ne_methods.applyfdr(double(mapv), 'F', ...
                    fdrt(:), vmpc.Map.DF1, vmpc.Map.DF2, true)];
        end
    end
    vmpc.Map.NrOfFDRThresholds = size(vmpc.Map.FDRThresholds, 1);
end
vmpc.Map.VMPData = smpd;
vmp.C = vmpc;

% try saving
try
    if ~isempty(vmpfile)
        aft_SaveAs(vmp, vmpfile);
    end
catch xfferror
    warning('neuroelf:xff:internalError', 'Error saving VMP file: ''%s''.', xfferror.message);
end

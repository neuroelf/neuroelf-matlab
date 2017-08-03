function [no1, no2, no3, no4] = srf_Combine(xo, xo2, cbopt)
% SRF::Combine  - combines two SRFs into one
%
% FORMAT:       [csrf [, ...]] = srf1.Combine(srf2, cbopt);
%
% Input fields:
%
%       srf2        surface to use to combine given SRF with
%       cbopt       struct with optional settings
%       .type       one of
%                   'backtoback' - rotate one mesh 180 in XY plane
%                   'custom'     - build custom scene, see below
%                   'gapped'     - join contents with a gap
%                   'wholebrain' - simply join contents (default)
%       .color1     1x4 double for .Color1 field in SRF options
%       .color2     1x4 double for .Color2 field in SRF options
%       .filename   store combined under new filename
%       .gap        1x1 double, mm to insert between two meshes
%                   (applied for all types accordingly, defaults:
%                   backtoback: 25, gapped: 100, outandin: 20,
%                   outintb: 25, patched: 25, spm2: 25, wholebrain: 0)
%       .linkedsrf  filename of linked SRF, set to empty if not given
%       .mtc1       xff MTC object for the first SRF
%       .mtc2       xff MTC object for the second SRF
%       .smp1       xff SMP object for the first SRF
%       .smp2       xff SMP object for the second SRF
%       .ssm1       xff SSM object for the first SRF
%       .ssm2       xff SSM object for the second SRF
%       .transform  1x2 cell array with 4x4 double transformation matrices,
%                   needed for custom scenaries
%
% Output fields:
%
%       csrf        combined SRF
%       ...         combined MTC/SMP/SSM (if given, in that order!)
%
% Using: tfmatrix.

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:52 PM EST
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

% check arguments
if nargin < 2 || numel(xo) ~= 1 || numel(xo2) ~= 1 || ...
   ~xffisobject(xo, true, {'fsbf', 'srf'}) || ~xffisobject(xo2, true, {'fsbf', 'srf'}) || ...
   ~strcmp(xo.S.Extensions{1}, xo2.S.Extensions{1})
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
bc1 = xo.C;
bc2 = xo2.C;
no1 = aft_CopyObject(xo);
no1.F = '';
no2 = [];
no3 = [];
no4 = [];

% argument check
if nargin < 3 || ~isstruct(cbopt) || numel(cbopt) ~= 1
    cbopt = struct;
end
if ~isfield(cbopt, 'type') || ~ischar(cbopt.type) || ...
   ~any(strcmpi(cbopt.type, {'backtoback', 'custom', 'gapped', 'wholebrain'}))
    cbopt.type = 'wholebrain';
end
if isfield(cbopt, 'mtc1') && isfield(cbopt, 'mtc2') && ...
    numel(cbopt.mtc1) == 1 && numel(cbopt.mtc2) == 1 && ...
    xffisobject(cbopt.mtc1, true, 'mtc') && xffisobject(cbopt.mtc2, true, 'mtc')
    mtcin = true;
    mbc1 = cbopt.mtc1.C;
    mbc2 = cbopt.mtc2.C;
    if mbc1.NrOfVertices ~= bc1.NrOfVertices || mbc2.NrOfVertices ~= bc2.NrOfVertices || ...
        size(mbc1.MTCData, 1) ~= size(mbc2.MTCData, 1)
        delete(no1);
        error('neuroelf:xff:badArgument', 'NrOfVertices and NrOfTimePoints must match in MTCs.');
    end
    no2 = aft_CopyObject(cbopt.mtc1);
    no2.F = '';
else
    mtcin = false;
end
if isfield(cbopt, 'smp1') && isfield(cbopt, 'smp2') && ...
    numel(cbopt.smp1) == 1 && numel(cbopt.smp2) == 1 && ...
    xffisobject(cbopt.smp1, true, 'smp') && xffisobject(cbopt.smp2, true, 'smp')
    smpin = true;
    sbc1 = cbopt.smp1.C;
    sbc2 = cbopt.smp2.C;
    if sbc1.NrOfVertices ~= bc1.NrOfVertices || sbc2.NrOfVertices ~= bc2.NrOfVertices || ...
        numel(sbc1.Map) ~= numel(sbc2.Map)
        delete(no1);
        if ~isempty(no2)
            delete(no2);
        end
        error('neuroelf:xff:badArgument', 'NrOfVertices and NrOfMaps must match in SMPs.');
    end
    no3 = aft_CopyObject(cbopt.smp1);
    no3.F = '';
else
    smpin = false;
end
if isfield(cbopt, 'ssm1') && isfield(cbopt, 'ssm2') && ...
    numel(cbopt.ssm1) == 1 && numel(cbopt.ssm2) == 1 && ...
    xffisobject(cbopt.ssm1, true, 'ssm') && xffisobject(cbopt.ssm2, true, 'ssm')
    ssmin = true;
    ssm1 = cbopt.ssm1.C;
    ssm2 = cbopt.ssm2.C;
    if ssm1.NrOfTargetVertices ~= bc1.NrOfVertices || ssm2.NrOfTargetVertices ~= bc2.NrOfVertices
        delete(no1);
        if ~isempty(no2)
            delete(no2);
        end
        if ~isempty(no3)
            delete(no3);
        end
        error('neuroelf:xff:badArgument', 'NrOfVertices and NrOfTimePoints must match in MTCs.');
    end
    no4 = aft_CopyObject(cbopt.ssm1);
    no4.F = '';
else
    ssmin = false;
end

% check some more things for custom type
cbtype = lower(cbopt.type);
if strcmp(cbtype, 'custom')
    if ~isfield(cbopt, 'transform') || ~iscell(cbopt.transform) || ...
        numel(cbopt.transform) ~= 2
        delete(no1);
        if ~isempty(no2)
            delete(no2);
        end
        if ~isempty(no3)
            delete(no3);
        end
        if ~isempty(no4)
            delete(no4);
        end
        error('neuroelf:xff:badArgument', 'Missing or bad transformation matrices for input SRF''s.');
    end
    try
        if ~all(size(cbopt.transform{1}) == 4)
            srf_t1 = ne_methods.tfmatrix(cbopt.transform{1});
        else
            srf_t1 = cbopt.transform{1}(:, :, 1);
        end
        if ~all(size(cbopt.transform{2}) == 4)
            srf_t2 = ne_methods.tfmatrix(cbopt.transform{2});
        else
            srf_t2 = cbopt.transform{2}(:, :, 1);
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
        delete(no1);
        if ~isempty(no2)
            delete(no2);
        end
        if ~isempty(no3)
            delete(no3);
        end
        if ~isempty(no4)
            delete(no4);
        end
        error('neuroelf:xff:badArgument', 'Bad transformation matrix for input SRF %d.', tc);
    end
end
if isfield(cbopt, 'color1') && isa(cbopt.color1, 'double') && numel(cbopt.color1) == 4 && ...
   ~any(isnan(cbopt.color1(:)) | isinf(cbopt.color1(:)) | ...
    cbopt.color1(:) < 0 | cbopt.color1(:) > 1)
    no1.C.ConvexColor = cbopt.color1(:)';
end
if isfield(cbopt, 'color2') && isa(cbopt.color2, 'double') && numel(cbopt.color2) == 4 && ...
   ~any(isnan(cbopt.color2(:)) | isinf(cbopt.color2(:) | ...
    cbopt.color2(:) < 0 | cbopt.color2(:) > 1))
    no1.C.ConcaveColor = cbopt.color2(:)';
end
if ~isfield(cbopt, 'gap') || ~isa(cbopt.gap, 'double') || numel(cbopt.gap) ~= 1 || ...
    isnan(cbopt.gap) || isinf(cbopt.gap) || cbopt.gap < -512 || cbopt.gap > 512
    cbopt.gap = [];
end
if isfield(cbopt, 'linkedsrf') && ischar(cbopt.linkedsrf) && ~isempty(cbopt.linkedsrf)
    no1.C.AutoLinkedSRF = cbopt.linkedsrf(:)';
end

if mtcin
    no2.C.MTCData = [mbc1.MTCData, mbc2.MTCData];
    no2.C.NrOfVertices = size(no2.C.MTCData, 2);
end
if smpin
    for mc = 1:numel(sbc1.Map)
        no3.C.Map(mc).SMPData = [sbc1.Map(mc).SMPData; sbc2.Map(mc).SMPData];
    end
    no3.C.NrOfVertices = size(no3.C.Map(1).SMPData, 1);
end
if ssmin
    no4.C.SourceOfTarget = [ssm1.SourceOfTarget(:); ...
        ssm2.SourceOfTarget(:) + ssm1.NrOfSourceVertices];
    no4.C.NrOfSourceVertices = ssm1.NrOfSourceVertices + ssm2.NrOfSourceVertices;
    no4.C.NrOfTargetVertices = numel(no4.C.SourceOfTarget);
end

% get coordinates, etc.
srf_p1 = bc1.VertexCoordinate;
srf_n1 = bc1.VertexNormal;
srf_np1 = size(srf_p1, 1);
srf_nt1 = size(bc1.TriangleVertex, 1);
srf_p2 = bc2.VertexCoordinate;
srf_n2 = bc2.VertexNormal;
srf_np2 = size(srf_p2, 1);
srf_nt2 = size(bc2.TriangleVertex, 1);

% summary
srf_np = srf_np1 + srf_np2;
srf_nt = srf_nt1 + srf_nt2;

% wholebrain is special (zero) gap case
if strcmp(cbtype, 'wholebrain')
    if isempty(cbopt.gap), cbopt.gap = 0; end
    cbtype = 'gapped';
end

% for new center
cpos = (bc1.MeshCenter + bc2.MeshCenter) / 2;

% which type
switch cbtype
    case 'backtoback'
        if isempty(cbopt.gap)
            cbopt.gap = 25;
        end
        srf_p1(:, 2) = srf_p1(:, 2) - min(srf_p1(:, 2)) + cbopt.gap / 2;
        srf_p2(:, 2) = srf_p2(:, 2) - min(srf_p2(:, 2)) + cbopt.gap / 2;
        srf_p2(:, 1:2) = -srf_p2(:, 1:2);
        srf_n2(:, 1:2) = -srf_n2(:, 1:2);
        cpos = [256, 128, 128];

    case 'custom'

        % prepare points and normals matrices
        srf_p1(:, 4) = 1;
        srf_p2(:, 4) = 1;

        % apply transformation matrices to points and normals
        srf_p1 = (srf_t1 * srf_p1')';
        srf_n1 = (srf_t1(1:3, 1:3) * srf_n1')';
        srf_p2 = (srf_t2 * srf_p2')';
        srf_n2 = (srf_t2(1:3, 1:3) * srf_n2')';

        % remove fourth coordinate again
        srf_p1(:, 4) = [];
        srf_p2(:, 4) = [];

    case {'gapped'}

        if isempty(cbopt.gap)
            cbopt.gap = 100;
        end
        srf_p1(:, 3) = srf_p1(:, 3) - cbopt.gap / 2;
        srf_p2(:, 3) = srf_p2(:, 3) + cbopt.gap / 2;
end

% build neighbors and triangle strip lists
neigh = [bc1.Neighbors; bc2.Neighbors];
for nc = (srf_np1 + 1):srf_np
    neigh{nc, 2} = neigh{nc, 2} + srf_np1;
end
trist = bc2.TriangleStripSequence;
rtrist = (trist > 0);
trist(rtrist) = trist(rtrist) + srf_np1;

% set fields correctly
no1.C.NrOfVertices = srf_np;
no1.C.NrOfTriangles = srf_nt;
no1.C.VertexCoordinate = [srf_p1; srf_p2];
no1.C.VertexNormal = [srf_n1; srf_n2];
no1.C.MeshCenter = cpos;
no1.C.VertexColor = [bc1.VertexColor; bc2.VertexColor];
no1.C.Neighbors = neigh;
no1.C.TriangleVertex = [bc1.TriangleVertex; (bc2.TriangleVertex + srf_np1)];
no1.C.TriangleStripSequence = [bc1.TriangleStripSequence; trist];
no1.C.NrOfTriangleStrips = numel(no1.C.TriangleStripSequence);

% write output?
if isfield(cbopt, 'filename') && ischar(cbopt.filename) && ~isempty(cbopt.filename)
    try
        aft_SaveAs(no1, cbopt.filename);
        if mtcin && numel(cbopt.filename) > 4 && ...
            strcmpi(cbopt.filename(end-3:end), '.srf')
            aft_SaveAs(no2, [cbopt.filename(1:end-4) '.mtc']);
        end
        if smpin && numel(cbopt.filename) > 4 && ...
            strcmpi(cbopt.filename(end-3:end), '.srf')
            aft_SaveAs(no3, [cbopt.filename(1:end-4) '.smp']);
        end
        if ssmin && numel(cbopt.filename) > 4 && ...
            strcmpi(cbopt.filename(end-3:end), '.srf')
            aft_SaveAs(no4, [cbopt.filename(1:end-4) '.ssm']);
        end
    catch xfferror
        warning('neuroelf:xff:errorWritingFile', ...
            'Couldn''t write new SRF/MTC/SMP/SSM file: ''%s''. (%s)', cbopt.filename, xfferror.message);
    end
end

% put into correct outputs
if isempty(no3)
    no3 = no4;
    no4 = [];
end
if isempty(no2)
    no2 = no3;
    no3 = no4;
    no4 = [];
end

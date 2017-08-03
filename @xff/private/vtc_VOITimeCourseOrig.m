function [voitc, weightv, weightr] = vtc_VOITimeCourseOrig(xo, voi, weight, fliplr, vmrbbox)
% VTC::VOITimeCourseOrig  - extract VOI time course data
%
% FORMAT:       voitc [, uvec, uvecr] = vtc.VOITimeCourseOrig(voi [, weight, fliplr, vmrbbox])
%
% Input fields:
%
%       voi         VOI file or coordinates (e.g. from VOI::BVCoords)
%       weight      (cell array of) Nx1 vector(s) with voxel weights
%                   give scalar 0 for unique, scalar 1 for no weighting,
%                   a scalar -1 for SVD after z-transform, or
%                   a scalar [Inf] to get a cell array of TxV arrays
%       fliplr      flip left/right (Z axes) for radiological convention
%       vmrbbox     VMR object's BoundingBox return value
%
% Output fields:
%
%       voitc       TxV time course of voi(s)
%       uvec        unique VTC voxel indices within VOI, so that
%                   voi.VOI(i).Voxels(uvec{i}) leads to those coordinates
%       uvecr       reverse indexing (to find out which anatomical voxels
%                   fall into which functional voxels)
%
% Note: this method is the original implementation for VTCs and obsolete!
%
% Using: ztrans.

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:34 PM EST
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
ztrans = ne_methods.ztrans;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ...
   (~all(xffisobject(voi(:), true, 'voi')) && ~isa(voi, 'double')) || isempty(voi)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 5 || ~isstruct(vmrbbox) || numel(vmrbbox) ~= 1 || ~isfield(vmrbbox, 'DimXYZ')
    vmrbbox = struct('DimXYZ', [256, 256, 256]);
end
bc = xo.C;
if numel(voi) == 1 && xffisobject(voi, true, 'voi')
    voic = voi.C;
    if isempty(voic.VOI)
        error('neuroelf:xff:badArgument', 'Invalid VOI object given.');
    end
    numvois = numel(voic.VOI);
    vois = cell(1, numvois);
    for vc = 1:numvois
        vois{vc} = voi_BVCoords(voi, vc, vmrbbox);
    end
else
    if ~any(size(voi)) == 3 || any(isinf(voi(:)) | isnan(voi(:)) | voi(:) < 0 | voi(:) > 255)
        error('neuroelf:xff:badArgument', 'Invalid VOI coordinates given.');
    end
    if size(voi, 1) == 3 && size(voi, 2) ~= 3
        voi = voi';
    end
    vois = {round(voi)};
end

% get VTC info
vres = bc.Resolution;
xstr = bc.XStart;
xend = bc.XEnd;
ystr = bc.YStart;
yend = bc.YEnd;
zstr = bc.ZStart;
zend = bc.ZEnd;
vtcsz = size(bc.VTCData);
numtp = vtcsz(1);
if vres * vtcsz(2) ~= (xend - xstr) || vres * vtcsz(3) ~= (yend - ystr) || vres * vtcsz(4) ~= (zend - zstr)
    error('neuroelf:xff:badObject', 'Invalid Resolution/Start/End/Size combination.');
end
voff = [xstr, ystr, zstr];

% check vois / weights
numvois = numel(vois);
if nargin < 3 || isempty(weight) || (~iscell(weight) && ~isa(weight, 'double'))
    weight = 1;
end
if ~iscell(weight)
    weight = {weight(:)};
end
if numel(weight) ~= numvois && numel(weight) == 1
    weight = weight(ones(1, numvois));
end
if numel(weight) ~= numvois
    error('neuroelf:xff:badArgument', 'Number of weighting vectors mismatches VOIs.');
end

% flipping
if nargin > 3 && (islogical(fliplr) || isa(fliplr, 'double')) && ~isempty(fliplr) && fliplr(1)
    fliplr = true;
else
    fliplr = false;
end

% make dimension check
cellout = false;
for vc = 1:numvois
    if isempty(weight{vc})
        weight{vc} = ones(size(vois{vc}, 1), 1);
    elseif numel(weight{vc}) ~= size(vois{vc}, 1)
        weight{vc} = weight{vc}(ones(1, size(vois{vc}, 1)), 1);
    end
    if any(isinf(weight{vc}))
        cellout = true;
    end
    if size(vois{vc}, 2) ~= 3 || size(weight{vc}, 2) ~= 1 || size(vois{vc}, 1) ~= size(weight{vc}, 1)
        error('neuroelf:xff:invalidArgument', 'Invalid VOI/weight combination (dim mismatch).');
    end

    % flip
    if fliplr
        vois{vc}(:, 3) = vmrbbox.DimXYZ(3) - vois{vc}(:, 3);
    end

    % prepare coord lists
    vois{vc} = floor(1 + (vois{vc} - repmat(voff, [size(vois{vc}, 1), 1])) / vres);

    % remove bad entries
    be = (vois{vc}(:, 1) < 1 | vois{vc}(:, 1) > vtcsz(2) | ...
        vois{vc}(:, 2) < 1 | vois{vc}(:, 2) > vtcsz(3) | ...
        vois{vc}(:, 3) < 1 | vois{vc}(:, 3) > vtcsz(4));
    vois{vc}(be, :) = [];
    vois{vc} = sub2ind(vtcsz(2:4), vois{vc}(:, 1), vois{vc}(:, 2), vois{vc}(:, 3));
    if numel(weight{vc}) > 1
        weight{vc}(be) = [];
    end
end

% initialize voitc
weightv = {};
weightr = {};
if cellout
    voitc = {zeros(numtp, 1)};
    if numvois > 1
        voitc(2:numvois) = voitc(1);
    end
    weightv = cell(1, numvois);
    weightr = cell(1, numvois);
else
    voitc = zeros(numtp, numvois);
end

% pre-load VTC data for more than 25% of the voxels
if istransio(bc.VTCData) && numel(cat(1, vois{:})) >= (0.25 * prod(vtcsz(2:4)))
    try
        bc.VTCData = resolve(bc.VTCData);
    catch xfferror
        rethrow(xfferror);
    end
end

% iterate over vois
for vc = 1:numvois

    % extract voxel time courses
    voxs = vois{vc};
    dosvd = false;
    if all(weight{vc} == weight{vc}(1))
        if weight{vc} == 0
            doweight = false;
            voxs = unique(voxs);
        elseif weight{vc}(1) == -1
            doweight = false;
            dosvd = true;
        elseif isinf(weight{vc}(1))
            doweight = false;
            dosvd = false;
        else
            doweight = true;
        end
    else
        doweight = true;
    end
    if istransio(bc.VTCData)
        [voxsu, voxsi, voxsj] = unique(voxs);
        vvtc = bc.VTCData(:, voxsu);
        vvtc = double(vvtc(:, voxsj));
    else
        vvtc = double(bc.VTCData(:, voxs));
    end

    % which method (weighting or unique)
    if doweight
        vvtc = vvtc .* repmat(weight{vc}', [numtp, 1]);
        voitc(:, vc) = sum(vvtc, 2) / sum(weight{vc});
    elseif dosvd
        [u{1:3}] = svd(ztrans(vvtc));
        u = u{1}(:, 1);
        u = u ./ std(u);
        cv = cov([mean(vvtc, 2), u]);
        if cv(1, 2) < 0
            u = -u;
        end
        voitc(:, vc) = u;
    else
        if cellout
            voitc{vc} = vvtc;
            [uvec{1:3}] = unique(voxs, 'rows');
            weightv{vc} = uvec{2};
            weightr{vc} = uvec{3};
        else
            voitc(:, vc) = mean(vvtc, 2);
        end
    end
end

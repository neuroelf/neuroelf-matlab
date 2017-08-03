% PUBLIC FUNCTION ne_markclusters: mark clusters in VMR
function varargout = ne_markclusters(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:35 AM EST
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% check if slicevar is VMR object
if numel(cc.SliceVar) ~= 1 || ...
   ~isxff(cc.SliceVar, 'vmr') || ...
   ~isxff(ne_gcfg.voi, 'voi')
    return;
end
slvar = cc.SliceVar;

% get list of clusters
clidx = ch.Clusters.Value(:);
if isempty(clidx)
    return;
end

% get matching VOIs
try
    vois = ne_gcfg.voi.VOI(clidx);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% ask for code to use
code = inputdlg({'Color code:', 'Resolution:'}, ...
    'NeuroElf GUI - draw VOI to VMR', 1, {sprintf('  %d', cc.paint.code), '  1'});
if numel(code) ~= 2 || ...
   ~iscell(code) || ...
    isempty(code{1}) || ...
   ~ischar(code{1}) || ...
   ~ischar(code{2})
    return;
end
code{1}(code{1} == ' ') = [];
code{2}(code{2} == ' ') = [];
if ~isempty(code{2})
    rad = str2double(code{2});
else
    rad = 1;
end
if numel(rad) ~= 1 || ...
    isinf(rad) || ...
    isnan(rad) || ...
    rad == 0
    rad = 1;
else
    rad = max(-7, min(7, fix(rad)));
end
code = str2double(code{1});
if numel(code) ~= 1 || ...
    isinf(code) || ...
    isnan(code) || ...
    code < 0 || ...
    code > 255
    return;
end
ne_gcfg.fcfg.paint.code = code;
code = fix(code);

% special case -> resolution < 0
if rad < 0

    % prepare
    rad = -0.5 - rad;
    from = ceil(-0.5 * rad);
    to = ceil(0.5 * (rad - 1));
    nrep = numel(from:to);
    nrep = nrep * nrep * nrep;

    % delete old RGBImage if any
    voih = handles(ne_gcfg.voi);
    if isfield(voih, 'RGBImage') && ...
        numel(voih.RGBImage) == 1 && ...
        isxff(voih.RGBImage, 'hdr')
        voih.RGBImage.ClearObject;
    end

    % create new RGBImage
    rgbd = uint8(0);
    rgbd(256, 256, 256, 1, 3) = 0;
    trf = slvar.BoundingBox.QuatB2T;
    ne_gcfg.voi.SetHandle('RGBImage', newhdr([256, 256, 256], trf, ...
        struct('data', rgbd(:, :, :, 1))));

    % get image object
    voih = handles(ne_gcfg.voi);
    rgbo = voih.RGBImage;

    % and "upgrade" to RGB
    rgbo.ImgDim.DataType = 128;
    rgbo.ImgDim.BitsPerPixel = 24;
    rgbo.VoxelData = [];

    % then fill with selected VOIs
    trf = inv(trf)';
    for vc = 1:numel(vois)
        voicol = vois(vc).Color;
        if nrep == 1
            voivox = round([vois(vc).Voxels, ones(size(vois(vc).Voxels, 1), 1)] * trf);
        else
            voivxo = vois(vc).Voxels;
            voivxn = size(voivxo, 1);
            voivx1 = ones(voivxn, 1);
            voivxo = [voivxo, voivx1];
            voivox = zeros(nrep * voivxn, 4);
            voivxi = 1;
            for xc = from:to
                for yc = from:to
                    for zc = from:to
                        voivox(voivxi:voivxi+voivxn-1, :) = ...
                            round((voivxo + voivx1 * [xc, yc, zc, 0]) * trf);
                        voivxi = voivxi + voivxn;
                    end
                end
            end
        end
        voivox(:, 4) = [];
        voivox(any(voivox < 1, 2) | any(voivox > 256, 2), :) = [];
        voivox = unique(voivox, 'rows');
        voivox = 1 + sum((voivox - 1) .* (ones(size(voivox, 1), 1) * [1, 256, 65536]), 2);
        curcol = double([lsqueeze(rgbd(voivox)), ...
            lsqueeze(rgbd(voivox + 16777216)), ...
            lsqueeze(rgbd(voivox + 33554432))]);
        ovrcol = all(curcol == 0, 2);
        curcol(ovrcol, :) = voicol(ones(sum(ovrcol), 1), :);
        curcol(~ovrcol, :) = ...
            round(0.5 .* (curcol(~ovrcol, :) + voicol(ones(sum(~ovrcol), 1), :)));
        rgbd(voivox) = curcol(:, 1);
        rgbd(voivox + 16777216) = curcol(:, 2);
        rgbd(voivox + 33554432) = curcol(:, 3);
    end

    % set new data
    rgbo.VoxelDataRGBA = rgbd;

    % and overlay
    slvar.SetHandle('Underlay', rgbo);

    % and mode
    ne_setoption(0, 0, 'joinulay', 4);
    return;
end

% put voxels into VMR
if rad ~= 1
    rad = rad - 0.5;
end
for vc = 1:numel(vois)
    if rad == 1
        slvar.VMRData(bvcoordconv(vois(vc).Voxels, 'tal2bvx', slvar.BoundingBox)) = code;
    else
        voxi = bvcoordconv(vois(vc).Voxels, 'tal2bvc', slvar.BoundingBox);
        from = ceil(-0.5 * rad);
        to = ceil(0.5 * (rad - 1));
        szvmr = size(slvar.VMRData);
        for xc = from:to
            for yc = from:to
                for zc = from:to
                    voxt = max(1, [min(szvmr(1), voxi(:, 1) + xc), ...
                        min(szvmr(2), voxi(:, 2) + yc), ...
                        min(szvmr(3), voxi(:, 3) + zc)]);
                    slvar.VMRData(sub2ind(szvmr, voxt(:, 1), voxt(:, 2), voxt(:, 3))) = code;
                end
            end
        end
    end
end

% redraw
ne_setslicepos;

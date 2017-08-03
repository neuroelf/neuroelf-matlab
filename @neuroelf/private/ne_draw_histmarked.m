% FUNCTION ne_draw_histmarked: histogram of marked voxels
function varargout = ne_draw_histmarked(varargin)

% Version:  v0.9d
% Build:    14071115
% Date:     Jul-11 2014, 3:33 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% check SliceVar
svar = cc.SliceVar;
if numel(svar) ~= 1 || ...
   ~isxff(svar, {'hdr', 'vmr'}) || ...
   ~isfield(svar.RunTimeVars, 'UndoBuffer')
    isempty(svar.RunTimeVars.UndoBuffer)
    return;
end
svartyp = svar.Filetype;
ub = svar.RunTimeVars.UndoBuffer;
[svarpath, svarfile] = fileparts(svar.FilenameOnDisk);
if isempty(svarfile)
    svarfile = sprintf('xff (%s) object ID=%s', svartyp, svar.RunTimeVars.xffID);
else
    svarfile(svarfile == '_') = ' ';
end

% for VMRs
if strcmpi(svartyp, 'vmr')

    % marked data
    if cc.showv16 && ...
       ~isempty(svar.VMRData16)
        ub = svar.RunTimeVars.UndoBuffer16;
        md = lsqueeze(ub(svar.VMRData16 == cc.paint.code(1)));
    else
        md = lsqueeze(ub(svar.VMRData == cc.paint.code(1)));
    end

% for HDRs
else
    if ~isempty(svar.VoxelData)
        vd = svar.VoxelData;
    else
        vd = svar.VoxelDataRGBA;
    end
    md = (vd(:, :, :, 1) == cc.paint.code(1));
    for vc = 2:min(size(vd, 5), numel(cc.paint.code))
        md = md & (vd(:, :, :, vc) == cc.paint.code(vc));
    end

    % marked data
    if size(vd, 5) == 1
        md = ub(md);
    else
        ub = reshape(ub, numel(md), round(numel(ub) / numel(md)));
        md = ub(md, :);
    end
end

% min and max
mm = minmaxmean(md);

% histogram
if (mm(2) - mm(1)) >= 60
    mhs = 1;
elseif (mm(2) - mm(1)) > 30
    mhs = 0.5;
elseif (mm(2) - mm(1)) > 12
    mhs = 0.2;
elseif (mm(2) - mm(1)) > 6
    mhs = 0.1;
elseif (mm(2) - mm(1)) > 2.5
    mhs = 0.05;
else
    mhs = 0.01;
end
mh = histcount(md(:, 1), floor(mm(1)), ceil(mm(2)), mhs);
mh = (1 / sum(mh)) .* mh(:)';

% show histogram
f = figure;
set(f, 'Name', 'NeuroElf - histogram');
ax = axes;
set(ax, 'FontSize', 12);
h = plot(ax, floor(mm(1)):mhs:ceil(mm(2)), mh);
set(h, 'LineWidth', 3, 'Color', [0, 0, 0]);
xlabel(ax, 'Intensities (a.u.)', 'FontSize', 16);
ylabel(ax, 'Relative frequency', 'FontSize', 16);
title(ax, sprintf('Histogram of ''%s'' (%d voxels)', svarfile, size(md, 1)), 'FontSize', 18);

% return handles
if nargout > 0
    varargout{1} = f;
    if nargout > 1
        varargout{2} = ax;
        if nargout > 2
            varargout{3} = h;
        end
    end
end

function hsrf = srf_Show(xo, opts)
% SRF::Show  - show a SRF in a MATLAB axes
%
% FORMAT:       srf.Show([opts]);
%
% Input fields:
%
%      opts         optional settings (struct)
%      .coldata     alternative coloring (either SMP object or Vx1 data)
%
% No output fields.

% Version:  v1.1
% Build:    16031616
% Date:     Mar-16 2016, 4:44 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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

% global config
global xffsngl;

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bgcolor') || ~isa(opts.bgcolor, 'double') || numel(opts.bgcolor) ~= 3
    bgc = [0, 0, 0];
else
    bgc = opts.bgcolor(:)';
end
if ~isfield(opts, 'coldata')
    coldata = [];
else
    coldata = opts.coldata;
end
[sbf{1:3}] = fileparts(xo.F);
bc = xo.C;

% get shortcuts
numv = bc.NrOfVertices;
vcol = bc.VertexColor(:, 1);
tvs = bc.TriangleVertex;
p = bc.VertexCoordinate - repmat(bc.MeshCenter, [numv, 1]);
pn = bc.VertexNormal;

% check coldata
cdvalues = [];
if ~isempty(coldata)

    % xff object?
    if numel(coldata) == 1 && xffisobject(coldata, true, 'smp')
        smpc = coldata.C;
        if smpc.NrOfVertices == numv && ~isempty(smpc.Map)
            coldata = smpc.Map(1);
        else
            error('neuroelf:xff:invalidObject', 'Invalid SMP object for coloring.');
        end

    % Map struct
    elseif isstruct(coldata) && numel(coldata) == 1 && ...
        isfield(coldata, 'SMPData') && size(coldata.SMPData, 1) == numv

    % double data
    elseif isa(coldata, 'double') && size(coldata, 1) == numv
        cda.SMPData = coldata(:);

        coldata = cda;

    % reject other options for the moment
    else
        coldata = [];
    end
end

% check coldata again
if ~isempty(coldata)
    cdvalues = coldata.SMPData(:);

    % check for LowerThreshold and UpperThreshold
    if isfield(coldata, 'LowerThreshold') && isfield(coldata, 'UpperThreshold') && ...
        isfield(coldata, 'RGBLowerThreshPos') && isfield(coldata, 'RGBUpperThreshPos') && ...
        isfield(coldata, 'RGBLowerThreshNeg') && isfield(coldata, 'RGBUpperThreshNeg')

        % generate additional colormap
        hsvposmin = rgb2hsv(coldata.RGBLowerThreshPos / 255);
        hsvposmax = rgb2hsv(coldata.RGBUpperThreshPos / 255);
        hsvnegmin = rgb2hsv(coldata.RGBLowerThreshNeg / 255);
        hsvnegmax = rgb2hsv(coldata.RGBUpperThreshNeg / 255);
        posmap = zeros(11, 3);
        negmap = zeros(11, 3);

        % fill colormaps
        for sc = 1:11
            scf = (sc - 1) / 10;
            posmap(sc, :) = hsv2rgb(hsvposmin * (1 - scf) + hsvposmax * scf);
            negmap(sc, :) = hsv2rgb(hsvnegmin * (1 - scf) + hsvnegmax * scf);
        end

        % set coldata.SMPData to values -1 ... 1
        cdvalues(isinf(cdvalues) | (abs(cdvalues) < coldata.LowerThreshold)) = NaN;
        cdvalues = cdvalues / coldata.UpperThreshold;
        cdvalmax = find(~isnan(cdvalues) & (abs(cdvalues) > 1));
        cdvalues(cdvalmax) = sign(cdvalues(cdvalmax));

    % otherwise
    else

        % scale between -1 and 1 for each side
        cdvalues(isinf(cdvalues)) = NaN;
        cdvok = ~isnan(cdvalues);
        cdvalues(cdvok) = cdvalues(cdvok) / max(abs(cdvalues(cdvok)));

        % build some standard maps
        posmap = hot(11);
        negmap = cool(11);
    end
end

% build palette entries
palette = ones(numv, 1);

% find statistical and base color vertexes
sc = find(isnan(vcol));
ic = find(~isnan(vcol) & (vcol < 2));

% find vertices to hide
hv = ~isnan(vcol) & (vcol > 1);

% remove hidden vertices triangles
keep = ~any(hv(tvs), 2);
tvs(~keep, :) = [];

% set palette of base colors to value
palette(ic) = vcol(ic, 1) + 1;

% find unique colors in all colors
rgbcol = bc.VertexColor(sc, 2:4);
[rgbuni{1:3}] = unique(rgbcol, 'rows');

% extend palette with non-all-zero colors
palette(sc) = rgbuni{3} + 2;

% build colormap
cmapl = 2 + size(rgbuni{1}, 1);
cmap = zeros(cmapl, 3);
cmap(1, :) = bc.ConvexRGBA(1:3);
cmap(2, :) = bc.ConcaveRGBA(1:3);
cmap(3:end, :) = rgbuni{1} / 255;

% add cdvalues colors now
if ~isempty(cdvalues)
    cmap((cmapl + 1):(cmapl + 11), :) = posmap;
    cmap((cmapl + 12):(cmapl + 22), :) = negmap;
    cpal = ~isnan(cdvalues);
    palette(cpal & (cdvalues >= 0)) = (cmapl + 1) + round(cdvalues(cpal & (cdvalues >= 0)) * 10);
    palette(cpal & (cdvalues < 0)) = (cmapl + 12) - round(cdvalues(cpal & (cdvalues < 0)) * 10);
end

% some calculations
mnp = min(p, [], 1);
mxp = max(p, [], 1);

% already shown
if isfield(xo.H, 'Surface') && (isa(xo.H.Surface, 'double') || ...
    isa(xo.H.Surface, 'matlab.graphics.primitive.Patch')) && ...
    numel(xo.H.Surface) == 1 && ishandle(xo.H.Surface) && ...
    strcmpi(get(xo.H.Surface, 'type'), 'patch')

    % then just get handle
    hsrf = xo.H.Surface;

    % and update position
    set(hsrf, 'Vertices', [p(:, 3), -p(:,[1, 2])]);

% otherwise
else

    % turn off HW acceleration ?
    if ispc && ~xffsngl.CONF.settings.OpenGL.HardwareAccelOnWindows
        try
            opengl('software', true);
        catch xfferror
            neuroelf_lasterr(xfferror);
            warning('neuroelf:xff:openGLError', ...
                'Error setting OpenGL to SoftwareAcceleration for PCs.');
        end
    end

    % create figure
    hf = figure;
    set(hf, 'Units', 'pixels', 'Color', bgc, 'Tag', 'SURFACE');
    set(hf, 'Position', [200, 200, 400, 400]);
    set(hf, 'Units', 'normalized');

    % Create 2D plot with no relevant content
    ha = axes('position', [0 0 1 1]);
    set(ha, 'Units', 'normalized');
    hold(ha, 'on');

    % create surface
    hsrf = trisurf(tvs, p(:,3), -p(:,1), -p(:,2), palette);

    % set coloring and limits
    set(hsrf, 'CData', palette);
    set(hsrf, 'CDataMapping', 'direct');
    set(hsrf, 'EdgeColor', 'none');
    set(hsrf, 'Faces', tvs(:, [1, 3, 2]))
    set(hsrf, 'LineStyle', 'none');
    set(hsrf, 'FaceColor', 'interp');
    set(hsrf, 'FaceLighting', 'gouraud');
    set(ha, 'YDir', 'reverse');
    if ~isempty(strfind(lower(sbf{2}), '_rh'))
        set(ha, 'View', [-90, 0]);
    else
        set(ha, 'View', [90, 0]);
    end
    light('Parent', ha, 'Position', [ 0,  1, -0.5], 'Color', [0.5, 0.5, 0.5]);
    light('Parent', ha, 'Position', [ 0, -1, -0.5], 'Color', [0.5, 0.5, 0.5]);
    light('Parent', ha, 'Position', [ 1,  0,  0.5], 'Color', [0.7, 0.7, 0.7]);
    light('Parent', ha, 'Position', [-1,  0,  0.5], 'Color', [0.7, 0.7, 0.7]);
    set(hf, 'Name', sprintf('Surface rendering: %s', sbf{2}));
    set(hf, 'NumberTitle', 'off');
    set(ha, 'Color', bgc, 'XColor', bgc, 'YColor', bgc, 'ZColor', bgc);
    set(hf, 'Color', bgc);
    colormap(ha, cmap);
    set(ha, 'XLim', [min(mnp(1), -144) , max(mxp(1), 144)]);
    set(ha, 'YLim', [min(mnp(2), -144) , max(mxp(2), 144)]);
    set(ha, 'ZLim', [min(mnp(3), -144) , max(mxp(3), 144)]);
    set(ha, 'XTick', []);
    set(ha, 'YTick', []);
    set(ha, 'ZTick', []);

    % disable hold
    hold(ha, 'off');
end

% (re-)set coloring
if ~isempty(coldata)
    set(hsrf, 'CData', palette);
    colormap(get(hsrf, 'Parent'), cmap);
end

% set/update normals
set(hsrf, 'VertexNormals', [-pn(:, 3), pn(:,[1, 2])]);

% make sure the update is done
pause(0.001);

% store handle
xo.H.Surface = hsrf;

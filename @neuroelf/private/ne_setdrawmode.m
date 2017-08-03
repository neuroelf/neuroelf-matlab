function varargout = ne_setdrawmode(varargin)
% ne_setdrawmode  - set drawing mode and perform actions
%
% FORMAT:       ne_setdrawmode(SRC, EVT, action [, subaction [, config]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       action      0x0 or 1x1 double value as action code:
%                   empty: follow subaction (see below)
%                    1: set to browsing mode (no drawing)
%                    2: set to 2D drawing mode
%                    3: set to 3D drawing mode
%                   -1: set undo mode (read from button state)
%       subaction   either a 1x1 double (empty action) or further config
%                   with empty action (subaction):
%                    0: mask (Buffer) with Data == ne_gcfg.fcfg.paint.code
%                    0.5: inverse mask with Data == ne_gcfg.fcfg.paint.code
%                    1: accept changes (copy data to UndoBuffer)
%                   -1: revert changes (copy data from UndoBuffer)
%                   -2: reload data from disk (and set UndoBuffer to Data)
%                   -2.5: sample UndoBuffer from Underlay object
%                   -3: flood-fill from ne_gcfg.fcfg.cpos coordinate
%                   -4: expand marked color
%                   -5: smooth border
%                   -6: smooth data
%                   with action 2 or 3 (configure drawing mode):
%                   1x5 cell array (all must be given) with strings for (!)
%                   {RADIUS, COLOR_CODES, DRAW_SHAPE, OVER_RANGE, SMOOTHNESS}
%                   whereas RADIUS can be 1x1 or 1x3 in millimeters
%                   COLOR_CODES can be 1x1 with >= 0 being written into
%                   the data, or < -4, in which case they will be added
%                   (i.e. the absolute value subtracted) from the data
%                   for [-2 .. 0[ multiply current data with (positive) value
%                   for [-4 .. -2[ multiply original data with (positive) value
%       config      optionally given cell array with string (!) values for
%                   flood-fill: 1x3 {COLOR_CODE, FILL_OVER_RANGE, CONNTYPE}
%                   expand: 1x3 {COLOR_CODE, FILL_INTO_RANGE, CONNTYPE}
%                   smooth-border: 1x3 {COLOR_CODE, SMOOTH_KERNEL, THRESH}
%                   smooth: 1x3 {SMOOTH_KERNEL, BOUNDING_BOX, VALUE_RANGE}
%
% Examples:
%
% % create an empty VMR
% vmr = xff('new:vmr');
% vmr.Browse;
%
% % draw a pyramid
% for ds = 30:-1:1
%     ne_setdrawmode(0, 0, 2, {sprintf('%d', ds), '240', 'c', '0 255', '0'});
%     neuroelf_gui('draw', [0, 0, -2 * ds + 24], [1, 3]);
%     neuroelf_gui('draw', [0, 0, -2 * ds + 25], [1, 3]);
% end
%
% % disable drawing
% ne_setdrawmode(0, 0, 1);

% Version:  v1.0
% Build:    16010821
% Date:     Jan-08 2016, 9:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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

% get current VMR and return for non-VMR objects
svar = cc.SliceVar;
if numel(svar) ~= 1 || ...
   (~isxff(svar, 'vmr') && ...
    ~isxff(svar, 'hdr')) || ...
    nargin < 3
    return;
end
svarh = handles(svar);
if isfield(svarh, 'Underlay') && ...
    numel(svarh.Underlay) == 1 && ...
    isxff(svarh.Underlay, {'hdr', 'head', 'msk', 'vmr', 'vtc'})
    uvar = svarh.Underlay;
else
    uvar = [];
end
svart = lower(svar.Filetype);

% for VMRs
if strcmp(svart, 'vmr')

    % and transformation values (if not == eye(4))
    svar.RunTimeVars.Trf = ...
        bvcoordconv(zeros(0, 3), 'tal2bvc', svar.BoundingBox)';
    if ~isequal(cc.strans, eye(4))
        trans = cc.strans;
    else
        trans = [];
    end

    % use Trf to get voxel position
    rtv = svar.RunTimeVars;
    if isfield(rtv, 'Trf')
        ttrf = rtv.Trf;
    else
        ttrf = eye(4);
    end
    if isfield(rtv, 'TrfPlus')
        tplus = inv(rtv.TrfPlus)';
    else
        tplus = eye(4);
    end
    if ~isempty(trans)
        bvpos = [cc.cpos, 1] * tplus * trans' * ttrf;
    else
        bvpos = [cc.cpos, 1] * tplus * ttrf;
    end
    bvpos(end) = [];
end

% make sure the Object has an UndoBuffer
if ~isfield(svar.RunTimeVars, 'UndoBuffer')
    if ~strcmp(svart, 'hdr')
        svar.RunTimeVars.UndoBuffer = svar.GetVolume(1);
    else
        if ~isempty(svar.VoxelData)
            svar.RunTimeVars.UndoBuffer = svar.VoxelData(:, :, :, :);
        elseif ~isempty(svar.VoxelDataRGBA)
            svar.RunTimeVars.UndoBuffer = svar.VoxelDataRGBA(:, :, :, :, :);
        end
    end
end

% empty input
if isempty(varargin{3})

    % accept (1)
    if varargin{4} == 1
        if ~strcmp(svart, 'hdr')
            svar.RunTimeVars.UndoBuffer = svar.GetVolume(1);
        else
            if ~isempty(svar.VoxelData)
                svar.RunTimeVars.UndoBuffer = svar.VoxelData(:, :, :, :);
            elseif ~isempty(svar.VoxelDataRGBA)
                svar.RunTimeVars.UndoBuffer = svar.VoxelDataRGBA(:, :, :, :, :);
            end
        end

    % revert changes (-1)
    elseif varargin{4} == -1
        if strcmp(svart, 'hdr')
            svdhn = svar.DataHist.NIftI1;
        end
        rtv = svar.RunTimeVars;
        try
            ch.MainFig.Pointer = 'watch';
            drawnow;
            svar.ReloadFromDisk;
            svar.RunTimeVars = rtv;
            if ~strcmp(svart, 'hdr')
                ub = svar.GetVolume(1);
                svar.VMRData = rtv.UndoBuffer;
            else
                if ~isempty(svar.VoxelData)
                    ub = svar.VoxelData(:, :, :, :);
                    svar.VoxelData = rtv.UndoBuffer;
                elseif ~isempty(svar.VoxelDataRGBA)
                    ub = svar.VoxelDataRGBA(:, :, :, :, :);
                    svar.VoxelDataRGBA = rtv.UndoBuffer;
                else
                    ub = [];
                end
                svar.DataHist.NIftI1 = svdhn;
            end
            rtv.UndoBuffer = ub;
            svar.RunTimeVars = rtv;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            if ~strcmp(svart, 'hdr')
                svar.VMRData = rtv.UndoBuffer;
            else
                if ~isempty(svar.VoxelData)
                    svar.VoxelData = rtv.UndoBuffer;
                elseif ~isempty(svar.VoxelDataRGBA)
                    svar.VoxelDataRGBA = rtv.UndoBuffer;
                end
            end
        end

    % reload all (-2)
    elseif varargin{4} == -2 && ...
       ~isempty(svar.FilenameOnDisk)
        rtv = svar.RunTimeVars;
        ch.MainFig.Pointer = 'watch';
        drawnow;
        svar.ReloadFromDisk;
        svar.RunTimeVars = rtv;
        if ~strcmp(svart, 'hdr')
            svar.RunTimeVars.UndoBuffer = svar.GetVolume(1);
        else
            if ~isempty(svar.VoxelData)
                svar.RunTimeVars.UndoBuffer = svar.VoxelData(:, :, :, :);
            elseif ~isempty(svar.VoxelDataRGBA)
                svar.RunTimeVars.UndoBuffer = svar.VoxelDataRGBA(:, :, :, :, :);
            end
        end

    % sample UndoBuffer from Underlay object (-2.5)
    elseif varargin{4} == -2.5 && ...
       ~isempty(uvar)

        % depending on object type
        ch.MainFig.Pointer = 'watch';
        drawnow;
        if ~strcmp(svart, 'hdr')
            bbox = svar.BoundingBox;
            if any(bbox.ResXYZ == 1)
                bbox.BBox(2, bbox.ResXYZ == 1) = bbox.BBox(2, bbox.ResXYZ == 1) + 1;
            end
            ubuff = uvar.SampleBVBox(bbox);
            if isxff(uvar, 'vmr')
                ubuff = uint8(round(ubuff));
            else
                ubuffmm = minmaxmean(ubuff);
                ubuffmd = (ubuffmm(2) - ubuffmm(1)) + sqrt(eps);
                ubuff = uint8(round((225 / ubuffmd) .* (ubuff - ubuffmm(1))));
            end
            svar.RunTimeVars.UndoBuffer = ubuff;
        else
        end

    % fill with code (floodfill, -3)
    elseif varargin{4} == -3
        if nargin < 5 || ...
           ~iscell(varargin{5}) || ...
            numel(varargin{5}) ~= 3
            pconf = inputdlg( ...
                {'Fill with color code:', 'Fill over color code range:', ...
                 'Connectivity type: (f)ace, (e)dge, or (v)ertex'}, ...
                'NeuroElf GUI - floodfill3 configuration', 1, ...
                {sprintf('  %g', cc.paint.code), ...
                 sprintf('  %g', cc.paint.over), '  f'});
        else
            pconf = varargin{5};
        end
        try
            pconf{1}(pconf{1} == ' ') = [];
            pconf{3}(pconf{3} == ' ') = [];
            if numel(pconf{3}) ~= 1 || ...
               ~any('efv' == lower(pconf{3}))
                return;
            end
            pconf{1} = str2double(pconf{1});
            pconf{2} = eval(['[' pconf{2} ']']);
            pconf{3} = lower(pconf{3});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        if strcmp(svart, 'vmr')
            flpos = round(bvcoordconv(cc.cpos, 'tal2bvc', svar.BoundingBox));
            if all(flpos > 0 & flpos <= size(svar.VMRData)) && ...
                pconf{1} >= 0 && ...
                pconf{1} <= 255 && ...
               ~isempty(pconf{2})
                ch.MainFig.Pointer = 'watch';
                drawnow;
                vmrmsk = (svar.VMRData >= min(pconf{2}) & svar.VMRData <= max(pconf{2}));
                [vmrmsk, flooded] = ...
                    floodfill3(vmrmsk, flpos(1), flpos(2), flpos(3), pconf{3});
                if flooded > 0
                    svar.VMRData(vmrmsk) = round(pconf{1});
                end
            end
        elseif strcmp(svart, 'hdr')
            flpos = round([cc.cpos, 1] * svar.RunTimeVars.Trf);
            flpos(end) = [];
            if all(flpos > 0 & flpos <= svar.ImgDim.Dim(2:4)) && ...
               ~isempty(pconf{2})
                ch.MainFig.Pointer = 'watch';
                drawnow;
                if ~isempty(svar.VoxelData)
                    vmrmsk = (svar.VoxelData(:, :, :) >= min(pconf{2}) & ...
                        svar.VoxelData(:, :, :) <= max(pconf{2}));
                    [vmrmsk, flooded] = ...
                        floodfill3(vmrmsk, flpos(1), flpos(2), flpos(3), pconf{3});
                    if flooded > 0
                        svar.VoxelData(vmrmsk) = round(pconf{1});
                    end
                end
            end
        else
            return;
        end
        ne_gcfg.fcfg.paint.code = pconf{1};
        ne_gcfg.fcfg.paint.over = [min(pconf{2}), max(pconf{2})];

    % expand marking (dilate3d, -4)
    elseif varargin{4} == -4
        if nargin < 5 || ...
           ~iscell(varargin{5}) || ...
            numel(varargin{5}) ~= 3
            pconf = inputdlg( ...
                {'Expand color code:', 'Fill over color code range:', ...
                 'Connectivity type: (f)ace, (e)dge, or (v)ertex'}, ...
                'NeuroElf GUI - dilate3d configuration', 1, ...
                {sprintf('  %g', cc.paint.code), ...
                 sprintf('  %g', cc.paint.over), '  f'});
        else
            pconf = varargin{5};
        end
        try
            pconf{1}(pconf{1} == ' ') = [];
            pconf{3}(pconf{3} == ' ') = [];
            if numel(pconf{3}) ~= 1 || ...
               ~any('efv' == lower(pconf{3}))
                return;
            end
            pconf{1} = str2double(pconf{1});
            pconf{2} = eval(['[' pconf{2} ']']);
            pconf{3} = lower(pconf{3});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        if pconf{1} >= 0 && ...
            pconf{1} <= 255 && ...
            ~isempty(pconf{2})
            dop = false(3, 3, 3);
            if pconf{3} == 'f'
                dop(2, 2, :) = true;
                dop(2, :, 2) = true;
                dop(:, 2, 2) = true;
            elseif pconf{3} == 'e'
                dop(2, :, :) = true;
                dop(:, 2, :) = true;
                dop(:, :, 2) = true;
            else
                dop(:) = true;
            end
            if strcmp(svart, 'vmr')
                ch.MainFig.Pointer = 'watch';
                drawnow;
                vmrmsk = dilate3d(svar.VMRData == round(pconf{1}), dop);
                vmrmsk = vmrmsk & (svar.VMRData >= min(pconf{2})) & ...
                    (svar.VMRData <= max(pconf{2})) & ...
                    (svar.VMRData ~= round(pconf{1}));
                svar.VMRData(vmrmsk) = round(pconf{1});
            elseif strcmp(svart, 'hdr')
                ch.MainFig.Pointer = 'watch';
                drawnow;
            else
                return;
            end
            ne_gcfg.fcfg.paint.code = pconf{1};
            ne_gcfg.fcfg.paint.over = [min(pconf{2}), max(pconf{2})];
        end

    % smooth border (-5)
    elseif varargin{4} == -5
        if nargin < 5 || ...
           ~iscell(varargin{5}) || ...
            numel(varargin{5}) ~= 3
            pconf = inputdlg( ...
                {'Smooth border of color code:', 'Smoothing kernel:', ...
                 'Relative threshold:'}, ...
                'NeuroElf GUI - smooth3 configuration', 1, ...
                {sprintf('  %g', cc.paint.code), '2', '0.5'});
        else
            pconf = varargin{5};
        end
        try
            pconf{1}(pconf{1} == ' ') = [];
            pconf{2}(pconf{2} == ' ') = [];
            pconf{3}(pconf{3} == ' ') = [];
            pconf{1} = str2double(pconf{1});
            pconf{2} = str2double(pconf{2});
            pconf{3} = str2double(pconf{3});
            if pconf{1} < 0 || ...
                pconf{1} > 32767 || ...
                pconf{2} < 0.5 || ...
                pconf{2} > 5 || ...
                pconf{3} < 0 || ...
                pconf{3} > 1
                return;
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        if strcmp(svart, 'vmr')
            ch.MainFig.Pointer = 'watch';
            drawnow;
            vmrmsk = (svar.VMRData == round(pconf{1}));
            newmsk = (smoothdata3(double(vmrmsk), pconf{2} .* ones(1, 3)) >= pconf{3});
            svar.VMRData(newmsk) = round(pconf{1});
            vmrmsk = vmrmsk & ~newmsk;
            svar.VMRData(vmrmsk) = svar.RunTimeVars.UndoBuffer(vmrmsk);
        elseif strcmp(svart, 'hdr')
            ch.MainFig.Pointer = 'watch';
            drawnow;
        else
            return;
        end
        ne_gcfg.fcfg.paint.code = pconf{1};

    % smooth data (-6)
    elseif varargin{4} == -6
        if nargin < 5 || ...
           ~iscell(varargin{5}) || ...
            numel(varargin{5}) ~= 3
            pconf = inputdlg( ...
                {'Smoothing kernel:', 'Bounding box', 'Value range:'}, ...
                'NeuroElf GUI - smooth VMR configuration', 1, {'  2', ...
                 sprintf('  [%d %d %d; %d %d %d]', lsqueeze(cc.paint.bbox')), ...
                 sprintf('  %g', cc.paint.over)});
        else
            pconf = varargin{5};
        end
        try
            pconf{1}(pconf{1} == ' ') = [];
            pconf{1} = str2double(pconf{1});
            pconf{2} = eval(pconf{2});
            pconf{3} = eval(['[' pconf{3} ']']);
            if pconf{1} <= 0 || ...
                pconf{1} > 8 || ...
               ~isa(pconf{2}, 'double') || ...
                numel(pconf{2}) ~= 6 || ...
                any(pconf{2}(:) < -256 | pconf{2}(:) > 256) || ...
               ~isa(pconf{3}, 'double') || ...
                numel(pconf{3}) ~= 2 || ...
                any(pconf{3} < 0 | pconf{3} > 32767) || ...
                pconf{3}(2) < pconf{3}(1)
                return;
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        pconf{2} = reshape(round(pconf{2}), 2, 3);
        smopt = struct('range', pconf{3});
        if strcmp(svart, 'vmr') && ...
            isa(svar.VMRData, 'uint8')
            smopt.nanrange = [226, 255];
        end
        ch.MainFig.Pointer = 'watch';
        drawnow;
        svar.SmoothData3D(pconf{1}, [Inf, Inf, Inf; pconf{2}(1, :); ...
            1, 1, 1; pconf{2}(2, :)], smopt);
        ne_gcfg.fcfg.paint.bbox = pconf{2};
        ne_gcfg.fcfg.paint.over = pconf{3};

    % mask (0)
    elseif varargin{4} == 0
        vd = svar.GetVolume(1);
        if size(vd, 5) == 1
            vdm = (vd ~= cc.paint.code);
        else
            codes = cc.paint.code(:)';
            if numel(codes) < size(size(vd, 5))
                codes = repmat(codes, 1, ceil(size(vd, 5) / numel(codes)));
            end
            vdm = (vd(:, :, :, 1, 1) ~= codes(1));
            for s5c = 2:size(vd, 5)
                vdm = vdm | (vd(:, :, :, 1, s5c) ~= codes(s5c));
            end
            vdm = repmat(vdm, [1, 1, 1, 1, s5c]);
        end
        rtv = svar.RunTimeVars;
        ch.MainFig.Pointer = 'watch';
        drawnow;
        svar.ReloadFromDisk;
        svar.LoadTransIOData;
        svar.RunTimeVars = rtv;
        if strcmp(svart, 'vmr')
            svar.VMRData(vdm) = 0;
            svar.RunTimeVars.UndoBuffer = svar.VMRData(:, :, :);
        elseif strcmp(svart, 'hdr')
            if ~isempty(svar.VoxelData)
                svar.VoxelData(vdm) = 0;
                svar.RunTimeVars.UndoBuffer = svar.VoxelData(:, :, :, :);
            elseif ~isempty(svar.VoxelDataRGBA)
                svar.VoxelDataRGBA(vdm) = 0;
                svar.RunTimeVars.UndoBuffer = svar.VoxelDataRGBA(:, :, :, :, :);
            end
        else
            ch.MainFig.Pointer = 'arrow';
            drawnow;
            return;
        end

    % inverse mask (0.5)
    elseif varargin{4} == 0.5
        vd = svar.GetVolume(1);
        if size(vd, 5) == 1
            vdm = (vd == cc.paint.code);
        else
            codes = cc.paint.code(:)';
            if numel(codes) < size(size(vd, 5))
                codes = repmat(codes, 1, ceil(size(vd, 5) / numel(codes)));
            end
            vdm = (vd(:, :, :, 1, 1) == codes(1));
            for s5c = 2:size(vd, 5)
                vdm = vdm & (vd(:, :, :, 1, s5c) == codes(s5c));
            end
            vdm = repmat(vdm, [1, 1, 1, 1, s5c]);
        end
        rtv = svar.RunTimeVars;
        ch.MainFig.Pointer = 'watch';
        drawnow;
        svar.ReloadFromDisk;
        svar.RunTimeVars = rtv;
        if strcmp(svart, 'vmr')
            svar.VMRData(vdm) = 0;
            svar.RunTimeVars.UndoBuffer = svar.VMRData(:, :, :);
        elseif strcmp(svart, 'hdr')
            if ~isempty(svar.VoxelData)
                svar.VoxelData(vdm) = 0;
                svar.RunTimeVars.UndoBuffer = svar.VoxelData(:, :, :, :);
            elseif ~isempty(svar.VoxelDataRGBA)
                svar.VoxelDataRGBA(vdm) = 0;
                svar.RunTimeVars.UndoBuffer = svar.VoxelDataRGBA(:, :, :, :, :);
            end
        else
            ch.MainFig.Pointer = 'arrow';
            drawnow;
            return;
        end
    end

    % update display (without painting again!)
    dm = ne_gcfg.fcfg.paint.mode;
    ne_gcfg.fcfg.paint.mode = 0;
    ch.MainFig.Pointer = 'arrow';
    ne_setslicepos;
    ne_gcfg.fcfg.paint.mode = dm;

% set drawing mode
elseif varargin{3} > 0

    % ask for radius?
    if varargin{3} > 1
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 5
            pconf = inputdlg({ ...
                'Please give a radius size:', ...
                'Drawing color code:', ...
                'Use (c)uboid or (s)pheroid shape:', ...
                'Draw over color code range:', ...
                'Drawing smoothness (0 - 1, spheroid only)'}, ...
                'NeuroElf GUI - paint configuration', ...
                1, { ...
                sprintf('  %g', cc.paint.rad), ...
                sprintf('  %g', cc.paint.code), ...
                sprintf('  %c', cc.paint.shape), ...
                sprintf('  %g', cc.paint.over), ...
                sprintf('  %0.3f', cc.paint.smooth)});
            if isempty(pconf) || ...
               ~iscell(pconf)
                return;
            end
        else
            pconf = varargin{4};
        end
        try
            pconf{1} = ddeblank(pconf{1});
            pconf{3}(pconf{3} == ' ') = [];
            if numel(pconf{3}) ~= 1 || ...
               ~any('cs' == lower(pconf{3}))
                return;
            end
            pconf{1} = eval(['[' pconf{1} ']']);
            if ~isa(pconf{1}, 'double') || ...
               ~any(numel(pconf{1}) == [1, 3]) || ...
                any(isinf(pconf{1}) | isnan(pconf{1}) | pconf{1} < 0)
                return;
            end
            pconf{2} = eval(['[' pconf{2} ']']);
            if ~isa(pconf{2}, 'double') || ...
                isempty(pconf{2})
                return;
            end
            pconf{2} = pconf{2}(:)';
            pconf{3} = lower(pconf{3});
            pconf{4} = eval(['[' pconf{4} ']']);
            if numel(pconf{4}) ~= 2 || ...
                any(isinf(pconf{4}) | isnan(pconf{4}) | pconf{4} < 0 | pconf{4} > 32767) || ...
                pconf{4}(2) < pconf{4}(1)
                return;
            end
            pconf{5} = str2double(pconf{5});
            if isinf(pconf{5}) || ...
                isnan(pconf{5}) || ...
                pconf{5} < 0 || ...
                pconf{5} > 1
                return;
            end
            pconf{5} = 0.001 * round(1000 * pconf{5});
            if strcmp(svart, 'vmr')
                pconf{4} = fix(pconf{4}(:)');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % set unto toggle button state
        ch.DrawUndo.Value = double(cc.paint.mode < 0);

    % show undo toggle button as un-set
    else
        ch.DrawUndo.Value = 0;

        % get radius/shape config for setting
        pconf = {cc.paint.rad, cc.paint.code, cc.paint.shape, cc.paint.over, cc.paint.smooth};
    end

    % set mode
    rad = pconf{1};
    if ~any(numel(rad) == [1, 3]) || ...
        any(isinf(rad) | isnan(rad) | rad < 0)
        varargin{3} = 1;
        rad = 0;
    end
    cod = pconf{2};
    if any(isinf(cod) | isnan(cod)) || ...
       (strcmp(svart, 'vmr') && ...
        numel(cod) ~= 1)
        cod = 240;
    elseif strcmp(svart, 'vmr') && ...
       (numel(cod) ~= 1 || ...
        cod > 0)
        cod = fix(cod(1));
    elseif strcmp(svart, 'vmr')
        cod = cod(1);
    end
    smv = pconf{5};
    if isinf(smv) || ...
        isnan(smv) || ...
        smv < 0 || ...
        smv > 1
        smv = 0;
    end
    ne_gcfg.fcfg.paint.mode = sign(ne_gcfg.fcfg.paint.mode) * varargin{3};
    ne_gcfg.fcfg.paint.rad = rad;
    ne_gcfg.fcfg.paint.code = cod;
    ne_gcfg.fcfg.paint.shape = pconf{3};
    ne_gcfg.fcfg.paint.over = pconf{4};
    ne_gcfg.fcfg.paint.smooth = smv;

    % smooth drawing
    if pconf{3} == 's' && ...
        smv > 0

        % large smoothing kernel
        smk = smoothkern(1000);

        % half
        smk(1:floor(0.5 * numel(smk))) = [];

        % fix center to 1
        smk = smk ./ smk(1);
        smk(1) = 1;

        % interpolate based on requirements
        if smv < 1
            smk = [ones(round(1000 * (1 - smv)), 1); ...
                flexinterpn(smk, [Inf; 1; 1 / smv; 1001])];
        else
            smk(1002:end) = [];
        end
        ne_gcfg.fcfg.paint.smootk = smk;
    else
        ne_gcfg.fcfg.paint.smootk = ones(1001, 1);
    end

    % generate shape
    if varargin{3} == 2
        rad = max(ceil(rad));
        [c1, c2] = ndgrid(-rad:rad, -rad:rad);
        c1 = [c1(:), c2(:)];
        if pconf{3} == 's'
            c1d = sqrt(sum(c1 .* c1, 2));
            c1(c1d > pconf{1}, :) = [];
        end
        ne_gcfg.fcfg.paint.shap2 = c1;
        if smv == 0 || ...
            pconf{3} ~= 's'
            ne_gcfg.fcfg.paint.shap2w = 1;
        else
            c1d(c1d > pconf{1}) = [];
            ne_gcfg.fcfg.paint.shap2w = flexinterpn(smk, 1 + 1000 .* (c1d(:) ./ rad));
        end
    elseif varargin{3} == 3
        if strcmp(svart, 'vmr')
            if numel(rad) == 1
                rad = rad ./ svar.BoundingBox.ResXYZ;
            else
                rad = lsqueeze(rad([2, 3, 1]))' ./ svar.BoundingBox.ResXYZ;
            end
        end
        crad = ceil(rad);
        if numel(rad) == 1
            [c1, c2, c3] = ndgrid(-crad:crad, -crad:crad, -crad:crad);
        else
            [c1, c2, c3] = ndgrid(-crad(1):crad(1), -crad(2):crad(2), -crad(3):crad(3));
        end
        c1 = [c1(:), c2(:), c3(:)];
        if pconf{3} == 's'
            if numel(rad) == 1
                c1d = sqrt(sum(c1 .* c1, 2));
                c1(c1d > pconf{1}, :) = [];
            else
                c1r = c1 * diag(1 ./ rad(:));
                c1d = sqrt(sum(c1r .* c1r, 2));
                c1(c1d > 1, :) = [];
            end
        end
        ne_gcfg.fcfg.paint.shap3 = c1;
        if smv == 0 || ...
            pconf{3} ~= 's'
            ne_gcfg.fcfg.paint.shap3w = 1;
        elseif numel(rad) == 1
            c1d(c1d > pconf{1}) = [];
            ne_gcfg.fcfg.paint.shap3w = flexinterpn(smk, 1 + 1000 .* (c1d(:) ./ rad));
        else
            c1d(c1d > 1) = [];
            ne_gcfg.fcfg.paint.shap3w = flexinterpn(smk, 1 + 1000 .* c1d(:));
        end
    end

% set undo mode
elseif varargin{3} < 0
    s = 1 - 2 * ch.DrawUndo.Value;
    ne_gcfg.fcfg.paint.mode = s * abs(ne_gcfg.fcfg.paint.mode);
end

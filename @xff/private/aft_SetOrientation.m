function mat44 = aft_SetOrientation(xo, opts)
% AFT::SetOrientation  - set orientation/origin of file
%
% FORMAT:       [mat44 = ] obj.SetOrientation([opts])
%
% Input fields:
%
%       opts        1x1 struct with options (or 1x1 TRF object of Type 2)
%        .updhdr    flag, update HDR/NII file with origin (default: false)
%
% Output fields:
%
%       mat44       4x4 quaternion, result of re-orientation
%
% TYPES: AVA, CMP, GLM, HDR, HEAD, MGH, MSK, NLF, SRF, TVL, VMP, VMR
%
% Note: this updates the .RunTimeVars.Trf, and if updhdr is true saves
%       changes to disk
%
% Warning: if the .updhdr flag is true and the file is an Analyze 7.5 file
%          (without 'ni1' token in the header, which does not support
%          storing a 4x4 quaternion matrix with offset values, the file
%          will be converted into a compatible 'ni1' magic file !!!
%
% Using: mtimesnd, spmitrf, spmtrf.

% Version:  v1.1
% Build:    16031615
% Date:     Mar-16 2016, 3:45 PM EST
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

% global variable (used to access in sub-functions) and library
global so_cfg ne_methods;
if numel(so_cfg) ~= 1 || ~isstruct(so_cfg)

    % load figure
    try
        hFig = xfigure([neuroelf_path('tfg') '/setorient.tfg']);
    catch xfferror
        error('neuroelf:xff:xfigureError', 'Error opening the required figure: %s.', xfferror.message);
    end
    so_cfg = struct;
    so_cfg.hFig = hFig;
    so_cfg.tags = hFig.TagStruct;
end

% check figure
hFig = so_cfg.hFig;
if ~isxfigure(hFig, true)

    % reload
    try
        hFig = xfigure([neuroelf_path('tfg') '/setorient.tfg']);
    catch xfferror
        error('neuroelf:xff:xfigureError', 'Error re-opening the required figure: %s.', xfferror.message);
    end
    so_cfg = struct;
    so_cfg.hFig = hFig;
    so_cfg.tags = hFig.TagStruct;
end
tags = so_cfg.tags;

% preset output
mat44 = eye(4);

% argument check
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, ...
    {'ava', 'cmp', 'glm', 'hdr', 'head', 'mgh', 'msk', 'nlf', 'srf', 'vmp', 'vmr'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% special case: TRF object
if nargin == 2 && numel(opts) == 1 && xffisobject(opts, true, 'trf')
    trfc = opts.C;

    % check content
    if ischar(trfc.DataFormat) && strcmpi(trfc.DataFormat, 'matrix') && ...
        isequal(size(trfc.TFMatrix), [4, 4]) && isequal(trfc.TransformationType, 2) && ...
        isequal(trfc.CoordinateSystem, 0) && isequal(trfc.ExtraVMRTransf, 0)

        % get content
        mat44 = trfc.TFMatrix;

        % adapt
        mat44 = mat44([3, 1, 2, 4], :);
        mat44 = mat44(:, [3, 1, 2, 4]);
        mat44(1:3, 1:3) = pinv(mat44(1:3, 1:3));

        % make sure fields are not minimally different from 0
        mat44(abs(mat44) < 1e-12) = 0;

        % final check
        if ~any(isinf(mat44(:)) | isnan(mat44(:))) && all(mat44(4, :) == [0, 0, 0, 1])

            % replace in RunTimeVars
            bc.RunTimeVars.TrfPlus = mat44;
            xo.C = bc;
            return;
        end
    end
end

% regular inputs
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'updhdr') || ~islogical(opts.updhdr) || numel(opts.updhdr) ~= 1
    opts.updhdr = false;
end

% de-compose current matrix
cm = ne_methods.spmitrf(bc.RunTimeVars.TrfPlus);

% don't allow shears at this time
if any(abs(cm{4}) > sqrt(eps))
    error('neuroelf:xff:badArgument', 'Shears are currently not supported. Please manually correct.');
end
cm{2} = (180 / pi) .* cm{2};

% cross hair display
so_cfg.chair = true;

% keep track of current settings
so_cfg.cm = cm;
so_cfg.cmb = cm;

% order in which to toggle directions
so_cfg.dirorder = {'sag', 'cor', 'tra'};

% sampling frame
so_cfg.sframe = [128, 128, 128; -127.999, -127.999, -127.999];

% key modifiers
so_cfg.mods = {};

% currently displayed page (xfigure property of figure object)
so_cfg.page = 1;

% position of 3-slice images
so_cfg.slicepos = [tags.IM_setorient_Slice_SAG.Position; ...
    tags.IM_setorient_Slice_COR.Position; tags.IM_setorient_Slice_TRA.Position];
so_cfg.slicepos(:, 3:4) = so_cfg.slicepos(:, 1:2) + so_cfg.slicepos(:, 3:4);

% sampling stepsize (default: 1mm)
so_cfg.sstep = 1;

% update flag
so_cfg.update = false;
so_cfg.updhdr = opts.updhdr;
tags.CB_setorient_updhdr.Value = double(so_cfg.updhdr);

% zoom flag (1, 2, 3 as different zoomed slice)
so_cfg.zoom = 0;

% position of zoomed slices (times 3, for three "objects")
so_cfg.zslicepos = tags.IM_setorient_Slice_Zoom.Position;
so_cfg.zslicepos(3:4) = so_cfg.zslicepos(1:2) + so_cfg.zslicepos(3:4);
so_cfg.zslicepos = so_cfg.zslicepos([1, 1, 1], :);

% initialiaze SliceVar and StatsVar to empty
so_cfg.SliceVar = xo;

% create layered image objects
so_cfg.tio = struct('imSag', transimg(256, 256), 'imCor', transimg(256, 256), ...
    'imTra', transimg(256, 256), 'imSlZ', transimg(512, 512));

% and set handles of images into transimg object (to allow use of display)
sethandle(so_cfg.tio.imSag, get(tags.IM_setorient_Slice_SAG.MLHandle, 'Children'));
sethandle(so_cfg.tio.imCor, get(tags.IM_setorient_Slice_COR.MLHandle, 'Children'));
sethandle(so_cfg.tio.imTra, get(tags.IM_setorient_Slice_TRA.MLHandle, 'Children'));
sethandle(so_cfg.tio.imSlZ, get(tags.IM_setorient_Slice_Zoom.MLHandle, 'Children'));

% main figure (xfigure and MLHandle)
so_cfg.hFig = hFig;
so_cfg.hFigMLH = hFig.MLHandle;
so_cfg.tags = tags;

% add crosshair lines to images axes objects -> SAG
chax = tags.AX_setorient_Slice_SAG.MLHandle;
set(chax, 'Units', 'pixels');
so_cfg.SagLineX = line([0; 0.999], [0.5; 0.5], 'Color', [0, 0, 1], 'Parent', chax);
so_cfg.SagLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', [0, 0, 1], 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');

% -> COR
chax = tags.AX_setorient_Slice_COR.MLHandle;
set(chax, 'Units', 'pixels');
so_cfg.CorLineX = line([0; 0.999], [0.5; 0.5], 'Color', [0, 0, 1], 'Parent', chax);
so_cfg.CorLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', [0, 0, 1], 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');

% -> TRA
chax = tags.AX_setorient_Slice_TRA.MLHandle;
set(chax, 'Units', 'pixels');
so_cfg.TraLineX = line([0; 0.999], [0.5; 0.5], 'Color', [0, 0, 1], 'Parent', chax);
so_cfg.TraLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', [0, 0, 1], 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');

% -> Zoom
chax = tags.AX_setorient_Slice_Zoom.MLHandle;
set(chax, 'Units', 'pixels');
so_cfg.ZoomLineX = line([0; 0.999], [0.5; 0.5], 'Color', [0, 0, 1], 'Parent', chax);
so_cfg.ZoomLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', [0, 0, 1], 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');

% set callbacks
tags.ED_setorient_XT.Callback = @so_fromtext;
tags.ED_setorient_YT.Callback = @so_fromtext;
tags.ED_setorient_ZT.Callback = @so_fromtext;
tags.ED_setorient_XR.Callback = @so_fromtext;
tags.ED_setorient_YR.Callback = @so_fromtext;
tags.ED_setorient_ZR.Callback = @so_fromtext;
tags.ED_setorient_XS.Callback = @so_fromtext;
tags.ED_setorient_YS.Callback = @so_fromtext;
tags.ED_setorient_ZS.Callback = @so_fromtext;
tags.BT_setorient_XTDD.Callback = {@so_frombtn, 'xtdd'};
tags.BT_setorient_XTHD.Callback = {@so_frombtn, 'xthd'};
tags.BT_setorient_XTHU.Callback = {@so_frombtn, 'xthu'};
tags.BT_setorient_XTDU.Callback = {@so_frombtn, 'xtdu'};
tags.BT_setorient_YTDD.Callback = {@so_frombtn, 'ytdd'};
tags.BT_setorient_YTHD.Callback = {@so_frombtn, 'ythd'};
tags.BT_setorient_YTHU.Callback = {@so_frombtn, 'ythu'};
tags.BT_setorient_YTDU.Callback = {@so_frombtn, 'ytdu'};
tags.BT_setorient_ZTDD.Callback = {@so_frombtn, 'ztdd'};
tags.BT_setorient_ZTHD.Callback = {@so_frombtn, 'zthd'};
tags.BT_setorient_ZTHU.Callback = {@so_frombtn, 'zthu'};
tags.BT_setorient_ZTDU.Callback = {@so_frombtn, 'ztdu'};
tags.BT_setorient_XRDD.Callback = {@so_frombtn, 'xrdd'};
tags.BT_setorient_XRHD.Callback = {@so_frombtn, 'xrhd'};
tags.BT_setorient_XRHU.Callback = {@so_frombtn, 'xrhu'};
tags.BT_setorient_XRDU.Callback = {@so_frombtn, 'xrdu'};
tags.BT_setorient_YRDD.Callback = {@so_frombtn, 'yrdd'};
tags.BT_setorient_YRHD.Callback = {@so_frombtn, 'yrhd'};
tags.BT_setorient_YRHU.Callback = {@so_frombtn, 'yrhu'};
tags.BT_setorient_YRDU.Callback = {@so_frombtn, 'yrdu'};
tags.BT_setorient_ZRDD.Callback = {@so_frombtn, 'zrdd'};
tags.BT_setorient_ZRHD.Callback = {@so_frombtn, 'zrhd'};
tags.BT_setorient_ZRHU.Callback = {@so_frombtn, 'zrhu'};
tags.BT_setorient_ZRDU.Callback = {@so_frombtn, 'zrdu'};
tags.BT_setorient_XSDD.Callback = {@so_frombtn, 'xsdd'};
tags.BT_setorient_XSHD.Callback = {@so_frombtn, 'xshd'};
tags.BT_setorient_XSHU.Callback = {@so_frombtn, 'xshu'};
tags.BT_setorient_XSDU.Callback = {@so_frombtn, 'xsdu'};
tags.BT_setorient_YSDD.Callback = {@so_frombtn, 'ysdd'};
tags.BT_setorient_YSHD.Callback = {@so_frombtn, 'yshd'};
tags.BT_setorient_YSHU.Callback = {@so_frombtn, 'yshu'};
tags.BT_setorient_YSDU.Callback = {@so_frombtn, 'ysdu'};
tags.BT_setorient_ZSDD.Callback = {@so_frombtn, 'zsdd'};
tags.BT_setorient_ZSHD.Callback = {@so_frombtn, 'zshd'};
tags.BT_setorient_ZSHU.Callback = {@so_frombtn, 'zshu'};
tags.BT_setorient_ZSDU.Callback = {@so_frombtn, 'zsdu'};
tags.CB_setorient_updhdr.Callback = @so_updhdr;
tags.BT_setorient_accept.Callback = @so_accept;
hFig.CloseRequestFcn = @so_closewindow;
hFig.KeyPressFcn = @so_keypress;
hFig.KeyReleaseFcn = @so_keyrelease;

% show the first page (3-slices view)
so_showpage(0, 0, 1);

% set initial slicing position (0, 0, 0)
so_setslicepos;

% hdr?
if strcmpi(xo.S.Extensions{1}, 'hdr')
    hFig.SetGroupEnabled('UpdHdr', 'on');
    if opts.updhdr
        tags.CB_setorient_updhdr.Value = 1;
    end
end

% set correct renderer and make figure visible
hFig.HandleVisibility = 'callback';
hFig.Visible = 'on';
hFig.WindowStyle = 'modal';

% then we wait
waitfor(hFig.MLHandle, 'Visible', 'off');

% clear transimg objects
delete(so_cfg.tio.imSag);
delete(so_cfg.tio.imCor);
delete(so_cfg.tio.imTra);
delete(so_cfg.tio.imSlZ);

% no changes requested
if ~so_cfg.update

    % make sure TrfPlus is back to what it was
    xo.C = bc;
    return;
end

% compute output
cm = so_cfg.cm;
cm{2} = (pi / 180) .* cm{2};
mat44 = ne_methods.spmtrf(cm{:});

% make sure fields are not minimally different from 0
mat44(abs(mat44) < 1e-12) = 0;

% if update in HDR and type == HDR
if strcmpi(xo.S.Extensions{1}, 'hdr') && so_cfg.updhdr

    % compute new Trf matrix
    if size(bc.VoxelData, 4) == 1 || ~isfield(bc.RunTimeVars, 'Mat44') || ...
        size(bc.RunTimeVars.Mat44, 1) ~= 4 || size(bc.RunTimeVars.Mat44, 2) ~= 4
        cfr = hdr_CoordinateFrame(xo);
        trf = cfr.Trf;
        trf = mat44 * trf;
    else
        bc.RunTimeVars.Mat44 = ne_methods.mtimesnd(mat44(:, :, ones(1, size(bc.RunTimeVars.Mat44, 3))), ...
            bc.RunTimeVars.Mat44(:, :, :, 1));
        trf = bc.RunTimeVars.Mat44(:, :, 1);
    end

    % make sure file is at least ni1
    if bc.NIIFileType < 1
        bc.FileMagic = 'ni1';
        bc.NIIFileType = 1;
        bc.DataHist.NIftI1.QuaternionB = 0;
        bc.DataHist.NIftI1.QuaternionC = 1;
        bc.DataHist.NIftI1.QuaternionD = 0;
        bc.DataHist.NIftI1.QuatOffsetX = bc.ImgDim.PixSpacing(2) * bc.DataHist.OriginSPM(1);
        bc.DataHist.NIftI1.QuatOffsetY = bc.ImgDim.PixSpacing(3) * bc.DataHist.OriginSPM(2);
        bc.DataHist.NIftI1.QuatOffsetZ = bc.ImgDim.PixSpacing(4) * bc.DataHist.OriginSPM(3);
    end

    % set NIftI data
    if ~any([2, 4] == bc.DataHist.NIftI1.QFormCode)
        bc.DataHist.NIftI1.QFormCode = 2;
    end
    bc.DataHist.NIftI1.SFormCode = 2;
    trf(1:3, end) = trf(1:3, end) + trf(1:3, 1:3) * [1; 1; 1];
    bc.DataHist.NIftI1.AffineTransX = trf(1, :);
    bc.DataHist.NIftI1.AffineTransY = trf(2, :);
    bc.DataHist.NIftI1.AffineTransZ = trf(3, :);

    % re-set TrfPlus
    bc.RunTimeVars.TrfPlus = eye(4);

    % set to memory
    xo.C = bc;

    % and write HDR/NII header
    try
        hdr_UpdateHeader(xo);
        aft_SaveRunTimeVars(xo);
    catch xfferror
        warning('neuroelf:xff:internalError', 'Error updating header: ''%s''.', xfferror.message);
    end

% save all the save
elseif so_cfg.updhdr

    % save
    aft_SaveRunTimeVars(xo);
end


% % % UI functions (called internally)


function so_keypress(src, ke, varargin)
global so_cfg;

% only valid from current figure
if src ~= so_cfg.hFigMLH
    return;
end

% get Key and Modifier from keyboard event (see Matlab docu!)
kk = ke.Key;
mn = ke.Modifier;
so_cfg.mods = mn;

% determine which modifiers are pressed
km = false(1, 4);
if ~isempty(mn)
    try
        km = [any(strcmpi('alt', mn)), any(strcmpi('control', mn)), ...
            any(strcmpi('shift', mn)), any(strcmpi('command', mn))];
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

% get current position
cm = so_cfg.cm;
cpag = so_cfg.page;

% handle events without modifiers
if ~any(km)
    switch (lower(kk))

        % cursor movements -> update translation
        case 'downarrow'
            so_cfg.cm{1}(3) = min(256, max(-256, cm{1}(3) + 1));
            so_setslicepos;
        case 'leftarrow'
            so_cfg.cm{1}(1) = min(256, max(-256, cm{1}(1) - 1));
            so_setslicepos;
        case 'rightarrow'
            so_cfg.cm{1}(1) = min(256, max(-256, cm{1}(1) + 1));
            so_setslicepos;
        case 'uparrow'
            so_cfg.cm{1}(3) = min(256, max(-256, cm{1}(3) - 1));
            so_setslicepos;

        % toggle lines on/off for VMR Browser
        case 'a'
            sl = ~so_cfg.chair;
            so_cfg.chair = sl;
            if sl
                sl = 'on';
            else
                sl = 'off';
            end

            % depending on currently shown page
            if cpag == 1
                set([so_cfg.CorLineX, so_cfg.CorLineY, so_cfg.SagLineX, so_cfg.SagLineY, ...
                     so_cfg.TraLineX, so_cfg.TraLineY], 'Visible', sl);
                set([so_cfg.ZoomLineX, so_cfg.ZoomLineY], 'Visible', 'off');
            else
                set([so_cfg.CorLineX, so_cfg.CorLineY, so_cfg.SagLineX, so_cfg.SagLineY, ...
                     so_cfg.TraLineX, so_cfg.TraLineY], 'Visible', 'off');
                set([so_cfg.ZoomLineX, so_cfg.ZoomLineY], 'Visible', sl);
            end

        % reset
        case 'r'
            so_cfg.cm = so_cfg.cmb;
            so_setslicepos;

        % zoomed view
        case 't'

            % switch through the 4 zoom modes
            so_cfg.zoom = mod(so_cfg.zoom + 1, 4);

            % then either show 3-slice (page 1) or zoom (page 2)
            if so_cfg.zoom == 0
                so_showpage(0, 0, 1);
            else
                so_showpage(0, 0, 2);
            end
    end

% handle events with CTRL
elseif km(2)
    switch (lower(kk))

        % rotation
        case 'downarrow'
            so_cfg.cm{2}(1) = mod(cm{2}(1) + 359, 360);
            so_setslicepos;
        case 'leftarrow'
            so_cfg.cm{2}(3) = mod(cm{2}(3) + 359, 360);
            so_setslicepos;
        case 'rightarrow'
            so_cfg.cm{2}(3) = mod(cm{2}(3) + 1, 360);
            so_setslicepos;
        case 'uparrow'
            so_cfg.cm{2}(1) = mod(cm{2}(1) + 1, 360);
            so_setslicepos;
    end


% handle events with SHIFT
elseif km(3)
    switch (lower(kk))

        % cursor movements
        case 'downarrow'
            so_cfg.cm{1}(2) = min(256, max(-256, cm{1}(2) + 1));
            so_setslicepos;
        case 'leftarrow'
            so_cfg.cm{2}(2) = mod(cm{2}(2) + 359, 360);
            so_setslicepos;
        case 'rightarrow'
            so_cfg.cm{2}(2) = mod(cm{2}(2) + 1, 360);
            so_setslicepos;
        case 'uparrow'
            so_cfg.cm{1}(2) = min(256, max(-256, cm{1}(2) - 1));
            so_setslicepos;
    end
end

function so_keyrelease(varargin)
global so_cfg;

% get key
if nargin < 2
    return;
end
ke = varargin{2};

% remove from list
mo = so_cfg.mods;
if isempty(mo)
    return;
end
mo(strcmpi(mo, ke.Key)) = [];

% and definitely remove command
if ~isempty(mo)
    mo(strcmpi(mo, 'command')) = [];
end
so_cfg.mods = mo;

function so_fromtext(varargin)
global so_cfg;

% update each element in turn
try
    el = str2double(so_cfg.tags.ED_setorient_XT.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{1}(1) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end
try
    el = str2double(so_cfg.tags.ED_setorient_YT.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{1}(2) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end
try
    el = str2double(so_cfg.tags.ED_setorient_ZT.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{1}(3) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end
try
    el = str2double(so_cfg.tags.ED_setorient_XR.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{2}(1) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end
try
    el = str2double(so_cfg.tags.ED_setorient_YR.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{2}(2) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end
try
    el = str2double(so_cfg.tags.ED_setorient_ZR.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{2}(3) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end
try
    el = str2double(so_cfg.tags.ED_setorient_XS.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{3}(1) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end
try
    el = str2double(so_cfg.tags.ED_setorient_YS.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{3}(2) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end
try
    el = str2double(so_cfg.tags.ED_setorient_ZS.String);
    if numel(el) == 1 && ~isinf(el) && ~isnan(el)
        so_cfg.cm{3}(3) = el;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end

% update slicepos
so_setslicepos;

function so_frombtn(varargin)
global so_cfg;

% only if third input is char
if nargin < 3 || ~ischar(varargin{3})
    disp(nargin);
    return;
end

% get config
cm = so_cfg.cm;

% which axes
ax = double(lower(varargin{3}(1))) - 119;

% what type
switch (lower(varargin{3}(2)));
    case 't'
        ty = 1;
    case 'r'
        ty = 2;
    case 's'
        ty = 3;
end

% magnitude
if varargin{3}(3) == 'd'
    mg = 1;
else
    mg = 0.1;
end

% direction
if varargin{3}(4) == 'd'
    mg = -mg;
end

% depending on type
switch (ty)
    case 1
        so_cfg.cm{ty}(ax) = min(256, max(-256, cm{ty}(ax) + mg));
    case 2
        if mg < 0
            mg = 360 + mg;
        end
        so_cfg.cm{ty}(ax) = mod(cm{ty}(ax) + mg, 360);
    case 3
        so_cfg.cm{ty}(ax) = min(4, max(.2, cm{ty}(ax) + 0.1 * mg));
end

% update
so_setslicepos;

function so_updhdr(varargin)
global so_cfg;

% update flag
so_cfg.updhdr = (so_cfg.tags.CB_setorient_updhdr.Value > 0);

function so_accept(varargin)
global so_cfg;

% set update to true (accept new settings)
so_cfg.update = true;

% then close window
so_closewindow;

function so_closewindow(varargin)
global so_cfg;

% then hide the main figure
try
    so_cfg.hFig.Visible = 'off';
catch xfferror
    neuroelf_lasterr(xfferror);
end

function so_setslicepos(varargin)
global so_cfg ne_methods;

% currently configured position
cm = so_cfg.cm;

% update texts
so_cfg.tags.ED_setorient_XT.String = sprintf('%.1f', cm{1}(1));
so_cfg.tags.ED_setorient_YT.String = sprintf('%.1f', cm{1}(2));
so_cfg.tags.ED_setorient_ZT.String = sprintf('%.1f', cm{1}(3));
so_cfg.tags.ED_setorient_XR.String = sprintf('%.1f', cm{2}(1));
so_cfg.tags.ED_setorient_YR.String = sprintf('%.1f', cm{2}(2));
so_cfg.tags.ED_setorient_ZR.String = sprintf('%.1f', cm{2}(3));
so_cfg.tags.ED_setorient_XS.String = sprintf('%.3f', cm{3}(1));
so_cfg.tags.ED_setorient_YS.String = sprintf('%.3f', cm{3}(2));
so_cfg.tags.ED_setorient_ZS.String = sprintf('%.3f', cm{3}(3));

% get shortcut with transimg objects
tio = so_cfg.tio;
if so_cfg.zoom == 0
    diropt = 'all';
    tio = [tio.imSag, tio.imCor, tio.imTra];
else
    diropt = so_cfg.dirorder{so_cfg.zoom};
    tio = tio.imSlZ;
end

% update svar
svar = so_cfg.SliceVar;
svarc = svar.C;
cm{2} = (pi / 180) .* cm{2};
mat44 = ne_methods.spmtrf(cm{:});
mat44(abs(mat44) < 1e-12) = 0;
svarc.RunTimeVars.TrfPlus = mat44;
svar.C = svarc;

% sample requested directions and get sampling voxel value and coord
aft_SliceToTransimg(svar, [0, 0, 0], tio, struct('dir', diropt, 'frame', so_cfg.sframe, ...
    'layers', 1, 'mapvol', 1, 'method', 'linear', 'trans', []));

% display
display(render(tio(1)));
if so_cfg.zoom == 0
    display(render(tio(2)));
    display(render(tio(3)));
end

function so_showpage(varargin)
global so_cfg;

% what about crosshair display
if so_cfg.chair
    sl = 'on';
else
    sl = 'off';
end

% show requested page
so_cfg.hFig.ShowPage(varargin{3}, 'norefresh');

% then get current page (allowing for string call)
so_cfg.page = so_cfg.hFig.ShowPage('cur');

% update
so_setslicepos;

% for first page
if so_cfg.page == 1

    % take show flag for 3-slice crosshairs
    sl1 = sl;

    % and set zoomed crosshairs invisible
    sl2 = 'off';

% for second page
else

    % vice versa
    sl1 = 'off';
    sl2 = sl;
end

% make the settings
set([so_cfg.CorLineX, so_cfg.CorLineY, so_cfg.SagLineX, so_cfg.SagLineY, ...
     so_cfg.TraLineX, so_cfg.TraLineY], 'Visible', sl1);
set([so_cfg.ZoomLineX, so_cfg.ZoomLineY], 'Visible', sl2);
set([get(so_cfg.CorLineX, 'Parent'), get(so_cfg.SagLineX, 'Parent'), ...
     get(so_cfg.TraLineX, 'Parent'), get(so_cfg.ZoomLineX, 'Parent')], 'Visible', 'off');

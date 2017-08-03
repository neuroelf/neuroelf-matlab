function varargout = ne_setslicepos(varargin)
% ne_setslicepos  - set slicing position (update window output)
%
% FORMAT:       ne_setslicepos(SRC, EVT [, cpos [, dodraw]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       cpos        position at which to slice data (otherwise from config)
%       dodraw      override drawing flag (1x1 boolean)
%
% No output fields.
%
% Example:
%
%     ne_setslicepos(0, 0, [30, -20, 0], true);
%
%     perform a drawing at [30, -20, 0] and slice output to viewer.
%
% See also ne_draw.

% Version:  v1.1
% Build:    16061500
% Date:     Jun-15 2016, 12:20 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% global variable
global ne_gcfg;

% get handle shortcuts
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
ci = ne_gcfg.c.ini;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% currently configured position
cpos = cc.cpos;
cpsu = false;

% update with input
if nargin > 2 && isa(varargin{3}, 'double') && numel(varargin{3}) == 3 && ...
   ~any(isinf(varargin{3}) | isnan(varargin{3}) | abs(varargin{3}) > 256)
    cpos = 0.5 .* round(2 .* varargin{3}(:)');
    ne_gcfg.fcfg.cpos = cpos;
    cpsu = true;
end

% only draw for mouse based drawing
if nargin > 3 && ischar(varargin{4}) && ~isempty(varargin{4}) && ...
    isfield(ci.Drawing, varargin{4}(:)')
    dodraw = ci.Drawing.(varargin{4}(:)');
else
    dodraw = false;
end

% get currently set orientation
o = lower(cc.orient(1));

% position of controls on which clicks allows updating of position
fPos = [cc.slicepos; cc.zslicepos; cc.tcpos; cc.histpos];

% get mouse position
nPos = ch.MainFig.CurrentPoint;

% the size of the click area is fixed to 256 pixels (integer coordinates)
vsz = 256;

% get number of currently displayed page and, if required, zoomed slice
cpg = cc.page;
showz = cc.zoom;

% get variable that is to be sliced
svar = cc.SliceVar;
if isxff(svar, true)
    svartyp = lower(svar.FileType);
else
    svar = struct('FileType', '', 'RunTimeVars', struct);
    svartyp = '';
end
uvar = cc.SliceUnder;
if isxff(uvar, true)
    urtv = uvar.RunTimeVars;
end

% get handle of StatsVar and current selection
stvar = cc.StatsVar;
if isxff(stvar, true)
    stvartyp = lower(stvar.FileType);
    strtv = stvar.RunTimeVars;
else
    stvartyp = '';
    strtv = struct;
end
svi = cc.StatsVarIdx;

% and transformation values (if not == eye(4))
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
if isfield(rtv, 'SPMsn') && ...
    isstruct(rtv.SPMsn) && ...
    numel(rtv.SPMsn) == 1
    spmsn = {'snmat', rtv.SPMsn};
else
    spmsn = {};
end
if ~isempty(trans)
    bvpos = [cpos, 1] * tplus * trans' * ttrf;
else
    bvpos = [cpos, 1] * tplus * ttrf;
end
bvpos(end) = [];

% update TCPlot handles
ccf = fieldnames(ne_gcfg.cc);
if ~isempty(ccf)
    ccf(cellfun('isempty', regexp(ccf, '^TC'))) = [];
end
if ~isempty(ccf) && numel(ne_gcfg.cc.(ccf{1})) == 1 && isstruct(ne_gcfg.cc.(ccf{1})) && ...
    isfield(ne_gcfg.cc.(ccf{1}), 'Satellite') && isxfigure(ne_gcfg.cc.(ccf{1}).Satellite, true)
    tcpct = ccf{1};
    tcpch = ne_gcfg.cc.(ccf{1});
    ch.TCPlot = tcpch.TCPlot;
    ch.TCPlotChild = tcpch.TCPlotChild;
    ch.TCPlotChildren = tcpch.TCPlotChildren;
    ch.TCPlotDiscards = tcpch.TCPlotDiscards;
else
    tcpct = '';
end

% for first page (3-slice view)
if cpg == 1

    % set object positions of SAG/COR/TRA image to -1 (cannot hit)
    fPos(4:6, :) = -1;

    % make hit-test
    cobj = findfirst( ...
        fPos(:, 1) <= nPos(1) & fPos(:, 2) <= nPos(2) & ...
        fPos(:, 3) >  nPos(1) & fPos(:, 4) >  nPos(2));

% for page 2
elseif cpg == 2

    % do the reverse (set 3-slice objects' positions to -1)
    fPos(1:3, :) = -1;
    cobj = findfirst( ...
        fPos(:, 1) <= nPos(1) & fPos(:, 2) <= nPos(2) & ...
        fPos(:, 3) >  nPos(1) & fPos(:, 4) >  nPos(2));

    % but as the three slices are the same
    if ~isempty(cobj) && ...
        cobj < 7

        % we must take the current zoom selection into account
        cobj = cobj + showz - 1;
    end

% otherwise (page ~= 1 and ~= 2)
else

    % no slices displayed
    fPos([1:6, 8], :) = -1;
    cobj = findfirst( ...
        fPos(:, 1) <= nPos(1) & fPos(:, 2) <= nPos(2) & ...
        fPos(:, 3) >  nPos(1) & fPos(:, 4) >  nPos(2));
end

% the click occurred on one of those objects (and within the main window)
if ~isempty(cobj) && any(ne_gcfg.c.btdown == ch.MainFigMLH) && ~cpsu

    % compute relative position
    cp = nPos - fPos(cobj, 1:2);
    slxsz = ch.MainFigTags.IM_NeuroElf_Slice_COR.Position(3:4);
    slzsz = ch.MainFigTags.IM_NeuroElf_Slice_Zoom.Position(3:4);
    tcpsz = ch.MainFigTags.AX_NeuroElf_TC_Plot.Position(3:4);

    % and depending on which object was hit
    switch (cobj)

        % SAG slice
        case {1}

            % update TAL-Y and -Z and set drawing directions
            cpos(2) = 128 - round(cp(1, 1)) * (vsz / slxsz(1));
            cpos(3) = -128 + round(cp(1, 2)) * (vsz / slxsz(2));
            ne_gcfg.fcfg.ddir = [1, 2];

        % COR slice
        case {2}

            % update TAL-X and -Z
            if o == 'r'
                cpos(1) = 128 - round(cp(1, 1)) * (vsz / slxsz(1));
            else
                cpos(1) = -128 + round(cp(1, 1)) * (vsz / slxsz(1));
            end
            cpos(3) = -128 + round(cp(1, 2)) * (vsz / slxsz(2));
            ne_gcfg.fcfg.ddir = [2, 3];

        % TRA slice
        case {3}

            % update TAL-X and -Y
            if o == 'r'
                cpos(1) = 128 - round(cp(1, 1)) * (vsz / slxsz(1));
            else
                cpos(1) = -128 + round(cp(1, 1)) * (vsz / slxsz(1));
            end
            cpos(2) = 128 - round(slxsz(2) - (cp(1, 2) + 1)) * (vsz / slxsz(2));
            ne_gcfg.fcfg.ddir = [1, 3];

        % and the same for the zoomed slices
        case {4}
            cpos(2) = 128 - round(cp(1, 1)) * (vsz / slzsz(1));
            cpos(3) = -128 + round(cp(1, 2)) * (vsz / slzsz(2));
            ne_gcfg.fcfg.ddir = [1, 2];
        case {5}
            if o == 'r'
                cpos(1) = 128 - round(cp(1, 1)) * (vsz / slzsz(1));
            else
                cpos(1) = -128 + round(cp(1, 1)) * (vsz / slzsz(1));
            end
            cpos(3) = -128 + round(cp(1, 2)) * (vsz / slzsz(2));
            ne_gcfg.fcfg.ddir = [2, 3];
        case {6}
            if o == 'r'
                cpos(1) = 128 - round(cp(1, 1)) * (vsz / slzsz(1));
            else
                cpos(1) = -128 + round(cp(1, 1)) * (vsz / slzsz(1));
            end
            cpos(2) = 128 - round(slzsz(2) - (cp(1, 2) + 1)) * (vsz / slzsz(2));
            ne_gcfg.fcfg.ddir = [1, 3];

        % for the TC-plot
        case {7}

            % check that the plot is enabled
            if cc.tcplot || isxff(stvar, 'vtc')

                % then get the number of volumes and adapt the position
                if isfield(cc.SliceVar, 'NrOfVolumes')
                    nvol = cc.SliceVar.NrOfVolumes;
                    tsvalue = min(nvol, max(1, round((nvol + 2) * cp(1) / tcpsz(1)) - 1));
                    ch.Coord.TempSlider.Value = tsvalue;
                elseif isfield(cc.SliceVar, 'VoxelData') && size(cc.SliceVar.VoxelData, 4) > 1
                    nvol = size(cc.SliceVar.VoxelData, 4);
                    tsvalue = min(nvol, max(1, round((nvol + 2) * cp(1) / tcpsz(1)) - 1));
                    ch.Coord.TempSlider.Value = tsvalue;
                elseif isxff(stvar, 'vtc')
                    nvol = stvar.RunTimeVars.NrOfVolumesPerTC;
                    tsvalue = min(nvol, max(1, ((nvol + 2) * cp(1) / (1.075 * tcpsz(1))) + 1));
                    stvar.RunTimeVars.SubMapVol = tsvalue;
                end
            end

        % for the histogram (min/max) sliders
        case {8}

            % only if a valid object is shown
            if any(strcmp(svartyp, {'fmr', 'hdr', 'head', 'mgh', 'vmr', 'vtc'}))

                % y axes and RunTimeVars
                hpos = cp(2) / 255;
                if strcmp(svartyp, 'vmr') && ...
                    rtv.ShowV16
                    swl = rtv.ScalingWindowLim16;
                else
                    swl = rtv.ScalingWindowLim;
                end

                % determine which boundary
                if cc.histset == 0

                    % the closer one
                    if abs(hpos - cc.histval(1)) < abs(hpos - cc.histval(2))
                        ne_gcfg.fcfg.histset = 1;
                    else
                        ne_gcfg.fcfg.histset = 2;
                    end
                end
                hset = ne_gcfg.fcfg.histset;

                % which histogram line and limit to update?
                if hset == 1

                    % compute new value
                    hpos = min(hpos, cc.histval(2) - 1/256);
                    ne_gcfg.fcfg.histval(1) = hpos;

                    % and set
                    set(ch.HistLine1, 'YData', 0.002 + 0.996 * [hpos; hpos]);
                else
                    hpos = max(hpos, cc.histval(1) + 1/256);
                    ne_gcfg.fcfg.histval(2) = hpos;
                    set(ch.HistLine2, 'YData', 0.002 + 0.996 * [hpos; hpos]);
                end

                % compute actual value
                hval = hpos * swl(2) + (1 - hpos) * swl(1);

                % and update
                if strcmp(svartyp, 'vmr') && ...
                    rtv.ShowV16
                    svar.RunTimeVars.ScalingWindow16(hset) = hval;
                    swn = svar.RunTimeVars.ScalingWindow16;
                else
                    svar.RunTimeVars.ScalingWindow(hset) = hval;
                    swn = svar.RunTimeVars.ScalingWindow;
                end
                rtv = svar.RunTimeVars;

                % as well as image
                hid = (0:255)';
                swld = abs(swl(2) - swl(1));
                swnd = abs(swn(2) - swn(1));
                hid = uint8(min(255, max(0, round( ...
                    (hid - (256 / swld) * (swn(1) - swl(1))) * (swld / swnd)))));
                set(ch.HistImage, 'CData', repmat(hid, [1, 16, 3]));
            end
    end

    % update ?
    if cobj > 0
        if cobj < 7 || (strcmp(stvartyp, 'vtc') && cobj == 7)
            cpsu = true;
        end
    end

    % make sure cpos is a valid coordinate
    if cc.cstep >= 1
        cpos = ceil(cpos);
    else
        cpos = cc.cstep * round(cpos ./ cc.cstep);
    end

    % use Trf to recompute position
    if ~isempty(trans)
        bvpos = [cpos, 1] * tplus * trans' * ttrf;
    else
        bvpos = [cpos, 1] * tplus * ttrf;
    end
    bvpos(end) = [];

    % and set into configuration
    ne_gcfg.fcfg.cpos = cpos;
    cc = ne_gcfg.fcfg;

% button pressed yet
elseif any(ne_gcfg.c.btdown == ch.MainFigMLH) && ...
   ~cpsu

    % then don't do anything
    return;
end

% return early
if cc.noupdate
    return;
end

% we updated!
ne_gcfg.c.lastupd = now;

% drawing ?
if abs(cc.paint.mode) > 1 && ...
    numel(ne_gcfg.fcfg.ddir) == 2 && ...
    dodraw

    % pass on
    ne_draw;
end

% update satellites
if cpsu && ne_gcfg.c.linked

    % get list of satellites
    sats = fieldnames(ne_gcfg.cc);
    for satc = 1:numel(sats)
        try
            if strcmp(ne_gcfg.cc.(sats{satc}).Config.sattype, 'slice')
                ne_gcfg.cc.(sats{satc}).Config.cpos = cpos;
                ne_setsatslicepos(0, 0, sats{satc}, cpos);
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
end

% truncate precision
cpos = 0.1 .* round(10 .* cpos);
bvpos = 0.1 .* round(10 .* bvpos);

% instantaneous seed-correlation
if cc.instscorr && isxff(cc.instscvar, {'hdr', 'head', 'vtc'}) && ...
    cc.instscvar.NrOfVolumes >= 20
    ne_instscorr;
    cc = ne_gcfg.fcfg;
    stvar = cc.StatsVar;
    svi = cc.StatsVarIdx;
end

% if first page is shown
if cpg == 1

    % sampled values doesn't require position
    cursorpos = '';

    % update position of slice lines (crosshair)
    if o == 'r'
        set(ch.CorLineY, 'XData', (128 - [cpos(1), cpos(1)]) / vsz);
        set(ch.TraLineY, 'XData', (128 - [cpos(1), cpos(1)]) / vsz);
    else
        set(ch.CorLineY, 'XData', (128 + [cpos(1), cpos(1)]) / vsz);
        set(ch.TraLineY, 'XData', (128 + [cpos(1), cpos(1)]) / vsz);
    end
    set(ch.SagLineY, 'XData', (128 - [cpos(2), cpos(2)]) / vsz);
    set(ch.TraLineX, 'YData', (128 + [cpos(2), cpos(2)]) / vsz);
    set(ch.SagLineX, 'YData', (128 + [cpos(3), cpos(3)]) / vsz);
    set(ch.CorLineX, 'YData', (128 + [cpos(3), cpos(3)]) / vsz);

    % and set text controls string
    set(ch.Coord.TEdX, 'String', sprintf('%g', cpos(1)));
    set(ch.Coord.TEdY, 'String', sprintf('%g', cpos(2)));
    set(ch.Coord.TEdZ, 'String', sprintf('%g', cpos(3)));
    set(ch.Coord.VEdX, 'String', sprintf('%g', bvpos(1)));
    set(ch.Coord.VEdY, 'String', sprintf('%g', bvpos(2)));
    set(ch.Coord.VEdZ, 'String', sprintf('%g', bvpos(3)));

% for second page
elseif cpg == 2

    % prepare cursor text (for sampled values info)
    cursorpos = sprintf(' (%g %g %g)', cpos);

    % this depends on which slice direction is shown
    switch (showz)

        % SAG slice
        case {1}
            set(ch.ZoomLineX, 'YData', (128 + [cpos(3), cpos(3)]) / vsz);
            set(ch.ZoomLineY, 'XData', (128 - [cpos(2), cpos(2)]) / vsz);

        % COR slice
        case {2}
            set(ch.ZoomLineX, 'YData', (128 + [cpos(3), cpos(3)]) / vsz);
            if o == 'r'
                set(ch.ZoomLineY, 'XData', (128 - [cpos(1), cpos(1)]) / vsz);
            else
                set(ch.ZoomLineY, 'XData', (128 + [cpos(1), cpos(1)]) / vsz);
            end

        % TRA slice
        case {3}
            set(ch.ZoomLineX, 'YData', (128 + [cpos(2), cpos(2)]) / vsz);
            if o == 'r'
                set(ch.ZoomLineY, 'XData', (128 - [cpos(1), cpos(1)]) / vsz);
            else
                set(ch.ZoomLineY, 'XData', (128 + [cpos(1), cpos(1)]) / vsz);
            end
    end
end

% 3-slice page?
show3 = (showz == 0);

% initialize array of sampled values
svals = [];

% get the to-be-sampled frame
if cc.szoom
    frame = cc.sframez;
else
    frame = cc.sframe;
end

% get shortcut with transimg objects
tio = ne_gcfg.tio;
if show3
    diropt = 'all';
    tio = [tio.imSag, tio.imCor, tio.imTra];
else
    diropt = cc.dirorder{showz};
    tio = tio.imSlZ;
end

% sanity check
cpos = limitrangec(cpos, -255.5, 255.5, 0);

% gray-scale lut
if ~isempty(cc.graylut)
    graylut = {'gcolblend', 'lut', 'gcollut', cc.graylut};
else
    graylut = {};
end

% and do the work if either page 1 or 2 is shown and SliceVar is valid
if any([1, 2] == cpg) && ...
    numel(svar) == 1 && ...
    isxff(svar, true)

    % get filetype
    ftype = lower(svar.Filetype);

    % and (for temporal vars) the volume number
    volnum = round(ch.Coord.TempSlider.Value);

    % make sure this is also in the text control
    ch.Coord.Temp.String = sprintf('%d', volnum);

    % LUT
    if isfield(rtv, 'GrayScaleLUT') && ...
        isequal(size(rtv.GrayScaleLUT), [256, 3])
        graylut = {'gcolblend', 'lut', 'gcollut', rtv.GrayScaleLUT};
    end

    % what interpolation method
    if ch.Interpolate.Value > 0 && ...
       ~strcmp(ftype, 'msk')
        smeth = 'linear';
    else
        smeth = 'nearest';
    end

    % show V16 data
    showv16 = {'v16', false};
    if strcmp(svartyp, 'vmr')
        showv16{2} = rtv.ShowV16;
    end

    % x-scaling
    if strcmp(ftype, 'hdr') && ...
        size(svar.VoxelData, 5) > 1 && ...
        isa(svar.VoxelData, 'single')
        xscale = {'xscale', max(abs(rtv.ScalingWindow))};
    else
        xscale = {};
    end

    % no underlay
    if ~isxff(uvar, {'fmr', 'hdr', 'head', 'mgh', 'vmr', 'vtc'})

        % sample requested directions and get sampling voxel value and coord
        [svals, voxc] = svar.SliceToTransimg(cpos, tio, ...
            struct('dir', diropt, 'frame', frame, 'layers', 1, ...
            'mapvol', volnum, 'method', smeth, 'orient', o, 'trans', trans, ...
            spmsn{:}, graylut{:}, showv16{:}, xscale{:}));
    else

        % normalization of underlay
        if isfield(urtv, 'SPMsn') && ...
            isstruct(urtv.SPMsn) && ...
            numel(urtv.SPMsn) == 1
            spmsnu = {'snmat', urtv.SPMsn};
        else
            spmsnu = {};
        end
        ushowv16 = {'v16', false};
        if isfield(urtv, 'ShowV16')
            ushowv16{2} = urtv.ShowV16;
        end

        % first, sample underlay into layer 1
        uvals = uvar.SliceToTransimg(cpos, tio, ...
            struct('dir', diropt, 'frame', frame, 'layers', 1, ...
            'mapvol', 1, 'method', smeth, 'orient', o, 'trans', trans, ...
            spmsnu{:}, graylut{:}, ushowv16{:}));

        % then sample requested stuff
        [svals, voxc] = svar.SliceToTransimg(cpos, tio, ...
            struct('dir', diropt, 'frame', frame, 'layers', 2, ...
            'mapvol', volnum, 'method', smeth, 'orient', o, 'trans', trans, ...
            spmsn{:}, graylut{:}, showv16{:}, xscale{:}));
        svals = [uvals(:)', svals(:)'];

        % then do a shorthand
        setlayerpixel(tio(1), 1, montagemix( ...
            tio(1).Layer(1).Pixel, tio(1).Layer(2).Pixel, cc.joinulay));
        dellayer(tio(1), 2);
        if show3
            setlayerpixel(tio(2), 1, montagemix( ...
                tio(2).Layer(1).Pixel, tio(2).Layer(2).Pixel, cc.joinulay));
            dellayer(tio(2), 2);
            setlayerpixel(tio(3), 1, montagemix( ...
                tio(3).Layer(1).Pixel, tio(3).Layer(2).Pixel, cc.joinulay));
            dellayer(tio(3), 2);
        end
    end
    
    % additional values
    if volnum == 1 && ...
        strcmp(ftype, 'hdr')
        svol = size(svar.VoxelData);
        if numel(svals) == 1 && ...
            numel(svol) > 3 && ...
            svol(4) > 1 && ...
            svol(4) < 5 && ...
            all(voxc >= 0.5 & voxc <= (svol(1:3) + 0.49))
            svals(2:svol(4)) = svar.VoxelData(1 + sum((round(voxc) - 1) .* cumprod([1, svol(1:2)])) + prod(svol(1:3)) .* (1:svol(4)-1));
        end
    end

    % udpate with actual position
    set(ch.Coord.VEdX, 'String', sprintf('%g', 0.1 * round(10 * voxc(1))));
    set(ch.Coord.VEdY, 'String', sprintf('%g', 0.1 * round(10 * voxc(2))));
    set(ch.Coord.VEdZ, 'String', sprintf('%g', 0.1 * round(10 * voxc(3))));

    % -> for functional files
    timc = [];
    timcx = [];
    srvaru = ch.StatsVarRefs.UserData;
    srvarh = struct;
    if any(strcmp(ftype, {'vtc', 'fmr'})) || ...
       (strcmp(ftype, 'hdr') && ...
        svar.NrOfVolumes > 1) || ...
       (strcmp(ftype, 'head') && ...
        svar.NrOfVolumes > 1)

        % get sampling voxel
        voxc = round(voxc);

        % for VTCs
        if strcmp(ftype, 'vtc')

            % the maximum number comes from VTCData
            dsize = size(svar.VTCData);
            dsize(1) = [];

            % check that the voxel position is valid
            if all(voxc > 0 & voxc <= dsize)

                % then get timecourse
                timc = double(svar.VTCData(:, voxc(1), voxc(2), voxc(3)));
                srvarh = handles(svar);

            % otherwise
            else

                % just a bunch of zeros to plot
                timc = zeros(svar.NrOfVolumes, 1);
            end

        % for FMRs
        elseif strcmp(ftype, 'fmr')

            % the z-dim is potentially defined by dim of Slice field
            dsize = [size(svar.Slice(1).STCData), numel(svar.Slice)];

            % check that the position is valid
            if all(voxc > 0 & voxc <= dsize([1, 2, 4]))

                % for the new fileversion
                if numel(svar.Slice) == 1

                    % simply access the correct data
                    timc = double(squeeze( ...
                        svar.Slice.STCData(voxc(1), voxc(2), :, voxc(3))));

                % for the old fileversion
                else

                    % access data in Slice(Z)
                    timc = double(squeeze( ...
                        svar.Slice(voxc(3)).STCData(voxc(1), voxc(2), :)));
                end

            % for invalid positions
            else

                % plot a series of 0s
                timc = zeros(svar.NrOfVolumes, 1);
            end

        % for HDR
        elseif strcmp(ftype, 'hdr')

            % the maximum number comes from VoxelData
            dsize = size(svar.VoxelData);
            dsize(4:end) = [];

            % check that the initial voxel position is valid
            if all(voxc > 0 & voxc <= dsize)

                % if no Mat44 is defined or all the same
                if ~isfield(rtv, 'Mat44') || ...
                   ~isa(rtv.Mat44, 'double') || ...
                   ~isequal(size(rtv.Mat44), [4, 4, size(svar.VoxelData, 4)]) || ...
                    all(lsqueeze(diff(rtv.Mat44, 1, 3)) == 0)

                    % then get timecourse from that one voxel
                    timc = lsqueeze(double( ...
                        svar.VoxelData(voxc(1), voxc(2), voxc(3), :, 1)));

                % otherwise
                else

                    % get VoxelData handle
                    vxd = svar.VoxelData;
                    m44 = rtv.Mat44;
                    b44 = find(all(all(m44 == 0, 1), 2));
                    if ~isempty(b44)
                        m44(:, :, b44) = ...
                            repmat(svar.CoordinateFrame(1).Trf, [1, 1, numel(b44)]);
                    end

                    % recompute voxel coordinates based on Mat44
                    vc4 = squeeze(mtimesnd(invnd(m44), ...
                        reshape([cpos(:); 1] * ones(1, size(m44, 3)), ...
                        [4, 1, size(m44, 3)])));
                    vc4(1, :) = round(limitrangec(vc4(1, :), 1, dsize(1), 1));
                    vc4(2, :) = round(limitrangec(vc4(2, :), 1, dsize(2), 1));
                    vc4(3, :) = round(limitrangec(vc4(3, :), 1, dsize(3), 1));
                    vc4(4, :) = 1:size(m44, 3);
                    dsize = ([1, cumprod(dsize)])';
                    vc4 = 1 + sum((vc4 - 1) .* dsize(:, ones(1, size(vc4, 2))), 1);

                    % sample at 4D positions
                    timc = lsqueeze(vxd(vc4));
                end

                % scaling
                idim = svar.ImgDim;
                if any([2, 4, 8, 130, 132, 136, 256, 512, 768] == idim.DataType) && ...
                   (idim.ScalingIntercept ~= 0 || ...
                    (idim.ScalingSlope ~= 1 && ...
                     idim.ScalingSlope ~= 0))
                    if all([0, 1] ~= idim.ScalingSlope)
                        timc = idim.ScalingIntercept + idim.ScalingSlope .* timc;
                    else
                        timc = idim.ScalingIntercept + timc;
                    end
                end

            % otherwise
            else

                % just a bunch of zeros to plot
                timc = zeros(svar.NrOfVolumes, 1);
            end

        % for HEAD
        elseif strcmp(ftype, 'head')

            % the maximum number comes from VoxelData
            bdata = svar.Brick;
            dsize = size(bdata(1).Data);
            timc = zeros(numel(bdata), 1);

            % check that the initial voxel position is valid
            if all(voxc > 0 & voxc <= dsize)
                
                % loop over volumes
                voxcx = 1 + sum((voxc - 1) .* [1, cumprod(dsize(1:2))]);
                for bdatac = 1:numel(timc)
                    timc(bdatac) = bdata(bdatac).Data(voxcx);
                    if ~any([0, 1] == bdata(bdatac).ScalingFactor)
                        timc(bdatac) = timc(bdatac) * bdata(bdatac).ScalingFactor;
                    end
                end
            end
        end

    % reference object
    elseif ~isempty(srvaru) && ...
        iscell(srvaru) && ...
        numel(srvaru{1}) == 1 && ...
        isxff(srvaru{1}, 'vtc')

        % the maximum number comes from VTCData
        srvar = srvaru{1};

        voxx = bvcoordconv(cpos, 'tal2bvx', srvar.BoundingBox);

        % check that the voxel position is valid
        if ~isnan(voxx)

            % then get timecourse
            timc = double(srvar.VTCData(:, voxx));
            srvarh = handles(srvar);

        % otherwise
        else

            % just a bunch of zeros to plot
            timc = zeros(srvar.NrOfVolumes, 1);
        end

    % stats var is a VTC
    elseif isxff(stvar, 'vtc') && ...
       ~isempty(svi)

        % get sampling voxel
        voxx = bvcoordconv(cpos, 'tal2bvx', stvar.BoundingBox);

        % the maximum number comes from VTCData
        tsize = strtv.NrOfVolumesPerTC;
        tfrom = strtv.AvgWindowFrom;
        tstep = strtv.AvgWindowStep;
        xtpos = 0.001 * (tstep * (strtv.SubMapVol - 1) - tfrom);
        ntcpc = strtv.NrOfTCsPerCondition;

        % just a bunch of zeros to plot
        timc = zeros(tsize, numel(svi));
        timce = timc;

        % check that the voxel position is valid
        if ~isnan(voxx)

            % for each condition
            for pcc = 1:numel(svi)

                % what indices
                tidx = tsize * ntcpc * (svi(pcc) - 1);
                tidx = tidx+1:tidx+tsize;
                ncon = strtv.NrOfConditionOnsets(svi(pcc));

                % then get timecourse
                timc(:, pcc) = stvar.VTCData(tidx, voxx);
                if ntcpc == 2
                    timce(:, pcc) = abs(1.96 .* (timc(:, pcc) ./ stvar.VTCData(tidx + tsize, voxx)));
                else
                    timcw = 0.51 .* sqrt(ncon .* stvar.VTCData(tidx+2.*tsize, voxx));
                    timce(:, pcc) = stvar.VTCData(tidx+tsize, voxx) ./ timcw;
                end
            end
        end
        timcx = 0.001 .* (strtv.AvgWindowFrom:strtv.AvgWindowStep:strtv.AvgWindowTo);
        srvarh = struct;
    end

    % time course data available?
    if ~isempty(timc)

        % enhanced features?
        if isfield(srvarh, 'StudyData') && ...
            isstruct(srvarh.StudyData) && ...
            isfield(srvarh.StudyData, 'SDMMatrix') && ...
            size(srvarh.StudyData.SDMMatrix, 1) == size(timc, 1) && ...
            isfield(srvarh.StudyData, 'SDMMatrixInv') && ...
            isfield(srvarh.StudyData, 'NrOfConfounds') && ...
            size(timc, 2) == 1

            % regress out nuisance variables
            sdmm = srvarh.StudyData.SDMMatrix;
            if ch.StatsVarRefNuis.Value > 0
                sdmnc = srvarh.StudyData.NrOfConfounds;
                timb = srvarh.StudyData.SDMMatrixInv * sdmm' * timc;
                timbf = timb([1:(end-sdmnc), end]);
                timb([1:(end-sdmnc), end]) = 0;
                timc = timc - sdmm * timb;
            elseif ch.StatsVarRefRegs.Value > 0
                sdmnc = 1;
                timbf = srvarh.StudyData.SDMMatrixInv * sdmm' * timc;
            else
                timbf = [];
            end
        else
            timbf = [];
        end

        % make sure no inf/nan values are in time-course
        iin = isinfnan(timc);
        if any(iin(:))
            for pcc = 1:size(timc, 2)
                timc(iin(:, pcc), pcc) = meannoinfnan(timc(:, pcc));
            end
        end

        % get minimum and maximum
        tmmm = minmaxmean(timc, 4);
        tmmm(isinf(tmmm) | isnan(tmmm)) = 0;

        % and get the range
        tmd = (tmmm(2) - tmmm(1)) + 0.0001;

        % delete old extra children
        if ~isempty(ch.TCPlotChildren) && ishandle(ch.TCPlotChildren(1))
            delete(ch.TCPlotChildren);
        end
        if isempty(tcpct)
            ne_gcfg.h.TCPlotChildren = zeros(0, 1);
        else
            ne_gcfg.cc.(tcpct).TCPlotChildren = zeros(0, 1);
        end

        % then plot the data for regular time courses
        timx = ch.TCPlot.MLHandle;
        hold(timx, 'on');
        if isempty(timcx)
            timcx = 1:numel(timc);
            set(ch.TCPlotChild, 'XData', 1:numel(timc), 'YData', double(timc(:, 1)));
            ch.TCPlot.XLim = [0, numel(timc) + 1];
            ch.TCPlot.YLim = [floor(tmmm(1) - 0.05 * tmd), ceil(tmmm(2) + 0.05 * tmd)];
            if ~isempty(tcpct)
                ne_gcfg.cc.(tcpct).Config.tcvar = cc.SliceVar;
            end

        % for VTC data
        elseif isxff(stvar, 'vtc')

            % recompute (with error)
            tmmm = minmaxmean([timc, timc + timce, timc - timce], 4);
            tmd = (tmmm(2) - tmmm(1)) + 0.0001;

            % then set data
            set(ch.TCPlotChild, 'XData', timcx(1), 'YData', double(timc(1)));
            ch.TCPlot.XLim = [timcx(1), timcx(end)];
            
            % override YLim?
            tcylim = [tmmm(1) - 0.05 * tmd, tmmm(2) + 0.05 * tmd];
            if ~all(isinf(cc.tcplotylim))
                tcylim(~isinf(cc.tcplotylim)) = cc.tcplotylim(~isinf(cc.tcplotylim));
            end
            ch.TCPlot.YLim = tcylim;

            % and use tcplot
            for pcc = 1:size(timc, 2)
                scolor = (1 / 255) .* strtv.Map(svi(pcc)).RGBLowerThreshPos;
                [tcpl, tcpp] = tcplot(timx, timcx, double(timc(:, pcc)), ...
                    timce(:, pcc), timce(:, pcc), struct( ...
                    'color', scolor, 'lwidth', 1.5, 'spline', false, ...
                    'scolor', 0.5 + 0.5 .* scolor, 'spalpha', 0.5));
                if isempty(tcpct)
                    ne_gcfg.h.TCPlotChildren = [ne_gcfg.h.TCPlotChildren; ...
                        tcpl(:); tcpp(:)];
                else
                    ne_gcfg.cc.(tcpct).TCPlotChildren = ...
                        [ne_gcfg.cc.(tcpct).TCPlotChildren; tcpl(:); tcpp(:)];
                end
            end

            % add vertical line
            if isempty(tcpct)
                ne_gcfg.h.TCPlotChildren(end+1) = plot(timx, [xtpos; xtpos], ...
                    ch.TCPlot.YLim(:));
                set(ne_gcfg.h.TCPlotChildren(end), 'Color', [0, 0, 0]);
            else
                ne_gcfg.cc.(tcpct).TCPlotChildren(end+1) = ...
                    plot(timx, [xtpos; xtpos], ch.TCPlot.YLim(:));
                set(ne_gcfg.cc.(tcpct).TCPlotChildren(end), 'Color', [0, 0, 0]);
                ne_gcfg.cc.(tcpct).Config.tcvar = stvar;
            end

        % (or other data)
        else
            set(ch.TCPlotChild, 'XData', timcx, 'YData', double(timc(:, 1)));
            ch.TCPlot.XLim = [timcx(1), timcx(end)];
            ch.TCPlot.YLim = [tmmm(1) - 0.05 * tmd, tmmm(2) + 0.05 * tmd];
            if ~isempty(tcpct)
                ne_gcfg.cc.(tcpct).Config.tcvar = cc.SliceVar;
            end
        end
        ne_gcfg.fcfg.tcplotdata = double(timc);

        % add further regressors
        if ch.StatsVarRefRegs.Value > 0 && ~isempty(timbf)
            timr = sdmm(:, [1:(end-sdmnc), end]) * timbf;
            if isempty(tcpct)
                ne_gcfg.h.TCPlotChildren = [ne_gcfg.h.TCPlotChildren; ...
                    plot(ch.TCPlot.MLHandle, timcx, timr)];
                set(ne_gcfg.h.TCPlotChildren, 'Color', [0.5, 0.5, 0]);
            else
                ne_gcfg.cc.(tcpct).TCPlotChildren = ...
                    [ne_gcfg.cc.(tcpct).TCPlotChildren; ...
                     plot(ch.TCPlot.MLHandle, timcx, timr)];
                set(ne_gcfg.cc.(tcpct).TCPlotChildren, 'Color', [0.5, 0.5, 0]);
            end
        end
        hold(timx, 'off');

        % and finally, make sure it's visible
        ch.TCPlot.Visible = 'on';
        if isempty(tcpct)
            ch.TCPlotUndock.Visible = 'on';
            ch.MainFig.SetGroupEnabled('TCPlot', 'on');
        else
            set(ch.TCPlot.Parent, 'Visible', 'on');
            ch.TCPlotUndock.Visible = 'off';
            ch.MainFig.SetGroupEnabled('TCPlot', 'off');
        end
    else
        % disable timecourse display
        ch.TCPlot.Visible = 'off';
        ch.TCPlotUndock.Visible = 'off';
        ch.MainFig.SetGroupEnabled('TCPlot', 'off');
    end

% no good variable
else

    % set image data to zeros
    if cpg == 1
        sag = uint8(zeros(256, 256));
        setlayer(tio(1), 1, sag);
        setlayer(tio(2), 1, sag);
        setlayer(tio(3), 1, sag);
    elseif cpg == 2
        setlayer(tio, 1, uint8(zeros(512, 512)));
    end

    % disable timecourse display
    ch.TCPlot.Visible = 'off';
    ch.TCPlotUndock.Visible = 'off';
    ch.MainFig.SetGroupEnabled('TCPlot', 'off');
end

% if either page is shown
if any([1, 2] == cpg)

    % remove stats layers (of previous update)
    if numel(tio(1).Layer) > 1
        dellayer(tio(1), 2:numel(tio(1).Layer));
        if show3
            dellayer(tio(2), 2:numel(tio(2).Layer));
            dellayer(tio(3), 2:numel(tio(3).Layer));
        end
    end

    % if is valid xff
    if numel(stvar) == 1 && ~isempty(cc.StatsVarIdx)

        % normalization
        srtv = stvar.RunTimeVars;
        if isfield(srtv, 'SPMsn') && isstruct(srtv.SPMsn) && numel(srtv.SPMsn) == 1
            spmsns = {'snmat', srtv.SPMsn};
        else
            spmsns = {};
        end

        % and interpolation method of choice
        if ch.Interpolate.Value ~= 0
            imeth = ne_gcfg.fcfg.imethod;
        else
            imeth = 'nearest';
        end

        % get tails flag (for threshmapc)
        lthr = [];
        uthr = [];
        stvtype = lower(stvar.Filetype);
        if any(strcmp(stvtype, {'cmp', 'glm', 'hdr', 'head', 'vmp', 'vtc'}))
            lut = [];
            tls = [];
        else
            lut = ne_gcfg.lut.Colors;
            tls = double(ch.Stats.PosTail.Value ~= 0) + ...
                2 * double(ch.Stats.NegTail.Value ~= 0);
        end

        % bar position
        if ci.Statistics.ShowThreshBars
            rbp = ci.Statistics.ThreshBarPos;
        else
            rbp = [];
        end

        % sample data
        [svals(end+1:end+numel(svi)), svoxc, ovals] = stvar.SliceToTransimg(cpos, tio, ...
            struct('dir', diropt, 'frame', frame, 'layers', 2, ...
            'mapvol', svi, 'method', imeth, 'rgbcol', lut, ...
            'rgbctails', tls, 'rgblthr', lthr, 'rgbuthr', uthr, spmsns{:}, ...
            'orient', o, 'trans', trans, 'type', 'rgb', 'rgbbars', rbp));
        if strcmp(stvtype, 'vtc')
            strct = strtv.ConditionThresholds;
            for pcc = 1:numel(svi)
                ovals(pcc) = double(ovals(pcc) > 0) .* ...
                    (strct(svi(pcc), 2, 1) + ovals(pcc) * strct(svi(pcc), 2, 2));
            end
            svals(end+1:end+numel(svi)) = ovals;
        end

        % if requested and necessary
        if cc.join && ...
            numel(tio(1).Layer) > 2

            % special join for two stats layers
            if cc.joinmd2 && numel(tio(1).Layer) == 3

                % find overlay voxels
                t1l = tio(1).Layer;
                mxc = lsqueeze((t1l(2).Alpha > 0 & t1l(3).Alpha > 0));
                tp1 = reshape(t1l(2).Pixel, numel(mxc), 3);
                tp2 = reshape(t1l(3).Pixel, numel(mxc), 3);
                tp3 = tp1;
                ta3 = reshape(zeros(numel(mxc), 1), size(t1l(2).Alpha));
                if any(mxc)
                    tp3(mxc, :) = maxdistcol(tp1(mxc, :), tp2(mxc, :));
                    ta3(mxc) = 0.5 .* (t1l(2).Alpha(mxc) + t1l(3).Alpha(mxc));
                end
                tp3 = reshape(tp3, size(t1l(2).Pixel));
                tio(1) = addlayer(tio(1), tp3, ta3);
                if show3
                    t1l = tio(2).Layer;
                    mxc = lsqueeze((t1l(2).Alpha > 0 & t1l(3).Alpha > 0));
                    tp1 = reshape(t1l(2).Pixel, numel(mxc), 3);
                    tp2 = reshape(t1l(3).Pixel, numel(mxc), 3);
                    tp3 = tp1;
                    ta3 = reshape(zeros(numel(mxc), 1), size(t1l(2).Alpha));
                    if any(mxc)
                        tp3(mxc, :) = maxdistcol(tp1(mxc, :), tp2(mxc, :));
                        ta3(mxc) = 0.5 .* (t1l(2).Alpha(mxc) + t1l(3).Alpha(mxc));
                    end
                    tp3 = reshape(tp3, size(t1l(2).Pixel));
                    tio(2) = addlayer(tio(2), tp3, ta3);
                    t1l = tio(3).Layer;
                    mxc = lsqueeze((t1l(2).Alpha > 0 & t1l(3).Alpha > 0));
                    tp1 = reshape(t1l(2).Pixel, numel(mxc), 3);
                    tp2 = reshape(t1l(3).Pixel, numel(mxc), 3);
                    tp3 = tp1;
                    ta3 = reshape(zeros(numel(mxc), 1), size(t1l(2).Alpha));
                    if any(mxc)
                        tp3(mxc, :) = maxdistcol(tp1(mxc, :), tp2(mxc, :));
                        ta3(mxc) = 0.5 .* (t1l(2).Alpha(mxc) + t1l(3).Alpha(mxc));
                    end
                    tp3 = reshape(tp3, size(t1l(2).Pixel));
                    tio(3) = addlayer(tio(3), tp3, ta3);
                end

            % otherwise
            else
                % apply "stats-map join" (BV-like display) to the three images
                joinlayers(tio(1), 2:numel(tio(1).Layer));
                if show3
                    joinlayers(tio(2), 2:numel(tio(2).Layer));
                    joinlayers(tio(3), 2:numel(tio(3).Layer));
                end
            end
        end

    end

    % display
    display(render(tio(1)));
    if show3
        display(render(tio(2)));
        display(render(tio(3)));
    end

    % get talairach label
    if ch.Stats.TDClient.Value > 0
        if ch.Stats.ICBM2TAL.Value > 0
            tpos = icbm2tal(cpos);
        else
            tpos = cpos;
        end
        tallabel = tdlocal2(2, tpos(1), tpos(2), tpos(3));
        if tallabel(1) == '*'
            tallabel = '';
        else
            tallabel = sprintf('TAL-label: %s', tallabel);
        end
        if isempty(strfind(tallabel, 'Gray Matter'))
            [gmlabel, gmcrd] = tdlabel(round(tpos));
            if ~isempty(gmlabel) && ...
                iscell(gmlabel) && ...
               ~isempty(gmlabel{1})
                gmlabel = deblank(regexprep(gmlabel{1}, '^.H ([^\(]+).*$', '$1'));
                if isempty(strfind(tallabel, gmlabel))
                    tallabel = strrep( ...
                        strrep(strrep(tallabel, 'TAL-label: ', ''), ',*', ''), ...
                        ',White Matter', '');
                    if ch.Stats.ICBM2TAL.Value > 0
                        tallabel = [tallabel ' NGM: ' ...
                            gmlabel sprintf(' (%d, %d, %d)', round(tal2icbm(gmcrd)))];
                    else
                        tallabel = [tallabel ' NGM: ' ...
                            gmlabel sprintf(' (%d, %d, %d)', gmcrd)];
                    end
                end
            end
        end
        tallabel = sprintf('\n%s', tallabel);
    else
        tallabel = '';
    end

    % set sampled values
    if ~isempty(svals)
        ch.SampledValues.String = sprintf('Values at cursor%s:  [%g%s]%s', ...
            cursorpos, svals(1), sprintf('  %.4g', svals(2:end)), tallabel);
    else
        ch.SampledValues.String = ['Values at cursor' cursorpos ':  [  ]'];
    end

    % extract PLP information
    if numel(cc.plp) == 1 && ...
        isxff(cc.plp, 'plp') && ...
        ci.MKDA.LookupOnCursor
        cpopt = struct( ...
            'fringe', ci.MKDA.LookupFringe);
        prtv = cc.plp.RunTimeVars;
        if isfield(prtv, 'MKDAAnalyses') && ...
            isfield(ne_gcfg.h.MKDA, 'MKDAFig') && ...
            isxfigure(ne_gcfg.h.MKDA.MKDAFig, true)
            plph = ne_gcfg.h.MKDA.h;
            cpopt.cond = plph.Points.UserData{1};
            cpopt.unitcol = plph.StudyColumn.String{plph.StudyColumn.Value};
        end
        if numel(stvar) == 1 && ...
            isxff(stvar, 'vmp') && ...
            numel(svi) == 1
            stvmap = stvar.Map(svi);
            if isfield(stvmap, 'RunTimeVars') && ...
                isstruct(stvmap.RunTimeVars) && ...
                numel(stvmap.RunTimeVars) == 1 && ...
                isfield(stvmap.RunTimeVars, 'MKDAMapType')
                if isfield(stvmap.RunTimeVars, 'Condition') && ...
                   ~isempty(stvmap.RunTimeVars.Condition)
                    cpopt.cond = stvmap.RunTimeVars.Condition;
                    cpopt.unitcol = cc.plp.ColumnNames{stvmap.RunTimeVars.UnitColumn};
                    cpopt.uvmp = stvar;
                end
                cpopt.fringe = stvmap.RunTimeVars.FWHM;
            end
        end
        ne_gcfg.h.ClusterTable.String = ...
            cc.plp.CPointsAtCoord(cpos, cpopt);
    end
end

% update for remote?
if ne_gcfg.c.remote
    ifmt = lower(ci.Remote.ImageFormat);
    if strcmp(ifmt, 'jpg')
        iqual = {'Quality', ci.Remote.ImageJPGQuality};
    else
        iqual = {};
    end
    ipath = [neuroelf_path('remote') '/images'];
    try
        if show3
            imwrite(tio(1).Rendered, sprintf('%s/sagslice.%s', ipath, ifmt), iqual{:});
            imwrite(tio(2).Rendered, sprintf('%s/corslice.%s', ipath, ifmt), iqual{:});
            imwrite(tio(3).Rendered, sprintf('%s/traslice.%s', ipath, ifmt), iqual{:});
        else
            imwrite(tio(1).Rendered, sprintf('%s/zoomslice.%s', ipath, ifmt), iqual{:});
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% update GLM beta plot/s?
plotc = fieldnames(ne_gcfg.cc);
for pcc = 1:numel(plotc)
    plotcc = ne_gcfg.cc.(plotc{pcc});
    if isfield(plotcc, 'Config') && ...
        isstruct(plotcc.Config) && ...
        isfield(plotcc.Config, 'glm') && ...
        isxff(plotcc.Config.glm, 'glm') && ...
        isfield(plotcc.Config, 'upplot') && ...
        plotcc.Config.upplot
        try
            ne_glmplotbetasup(0, 0, plotcc.Config.glm, '', plotc{pcc});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
end

% update sampled values
ne_gcfg.c.svals = svals;

% FUNCTION ne_setsatslicepos: set slicing position (update sat window)
function varargout = ne_setsatslicepos(varargin)

% Version:  v1.1
% Build:    16052818
% Date:     May-28 2016, 6:04 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% input
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~isfield(ne_gcfg.cc, varargin{3})
    return;
end

% get handle shortcuts
iSat = varargin{3}(:)';
ch = ne_gcfg.cc.(iSat);
cc = ch.Config;
ci = ne_gcfg.c.ini;
dodraw = ci.Drawing.OnSatellite;
    
% currently configured position
cpos = cc.cpos;

% update with input
if nargin > 3 && ...
    isa(varargin{4}, 'double') && ...
    numel(varargin{4}) == 3 && ...
   ~any(isinf(varargin{4}) | isnan(varargin{4}) | abs(varargin{4}) > 128)
    cpos = round(varargin{4});
end

% get currently set orientation
o = lower(cc.orient(1));

% position of controls on which clicks allows updating of position
fPos = [cc.slicepos; cc.zslicepos];

% get mouse position
nPos = ch.Satellite.CurrentPoint;

% the size of the click area is fixed to 256 pixels (integer coordinates)
tSat = ch.Satellite.TagStruct;
isz = 1 / 256;
vsz = tSat.(sprintf('IM_%s_Slice_SAG', iSat)).Position(3);
vzz = tSat.(sprintf('IM_%s_Slice_Zoom', iSat)).Position(3);

% get number of currently displayed page and, if required, zoomed slice
cpg = cc.page;
showz = cc.zoom;

% get variable that is to be sliced
svar = cc.SliceVar;
if ~isxff(svar, true)
    svar = struct('Filetype', 'NONE', ...
        'RunTimeVars', struct('Trf', eye(4), 'TrfPlus', eye(4)));
end
svartyp = lower(svar.Filetype);
uvar = cc.SliceUnder;

% get handle of StatsVar and current selection
stvar = cc.StatsVar;
svi = cc.StatsVarIdx;

% and transformation values (if not == eye(4))
if ~isequal(cc.strans, eye(4))
    trans = cc.strans;
else
    trans = [];
end

% use Trf to get voxel position
tplus = inv(svar.RunTimeVars.TrfPlus)';
if ~isempty(trans)
    bvpos = [cpos, 1] * tplus * trans' * svar.RunTimeVars.Trf;
else
    bvpos = [cpos, 1] * tplus * svar.RunTimeVars.Trf;
end
bvpos(end) = [];

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
    if ~isempty(cobj)

        % we must take the current zoom selection into account
        cobj = cobj + showz - 1;
    end

% otherwise (page ~= 1 and ~= 2)
else

    % no slices displayed
    fPos(1:6, :) = -1;
    cobj = [];
end

% the click occurred on one of those objects (and within the main window)
if ~isempty(cobj) && ...
    any(ne_gcfg.c.btdown == ch.SatelliteMLH)

    % compute relative position
    cp = nPos - fPos(cobj, 1:2);

    % and depending on which object was hit
    switch (cobj)

        % SAG slice
        case {1}

            % update TAL-Y and -Z and set drawing directions
            cpos(2) = 128 - round(cp(1, 1)) * (256 / vsz);
            cpos(3) = -128 + round(cp(1, 2)) * (256 / vsz);
            ne_gcfg.fcfg.ddir = [1, 2];

        % COR slice
        case {2}

            % update TAL-X and -Z
            if o == 'r'
                cpos(1) = 128 - round(cp(1, 1)) * (256 / vsz);
            else
                cpos(1) = -128 + round(cp(1, 1)) * (256 / vsz);
            end
            cpos(3) = -128 + round(cp(1, 2)) * (256 / vsz);
            ne_gcfg.fcfg.ddir = [2, 3];

        % TRA slice
        case {3}

            % update TAL-X and -Y
            if o == 'r'
                cpos(1) = 128 - round(cp(1, 1)) * (256 / vsz);
            else
                cpos(1) = -128 + round(cp(1, 1)) * (256 / vsz);
            end
            cpos(2) = 128 - round(vsz - (cp(1, 2) + 1)) * (256 / vsz);
            ne_gcfg.fcfg.ddir = [1, 3];

        % and the same for the zoomed slices
        case {4}
            cpos(2) = 128 - round(cp(1, 1)) * (256 / vzz);
            cpos(3) = -128 + round(cp(1, 2)) * (256 / vzz);
            ne_gcfg.fcfg.ddir = [1, 2];
        case {5}
            if o == 'r'
                cpos(1) = 128 - round(cp(1, 1)) * (256 / vzz);
            else
                cpos(1) = -128 + round(cp(1, 1)) * (256 / vzz);
            end
            cpos(3) = -128 + round(cp(1, 2)) * (256 / vzz);
            ne_gcfg.fcfg.ddir = [2, 3];
        case {6}
            if o == 'r'
                cpos(1) = 128 - round(cp(1, 1)) * (256 / vzz);
            else
                cpos(1) = -128 + round(cp(1, 1)) * (256 / vzz);
            end
            cpos(2) = 128 - round(2 * vsz - (cp(1, 2) + 1)) * (256 / vzz);
            ne_gcfg.fcfg.ddir = [1, 3];
    end

    % make sure cpos is a valid coordinate
    if any(cc.cstep >= 1)
        cpos = ceil(cpos);
    else
        cpos = cc.cstep .* round(cpos ./ cc.cstep);
    end

    % use Trf to recompute position
    if ~isempty(trans)
        bvpos = [cpos, 1] * tplus * trans' * svar.RunTimeVars.Trf;
    else
        bvpos = [cpos, 1] * tplus * svar.RunTimeVars.Trf;
    end
    bvpos(end) = [];

    % and set into configuration
    ne_gcfg.cc.(varargin{3}).Config.cpos = cpos;
    cc = ne_gcfg.cc.(varargin{3}).Config;

    % if linked, do differently!
    if ne_gcfg.c.linked && ...
        nargin < 4
        ne_setslicepos(0, 0, cpos, 'OnLinked');
        return;
    end

% button pressed yet
elseif any(ne_gcfg.c.btdown == ch.SatelliteMLH)

    % then don't do anything
    return;
end

% truncate precision
if all(cc.cstep == cc.cstep(1)) && all(cc.cstep >= 0.5)
    cpos = 0.1 .* round(10 .* cpos);
    bvpos = 0.1 .* round(10 .* bvpos);
end

% drawing ?
if abs(ne_gcfg.fcfg.paint.mode) > 1 && numel(ne_gcfg.fcfg.ddir) == 2 && dodraw

    % pass on
    ne_draw(0, 0, cpos);
end

% if first page is shown
if cpg == 1

    % update position of slice lines (crosshair)
    if o == 'r'
        set(ch.CorLineY, 'XData', isz .* (128.25 - [cpos(1), cpos(1)]));
        set(ch.TraLineY, 'XData', isz .* (128.25 - [cpos(1), cpos(1)]));
    else
        set(ch.CorLineY, 'XData', isz .* (128.25 + [cpos(1), cpos(1)]));
        set(ch.TraLineY, 'XData', isz .* (128.25 + [cpos(1), cpos(1)]));
    end
    set(ch.SagLineY, 'XData', isz .* (128.25 - [cpos(2), cpos(2)]));
    set(ch.TraLineX, 'YData', isz .* (128.25 + [cpos(2), cpos(2)]));
    set(ch.SagLineX, 'YData', isz .* (128.25 + [cpos(3), cpos(3)]));
    set(ch.CorLineX, 'YData', isz .* (128.25 + [cpos(3), cpos(3)]));

    % update texts
    ch.Text.BVSX.String = sprintf('%g', bvpos(1));
    ch.Text.BVSY.String = sprintf('%g', bvpos(2));
    ch.Text.BVSZ.String = sprintf('%g', bvpos(3));
    ch.Text.TALX.String = sprintf('%g', cpos(1));
    ch.Text.TALY.String = sprintf('%g', cpos(2));
    ch.Text.TALZ.String = sprintf('%g', cpos(3));

% for second page
elseif cpg == 2

    % this depends on which slice direction is shown
    switch (showz)

        % SAG slice
        case {1}
            set(ch.ZoomLineX, 'YData', isz .* (128.25 + [cpos(3), cpos(3)]));
            set(ch.ZoomLineY, 'XData', isz .* (128.25 - [cpos(2), cpos(2)]));

        % COR slice
        case {2}
            set(ch.ZoomLineX, 'YData', isz .* (128.25 + [cpos(3), cpos(3)]));
            if o == 'r'
                set(ch.ZoomLineY, 'XData', isz .* (128.25 - [cpos(1), cpos(1)]));
            else
                set(ch.ZoomLineY, 'XData', isz .* (128.25 + [cpos(1), cpos(1)]));
            end

        % TRA slice
        case {3}
            set(ch.ZoomLineX, 'YData', isz .* (128.25 + [cpos(2), cpos(2)]));
            if o == 'r'
                set(ch.ZoomLineY, 'XData', isz .* (128.25 - [cpos(1), cpos(1)]));
            else
                set(ch.ZoomLineY, 'XData', isz .* (128.25 + [cpos(1), cpos(1)]));
            end
    end
end

% 3-slice page?
show3 = (showz == 0);

% get the to-be-sampled frame
if cc.szoom
    frame = cc.sframez;
else
    frame = cc.sframe;
end

% and transformation values (if not == eye(4))
if ~isequal(ne_gcfg.fcfg.strans, eye(4))
    trans = ne_gcfg.fcfg.strans;
else
    trans = [];
end

% delete plothandles (if any)
if ~isempty(cc.plpph)
    try
        delete(cc.plpph);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    ne_gcfg.cc.(iSat).Config.plpph = [];
end

% get shortcut with transimg objects
tio = cc.tio;
if show3
    diropt = 'all';
    tio = [tio.imSag, tio.imCor, tio.imTra];
else
    diropt = cc.dirorder{showz};
    tio = tio.imSlZ;
end

% preset sampled values
svals = [];

% gray-scale lut
if ~isempty(cc.graylut)
    graylut = {'gcolblend', 'lut', 'gcollut', cc.graylut};
else
    graylut = {};
end

% and do the work if either page 1 or 2 is shown and SliceVar is valid
if any([1, 2] == cpg) && ...
    numel(svar) == 1 && ...
    ~strcmp(svartyp, 'none')

    % spatial normalization
    rtv = svar.RunTimeVars;
    if isfield(rtv, 'SPMsn') && ...
        isstruct(rtv.SPMsn) && ...
        numel(rtv.SPMsn) == 1
        spmsn = {'snmat', rtv.SPMsn};
    else
        spmsn = {};
    end
    if isxff(uvar, true)
        urtv = uvar.RunTimeVars;
    end

    % and (for temporal vars) the volume number
    volnum = round(ch.Coord.TempSlider.Value);

    % LUT
    if isfield(rtv, 'GrayScaleLUT') && ...
        isequal(size(rtv.GrayScaleLUT), [256, 3])
        graylut = {'gcolblend', 'lut', 'gcollut', rtv.GrayScaleLUT};
    end

    % what interpolation method
    if ch.Interpolate.Value > 0 && ~strcmp(svartyp, 'msk')
        smeth = cc.imethod;
    else
        smeth = 'nearest';
    end

    % show V16 data
    showv16 = {'v16', false};
    if strcmp(svartyp, 'vmr')
        showv16{2} = rtv.ShowV16;
    end

    % x-scaling
    if strcmp(svartyp, 'hdr') && ...
        size(svar.VoxelData, 5) > 1 && ...
        isa(svar.VoxelData, 'single')
        xscale = {'xscale', max(abs(rtv.ScalingWindow))};
    else
        xscale = {};
    end

    % no underlay
    if ~isxff(uvar, {'fmr', 'hdr', 'head', 'vmr', 'vtc'})

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
            gradtypeu{:}, spmsnu{:}, graylut{:}, ushowv16{:}));

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
        strcmp(svartyp, 'hdr')
        svol = size(svar.VoxelData);
        if numel(svals) == 1 && ...
            numel(svol) > 3 && ...
            svol(4) > 1 && ...
            svol(4) < 5 && ...
            all(voxc >= 0.5 & voxc <= (svol(1:3) + 0.49))
            svals(2:svol(4)) = svar.VoxelData(1 + sum((round(voxc) - 1) .* cumprod([1, svol(1:2)])) + prod(svol(1:3)) .* (1:svol(4)-1));
        end
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
    if numel(stvar) == 1 && ...
        ~isempty(cc.StatsVarIdx)

        % normalization
        srtv = stvar.RunTimeVars;
        if isfield(srtv, 'SPMsn') && ...
            isstruct(srtv.SPMsn) && ...
            numel(srtv.SPMsn) == 1
            spmsns = {'snmat', srtv.SPMsn};
        else
            spmsns = {};
        end

        % and interpolation method of choice
        if ch.Interpolate.Value ~= 0
            imeth = cc.imethod;
        else
            imeth = 'nearest';
        end

        % get tails flag (for threshmapc)
        lthr = [];
        uthr = [];
        if any(strcmpi(stvar.Filetype, {'cmp', 'glm', 'hdr', 'head', 'vmp', 'vtc'}))
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
        svals(end+1:end+numel(svi)) = stvar.SliceToTransimg(cpos, tio, ...
            struct('dir', diropt, 'frame', frame, 'layers', 2, ...
            'mapvol', svi, 'method', imeth, 'rgbcol', lut, ...
            'rgbctails', tls, 'rgblthr', lthr, 'rgbuthr', uthr, spmsns{:}, ...
            'orient', o, 'trans', trans, 'type', 'rgb', 'rgbbars', rbp));

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

    % update text
    if numel(svals) < 2
        ch.Text.Values.String = sprintf('Values at cursor:  %g', svals);
    else
        sval1 = sprintf('  %g', svals(1));
        if mod(numel(svals), 2) == 0
            sval = sprintf('%10.5g', svals(end));
            svals(end) = [];
        else
            sval = '';
        end
        if numel(svals) > 2
            svals = sprintf('%10.5g\t%10.5g\n', svals(2:end));
        else
            svals = '';
        end
        ch.Text.Values.String = sprintf('Values at cursor:%s\n%s%s', ...
            sval1, svals, sval);
    end
end

% plot PLP points ?
if numel(cc.plp) == 1 && ...
    isxff(cc.plp, 'plp')

    % prepare handles array and config
    ph1 = [];
    plppcfg = struct( ...
        'bcolor',  cc.plpcfg.bcolor, ...
        'color',   cc.plpcfg.color, ...
        'cond',    cc.plpcfg.cond, ...
        'dist',    cc.plpcfg.range, ...
        'slice',   cpos, ...
        'symsize', cc.plpcfg.symsize);

    % test color, symbol size
    if ~any(strcmpi(cc.plp.ColumnNames, plppcfg.color)) || ...
        numel(unique(cc.plp.(plppcfg.color))) == 1
        plppcfg.color = ne_gcfg.c.ini.PLPPlot.DefaultColor;
    end
    if ~any(strcmpi(cc.plp.ColumnNames, plppcfg.symsize)) || ...
        any(cc.plp.(plppcfg.symsize) < 1)
        plppcfg.symsize = -ne_gcfg.c.ini.PLPPlot.DefaultSymbolSize;
    end

    % label plot points
    if ~isempty(cc.plp.RunTimeVars.Config.LabelColumn)
        plppcfg.label = cc.plp.RunTimeVars.Config.LabelColumn;
        plppcfg.labcolor = cc.plp.RunTimeVars.Config.LabelColor;
    end

    % depending on page (all three)
    if cpg == 1

        % plot three different slices
        ph1 = [];
        th1 = [];
        ph2 = [];
        th2 = [];
        ph3 = [];
        th3 = [];
        try
            plppcfg.axes = get(cc.tio.imSag.Handle, 'Parent');
            plppcfg.view = 'sa';
            [ph1, th1] = cc.plp.PlotOnSlice(plppcfg);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        try
            plppcfg.axes = get(cc.tio.imCor.Handle, 'Parent');
            plppcfg.view = 'co';
            [ph2, th2] = cc.plp.PlotOnSlice(plppcfg);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        try
            plppcfg.axes = get(cc.tio.imTra.Handle, 'Parent');
            plppcfg.view = 'ax';
            [ph3, th3] = cc.plp.PlotOnSlice(plppcfg);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

        % catenate handles
        ph1 = [ph1(:); ph2(:); ph3(:); th1(:); th2(:); th3(:)];

    % single slice view
    elseif cpg == 2

        % plot on single slice
        plppcfg.axes = get(cc.tio.imSlZ.Handle, 'Parent');
        switch showz
            case {1}
                plppcfg.view = 'sa';
            case {2}
                plppcfg.view = 'co';
            case {3}
                plppcfg.view = 'ax';
        end
        ph1 = [];
        th1 = [];
        try
            [ph1, th1] = cc.plp.PlotOnSlice(plppcfg);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ph1 = [ph1(:); th1(:)];
    end

    % remove invalid handles and then set handles for later deletion
    if ~isempty(ph1)
        ph1(isnan(ph1)) = [];
    end
    ne_gcfg.cc.(iSat).Config.plpph = ph1;
end

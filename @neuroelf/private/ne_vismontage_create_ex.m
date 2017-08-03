function varargout = ne_vismontage_create_ex(varargin)
% ne_vismontage_create_ex  - create montage image (actual process)
%
% FORMAT:       [m, ax, ma] = ne_vismontage_create_ex(SRC, EVT, options)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       options     1x1 (mandatory) struct with (mandatory!) fields
%        .atrans    boolean flag, anatomical background (0) is transparent
%        .atranscol 1x3 RGB (int) color behind transparent anatomical data
%        .blx       1x2 blocking factors (e.g. [9, 6] for 9-by-6 slices)
%        .brds      border spacing (around image)
%        .drc       1x1 double direction (1 := sag, 2 := cor, 3 := tra)
%        .drs       1x1 slicing step or 1xP slice positions
%        .filename  output filename (written if not empty and ~.showinfig)
%        .flp       boolean flag, flip slicing direction 
%        .flx       boolean flag, flip output (in output) in X direction
%        .fontcolor slice-labeling font color
%        .fontname  slice-labeling font name
%        .fontsize  slice-labeling font size
%        .frame     2x3 bounding box (frame) from within to sample
%        .imeth     interpolation method (for stats)
%        .imetha    interpolation method for anatomical data
%        .join      boolean flag, join maps (otherwise overlay)
%        .ppv       output pixel-per-voxel ratio (scaling)
%        .showinfig boolean flag, show output in new figure
%        .slcoord   slice-labeling with coordinate (text on axes in figure)
%        .slvar     1x1 anatomical xff object (may have Underlay set)
%        .stalp     statistical alpha (boundary, use 1 for default)
%        .stthr     1x2 lower and upper threshold for non-typical stats
%        .stvar     1x1 statistical xff object (VMP, etc.)
%        .stvix     1xM index into stvar's Maps or 3D volumes
%        .sws       boolean flag, show-while-slicing
%        .tpvol     time-point for volume-slicing (e.g. for VTC, etc.)
%
% Output fields:
%
%       m           HxWx3 montage image (blocked slices)
%       ax          axes handle (if created, use .showinfig = true)
%       ma          HxWx1 alpha image (requires .atrans = true)
%
% Example:
%
%   [mi, mx] = ne_vismontage_create_ex(0, 0, struct( ...
%       'atrans', false, 'atranscol', [0, 0, 0], 'blx', [6, 9], 'brds', 4, ...
%       'drc', 1, 'drs', 3, 'filename', '', 'flp', false, 'flx', false, ...
%       'fontcolor', [1, 1, 1], 'fontname', 'Helvetica', 'fontsize', 16, ...
%       'frame', [78, 82, 80; -78, -110, -64], 'imeth', 'cubic', ...
%       'imetha', 'cubic', 'join', true, 'ppv', 3, 'showinfig', true, ...
%       'slcoord', true, 'slvar', VMR_OBJECT, 'stalp', 1, 'stthr', [3, 9], ...
%       'stvar', VMP_OBJECT, 'stvix', [1, 2, 3], 'sws', false, 'tpvol',1));
%   neuroelf_gui('screenshot', get(mx, 'Parent'), 'montage.png', 'high-q');
%   delete(get(mx, 'Parent'));


% Version:  v1.1
% Build:    17050312
% Date:     May-03 2017, 12:56 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, 2017, Jochen Weber
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

% pre-set output
varargout = cell(1, nargout);

% check input
if nargin < 3 || ...
   ~isstruct(varargin{3}) || ...
    numel(varargin{3}) ~= 1
    return;
end
o = varargin{3};
if ~isfield(o, 'atrans') || ...
   ~isfield(o, 'atranscol') || ...
   ~isfield(o, 'blx') || ...
   ~isfield(o, 'brds') || ...
   ~isfield(o, 'drc') || ...
   ~isfield(o, 'drs') || ...
   ~isfield(o, 'filename') || ...
   ~isfield(o, 'flp') || ...
   ~isfield(o, 'flx') || ...
   ~isfield(o, 'fontcolor') || ...
   ~isfield(o, 'fontname') || ...
   ~isfield(o, 'fontsize') || ...
   ~isfield(o, 'frame') || ...
   ~isfield(o, 'imeth') || ...
   ~isfield(o, 'imetha') || ...
   ~isfield(o, 'join') || ...
   ~isfield(o, 'ppv') || ...
   ~isfield(o, 'showinfig') || ...
   ~isfield(o, 'slcoord') || ...
   ~isfield(o, 'slvar') || ...
   ~isfield(o, 'stalp') || ...
   ~isfield(o, 'stthr') || ...
   ~isfield(o, 'stvar') || ...
   ~isfield(o, 'stvix') || ...
   ~isfield(o, 'sws') || ...
   ~isfield(o, 'tpvol')
    return;
end

% get handle shortcuts
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
ci = ne_gcfg.c.ini;

% get frame construction
switch (o.drc)
    case {1}
        dr = 'sag';
        frc = [3, 2];
    case {2}
        dr = 'cor';
        frc = [3, 1];
    case {3}
        dr = 'tra';
        frc = [2, 1];
    otherwise
        return;
end

% and round the numbers
fr = o.frame;
ir = fr;
rf = round(fr);

% get size of requested box
fri = abs(diff(rf(:, :)));

% find cube that fits
frs = max(fri);

% make sure that for skipped voxels, this is an integer
if o.ppv < 1
    frs = round(1 / o.ppv) * ceil(o.ppv * frs);
end

% get required transimg size
tis = round(o.ppv * frs);

% patch frame
ir(2, :) = ir(1, :) - (frs - 0.01);

% then find the offset for the image
ofs = tis - ceil(o.ppv * (round(fr(2, frc)) - ir(2, frc)));
ofo = (1 + tis) - ofs;

% get slicing positions
if numel(o.drs) == 1
    slp = rf(2, o.drc):o.drs:rf(1, o.drc);
else
    slp = o.drs(:)';
end
slp(slp < -127 | slp > 128) = [];
nslp = numel(slp);

% flipping
if o.flp
    slp = slp(end:-1:1);
end

% get image size
ims = o.ppv * fri(frc) + o.brds;

% create transimg for slicing
tio = transimg(tis, tis);

% initialize image variable
mnt = uint8(0);
mnt(o.blx(2) * ims(1) + o.brds, o.blx(1) * ims(2) + o.brds, 3) = 0;

% initialize alpha map
if ~isempty(regexpi(o.filename, '\.png$')) && ...
    o.atrans
    mnta = single(0);
    mnta(o.blx(2) * ims(1) + o.brds, o.blx(1) * ims(2) + o.brds) = 0;
else
    mnta = [];
end

% transparency
if o.atrans
    color = o.atranscol;
    setbackground(tio, color);

    % depending on filetype
    if ~isempty(o.slvar) && ...
        isxff(o.slvar, 'vmr')
        tmax = 1 / 125;
    else
        tmax = 1 / 142;
    end

    % fill image with color
    mnt(:, :, 1) = color(1);
    mnt(:, :, 2) = color(2);
    mnt(:, :, 3) = color(3);
end

% hide montage figure
o.hFig.Visible = 'off';

% make progress bar visible
cprog = ne_progress(0, 0, {true, 0, 'Creating slice-montage'});
hPrg = ne_gcfg.h.Progress;

% initiate position for slicing and target image sub-pixel offset
cpos = [0, 0, 0];
xp = 1 + o.brds;
yp = 1 + o.brds;

% underlay
ulvar = [];
if isxff(o.slvar, true)
    slhnd = handles(o.slvar);
    if isfield(slhnd, 'Underlay') && ...
        numel(slhnd.Underlay) == 1 && ...
        isxff(slhnd.Underlay, {'fmr', 'hdr', 'head', 'vmr', 'vtc'})
        ulvar = slhnd.Underlay;
        urtv = ulvar.RunTimeVars;
        if isfield(urtv, 'SPMsn') && ...
            isstruct(urtv.SPMsn) && ...
            numel(urtv.SPMsn) == 1
            spmsnu = {'snmat', urtv.SPMsn};
        else
            spmsnu = {};
        end
    end
    rtv = o.slvar.RunTimeVars;
    if isfield(rtv, 'SPMsn') && ...
        isstruct(rtv.SPMsn) && ...
        numel(rtv.SPMsn) == 1
        spmsn = {'snmat', rtv.SPMsn};
    else
        spmsn = {};
    end
else
    rtv = struct;
end

% gray-scale lut
if isfield(rtv, 'GrayScaleLUT') && ...
    isequal(size(rtv.GrayScaleLUT), [256, 3])
    graylut = {'gcolblend', 'lut', 'gcollut', rtv.GrayScaleLUT};
elseif ~isempty(cc.graylut)
    graylut = {'gcolblend', 'lut', 'gcollut', cc.graylut};
else
    graylut = {};
end

% stats object OK
if ~isxff(o.stvar, true)
    o.stvix = [];
else
    srtv = o.stvar.RunTimeVars;
    if isfield(srtv, 'SPMsn') && ...
        isstruct(srtv.SPMsn) && ...
        numel(srtv.SPMsn) == 1
        spmsns = {'snmat', srtv.SPMsn};
    else
        spmsns = {};
    end
end

% iterate over slice positions
for slc = 1:nslp

    % set position
    cpos(o.drc) = slp(slc);

    % update progress bar
    hPrg.Progress((slc - 1) / nslp, ...
        sprintf('Creating montage... (slice %d/%d)', slc, numel(slp)));

    % slice SliceVar
    if isxff(o.slvar, true)

        % no underlay
        if isempty(ulvar)
            o.slvar.SliceToTransimg(cpos, tio, struct( 'dir', dr, 'frame', ir, ...
                'layers', 1, 'mapvol', o.tpvol, 'method', o.imetha, ...
                'trans', cc.strans, spmsn{:}, graylut{:}));

        % with underlay
        else

            % first slice underlay
            ulvar.SliceToTransimg(cpos, tio, struct( 'dir', dr, 'frame', ir, ...
                'layers', 1, 'mapvol', 1, 'method', o.imetha, ...
                'trans', cc.strans, spmsnu{:}));

            % then slice overlay
            o.slvar.SliceToTransimg(cpos, tio, struct( 'dir', dr, 'frame', ir, ...
                'layers', 2, 'mapvol', o.tpvol, 'method', o.imetha, ...
                'trans', cc.strans, spmsn{:}));

            % code copied from ne_setslicepos
            setlayerpixel(tio, 1, ...
                montagemix(tio.Layer(1).Pixel, tio.Layer(2).Pixel, cc.joinulay));
            dellayer(tio, 2);
        end

        % transparency ?
        if o.atrans

            % get data
            sdata = tio.Layer(1).Pixel;

            % make outer space mask
            smask = all(sdata < (0.166 / tmax), 3);

            % limit to shape outside
            msz = size(smask, 2);
            smask = ( ...
                floodfill3(smask,  1 ,  1 , 1, 'xyface') | ...
                floodfill3(smask, msz,  1 , 1, 'xyface') | ...
                floodfill3(smask,  1 , msz, 1, 'xyface') | ...
                floodfill3(smask, msz, msz, 1, 'xyface'));

            % and make a slightly smaller/fringe version
            smsks = squeeze(~dilate3d(shiftdim(~smask, -1)));
            smskf = smask;
            if o.ppv >= 1
                for frc = 1:ceil(sqrt(o.ppv + 1.5))
                    smskf = squeeze(dilate3d(shiftdim(smskf, -1)));
                end
            end
            smskf(smsks) = false;

            % create alpha map
            alpha = ones(size(smsks));

            % smooth other data
            [frx, fry] = find(smskf);
            if ~isempty(frx)
                for lc = 1:size(sdata, 3)
                    if o.ppv < 2
                        sdata(sub2ind(size(sdata), frx, fry, lc .* ones(numel(frx), 1))) = ...
                            flexinterpn_method(sdata(:, :, lc), [frx, fry], ...
                            'gauss', sqrt(o.ppv + 1.5));
                    else
                        sdata(sub2ind(size(sdata), frx, fry, lc .* ones(numel(frx), 1))) = ...
                            flexinterpn_method(sdata(:, :, lc), [frx, fry], ...
                            'gauss', (o.ppv + 1) ^ 0.666);
                    end
                end
            end

            % join masks again
            smskf = (smsks | smskf);

            % and voxels in fringe to a scaled gray value
            alpha(smskf) = tmax .* double(sdata(smskf));

            % set into transimg
            setlayeralpha(tio, 1, alpha);
        end

        % get filename
        [slpath, slname, slext] = fileparts(o.slvar.FilenameOnDisk);
    else
        slname = '';
        slext = 'blank';
    end

    % if valid StarsVar
    if ~isempty(o.stvix)

        % and create options struct
        slopt = struct('dir', dr, 'frame', ir, 'layers', 2, 'mapvol', o.stvix, ...
            'method', o.imeth, 'rgbalpha', o.stalp, 'trans', cc.strans, ...
            'type', 'rgb', spmsns{:});
        stvtyp = lower(o.stvar.Filetype);

        [slpath, stname, stext] = fileparts(o.stvar.FilenameOnDisk);

        % fill in further values for non-VMP vars
        if ~any(strcmp(stvtyp, {'cmp', 'glm', 'hdr', 'head', 'vmp'}))
            slopt.rgbcol = ne_gcfg.lut.Colors;
            slopt.rgbctails = double(ch.Stats.PosTail.Value ~= 0) + ...
                2 * double(ch.Stats.NegTail.Value ~= 0);
            slopt.rgblthr = o.stthr(1);
            slopt.rgbuthr = o.stthr(2);

            % name part
            if numel(o.stvix) ~= 1
                plur = 's';
            else
                plur = '';
            end

            % extend name
            slname = sprintf('%s%s (+ %d map%s from %s%s)', ...
                slname, slext, numel(o.stvix), plur, stname, stext);

        % at least extend name
        else

            % also Map.Name?
            if numel(o.stvix) == 1 && ...
                any(strcmp(stvtyp, {'cmp', 'vmp'}))
                slname = sprintf('%s%s (+ Map ''%s'' from %s%s)', ...
                    slname, stext, o.stvar.Map(o.stvix).Name, stname, stext);
            else
                if numel(o.stvix) ~= 1
                    plur = 's';
                else
                    plur = '';
                end
                slname = sprintf('%s%s (+ %d map%s from %s%s)', ...
                    slname, stext, numel(o.stvix), plur, stname, stext);
            end
        end

        % then slice to transimg object
        o.stvar.SliceToTransimg(cpos, tio, slopt);

        % if requested and necessary
        if o.join && ...
            numel(tio.Layer) > 2

            % apply "stats-map join" (BV-like display) to the three images
            joinlayers(tio, 2:numel(tio.Layer));
        end
    end

    % show output (rendered)
    if o.sws
        display(render(tio));

    % or render only
    else
        render(tio);
    end

    % then copy from rendered to target image
    if o.flx
        mnt(yp:yp+ofs(1)-1, xp+ofs(2)-1:-1:xp, :) = ...
            tio.Rendered(1:ofs(1), ofo(2):tis, :);
        if ~isempty(mnta)
            mnta(yp:yp+ofs(1)-1, xp+ofs(2)-1:-1:xp) = ...
                tio.Layer(1).Alpha(1:ofs(1), ofo(2):tis);
        end
    else
        mnt(yp:yp+ofs(1)-1, xp:xp+ofs(2)-1, :) = ...
            tio.Rendered(1:ofs(1), ofo(2):tis, :);
        if ~isempty(mnta)
            mnta(yp:yp+ofs(1)-1, xp:xp+ofs(2)-1) = ...
                tio.Layer(1).Alpha(1:ofs(1), ofo(2):tis);
        end
    end

    % and increase target offset counter(s)
    xp = xp + ims(2);
    if xp > (size(mnt, 2) - o.brds)
        xp = 1;
        yp = yp + ims(1);
    end
end

% close TIO display and delete object from memory
delete(tio);

% make figure visible again
o.hFig.Visible = 'on';

% show in figure ?
rs = get(0, 'ScreenSize');
ax = [];
if (o.showinfig || ...
    o.slcoord) && ...
    rs(3) > 16

    % create figure
    nf = figure;
    figure(nf);

    % create axes
    ax = axes('Parent', nf);

    % compute figure size from image
    mims = size(mnt);
    mims = mims([2, 1]);

    % compare to screen size
    rc = floor(0.5 * rs(3:4));
    rs = 2 * floor(0.45 * rs(3:4));

    % if image (either width/height) larger than screen
    if any(mims > rs)

        % reduce size to match available space
        di = max(mims ./ rs);
        np = [rc - ceil(0.5 * (mims / di)), ceil(mims / di)];

    % or simply use image size
    else
        np = [rc - ceil(0.5 * mims), mims];
    end

    % figure settings
    set(nf, 'Units', 'pixels');
    set(nf, 'Position', np);
    set(nf, 'Units', 'normalized', 'NumberTitle', 'off', 'Name', ...
        sprintf('Montage: %s', slname));
    drawnow;
    pause(0.1);
    drawnow;

    % create image
    set(0, 'CurrentFigure', nf);
    set(nf, 'CurrentAxes', ax);
    image(mnt);

    % labels
    if o.slcoord
        slt = zeros(1, nslp);
        xp = 1;
        yp = 1;
        for slc = 1:nslp
            slt(slc) = text(xp + 0.04 * ims(2), yp + 0.1 * ims(1), ...
                sprintf('%g', slp(slc)));
            xp = xp + ims(2);
            if xp > (size(mnt, 2) - o.brds)
                xp = 1;
                yp = yp + ims(1);
            end
        end
        set(slt, 'Color', o.fontcolor, 'FontName', o.fontname, ...
            'FontSize', o.fontsize);
    end

    % axes settings
    set(ax, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
    set(ax, 'Visible', 'off');

    % set up screen-shot keypress
    set(nf, 'WindowKeyPressFcn', @vmcex_shot);

    % make sure to update screen
    drawnow;

% or write to file
elseif ~isempty(o.filename)
    try
        q = {};
        if ~isempty(regexpi(o.filename, '\.jpe?g$'))
            q = {'Quality', 90};
        elseif ~isempty(regexpi(o.filename, '\.png$')) && ...
            o.atrans
            q = {'Alpha', double(mnta)};
        end
        imwrite(mnt, o.filename, q{:});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(warndlg('Error writing image file.', 'NeuroElf GUI - warning', 'modal'));
    end
end

% hide progress bar if indicated
ne_progress(0, 0, cprog);

% return if requested
if nargout > 0
    varargout{1} = mnt;
    if nargout > 1
        varargout{2} = ax;
        if nargout > 2
            varargout{3} = mnta;
        end
    end
end


% keypress screenshot
function vmcex_shot(src, ke, varargin)

% get Key and Modifier from keyboard event (see Matlab docu!)
kk = ke.Key;
mn = ke.Modifier;

% screenshot key combination
if numel(mn) == 1 && ...
    strcmpi(mn{1}, 'shift') && ...
    strcmpi(kk, 's')

    % create screenshot
    if strcmpi(get(src, 'Type'), 'figure')
        ne_screenshot(0, 0, src, '', 'high-q');
    else
        ne_screenshot(0, 0, src);
    end
end

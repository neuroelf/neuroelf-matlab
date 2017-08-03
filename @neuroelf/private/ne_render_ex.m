% FUNCTION ne_render_ex: create rendered image (actual process)
function varargout = ne_render_ex(varargin)

% Version:  v1.1
% Build:    16042608
% Date:     Apr-26 2016, 8:05 AM EST
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
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
% AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
ci = ne_gcfg.c.ini;

% pre-set output
varargout = cell(1, nargout);

% check input
if nargin < 3 || ...
   ~isstruct(varargin{3}) || ...
    numel(varargin{3}) ~= 1

    % check global fig
    if ~isstruct(ne_gcfg.fcfg.Render) || ...
        numel(ne_gcfg.fcfg.Render) ~= 1
        return;
    else
        o = ne_gcfg.fcfg.Render;
    end
else
    o = varargin{3};
end

% check options
if ~isfield(o, 'frame') || ...
    numel(o.frame) ~= 2 || ...
   ~isa(o.frame, 'double') || ...
    any(isinf(o.frame) | isnan(o.frame) | o.frame < -256 | o.frame > 256) || ...
    o.frame(2) < o.frame(1) || ...
   ~isfield(o, 'join') || ...
    numel(o.join) ~= 1 || ...
   ~islogical(o.join) || ...
   ~isfield(o, 'slfrom') || ...
    numel(o.slfrom) ~= 1 || ...
   ~isa(o.slfrom, 'double') || ...
    isinf(o.slfrom) || ...
    isnan(o.slfrom) || ...
   ~isfield(o, 'slstep') || ...
    numel(o.slstep) ~= 1 || ...
   ~isa(o.slstep, 'double') || ...
    isinf(o.slstep) || ...
    isnan(o.slstep) || ...
   ~isfield(o, 'slto') || ...
    numel(o.slto) ~= 1 || ...
   ~isa(o.slto, 'double') || ...
    isinf(o.slto) || ...
    isnan(o.slto) || ...
   ~any(sign(o.slto - o.slfrom) == [0, sign(o.slstep)]) || ...
   ~isfield(o, 'slvar') || ...
    numel(o.slvar) ~= 1 || ...
   ~isxff(o.slvar) || ...
   ~isfield(o, 'showinfig') || ...
   ~isfield(o, 'stsmooth') || ...
   ~isfield(o, 'stvar') || ...
   ~isfield(o, 'stvix')
    return;
end
slfrom = o.slfrom;
slto = o.slto;
slstep = abs(o.slstep) * sign(o.slto - o.slfrom) + double(o.slfrom == o.slto);
steps = slfrom:slstep:slto;
if diff(o.frame) < 2
    o.frame = [-128, 128];
    fframe = [128, 128, 128; -127.999, -127.999, -127.999];
else
    fframe = [max(o.frame); min(o.frame)] * ones(1, 3);
end
frame = fframe;
if ~isfield(o, 'actdepth') || ...
    numel(o.actdepth) ~= 1 || ...
   ~isa(o.actdepth, 'double') || ...
    isinf(o.actdepth) || ...
    isnan(o.actdepth) || ...
    o.actdepth < 0
    o.actdepth = ci.Render.ActivationDepth;
end
if o.actdepth > 16
    o.actdepth = 16;
end
if ~isfield(o, 'agalpha') || ...
    numel(o.agalpha) ~= 1 || ...
   ~isa(o.agalpha, 'double') || ...
    isinf(o.agalpha) || ...
    isnan(o.agalpha) || ...
    o.agalpha < 0 || ...
    o.agalpha > 1
    o.agalpha = ci.Render.GrayAlphaFactor;
end
if ~isfield(o, 'bgcol') || ...
    numel(o.bgcol) ~= 3 || ...
   ~isa(o.bgcol, 'double') || ...
    any(isinf(o.bgcol) | isnan(o.bgcol) | o.bgcol < 0)
    o.bgcol = ci.Render.BackColor;
end
if ~isfield(o, 'colblend') || ...
   ~ischar(o.colblend) || ...
   ~any(strcmpi(o.colblend(:)', {'hsv', 'lut', 'rgb'}))
    o.colblend = upper(ci.Render.ColorBlend);
else
    o.colblend = upper(o.colblend(:)');
end
if ~isfield(o, 'dist') || ...
    numel(o.dist) ~= 1 || ...
   ~isa(o.dist, 'double') || ...
    isnan(o.dist) || ...
    o.dist < max(o.slfrom, o.slto)
    o.dist = ci.Render.Distance;
end
if ~isfield(o, 'filename') || ...
   ~ischar(o.filename)
    o.filename = ci.Render.Filename;
else
    o.filename = ddeblank(o.filename(:)');
end
if isempty(o.filename)
    o.filewrt = false;
end
if ~isfield(o, 'filenmov') || ~ischar(o.filenmov)
    o.filenmov = ci.Render.FilenameMovie;
else
    o.filenmov = ddeblank(o.filenmov(:)');
end
if isempty(o.filenmov)
    o.filewmov = false;
end
if ~isfield(o, 'filewmov') || numel(o.filewmov) ~= 1 || ~islogical(o.filewmov)
    o.filewmov = false;
end
if ~isfield(o, 'filewrt') || numel(o.filewrt) ~= 1 || ~islogical(o.filewrt)
    o.filewrt = false;
end
o.filesrs = false;
if o.filewrt
    o.filesrs = ~isempty(regexpi(o.filename, '\%0\d+d'));
end
if ~isfield(o, 'galpha') || ...
    numel(o.galpha) ~= 1 || ...
   ~isa(o.galpha, 'double') || ...
    isinf(o.galpha) || ...
    isnan(o.galpha) || ...
    o.galpha < 0 || ...
    o.galpha > 1
    o.galpha = ci.Render.GrayAlpha;
end
if ~isfield(o, 'hicol') || ...
    numel(o.hicol) ~= 3 || ...
   ~isa(o.hicol, 'double') || ...
    any(isinf(o.hicol) | isnan(o.hicol) | o.hicol < 0)
    o.hicol = ci.Render.HighColor;
end
if ~isfield(o, 'imeth') || ...
   ~ischar(o.imeth) || ...
   ~any(strcmpi(o.imeth(:)', {'cubic', 'lanczos3', 'linear'}))
    o.imeth = ci.Render.StatsInterp;
else
    o.imeth = lower(o.imeth(:)');
end
if ~isfield(o, 'imetha') || ...
   ~ischar(o.imetha) || ...
   ~any(strcmpi(o.imetha(:)', {'cubic', 'lanczos3', 'linear'}))
    o.imetha = ci.Render.AnatInterp;
else
    o.imetha = lower(o.imetha(:)');
end
if ~isfield(o, 'joinulay') || ...
    numel(o.joinulay) ~= 1 || ...
   ~isa(o.joinulay, 'double') || ...
    isinf(o.joinulay) || ...
    isnan(o.joinulay) || ...
    o.joinulay < 0 || ...
    o.joinulay > 6
    o.joinulay = 5;
end
if ~isfield(o, 'locol') || ...
    numel(o.locol) ~= 3 || ...
   ~isa(o.locol, 'double') || ...
    any(isinf(o.locol) | isnan(o.locol) | o.locol < 0)
    o.locol = ci.Render.LowColor;
end
if ~isfield(o, 'movhdfrom') || ...
   ~isa(o.movhdfrom, 'double') || ...
    numel(o.movhdfrom) ~= 1 || ...
    isinf(o.movhdfrom) || ...
    isnan(o.movhdfrom) || ...
    o.movhdfrom < 1 || ...
    o.movhdfrom > 561
    o.movhdfrom = ci.Render.MovieHDFrom;
else
    o.movhdfrom = round(o.movhdfrom);
end
if ~isfield(o, 'proty') || ...
    numel(o.proty) ~= 1 || ...
   ~isa(o.proty, 'double') || ...
    isinf(o.proty) || ...
    isnan(o.proty) || ...
    o.proty < -60 || ...
    o.proty > 60
    o.proty = 0;
end
if ~isfield(o, 'protz') || ...
    numel(o.protz) ~= 1 || ...
   ~isa(o.protz, 'double') || ...
    isinf(o.protz) || ...
    isnan(o.protz) || ...
    o.protz < -60 || ...
    o.protz > 60
    o.protz = 0;
end
if o.proty == 0 && ...
    o.protz == 0
    prot = false;
else
    prot = true;
    o.dist = min(5000, max(250, o.dist));
end
if ~isfield(o, 'res') || ...
    numel(o.res) ~= 1 || ...
   ~isa(o.res, 'double') || ...
    isinf(o.res) || ...
    isnan(o.res) || ...
    o.res < 32 || ...
    o.res > 4096
    o.res = ci.Render.Resolution;
end
res = o.res;
if ~isfield(o, 'roty') || ...
    numel(o.roty) ~= 1 || ...
   ~isa(o.roty, 'double') || ...
    isinf(o.roty) || ...
    isnan(o.roty)
    o.roty = 0;
end
if ~isfield(o, 'rotz') || ...
    numel(o.rotz) ~= 1 || ...
   ~isa(o.rotz, 'double') || ...
    isinf(o.rotz) || ...
    isnan(o.rotz)
    o.rotz = 180;
end
if ~isfield(o, 'slavar') || ...
    numel(o.slavar) ~= 1 || ...
   ~isxff(o.slavar, {'hdr', 'head', 'nlf', 'vmr', 'vtc'})
    o.slavar = [];
end
if ~isfield(o, 'stalp') || ...
    numel(o.stalp) ~= 1 || ...
   ~isa(o.stalp, 'double') || ...
    isinf(o.stalp) || ...
    isnan(o.stalp) || ...
    o.stalp > 1
    o.stalp = ci.Render.StatsAlpha;
end
if ~isfield(o, 'stminaalp') || ...
    numel(o.stminaalp) ~= 1 || ...
   ~isa(o.stminaalp, 'double') || ...
    isinf(o.stminaalp) || ...
    isnan(o.stminaalp) || ...
    o.stminaalp < 0 || ...
    o.stminaalp > 1
    o.stminaalp = ci.Render.StatsMinAnaAlpha;
end
if ~isfield(o, 'stmulaalp') || ...
    numel(o.stmulaalp) ~= 1 || ...
   ~islogical(o.stmulaalp)
    o.stmulaalp = ci.Render.StatsMultWithAnaAlpha;
end

% this method requires a valid slicing var
slvar = o.slvar;
if numel(slvar) ~= 1 || ...
   ~isxff(slvar, {'hdr', 'head', 'nlf', 'vmr', 'vtc'})
    return;
end
rtv = slvar.RunTimeVars;
if isfield(rtv, 'GrayScaleLUT') && ...
    isequal(size(rtv.GrayScaleLUT), [256, 3])
    graylut = {'gcollut', rtv.GrayScaleLUT};
else
    graylut = {};
end
if isfield(rtv, 'SPMsn') && ...
    isstruct(rtv.SPMsn) && ...
    numel(rtv.SPMsn) == 1
    spmsn = {'snmat', rtv.SPMsn};
else
    spmsn = {};
end
slhnd = handles(slvar);
if isfield(slhnd, 'Underlay') && ...
    numel(slhnd.Underlay) == 1 && ...
    isxff(slhnd.Underlay, {'hdr', 'head', 'vmr', 'vtc'})
    ulvar = slhnd.Underlay;
    urtv = ulvar.RunTimeVars;
    if isfield(urtv, 'SPMsn') && ...
        isstruct(urtv.SPMsn) && ...
        numel(urtv.SPMsn) == 1
        spmsnu = {'snmat', urtv.SPMsn};
    else
        spmsnu = {};
    end
else
    ulvar = [];
end
stvar = o.stvar;
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, {'cmp', 'hdr', 'head', 'vmp'})
    stvar = [];
    stvix = [];
else
    stvtyp = lower(stvar.Filetype);
    stvix = o.stvix;
    if ~isa(stvix, 'double') || ...
        isempty(stvix) || ...
        any(isinf(stvix(:)) | isnan(stvix(:)) | stvix(:) < 1) || ...
        numel(unique(round(stvix(:)))) ~= numel(stvix)
        stvix = [];
        stvar = [];
    else
        stvix = round(stvix(:)');
    end
    srtv = o.stvar.RunTimeVars;
    if isfield(srtv, 'SPMsn') && ...
        isstruct(srtv.SPMsn) && ...
        numel(srtv.SPMsn) == 1
        spmsns = {'snmat', srtv.SPMsn};
    else
        spmsns = {};
    end
end

% check fields
if numel(o.showinfig) ~= 1 || ...
   ~isa(o.showinfig, 'logical')
    o.showinfig = false;
end
if numel(o.slavar) ~= 1 || ...
   ~isxff(o.slavar, true)
    o.slavar = [];
end
if slfrom < slto
    slfromb = slfrom;
    slfrom = slto;
    slto = slfromb;
end
if numel(o.stsmooth) ~= 1 || ...
   ~isa(o.stsmooth, 'double') || ...
    isinf(o.stsmooth) || ...
    isnan(o.stsmooth) || ...
    o.stsmooth <= 0
    o.stsmooth = 0;
    stsmooth = [];
else
    if ~isempty(stvar)
        switch (stvtyp)
            case {'cmp', 'vmp'}
                mres = stvar.Resolution;
            case {'hdr', 'head'}
                cfr = stvar.CoordinateFrame;
                mres = prod(cfr.Resolution) .^ (1 / 3);
        end
        stsmooth = smoothkern(min(3 * mres, o.stsmooth / mres), 0.004);
        if numel(stsmooth) < 3
            stsmooth = [];
        else
            stsmooth = {{stsmooth, stsmooth, [0;1;0]}, {1, 1, 1}};
        end
    end
end

% update config
ci.Render.ActivationDepth = o.actdepth;
ci.Render.AnatInterp = o.imetha;
ci.Render.BackColor = o.bgcol;
ci.Render.ColorBlend = upper(o.colblend);
ci.Render.Distance = o.dist;
ci.Render.Filename = o.filename;
ci.Render.FilenameMovie = o.filenmov;
ci.Render.FileWrite = o.filewrt;
ci.Render.FileWriteMovie = o.filewmov;
ci.Render.Frame = o.frame;
ci.Render.GrayAlpha = o.galpha;
ci.Render.GrayAlphaFactor = o.agalpha;
ci.Render.HighColor = o.hicol;
ci.Render.LowColor = o.locol;
ci.Render.MovieHDFrom = o.movhdfrom;
ci.Render.Resolution = o.res;
ci.Render.ShowWhileRendering = o.showinfig;
ci.Render.SliceFrom = slfrom;
ci.Render.SliceStep = slstep;
ci.Render.SliceTo = slto;
ci.Render.SmoothStats = o.stsmooth;
ci.Render.StatsAlpha = o.stalp;
ci.Render.StatsInterp = o.imeth;
ci.Render.StatsMinAnaAlpha = o.stminaalp;
ci.Render.StatsMultWithAnaAlpha = o.stmulaalp;

% hide dialog
if isstruct(ne_gcfg.h.Render) && ...
    isxfigure(ne_gcfg.h.Render.RendFig, true)
    ne_gcfg.h.Render.RendFig.Visible = 'off';
end

% echo
if ne_gcfg.c.echo
    ne_echo({'neuroelf_gui(''render_ex'', %s)', any2ascii(o)});
end

% open movie file
if o.filewmov
    [moviefp, moviefn, moviefx] = fileparts(o.filenmov);
    if isempty(moviefp)
        moviefp = pwd;
    end
    moviefn = [moviefp '/' moviefn moviefx];
    mobj = VideoWriter(moviefn, 'MPEG-4');
    open(mobj);
end

% parse rotation
if mod(o.roty, 360) ~= 0 || ...
    mod(o.rotz, 360) ~= 180
    rot = {'trans', ...
        spmtrf([0, 0, 0], [0, 0, -pi * mod(o.rotz + 180, 360) / 180]) * ...
        spmtrf([0, 0, 0], [0, pi * mod(o.roty, 360) / 180, 0])};
else
    rot = {};
end

% post rotational stuff
if prot
    frame = fframe;
    resh = 0.5 * res;
    res2 = res * res;
    res22 = [res, res];
    resfac = res / abs(diff(frame(:, 1)));
    o.pdist = o.dist * resfac;
    hfov = atan((resh - 0.5) / o.pdist);
    hfovs = 2 * hfov / (res - 1);
    px = o.pdist .* tan(-hfov:hfovs:(hfov + 0.5 * hfovs));
    [px, py] = ndgrid(px, px);
    px = [px(:), py(:), o.pdist .* ones(numel(px), 1)]';
    py = 1 / sqrt(sum(px(:, 1) .* px(:, 1)));
    px = py .* px;
    cpo = [0, 0, o.pdist, 1];
    prt = (spmtrf([0, 0, o.pdist]) * ...
        spmtrf([0, 0, 0], [0, pi * mod(360 - o.proty, 360) / 180, 0]) * ...
        spmtrf([0, 0, 0], [pi * mod(o.protz, 360) / 180, 0, 0]) * ...
        spmtrf([0, 0, -o.pdist]));
    prtt = prt';
end

% create new transimg object
ti = transimg(res, res, o.bgcol);
setlayer(ti, 1, zeros(res), zeros(res));
render(ti);
ltrd = ti.Rendered;
tia = transimg(res, res);
setlayer(tia, 1, zeros(res), 1);

% compute depth factor
dfac = min(1, single((0.5 ^ (abs(slstep) / o.actdepth)) / o.galpha));

% get handle shortcut
hPrg = ne_gcfg.h.Progress;

% initialize progress
cprog = ne_progress(0, 0, {true, 0, 'Rendering brain...'});

% loop along steps
frnum = 1;
for r = 1:numel(steps)

    % get step
    rs = steps(r);

    % progress
    hPrg.Progress((r - 1) / numel(steps), sprintf( ...
        'Rendering brain... dist=%g, slice %d/%d', ...
        o.dist + rs, r, numel(steps)));

    % allow to break to function
    if any(strcmp(ne_gcfg.c.breakcb, 'render_ex'))
        break;
    end

    % render ti, so we know which voxels are to be kept
    if o.showinfig
        display(render(ti));

        % on first pass
        if r == 1

            % set up break handler
            set(ti.Handle, ...
                'DeleteFcn', 'neuroelf_gui(''breakcb'', ''render_ex'');');

            % and if OutputPosition is valid
            tifh = get(get(ti.Handle, 'Parent'), 'Parent');
            set(tifh, 'Units', 'pixels');
            fpos = get(tifh, 'Position');
            if any(ci.Render.OutputPosition ~= -1)
                fpos(1:2) = ci.Render.OutputPosition;
                set(tifh, 'Position', fpos);
            end
        end

        % update screen
        drawnow;
    elseif dfac < 1
        render(ti);
    end
    if dfac < 1
        trd = hsvconv(ti.Rendered, 2);
        trd = 1 - dfac .* limitrangec(trd(:, :, 2), 0, 1, 0);
    end

    % for last slice, increase grayalpha!
    if r == numel(steps)
        o.galpha = 2 * o.galpha;
    end

    % distance?
    if ~isinf(o.dist) && ...
       ~prot

        % compute factor, so that coordinate 0 is full frame
        frame = ((o.dist + rs) / o.dist) .* fframe;
    end

    % no underlay
    if isempty(ulvar)

        % simply slice SliceVar
        slvar.SliceToTransimg([rs, 0, 0], ti, struct( ...
            'dir',       'sag', ...
            'gcolblend', o.colblend, ...
            'gcolhigh',  o.hicol, ...
            'gcollow',   o.locol, ...
            'frame',     frame, ...
            'layers',    2, ...
             rot{:},     ...
            'method',    o.imetha, ...
            graylut{:},  ...
            spmsn{:}));

    % with underlay
    else

        % first slice Underlay then SliceVar
        ulvar.SliceToTransimg([rs, 0, 0], ti, struct( ...
            'dir',       'sag', ...
            'gcolblend', o.colblend, ...
            'gcolhigh',  o.hicol, ...
            'gcollow',   o.locol, ...
            'frame',     frame, ...
            'layers',    2, ...
             rot{:},     ...
            'method',    o.imetha, ...
            spmsnu{:}));
        slvar.SliceToTransimg([rs, 0, 0], ti, struct( ...
            'dir',       'sag', ...
            'gcolblend', o.colblend, ...
            'gcolhigh',  o.hicol, ...
            'gcollow',   o.locol, ...
            'frame',     frame, ...
            'layers',    3, ...
            rot{:},      ...
            'method',    o.imetha, ...
            spmsn{:}));

        % taken from ne_setslicepos
        setlayerpixel(ti, 2, montagemix( ...
            ti.Layer(2).Pixel, ti.Layer(3).Pixel, o.joinulay));
        dellayer(ti, 3);
    end

    % replace alpha channel
    if isempty(o.slavar)
        galp = o.galpha;
        if isempty(ulvar)
            slvar.SliceToTransimg([rs, 0, 0], tia, struct( ...
                'dir',       'sag', ...
                'frame',     frame, ...
                'grayalpha', galp, ...
                'layers',    1, ...
                rot{:},      ...
                'method',    o.imetha, ...
                spmsn{:}));
        else
            ulvar.SliceToTransimg([rs, 0, 0], tia, struct( ...
                'dir',       'sag', ...
                'frame',     frame, ...
                'grayalpha', galp, ...
                'layers',    1, ...
                rot{:},      ...
                'method',    o.imetha, ...
                spmsnu{:}));
        end
    else
        galp = o.agalpha;
        o.slavar.SliceToTransimg([rs, 0, 0], tia, struct( ...
            'dir',       'sag', ...
            'frame',     frame, ...
            'grayalpha', galp, ...
            'layers',    1, ...
             rot{:},     ...
            'method',    o.imetha));
    end
    setlayeralpha(ti, 2, tia.Layer(1).Alpha);

    % set layer of formerly colored pixel to lower alpha
    if dfac < 1
        setlayeralpha(ti, 2, trd .* ti.Layer(2).Alpha);
    end

    % and if required
    if ~isempty(stvix)

        % alpha required for multiplication or restriction
        if o.stmulaalp || ...
            o.stminaalp > 0

            % get information to alpha multiply stats information
            stga = limitrangec(single(1 / galp) .* tia.Layer(1).Alpha, 0, 1, 0);

            % restrict
            if o.stminaalp > 0

                % and multiply
                if o.stmulaalp
                    stga(stga < o.stminaalp) = 0;

                % restrict only
                else
                    stga = single(stga >= o.stminaalp);
                end
            end

        % otherwise
        else

            % simply use 1
            stga = single(1);
        end
    end

    % slice selected maps
    for mc = 1:numel(stvix)
        stvar.SliceToTransimg([rs, 0, 0], ti, struct( ...
            'dir',      'sag', ...
            'frame',    frame, ...
            'mapvol',   stvix(mc), ...
            'layers',   2 + mc, ...
            'method',   o.imeth, ...
            'rgbalpha', o.stalp, ...
            rot{:},     ...
            'type',     'rgb', ...
            spmsns{:}));

        % and ensure that no values are shown beyond gray matter
        if isempty(stsmooth)
            setlayeralpha(ti, 2 + mc, stga .* ti.Layer(2 + mc).Alpha);
        else
            til = ti.Layer(2 + mc);
            tis = [Inf, Inf, Inf; 1, 1, 1; 1, 1, 1; size(til.Pixel)];
            setlayer(ti, 2 + mc, ...
                flexinterpn(til.Pixel, tis, stsmooth{:}, 0), ...
                flexinterpn(stga .* til.Alpha, tis(:, 1:2), stsmooth{1}{1}, 1, 0));
        end
    end

    % join layers of stats first
    if numel(stvix) > 1
        joinlayers(ti, 3:numel(ti.Layer), ~o.join);
    end

    % post-rotation
    if prot

        % compute sampling coordinates by raytracing -> central point
        cpo(3) = o.pdist + resfac * rs;
        cp = cpo * prtt;
        cp = cp(:);
        [sn, sx, sy] = isectlp(zeros(3, res2), px, cp(1:3, ones(1, res2)), ...
            prt(1:3, ones(1, res2)), prtt(1:3, 2 .* ones(1, res2)));
        sx = (resh + 0.5) + sx;
        sy = (resh + 0.5) + sy;

        % resample layers
        for lc = 2:numel(ti.Layer)
            alay = ti.Layer(lc);
            if size(alay.Pixel, 3) == 1
                setlayer(ti, lc, ...
                    reshape(flexinterpn_method(alay.Pixel, [sx, sy], 'cubic'), res22), ...
                    reshape(flexinterpn_method(alay.Alpha, [sx, sy], 'cubic'), res22));
            else
                setlayer(ti, lc, cat(3, ...
                    reshape(flexinterpn_method(alay.Pixel(:, :, 1), [sx, sy], 'cubic'), res22), ...
                    reshape(flexinterpn_method(alay.Pixel(:, :, 2), [sx, sy], 'cubic'), res22), ...
                    reshape(flexinterpn_method(alay.Pixel(:, :, 3), [sx, sy], 'cubic'), res22)), ...
                    reshape(flexinterpn_method(alay.Alpha, [sx, sy], 'cubic'), res22));
            end
        end

    end

    % then join layer 1, 2, 3
    joinlayers(ti, 1:3);

    % write movie
    trd = ti.Rendered;
    if (o.filewmov || o.filesrs) && ~isequal(trd, ltrd)
        if o.filewmov
            if size(trd, 1) == 1280
                writeVideo(mobj, trd(o.movhdfrom:o.movhdfrom+719, :, :));
            else
                writeVideo(mobj, trd);
            end
        end
        if o.filesrs
            if numel(o.filename) > 4 && strcmpi(o.filename(end-3:end), '.jpg')
                imwrite(trd, sprintf(o.filename, frnum), 'Quality', 90);
            else
                imwrite(trd, sprintf(o.filename, frnum));
            end
            frnum = frnum + 1;
        end
        ltrd = trd;
    end
end

% get final rendering
render(ti);
trd = ti.Rendered;
if o.filewmov
    if size(trd, 1) == 1280
        writeVideo(mobj, trd(o.movhdfrom:o.movhdfrom+719, :, :));
    else
        writeVideo(mobj, trd);
    end
    close(mobj);
end
varargout{1} = trd;

% or write to file
if o.filewrt && ~isempty(o.filename) && ~o.filesrs
    try
        if numel(o.filename) > 4 && strcmpi(o.filename(end-3:end), '.jpg')
            imq = {'Quality', 90};
        elseif numel(o.filename) > 4 && ...
            strcmpi(o.filename(end-3:end), '.png')
            imq = {'Alpha', double(ti.Layer(1).Alpha)};
        end
        imwrite(trd, o.filename, imq{:});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(warndlg('Error writing image file.', 'NeuroElf GUI - warning', 'modal'));
    end
end

% keep track of last position
if o.showinfig
    try
        set(tifh, 'Units', 'pixels');
        fpos = get(tifh, 'Position');
        ci.Render.OutputPosition = fpos(1:2);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% close ti display and delete object from memory
delete(ti);
delete(tia);

% re-set break feature
ne_gcfg.c.breakcb(strcmp(ne_gcfg.c.breakcb, 'render_ex')) = [];

% hide progress bar if indicated
ne_progress(0, 0, cprog);

% show again
if o.showinfig

    % create figure and show
    f = figure;
    figure(f);
    drawnow;
    image(trd);
    drawnow;
end

% unhide dialog
if isstruct(ne_gcfg.h.Render) && ...
    isxfigure(ne_gcfg.h.Render.RendFig, true)
    ne_gcfg.h.Render.RendFig.Visible = 'on';
end

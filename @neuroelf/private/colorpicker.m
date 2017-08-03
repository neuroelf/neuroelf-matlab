function [col, cancelled] = colorpicker(col, cnames)
% colorpicker  - pick a color from the RGB/HSV color space
%
% FORMAT:       color = colorpicker([icolor [, icnames]]);
%
% Input fields:
%
%       icolor      initial color or colors (default: [0, 0, 0]);
%       icnames     color name(s) (labels)
%
% Output fields:
%
%       color       Cx3 RGB color (values in the [0...255] range)

% Version:  v0.9d
% Build:    14062315
% Date:     Jun-23 2014, 3:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2012, 2014, Jochen Weber
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

% variable for UI stuff
global ne_ui;

% argument check
if nargin > 0 && ...
   (~isa(col, 'double') || ...
    isempty(col) || ...
    size(col, 2) ~= 3 || ... || ...
    size(col, 1) > 256 || ...
    any(isinf(col(:)) | isnan(col(:)) | col(:) < 0 | col(:) > 255))
    error( ...
        'neuroelf:InvalidCall', ...
        'Invalid call to colorpicker.' ...
    );
end
if nargin == 0
    col = [0, 0, 0];
end
ncol = size(col, 1);
gridbase = 1;
sgrid = [1, 0];
sprod = 1;
if ncol == 1
    if nargin < 2 || ...
       (~ischar(cnames) && ...
        (~iscell(cnames) || ...
         numel(cnames) ~= 1 || ...
         ~ischar(cnames{1})))
        cnames = {''};
    elseif ischar(cnames)
        cnames = {cnames(:)'};
    else
        cnames{1} = cnames{1}(:)';
    end
else
    if nargin < 2 || ...
       ~iscell(cnames) || ...
        numel(cnames) ~= ncol
        cnames = cell(ncol, 1);
    else
        cnames = cnames(:);
    end
    for cnc = 1:ncol
        if ~ischar(cnames{cnc})
            cnames{cnc} = sprintf('color %d', cnc);
        elseif numel(cnames{cnc}) > 64
            cnames{cnc} = sprintf('%s...%s', ...
                cnames{cnc}(1:31), cnames{cnc}(end-30:end));
        else
            cnames{cnc} = cnames{cnc}(:)';
        end
    end
end

% load figure
shho = get(0, 'ShowHiddenHandles');
try
    hFig = xfigure([neuroelf_path('tfg') '/colorpicker.tfg']);
    hTag = hFig.TagStruct;
    set(0, 'ShowHiddenHandles', 'off');
    hTag.IM_colorpicker_SV_im = hTag.IM_colorpicker_SV.Children(1);
    hTag.IM_colorpicker_H_im = hTag.IM_colorpicker_H.Children(1);
    hTag.IM_colorpicker_R_im = hTag.IM_colorpicker_R.Children(1);
    hTag.IM_colorpicker_G_im = hTag.IM_colorpicker_G.Children(1);
    hTag.IM_colorpicker_B_im = hTag.IM_colorpicker_B.Children(1);
    hTag.IM_colorpicker_Pick_im = hTag.IM_colorpicker_Pick.Children(1);
    set(0, 'ShowHiddenHandles', shho);
catch ne_eo;
    set(0, 'ShowHiddenHandles', shho);
    error( ...
        'neuroelf:xfigureError', ...
        'Error creating UI for colorpicker: %s.', ...
        ne_eo.message ...
    );
end

% initial color data for picked color(s)
icd = repmat(uint8(204 .* ones(sprod, sprod)), [1, 1, 3]);

% extend grid (for multiple colors)
if ncol > 1
    gridform = [ 1, 43, 28, 20, 16, 13, 11, 10, 9, 8, 7, 6, 6, 5, 5, 4];
    gridline = [ 0,  2,  2,  2,  2,  2,  2,  1, 1, 1, 1, 1, 1, 1, 1, 1];
    gridbase = ceil(sqrt(ncol));
    sgrid = [gridform(gridbase), gridline(gridbase)];
    sprod = gridbase * sgrid(1) + (gridbase - 1) * sgrid(2);
    icd = repmat(uint8(204 .* ones(sprod, sprod)), [1, 1, 3]);

    % initialize colors more sophisticatedly
    for cc = 1:ncol

        % get from/to arguments
        xpos = mod(cc - 1, gridbase);
        ypos = ((cc - 1) - xpos) / gridbase;
        xfrom = xpos * sum(sgrid) + 1;
        xto = xfrom + sgrid(1) - 1;
        yfrom = ypos * sum(sgrid) + 1;
        yto = yfrom + sgrid(1) - 1;

        % fill in color
        icd(yfrom:yto, xfrom:xto, 1) = col(cc, 1);
        icd(yfrom:yto, xfrom:xto, 2) = col(cc, 2);
        icd(yfrom:yto, xfrom:xto, 3) = col(cc, 3);
    end

    % add "s" to label
    hTag.TX_colorpicker_Pick.String = 'Picked colors';
else
    icd(:) = col(:);
end

% perform some computations already
h = (1/512):1/256:1;
[v, s] = ndgrid((1-0.002):-1/256:0, h);
hsv = [zeros(65536, 1), s(:), v(:)];

% initialize global variable
col = round(real(col));
chsv = hsvconv(uint8(col), 2);
ne_ui.colorpicker = struct( ...
    'btst',  false, ...
    'cncl',  true, ...
    'ccol',  col, ...
    'cname', {cnames}, ...
    'col',   col, ...
    'chsv',  chsv, ...
    'grad',  (0:255)', ...
    'grid',  gridbase, ...
    'hndo',  'rgbhsv', ...
    'hcol',  chsv, ...
    'hFig',  hFig, ...
    'hTag',  hTag, ...
    'hsv',   reshape(hsv, [256, 256, 3]), ...
    'huei',  hsvconv(cat(3, h(:), ones(256, 1, 2))), ...
    'ncol',  ncol, ...
    'Pos',   [6, 265; 270, 301; 310, 341; 346, 377; 382, 413], ...
    'Poso',  'xhrgb', ...
    'scol',  1, ...
    'sgrid', sgrid, ...
    'sPosX', hTag.IM_colorpicker_Pick.Position(1), ...
    'sPosY', hTag.IM_colorpicker_Pick.Position(2), ...
    'sprod', sprod);

% prepate SV image and H bar
set(hTag.IM_colorpicker_H_im, 'CData', ne_ui.colorpicker.huei);

% set initial color
set(hTag.IM_colorpicker_Pick_im, 'CData', icd);
if size(icd, 1) > 1 || ...
    size(icd, 2) > 1
    hTag.IM_colorpicker_Pick.XLim = 0.5 + [0, size(icd, 2)];
    hTag.IM_colorpicker_Pick.YLim = 0.5 + [0, size(icd, 1)];
end

% and label
hTag.TX_colorpicker_Name.String = sprintf('Pick a color: %s', cnames{1});

% set color to images
cp_setcolors;

% make sure axes are correctly order
hTag.IM_colorpicker_R.YDir = 'normal';
hTag.IM_colorpicker_G.YDir = 'normal';
hTag.IM_colorpicker_B.YDir = 'normal';
hTag.IM_colorpicker_H.YDir = 'normal';

% set callback hooks
hTag.ED_colorpicker_R.Callback = @cp_callback;
hTag.ED_colorpicker_G.Callback = @cp_callback;
hTag.ED_colorpicker_B.Callback = @cp_callback;
hTag.ED_colorpicker_H.Callback = @cp_callback;
hTag.ED_colorpicker_S.Callback = @cp_callback;
hTag.ED_colorpicker_V.Callback = @cp_callback;
hTag.BT_colorpicker_reset.Callback = @cp_reset;
hTag.BT_colorpicker_cancel.Callback = @cp_closeui;
hTag.BT_colorpicker_accept.Callback = @cp_accept;
hFig.WindowButtonDownFcn = @cp_buttondown;
hFig.WindowButtonMotionFcn = @cp_buttonmove;
hFig.WindowButtonUpFcn = @cp_buttonup;

% create list of matlab handles for fast access
ne_ui.colorpicker.hndl = [ ...
    hTag.ED_colorpicker_R.MLHandle, ...
    hTag.ED_colorpicker_G.MLHandle, ...
    hTag.ED_colorpicker_B.MLHandle, ...
    hTag.ED_colorpicker_H.MLHandle, ...
    hTag.ED_colorpicker_S.MLHandle, ...
    hTag.ED_colorpicker_V.MLHandle];

% make figure modal
hFig.HandleVisibility = 'callback';
hFig.Visible = 'on';
hFig.WindowStyle = 'modal';

% wait for figure to close
uiwait(hFig.mlhandle);

% assign outputs
col = ne_ui.colorpicker.col;
cancelled = ne_ui.colorpicker.cncl;

% and clean up
ne_ui.colorpicker = [];


% internal functions
function cp_setcolors(input, val)
global ne_ui;
cp = ne_ui.colorpicker;
hTag = ne_ui.colorpicker.hTag;
if nargin == 0
    input = 'c';
end
switch (lower(input(:)'))
    case {'r'}
        if nargin < 2
            try
                val = min(255, max(0, floor(str2double( ...
                    hTag.ED_colorpicker_R.String))));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                val = 0;
            end
        end
        ne_ui.colorpicker.ccol(cp.scol, 1) = val;
    case {'g'}
        if nargin < 2
            try
                val = min(255, max(0, floor(str2double( ...
                    hTag.ED_colorpicker_G.String))));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                val = 0;
            end
        end
        ne_ui.colorpicker.ccol(cp.scol, 2) = val;
    case {'b'}
        if nargin < 2
            try
                val = min(255, max(0, floor(str2double( ...
                    hTag.ED_colorpicker_B.String))));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                val = 0;
            end
        end
        ne_ui.colorpicker.ccol(cp.scol, 3) = val;
    case {'h'}
        if nargin < 2
            try
                val = min(255, max(0, floor(str2double( ...
                    hTag.ED_colorpicker_H.String))));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                val = 0;
            end
        end
        ne_ui.colorpicker.chsv(cp.scol, 1) = val / 255;
    case {'s'}
        if nargin < 2
            try
                val = min(255, max(0, floor(str2double( ...
                    hTag.ED_colorpicker_S.String))));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                val = 0;
            end
        end
        ne_ui.colorpicker.chsv(cp.scol, 2) = val / 255;
    case {'v'}
        if nargin < 2
            try
                val = min(255, max(0, floor(str2double( ...
                    hTag.ED_colorpicker_V.String))));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                val = 0;
            end
        end
        ne_ui.colorpicker.chsv(cp.scol, 3) = val / 255;
    case {'x'}
        if nargin < 2 || ...
            numel(val) ~= 2
            return;
        end
        ne_ui.colorpicker.chsv(cp.scol, 2:3) = val ./ 255;
    case {'t'}
        col = ne_ui.colorpicker.col;
        ne_ui.colorpicker.ccol = col;
        ne_ui.colorpicker.chsv = hsvconv(uint8(round(col)), 2);
        if cp.ncol > 1
            icd = repmat(uint8(204 .* ones(cp.sprod, cp.sprod)), [1, 1, 3]);
            for cc = 1:cp.ncol
                xpos = mod(cc - 1, cp.grid);
                ypos = ((cc - 1) - xpos) / cp.grid;
                xfrom = xpos * sum(cp.sgrid) + 1;
                xto = xfrom + cp.sgrid(1) - 1;
                yfrom = ypos * sum(cp.sgrid) + 1;
                yto = yfrom + cp.sgrid(1) - 1;
                icd(yfrom:yto, xfrom:xto, 1) = col(cc, 1);
                icd(yfrom:yto, xfrom:xto, 2) = col(cc, 2);
                icd(yfrom:yto, xfrom:xto, 3) = col(cc, 3);
            end
            set(hTag.IM_colorpicker_Pick_im, 'CData', icd);
        end
end

% get current color
if ~any('hsvx' == input)
    col = min(255, max(0, ne_ui.colorpicker.ccol(cp.scol, :)));
    hcol = hsvconv(uint8(round(col)), 2);
    ne_ui.colorpicker.chsv(cp.scol, :) = hcol;
else
    hcol = ne_ui.colorpicker.chsv(cp.scol, :);
    col = hsvconv(hcol);
    ne_ui.colorpicker.ccol(cp.scol, :) = col;
end
hcolo = floor(255.999 .* hcol);

% set picked color
if cp.ncol == 1
    set(hTag.IM_colorpicker_Pick_im, 'CData', uint8(shiftdim(col(:), -2)));
else
    icd = get(hTag.IM_colorpicker_Pick_im, 'CData');
    cc = cp.scol;
    xpos = mod(cc - 1, cp.grid);
    ypos = ((cc - 1) - xpos) / cp.grid;
    xfrom = xpos * sum(cp.sgrid) + 1;
    xto = xfrom + cp.sgrid(1) - 1;
    yfrom = ypos * sum(cp.sgrid) + 1;
    yto = yfrom + cp.sgrid(1) - 1;
    icd(yfrom:yto, xfrom:xto, 1) = col(1);
    icd(yfrom:yto, xfrom:xto, 2) = col(2);
    icd(yfrom:yto, xfrom:xto, 3) = col(3);
    set(hTag.IM_colorpicker_Pick_im, 'CData', icd);
end
hTag.TX_colorpicker_Name.String = sprintf('Pick a color: %s', cp.cname{cp.scol});

% put into edit fields
hTag.ED_colorpicker_R.String = sprintf('%d', col(1));
hTag.ED_colorpicker_G.String = sprintf('%d', col(2));
hTag.ED_colorpicker_B.String = sprintf('%d', col(3));
hTag.ED_colorpicker_H.String = sprintf('%d', hcolo(1));
hTag.ED_colorpicker_S.String = sprintf('%d', hcolo(2));
hTag.ED_colorpicker_V.String = sprintf('%d', hcolo(3));

% compute RGB/HSV gradient images
gcol = col(ones(1, 256), :);
grad = ne_ui.colorpicker.grad;
hsv = ne_ui.colorpicker.hsv;
hsv(:, :, 1) = hcol(1);
H_im = cp.huei;
SV_im = hsvconv(hsv);
R_im = uint8(cat(3, grad, gcol(:, 2), gcol(:, 3)));
G_im = uint8(cat(3, gcol(:, 1), grad, gcol(:, 3)));
B_im = uint8(cat(3, gcol(:, 1), gcol(:, 2), grad));

% mark current position
dcol = double(col) + 1;
dhcol = 1 + min(255, max(0, round(255 .* hcol)));
dhcol(3) = 257 - dhcol(3);
H_im(dhcol(1), :, :) = 0;
if sum(double(col)) < 386
    SV_im(max(1, dhcol(3)-4):min(256, dhcol(3)+4), dhcol(2), :) = 255;
    SV_im(dhcol(3), max(1, dhcol(2)-4):min(256, dhcol(2)+4), :) = 255;
    R_im(dcol(1), :, :) = 255;
    G_im(dcol(2), :, :) = 255;
    B_im(dcol(3), :, :) = 255;
else
    SV_im(max(1, dhcol(3)-4):min(256, dhcol(3)+4), dhcol(2), :) = 0;
    SV_im(dhcol(3), max(1, dhcol(2)-4):min(256, dhcol(2)+4), :) = 0;
    R_im(dcol(1), :, :) = 0;
    G_im(dcol(2), :, :) = 0;
    B_im(dcol(3), :, :) = 0;
end

% set images
set(hTag.IM_colorpicker_H_im, 'CData', H_im);
set(hTag.IM_colorpicker_SV_im, 'CData', SV_im);
set(hTag.IM_colorpicker_R_im, 'CData', R_im);
set(hTag.IM_colorpicker_G_im, 'CData', G_im);
set(hTag.IM_colorpicker_B_im, 'CData', B_im);

% UI functions
function cp_closeui(varargin)
global ne_ui;
ne_ui.colorpicker.hFig.Delete;

function cp_buttondown(varargin)
global ne_ui;
ne_ui.colorpicker.btst = true;
cp_buttonmove;

function cp_buttonup(varargin)
global ne_ui;
ne_ui.colorpicker.btst = false;

function cp_buttonmove(varargin)
global ne_ui;
cp = ne_ui.colorpicker;
if ~cp.btst
    return;
end
hFig = cp.hFig;
npos = hFig.CurrentPoint;
npos(2) = npos(2) - 40;
pos = cp.Pos;
fpos = find(pos(:, 1) <= npos(1) & pos(:, 2) >= npos(1));
if ~isempty(fpos)
    code = cp.Poso(fpos(1));
    if fpos(1) == 1
        val = [npos(1) - pos(fpos(1), 1), npos(2)];
    else
        val = npos(2);
    end
    val = min(255, max(0, val));
    cp_setcolors(code, val);
    return;
end
posx = cp.sPosX;
posy = cp.sPosY - 42;
posx(2) = posx + 88;
posy(2) = posy + 88;
if npos(1) > posx(1) && ...
    npos(1) <= posx(2) && ...
    npos(2) > posy(1) && ...
    npos(2) <= posy(2)
    npos = round((cp.sprod / 88) .* (npos - [posx(1), posy(1)]));
    npos(2) = min(cp.sprod - 1, max(0, cp.sprod - npos(2)));
    scol = 1 + floor(npos(1) / sum(cp.sgrid)) + ...
        cp.grid * floor(npos(2) / sum(cp.sgrid));
    if scol <= cp.ncol
        ne_ui.colorpicker.scol = scol;
        cp_setcolors('p');
    end
end

function cp_callback(varargin)
global ne_ui;
if nargin > 0 && ...
    ishandle(varargin{1})
    h = (ne_ui.colorpicker.hndl == varargin{1});
    if any(h)
        cp_setcolors(ne_ui.colorpicker.hndo(h));
    end
end

function cp_reset(varargin)
cp_setcolors('t');

function cp_accept(varargin)
global ne_ui;
ne_ui.colorpicker.cncl = false;
ne_ui.colorpicker.col = ne_ui.colorpicker.ccol;
ne_ui.colorpicker.hFig.Delete;

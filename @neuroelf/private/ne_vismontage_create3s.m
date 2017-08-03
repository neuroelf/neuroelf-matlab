% FUNCTION ne_vismontage_create3s: create 3 single-slice montage images
function [m, ax, malp] = ne_vismontage_create3s(varargin)

% Version:  v1.0
% Build:    15122814
% Date:     Dec-28 2015, 2:15 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, Jochen Weber
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
vm = cc.VisMontage;
tags = ch.VisMontage.h;

% only allow one concurrent call
if any(strcmp('vismontagecreate3s', ne_gcfg.c.blockcb))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'vismontagecreate3s';

% position
if nargin > 2 && ...
    isa(varargin{3}, 'double') && ...
    numel(varargin{3}) == 3 && ...
   ~any(isinf(varargin{3}) | isnan(varargin{3}) | varargin{3} < -128 | varargin{3} > 128)
    cpos = round(varargin{3}(:)');
else
    cpos = round(ne_gcfg.fcfg.cpos);
end

% original values
sldir = tags.DD_vismontage_dir.Value;
xfrom = tags.ED_vismontage_xfrom.String;
xstep = tags.DD_vismontage_xstep.Value;
yfrom = tags.ED_vismontage_yfrom.String;
ystep = tags.DD_vismontage_ystep.Value;
zfrom = tags.ED_vismontage_zfrom.String;
zstep = tags.DD_vismontage_zstep.Value;

% show in figure
showinfig = (tags.RB_vismontage_showinfig.Value > 0);
[ofilepath, ofilename, ofileext] = fileparts(ddeblank(tags.ED_vismontage_filename.String));
if isempty(ofilepath)
    ofilepath = '.';
end

% no slices
xax = [];
yax = [];
zax = [];

% try
try
    
    % X-slicing
    tags.DD_vismontage_dir.Value = 1;
    ne_vismontage_updateui(0, 0, 'dir');
    tags.ED_vismontage_xfrom.String = sprintf('%d', cpos(1));
    ne_vismontage_updateui(0, 0, 'xfrom');
    tags.DD_vismontage_xstep.Value = 11;
    ne_vismontage_updateui(0, 0, 'xstep');
    if ~showinfig
        tags.ED_vismontage_filename.String = sprintf('%s%s%s_X%d%s', ...
            ofilepath, filesep, ofilename, cpos(1), ofileext);
    end
    [xm, xax, xmalp] = ne_vismontage_create;
    tags.ED_vismontage_xfrom.String = xfrom;
    ne_vismontage_updateui(0, 0, 'xfrom');
    tags.DD_vismontage_xstep.Value = xstep;
    ne_vismontage_updateui(0, 0, 'xstep');

    % Y-slicing
    tags.DD_vismontage_dir.Value = 2;
    ne_vismontage_updateui(0, 0, 'dir');
    tags.ED_vismontage_yfrom.String = sprintf('%d', cpos(2));
    ne_vismontage_updateui(0, 0, 'yfrom');
    tags.DD_vismontage_ystep.Value = 11;
    ne_vismontage_updateui(0, 0, 'ystep');
    if ~showinfig
        tags.ED_vismontage_filename.String = sprintf('%s%s%s_Y%d%s', ...
            ofilepath, filesep, ofilename, cpos(2), ofileext);
    end
    [ym, yax, ymalp] = ne_vismontage_create;
    tags.ED_vismontage_yfrom.String = yfrom;
    ne_vismontage_updateui(0, 0, 'yfrom');
    tags.DD_vismontage_ystep.Value = ystep;
    ne_vismontage_updateui(0, 0, 'ystep');

    % Z-slicing
    tags.DD_vismontage_dir.Value = 3;
    ne_vismontage_updateui(0, 0, 'dir');
    tags.ED_vismontage_zfrom.String = sprintf('%d', cpos(3));
    ne_vismontage_updateui(0, 0, 'zfrom');
    tags.DD_vismontage_zstep.Value = 11;
    ne_vismontage_updateui(0, 0, 'zstep');
    if ~showinfig
        tags.ED_vismontage_filename.String = sprintf('%s%s%s_Z%d%s', ...
            ofilepath, filesep, ofilename, cpos(3), ofileext);
    end
    [zm, zax, zmalp] = ne_vismontage_create;
    tags.ED_vismontage_zfrom.String = zfrom;
    ne_vismontage_updateui(0, 0, 'zfrom');
    tags.DD_vismontage_zstep.Value = zstep;
    ne_vismontage_updateui(0, 0, 'zstep');
    tags.DD_vismontage_dir.Value = sldir;
    ne_vismontage_updateui(0, 0, 'dir');
catch ne_eo;
    tags.ED_vismontage_xfrom.String = xfrom;
    ne_vismontage_updateui(0, 0, 'xfrom');
    tags.DD_vismontage_xstep.Value = xstep;
    ne_vismontage_updateui(0, 0, 'xstep');
    tags.ED_vismontage_yfrom.String = yfrom;
    ne_vismontage_updateui(0, 0, 'yfrom');
    tags.DD_vismontage_ystep.Value = ystep;
    ne_vismontage_updateui(0, 0, 'ystep');
    tags.ED_vismontage_zfrom.String = zfrom;
    ne_vismontage_updateui(0, 0, 'zfrom');
    tags.DD_vismontage_zstep.Value = zstep;
    ne_vismontage_updateui(0, 0, 'zstep');
    tags.DD_vismontage_dir.Value = sldir;
    ne_vismontage_updateui(0, 0, 'dir');
    tags.ED_vismontage_filename.String = [ofilepath filesep ofilename ofileext];
    if ~isempty(xax) && ...
        ishandle(xax)
        delete(get(xax, 'Parent'));
    end
    if ~isempty(yax) && ...
        ishandle(yax)
        delete(get(yax, 'Parent'));
    end
    if ~isempty(zax) && ...
        ishandle(zax)
        delete(get(zax, 'Parent'));
    end
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'vismontagecreate3s')) = [];
    rethrow(ne_eo);
end

% combine
xsize = size(xm);
ysize = size(ym);
zsize = size(zm);
maxy = max([xsize(1), ysize(1), zsize(1)]);
m = uint8(0);
m(maxy, xsize(2) + ysize(2) + max(zsize(1:2)), 3) = 0;
if ~isempty(xmalp)
    malp = single(zeros(maxy, size(m, 2)));
else
    malp = single([]);
end
yoff = 1 + floor(0.5 * (maxy - xsize(1)));
m(yoff:yoff+xsize(1)-1, 1:xsize(2), :) = xm;
if ~isempty(malp)
    malp(yoff:yoff+xsize(1)-1, 1:xsize(2)) = xmalp;
end
xoff = xsize(2) + 1;
yoff = 1 + floor(0.5 * (maxy - ysize(1)));
m(yoff:yoff+ysize(1)-1, xoff:xoff+ysize(2)-1, :) = ym;
if ~isempty(malp)
    malp(yoff:yoff+ysize(1)-1, xoff:xoff+ysize(2)-1) = ymalp;
end
xoff = xoff + ysize(2) + floor(0.5 * (max(zsize(1:2)) - zsize(2)));
yoff = 1 + floor(0.5 * (maxy - zsize(1)));
m(yoff:yoff+zsize(1)-1, xoff:xoff+zsize(2)-1, :) = zm;
if ~isempty(malp)
    malp(yoff:yoff+zsize(1)-1, xoff:xoff+zsize(2)-1) = zmalp;
end
if showinfig
    xf = get(xax, 'Parent');
    yf = get(yax, 'Parent');
    zf = get(zax, 'Parent');
    set(xf, 'Units', 'pixels');
    set(yf, 'Units', 'pixels');
    set(zf, 'Units', 'pixels');
    set(xax, 'Units', 'pixels');
    set(yax, 'Units', 'pixels');
    set(zax, 'Units', 'pixels');
    xfpos = get(xf, 'Position');
    xfpos(3:4) = [size(m, 2), size(m, 1)];
    set(xf, 'Position', xfpos, 'Color', [0, 0, 0]);
    set(xax, 'Position', [0, floor(0.5 * (maxy - xsize(1))), xsize(2), xsize(1)]);
    set(yax, 'Parent', xf);
    delete(yf);
    set(yax, 'Position', [xsize(2), floor(0.5 * (maxy - ysize(1))), ysize(2), ysize(1)]);
    set(zax, 'Parent', xf);
    delete(zf);
    set(zax, 'Position', [xsize(2) + ysize(2), floor(0.5 * (maxy - zsize(1))), zsize(2), zsize(1)]);
end

% allow further montage creations
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'vismontagecreate3s')) = [];

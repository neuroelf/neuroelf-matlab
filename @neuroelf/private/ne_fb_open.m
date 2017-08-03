% FUNCTION ne_fb_open: loads the fractal browser UI
function varargout = ne_fb_open(varargin)

% Version:  v1.0
% Build:    14102818
% Date:     Oct-28 2014, 6:08 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% already open?
if isfield(ne_gcfg.h, 'FB') && ...
    isstruct(ne_gcfg.h.FB) && ...
    isfield(ne_gcfg.h.FB, 'FBFig') && ...
    numel(ne_gcfg.h.FB.FBFig) == 1 && ...
    isxfigure(ne_gcfg.h.FB.FBFig, true)
    ne_gcfg.h.FB.FBFig.Visible = 'on';
    figure(ne_gcfg.h.FB.FBFig.MLHandle);
    return;
end

% blocked?
if any(strcmp(ne_gcfg.c.blockcb, 'fb_open'))
    return;
end

% echo
if ne_gcfg.c.echo
    ne_echo('neuroelf_gui(''fb_open'');');
end

% load fractal browser
try
    hFig = xfigure([neuroelf_path('tfg') '/ne_fbrowse.tfg']);
    ne_gcfg.h.FB.FBFig = hFig;

    % block further callbacks
    ne_gcfg.c.blockcb{end+1} = 'fb_open';

    % get required controls
    tags = hFig.TagStruct;
    hTags = fieldnames(tags);
    hTag = hTags{end}(4:11);
    tags.Browser = tags.(['IM_' hTag '_Browse']);
    tags.Browser.Position(1:2) = 0;
    tags.Browser.YDir = 'normal';

    % set h struct
    ne_gcfg.h.FB.h = tags;
    
    % create config
    cc = struct;
    cc.down = struct('dPos', [0, 0], 'mPos', [0, 0], 'time', now, 'move', false);
    cc.tag = hTag;
    cc.xbounds = [-2.6, 1];
    cc.ybounds = [-1.8, 1.8];
    cc.orgres = size(tags.Browser.CData);
    cc.xyres = 2 .* cc.orgres(1:2);

    % gpuArray available
    try
        cc.maxiter = gpuArray(5000);
        cc.xspace = gpuArray.linspace(cc.xbounds(1), cc.xbounds(2), cc.xyres(1));
        cc.yspace = gpuArray.linspace(cc.ybounds(1), cc.ybounds(2), cc.xyres(2));

    % gpuArray unavailable
    catch
        cc.maxiter = 5000;
        cc.xspace = linspace(cc.xbounds(1), cc.xyres(1), cc.xbounds(2));
        cc.yspace = linspace(cc.ybounds(1), cc.xyres(2), cc.ybounds(2));
    end
    cc.colors = (1 / 255) .* [192, 192, 192; ...
        double(hsvconv([(0:1/1249:1)', ones(1250, 2)])); 64, 64, 64];
    cc.colbands = 100;
    cc.mincol = 1;
    cc.maxcol = cc.maxiter;

    % build grid and store
    [cc.xgrid, cc.ygrid] = ndgrid(cc.xspace, cc.yspace);
    ne_gcfg.h.FB.Config = cc;

    % update
    ne_fb_update;

    % set visible, modal and wait
    hFig.CloseRequestFcn = @ne_fb_closeui;
    hFig.WindowButtonDownFcn = @ne_fb_btdown;
    hFig.WindowButtonMotionFcn = @ne_fb_btmove;
    hFig.WindowButtonUpFcn = @ne_fb_btup;
    hFig.ResizeFcn = @ne_fb_resize;
    hFig.HandleVisibility = 'callback';
    hFig.Resize = 'on';
    hFig.Visible = 'on';
    drawnow;

% give warning
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error loading fractal browser.', 'NeuroElf GUI - error', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'fb_open')) = [];
end




% sub functions
function ne_fb_closeui(varargin)

% global storage
global ne_gcfg;

% delete figure and config
if isxfigure(ne_gcfg.h.FB.FBFig, true)
    ne_gcfg.h.FB.FBFig.Delete;
end
ne_gcfg.h.FB = [];
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'fb_open')) = [];



function ne_fb_resize(varargin)

% global storage
global ne_gcfg;

% adapt size of image
ne_gcfg.h.FB.h.Browser.Position(3:4) = ne_gcfg.h.FB.FBFig.Position(3:4);



function ne_fb_update(varargin)

% global storage
global ne_gcfg;
cc = ne_gcfg.h.FB.Config;

% on GPU
if isa(cc.xspace, 'gpuArray')
    count = gather(arrayfun(@ne_fb_mbcount, cc.xgrid, cc.ygrid, cc.maxiter)');

% on CPU
else
    count = ne_fb_mbcountarray(cc.xgrid, cc.ygrid, cc.maxiter)';
end

% no color bands
maxi = gather(cc.maxiter);
if cc.colbands == 1
    count = floor(max(1, size(cc.colors, 1) .* (log(count) ./ log(maxi))));
else
    scol = size(cc.colors, 1) - 2;
    mcount = (count >= gather(cc.maxiter));
    count = 1 + mod(round((cc.colbands * scol / gather(cc.maxiter)) .* count), scol);
    count(mcount) = scol + 2;
end

% look up colors
szc = size(count);
count = reshape(cc.colors(count(:), :), [szc, 3]);

% set image
ne_gcfg.h.FB.h.Browser.CData = image_resize(count, 0.5);



% mouse
function ne_fb_btdown(varargin)

% global storage
global ne_gcfg;
cc = ne_gcfg.h.FB.Config;
ch = ne_gcfg.h.FB.FBFig;
fSize = ch.Position(3:4);
mPos = ch.CurrentPoint;
dPos = min(fSize, max([0, 0], mPos)) ./ fSize;
dPos = [cc.xbounds(1), cc.ybounds(1)] + ...
    dPos .* [(cc.xbounds(2) - cc.xbounds(1)), (cc.ybounds(2) - cc.ybounds(1))];
ne_gcfg.h.FB.Config.down = struct('dPos', dPos, 'mPos', mPos, 'time', now, 'move', false);

function ne_fb_btmove(varargin)

function ne_fb_btup(varargin)

% global storage
global ne_gcfg;
cc = ne_gcfg.h.FB.Config;
if ~cc.down.move
    if (now - ne_gcfg.h.FB.Config.down.time) < 1e-5
        bsize = cc.xbounds(2) - cc.xbounds(1);
        bmid = cc.down.dPos;
        cc.xbounds = bmid(1) + bsize .* [-0.25, 0.25];
        cc.ybounds = bmid(2) + bsize .* [-0.25, 0.25];
        if isa(cc.xspace, 'gpuArray')
            cc.xspace = gpuArray.linspace(cc.xbounds(1), cc.xbounds(2), cc.xyres(1));
            cc.yspace = gpuArray.linspace(cc.ybounds(1), cc.ybounds(2), cc.xyres(2));
        else
            cc.xspace = linspace(cc.xbounds(1), cc.xbounds(2), cc.xyres(1));
            cc.yspace = linspace(cc.ybounds(1), cc.ybounds(2), cc.xyres(2));
        end
        [cc.xgrid, cc.ygrid] = ndgrid(cc.xspace, cc.yspace);
        ne_gcfg.h.FB.Config = cc;
        ne_fb_update;
    end
end


function count = ne_fb_mbcount(xr, xi, mi)
z0 = complex(xr, xi);
z = z0;
count = 0;
zs = real(z) * real(z) + imag(z) * imag(z);
while count <= mi && ...
    zs <= 4
    count = count + 1;
    zs = real(z) * real(z) + imag(z) * imag(z);
    z = z * z + z0;
end
count = count + 0.5 / (zs - 3.5);

function count = ne_fb_mbcountarray(xr, xi, mi)
z0 = complex(xr, xi);
z = z0;
count = zeros(size(z));
docount = 1:numel(z);
while ~isempty(docount) && ...
    count(docount(1)) <= mi
    count(docount) = count(docount) + 1;
    zd = z(docount);
    z0d = z0(docount);
    z(docount) = zd .* zd + z0d;
    docount(real(zd) .* real(zd) + imag(zd) .* imag(zd) > 4) = [];
end

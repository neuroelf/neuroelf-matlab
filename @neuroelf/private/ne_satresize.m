function varargout = ne_satresize(varargin)
% ne_satresize  - resize satellite window
%
% FORMAT:       ne_satresize(SRC, EVT, window, winsize)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       window      window specifier (used to target satellite windows)
%       winsize     width and height for satellite window
%
% No output fields.
%
% Example:
%
%       ne_satresize(0, 0, 'BS123456', [800, 600]);

% Version:  v1.1
% Build:    16052716
% Date:     May-27 2016, 4:34 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% if already resizing
if ne_gcfg.c.satresize
    return;
end

% size given
if nargin > 3 && isa(varargin{4}, 'double') && numel(varargin{4}) == 2 && ...
   ~any(isinf(varargin{4}) | isnan(varargin{4}) | varargin{4} < 32 | varargin{4} > 2880)
    reqsize = round(varargin{4}(:)');
    cfn = fieldnames(ne_gcfg.cc);
    if isempty(varargin{3}) && numel(cfn) == 1
        varargin{3} = cfn{1};
    end
else
    reqsize = [];
end

% with error handling
try

    % only for valid satellites
    if nargin < 3 || ~ischar(varargin{3}) || ~isfield(ne_gcfg.cc, varargin{3})
        return;
    end
    ne_gcfg.c.satresize = true;
    iSat = varargin{3}(:)';
    ch = ne_gcfg.cc.(iSat);

    % wait until position stable
    if isempty(reqsize)
        cp = [-1, -1, -1, -1];
        while ~isequal(cp, ch.Satellite.Position)
            cp = ch.Satellite.Position;
            pause(0.333);
        end
    end

    % slice display
    if strcmpi(ch.Config.sattype, 'slice')

        % get size
        if isempty(reqsize)
            drawnow;
            cp = ch.Satellite.Position;
            sz = max(528, min(cp(3), cp(4)));
        else
            cp = ch.Satellite.Position;
            sz = max(528, min(reqsize(:)));
        end
        hsz = round(0.5 * sz);
        hsz8 = max(1, hsz([1, 1]) - 8);
        sz16 = max(1, sz([1, 1]) - 16);

        % set new position
        ch.Satellite.Position = [cp(1:2), sz, sz];

        % get satellite tags
        tSat = ch.Satellite.TagStruct;

        % update position of controls
        tSat.(sprintf('AX_%s_Slice_SAG', iSat)).Position = [6, hsz + 2, hsz8];
        tSat.(sprintf('AX_%s_Slice_COR', iSat)).Position = [hsz + 2, hsz + 2, hsz8];
        tSat.(sprintf('AX_%s_Slice_TRA', iSat)).Position = [hsz + 2, 6, hsz8];
        tSat.(sprintf('AX_%s_Slice_Zoom', iSat)).Position = [8, 8, sz16];
        tSat.(sprintf('IM_%s_Slice_SAG', iSat)).Position = [6, hsz + 2, hsz8];
        tSat.(sprintf('IM_%s_Slice_COR', iSat)).Position = [hsz + 2, hsz + 2, hsz8];
        tSat.(sprintf('IM_%s_Slice_TRA', iSat)).Position = [hsz + 2, 6, hsz8];
        tSat.(sprintf('IM_%s_Slice_Zoom', iSat)).Position = [8, 8, sz16];

        % update config
        ch.Config.slicepos = [ ...
            tSat.(sprintf('IM_%s_Slice_SAG', iSat)).Position; ...
            tSat.(sprintf('IM_%s_Slice_COR', iSat)).Position; ...
            tSat.(sprintf('IM_%s_Slice_TRA', iSat)).Position];
        ch.Config.slicepos(:, 3:4) = ...
            ch.Config.slicepos(:, 1:2) + ch.Config.slicepos(:, 3:4);
        ch.Config.zslicepos = ...
            tSat.(sprintf('IM_%s_Slice_Zoom', iSat)).Position;
        ch.Config.zslicepos(3:4) = ...
            ch.Config.zslicepos(1:2) + ch.Config.zslicepos(3:4);
        ch.Config.zslicepos = ch.Config.zslicepos([1, 1, 1], :);

        % update transimg?
        if ne_gcfg.c.ini.Satellites.ResizeTransImg

            % delete transimg handles
            imh = ch.Config.tio.imCor.Handle;
            delete(ch.Config.tio.imCor);
            set(imh, 'CData', uint8(zeros([hsz8, 3])));
            ch.Config.tio.imCor = transimg(hsz8(1), hsz8(2));
            sethandle(ch.Config.tio.imCor, imh);

            % also for other transimg's
            imh = ch.Config.tio.imSag.Handle;
            delete(ch.Config.tio.imSag);
            set(imh, 'CData', uint8(zeros([hsz8, 3])));
            ch.Config.tio.imSag = transimg(hsz8(1), hsz8(2));
            sethandle(ch.Config.tio.imSag, imh);
            imh = ch.Config.tio.imTra.Handle;
            delete(ch.Config.tio.imTra);
            set(imh, 'CData', uint8(zeros([hsz8, 3])));
            ch.Config.tio.imTra = transimg(hsz8(1), hsz8(2));
            sethandle(ch.Config.tio.imTra, imh);
            imh = ch.Config.tio.imSlZ.Handle;
            delete(ch.Config.tio.imSlZ);
            set(imh, 'CData', uint8(zeros([sz16, 3])));
            ch.Config.tio.imSlZ = transimg(sz16(1), sz16(2));
            sethandle(ch.Config.tio.imSlZ, imh);

            % update config
            ne_gcfg.cc.(iSat) = ch;
            ne_setsatslicepos(0, 0, iSat);
        else
            ne_gcfg.cc.(iSat) = ch;
        end

        % draw
        ch.Satellite.Position = [cp(1:2), sz, sz];
        if isempty(reqsize)
            pause(0.01);
            drawnow;
        end

    % surface display
    elseif strcmpi(ch.Config.sattype, 'surf')

        % new limits
        xlim = [-128, 128];
        ylim = [-128, 128];

        % larger X than Y (or vice versa)
        cp = ch.Satellite.Position;
        if ~isempty(reqsize)
            cp(3:4) = reqsize;
            ch.Satellite.Position = cp;
        end
        if cp(3) > cp(4)
            xlim = xlim .* (cp(3) / cp(4));
        else
            ylim = ylim .* (cp(4) / cp(3));
        end

        % set new axes properties
        set(ch.Surface, 'Ylim', xlim, 'ZLim', ylim);
        if isempty(reqsize)
            drawnow;
        end

        % and allow for new hit-test
        ne_gcfg.cc.(varargin{3}).Config.surfpos(3:4) = cp(3:4);

    % render display
    elseif strcmpi(ch.Config.sattype, 'render')

        % get current figure size
        if isempty(reqsize)
            drawnow;
            cp = ch.Satellite.Position;
        else
            cp = ch.Satellite.Position;
            cp(3:4) = reqsize;
        end

        % get size
        sz = max([256, 256], cp(3:4));

        % set new position
        ch.Satellite.Position = [cp(1:2), sz];

        % get satellite tags
        tSat = ch.Satellite.TagStruct;

        % update position of controls
        tSat.(sprintf('AX_%s_Slice_Rend', iSat)).Position = [0, 0, sz];
        tSat.(sprintf('IM_%s_Slice_Rend', iSat)).Position = [0, 0, sz];

        % delete and recreate transimg
        imh = ch.Config.tio.imRnd.Handle;
        delete(ch.Config.tio.imRnd);
        set(imh, 'CData', uint8(zeros([sz, 3])));
        ch.Config.tio.imRnd = transimg(sz(1), sz(2));
        sethandle(ch.Config.tio.imRnd, imh);

        % update config
        ch.Config.surfpos(3:4) = sz;
        ne_gcfg.cc.(iSat) = ch;
        ne_render_setview(0, 0, iSat);

    % tcplot display
    elseif strcmpi(ch.Config.sattype, 'tcplot')

        % get size
        if isempty(reqsize)
            drawnow;
            cp = ch.Satellite.Position;
            sz = max([180, 36], cp(1, 3:4));
        else
            sz = max([180, 36], reqsize(:)');
        end
        ch.Satellite.Position(3:4) = sz;

        % re-set time course size
        ch.TCPlot.FontSize = 12;
        if prod(sz) < 50000
            ch.TCPlot.Position = [2, 2, sz - 2];
        elseif prod(sz) < 120000
            ch.TCPlot.Position = [32, 20, sz(1) - 40, sz(2) - 24];
        else
            ch.TCPlot.FontSize = 16;
            ch.TCPlot.Position = [44, 32, sz(1) - 56, sz(2) - 44];
        end
    end

% error?
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
ne_gcfg.c.satresize = false;

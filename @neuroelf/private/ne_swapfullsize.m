% FUNCTION ne_swapfullsize: switch between smaller and full size modes
function varargout = ne_swapfullsize(varargin)

% Version:  v1.1
% Build:    16052513
% Date:     May-25 2016, 1:09 PM EST
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

% only once
if ne_gcfg.c.resize
    return;
end
ne_gcfg.c.resize = true;

% shorthands
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
mf = ch.MainFig;
tags = mf.TagStruct;
uici = cc.fullsizex(:, 2) > 0;
uics = cc.fullsizexh(uici, :);
uicd = cc.fullsizex(uici, :);
mfp = mf.Position;

% error handling
try

    % initial config
    ch.Maximize.Visible = 'off';
    ch.Minimize.Visible = 'on';
    ch.Menu.CloseFile.Visible = 'off';
    ch.Menu.SelectFile.Visible = 'off';
    ch.Menu.Stats.Visible = 'off';
    ch.Menu.VOI.Visible = 'off';
    ch.Progress.Progress(408, 'newwidth');

    % for resize call
    if nargin < 3 || ~ischar(varargin{3}) || ~strcmpi(varargin{3}(:)', 'swap')

        % give it a little bit (natural delay to avoid resizing too often)
        pause(0.1);

        % get current size
        if nargin > 2 && isa(varargin{3}, 'double') && numel(varargin{3}) == 2 && ...
           ~any(isinf(varargin{3}) | isnan(varargin{3}) | varargin{3} ~= round(varargin{3}) | varargin{3} < cc.fullsize)
            cpos = varargin{3};
            mf.Position(3:4) = cpos;
            drawnow;
        end
        lcpos = [0, 0];
        cpos = mf.Position(3:4);
        if nargin < 3 || ~isa(varargin{3}, 'double') || numel(varargin{3}) ~= 2 || ...
            any(isinf(varargin{3}) | isnan(varargin{3}) | varargin{3} ~= round(varargin{3}) | varargin{3} < cc.fullsize)
            while any(cpos ~= lcpos)
                lcpos = cpos;
                cpos = mf.Position(3:4);
                pause(0.16);
            end
        end
        cpos = max(cc.fullsize, min(cpos, [2400, 1344]));
        aw = max(0, cpos(1) - cc.fullsize(1));
        ah = max(0, cpos(2) - cc.fullsize(2));
        if ah > aw
            cpos = cc.fullsize + [aw, aw];
            ah = aw;
        end
        mf.Position(3:4) = cpos;
        aw = aw - ah;
        hs = round(ah / 2);
        qs = round(ah / 4);

        % compute shift vectors
        cshift = [ ...
             aw, ah,  0,  0; ... % left-side buttons (shift up and possibly to the right)
             aw,  0,  0, ah; ... % dividing frame (possibly shift right and size up)
        ah + aw, ah,  0,  0; ... % right-side buttons (shift right, up)
             aw,  0, ah, ah; ...
        ah + aw,  0,  0,  0; ...
             aw,  0, ah,  0; ...
             aw, hs, hs, hs; ...
        hs + aw, hs, hs, hs; ...
        hs + aw,  0, hs, hs; ...
        qs + aw, qs,  0,  0; ... % central controls, shift a quarter up + a quarter right
              0, ah, aw,  0; ...
              0,  0, aw, ah; ...
             aw,  0,  0,  0];    % surface stats controls

        % update positions
        for uicc = 1:size(uics, 1)
            set(uics(uicc, 1), 'Position', ...
                uicd(uicc, 3:6) + cshift(uicd(uicc, 2), :));
        end

        % and finally patch progress bar
        if aw > 0
            ch.Progress.Progress(408 + aw, 'newwidth');
        end

    % swapping if currently fullsized
    elseif cc.fullsized && ...
        all(mf.Position(3:4) == cc.fullsize)

        % hide first
        mf.Visible = 'off';
        drawnow;

        % then set positions to small size
        for uicc = 1:size(uics, 1)
            set(uics(uicc, 1), 'Position', uicd(uicc, 7:10));
        end
    else

        % hide first
        mf.Visible = 'off';
        drawnow;
        cc.fullsized = false;

        % then set position to large size
        for uicc = 1:size(uics, 1)
            set(uics(uicc, 1), 'Position', uicd(uicc, 3:6));
        end
    end

    % record current positions (for hit tests)
    ne_gcfg.fcfg.histpos = tags.AX_NeuroElf_Slice_Hist.Position;
    ne_gcfg.fcfg.histpos(3:4) = ...
        ne_gcfg.fcfg.histpos(1:2) + ne_gcfg.fcfg.histpos(3:4);
    ne_gcfg.fcfg.slicepos = [ ...
        tags.IM_NeuroElf_Slice_SAG.Position; ...
        tags.IM_NeuroElf_Slice_COR.Position; ...
        tags.IM_NeuroElf_Slice_TRA.Position];
    ne_gcfg.fcfg.slicepos(:, 3:4) = ...
        ne_gcfg.fcfg.slicepos(:, 1:2) + ne_gcfg.fcfg.slicepos(:, 3:4);
    ne_gcfg.fcfg.surfpos = tags.AX_NeuroElf_Surface.Position;
    ne_gcfg.fcfg.surfpos(:, 3:4) = ...
        ne_gcfg.fcfg.surfpos(:, 1:2) + ne_gcfg.fcfg.surfpos(:, 3:4);
    ne_gcfg.fcfg.tcpos = tags.AX_NeuroElf_TC_Plot.Position;
    ne_gcfg.fcfg.tcpos(3:4) = ...
        ne_gcfg.fcfg.tcpos(1:2) + ne_gcfg.fcfg.tcpos(3:4);
    ne_gcfg.fcfg.zslicepos = tags.IM_NeuroElf_Slice_Zoom.Position;
    ne_gcfg.fcfg.zslicepos(3:4) = ...
        ne_gcfg.fcfg.zslicepos(1:2) + ne_gcfg.fcfg.zslicepos(3:4);
    ne_gcfg.fcfg.zslicepos = ne_gcfg.fcfg.zslicepos([1, 1, 1], :);

    % initial config
    ne_gcfg.fcfg.fullsized = true;
    set(ch.MainFigMLH, 'Resize', 'on');

    % nothing else to do
    if nargin > 2 && ...
        ischar(varargin{3}) && ...
        strcmpi(varargin{3}(:)', 'swap')
        if cc.fullsized
            ne_gcfg.fcfg.fullsized = false;
            mf.Position(3:4) = cc.fullsize - double(cc.fullsized) .* cc.fullsizes;
            set(ch.MainFigMLH, 'Resize', 'off');
            ch.Maximize.Visible = 'on';
            ch.Minimize.Visible = 'off';
            ch.Progress.Progress(530, 'newwidth');
            ch.Listener.Position(1:2) = [550, 12];
            ch.SampledValues.Position(2) = -40;
            ch.SrfStats.PThresh.Position(2) = -40;
            ch.SrfStats.ClusterTable.Position(2) = -40;
            ch.SrfStats.RGBLPos.Position(2) = -40;
            ch.SrfStats.RGBUPos.Position(2) = -40;
            ch.SrfStats.RGBLNeg.Position(2) = -40;
            ch.SrfStats.RGBUNeg.Position(2) = -40;
            ch.Menu.CloseFile.Visible = 'on';
            ch.Menu.SelectFile.Visible = 'on';
            ch.Menu.Stats.Visible = 'on';
            ch.Menu.VOI.Visible = 'on';
        else
            mf.Position(3:4) = cc.fullsize;
        end
    end

    % keep track
    ne_gcfg.c.ini.MainFig.FullSize = ne_gcfg.fcfg.fullsized;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    mf.Position = mfp;
end

% visualize
if nargin < 3 || ...
   ~ischar(varargin{3})
    mf.Position(3:4) = cpos;
else
    mf.Redraw;
    mf.Visible = 'on';
end
drawnow;

% update child positions
if ~isempty(cc.chpos{1})
    ne_gcfg.fcfg.chpos{3} = get(cc.chpos{1}, 'Position');
    ne_gcfg.fcfg.chpos{3} = cat(1, ne_gcfg.fcfg.chpos{3}{:});
    ne_gcfg.fcfg.chpos{3}(:, 3:4) = ne_gcfg.fcfg.chpos{3}(:, 1:2) + ne_gcfg.fcfg.chpos{3}(:, 3:4);
end

% release
ne_gcfg.c.resize = false;

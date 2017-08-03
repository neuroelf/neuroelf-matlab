% FUNCTION ne_render_updateui: update UI/struct based on callback
function ne_render_updateui(varargin)

% Version:  v0.9d
% Build:    14062416
% Date:     Jun-24 2014, 4:54 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% get tags and config
rf = ne_gcfg.h.Render;
tags = rf.h;
cc = ne_gcfg.fcfg.Render;
ci = ne_gcfg.c.ini;

% depending on what to do
switch (lower(varargin{3}))

    % update basic dropdown and checkboxes
    case {'check'}
        colblend = {'RGB', 'HSV'};
        interpm = {'linear', 'cubic', 'lanczos3'};
        ne_gcfg.fcfg.Render.colblend = colblend{tags.DD_render_colblend.Value};
        ne_gcfg.fcfg.Render.filewrt = (tags.CB_render_filewrt.Value > 0);
        ci.Render.FileWrite = ne_gcfg.fcfg.Render.filewrt;
        ne_gcfg.fcfg.Render.filewmov = (tags.CB_render_movfilewrt.Value > 0);
        ci.Render.FileWriteMovie = ne_gcfg.fcfg.Render.filewmov;
        ne_gcfg.fcfg.Render.imetha = interpm{tags.DD_render_avarinterp.Value};
        ci.Render.AnatInterp = ne_gcfg.fcfg.Render.imetha;
        ne_gcfg.fcfg.Render.imeth = interpm{tags.DD_render_stvarinterp.Value};
        ci.Render.StatsInterp = ne_gcfg.fcfg.Render.imeth;
        ne_gcfg.fcfg.Render.join = (tags.CB_render_join.Value > 0);
        ne_gcfg.fcfg.Render.showinfig = (tags.CB_render_showinfig.Value > 0);
        ci.Render.ShowWhileRendering = ne_gcfg.fcfg.Render.showinfig;
        ne_gcfg.fcfg.Render.smstatk = (tags.CB_render_smstatk.Value > 0);
        ci.Render.SmoothStatBorderKeep = ne_gcfg.fcfg.Render.smstatk;
        ne_gcfg.fcfg.Render.stmulaalp = (tags.CB_render_stmulaalp.Value > 0);
        ci.Render.StatsMultWithAnaAlpha = ne_gcfg.fcfg.Render.stmulaalp;

    % udpate dropdown
    case {'slavar'}

        % set slavar to value of UserData
        slavidx = tags.DD_render_slavarfile.Value;
        ne_gcfg.fcfg.Render.slavar = ...
            tags.DD_render_slavarfile.UserData{slavidx};

        % get slicing var
        slvar = ne_gcfg.fcfg.Render.slvar;

        % set enabled status of alpha edit field
        if slavidx > 1
            slavidx = 'on';
            slvar.RunTimeVars.AlphaVolume = ne_gcfg.fcfg.Render.slavar.RunTimeVars.xffID;
        else
            slavidx = 'off';
            slvar.RunTimeVars.AlphaVolume = '';
        end
        tags.ED_render_agalpha.Enable = slavidx;

        % update render
        ne_gcfg.fcfg.Render.hFig.Pointer = 'watch';
        drawnow;
        try
            ne_render_setview;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ne_gcfg.fcfg.Render.hFig.Pointer = 'arrow';
        drawnow;

    % activations min/max depth
    case {'actdepth'}
        try
            val = str2double(tags.ED_render_actdepth.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val < 0 || ...
                val > 30
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.actdepth;
        end
        tags.ED_render_actdepth.String = sprintf('%g', val);
        ne_gcfg.fcfg.Render.actdepth = val;
        ci.Render.ActivationDepth = val;

    % distance
    case {'dist'}

        % get text value (with error allowance)
        try
            val = str2double(tags.ED_render_dist.String);
            if numel(val) ~= 1 || ...
                isnan(val) || ...
                val < 160
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.dist;
        end

        % update in frame
        tags.ED_render_dist.String = sprintf('%d', val);
        ne_gcfg.fcfg.Render.dist = val;
        ci.Render.Distance = val;

    % resolution
    case {'res'}
        try
            val = str2double(tags.ED_render_res.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val < 16 || ...
                val > 4096
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.res;
        end
        tags.ED_render_res.String = sprintf('%d', val);
        ne_gcfg.fcfg.Render.res = val;
        ci.Render.Resolution = val;

    % filename
    case {'filename'}

        % simply copy
        ne_gcfg.fcfg.Render.filename = ddeblank(tags.ED_render_filename.String);
        ci.Render.Filename = ne_gcfg.fcfg.Render.filename;

    % movie filename
    case {'filenmov'}

        % simply copy
        ne_gcfg.fcfg.Render.filenmov = ddeblank(tags.ED_render_movfilename.String);
        ci.Render.FilenameMovie = ne_gcfg.fcfg.Render.filenmov;

    % alpha values
    case {'agalpha'}
        try
            val = str2double(tags.ED_render_agalpha.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val <= 0 || ...
                val > 1
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.agalpha;
        end
        tags.ED_render_agalpha.String = sprintf('%g', val);
        ne_gcfg.fcfg.Render.agalpha = val;
    case {'galpha'}
        try
            val = str2double(tags.ED_render_galpha.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val <= 0 || ...
                val > 1
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.galpha;
        end
        tags.ED_render_galpha.String = sprintf('%g', val);
        ne_gcfg.fcfg.Render.galpha = val;
    case {'stalp'}
        try
            val = str2double(tags.ED_render_stalp.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val > 1
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.stalp;
        end
        tags.ED_render_stalp.String = sprintf('%g', val);
        ne_gcfg.fcfg.Render.stalp = val;
    case {'stminaalp'}
        try
            val = str2double(tags.ED_render_stminaalp.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val < 0 || ...
                val > 1
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.stminaalp;
        end
        tags.ED_render_stminaalp.String = sprintf('%.2g', val);
        ne_gcfg.fcfg.Render.stminaalp = val;
    case {'stsmooth'}
        try
            val = str2double(tags.ED_render_stsmooth.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val < 0 || ...
                val > 10
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.stsmooth;
        end
        tags.ED_render_stsmooth.String = sprintf('%g', val);
        ne_gcfg.fcfg.Render.stsmooth = val;
    case {'smstat'}
        try
            val = str2double(tags.ED_render_smstat.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val < 0 || ...
                val > 16
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.smstat;
        end
        tags.ED_render_smstat.String = sprintf('%.3g', val);
        ne_gcfg.fcfg.Render.smstat = val;
        ci.Render.SmoothStatBorder = val;

    % update xfrom value
    case {'xfrom'}

        % get text value (with error allowance)
        try
            val = str2double(tags.ED_render_xfrom.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val)
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.slfrom;
        end

        % update in frame
        tags.ED_render_xfrom.String = sprintf('%d', val);
        ne_gcfg.fcfg.Render.slfrom = val;

    % update xstep value
    case {'xstep'}

        % get text value (with error allowance)
        try
            val = str2double(tags.ED_render_xstep.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val)
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.slstep;
        end

        % update in frame
        tags.ED_render_xstep.String = sprintf('%.1g', val);
        ne_gcfg.fcfg.Render.slstep = val;

    % update xto value
    case {'xto'}

        % get text value (with error allowance)
        try
            val = str2double(tags.ED_render_xto.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val)
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.slto;
        end

        % update in frame
        tags.ED_render_xto.String = sprintf('%d', val);
        ne_gcfg.fcfg.Render.slto = val;

    case {'yzfrom'}
        try
            val = str2double(tags.ED_render_yzfrom.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val < -255 || ...
                val > 255
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.frame(1);
        end
        if val > cc.frame(2)
            cc.frame = [cc.frame(2) val];
        else
            cc.frame(1) = val;
        end
        ne_gcfg.fcfg.Render.frame = cc.frame;
        tags.ED_render_yzfrom.String = sprintf('%d', cc.frame(1));
        tags.ED_render_yzto.String = sprintf('%d', cc.frame(2));
    case {'yzto'}
        try
            val = str2double(tags.ED_render_yzto.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val) || ...
                val < -255 || ...
                val > 255
                error('BAD_VAL');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.frame(2);
        end
        if val < cc.frame(1)
            cc.frame = [val cc.frame(2)];
        else
            cc.frame(2) = val;
        end
        ne_gcfg.fcfg.Render.frame = cc.frame;
        tags.ED_render_yzfrom.String = sprintf('%d', cc.frame(1));
        tags.ED_render_yzto.String = sprintf('%d', cc.frame(2));

    % rotations
    case {'yrot'}

        % get text value (with error allowance)
        try
            val = str2double(tags.ED_render_yzenith.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val)
                error('BAD_VAL');
            end
            val = mod(val, 360);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.roty;
        end

        % update in frame
        tags.ED_render_yzenith.String = sprintf('%d', val);
        ne_gcfg.fcfg.Render.roty = val;
    case {'zrot'}
        try
            val = str2double(tags.ED_render_zazimuth.String);
            if numel(val) ~= 1 || ...
                isinf(val) || ...
                isnan(val)
                error('BAD_VAL');
            end
            val = mod(val, 360);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = cc.rotz;
        end
        tags.ED_render_zazimuth.String = sprintf('%d', val);
        ne_gcfg.fcfg.Render.rotz = val;

    % update color buttons
    case {'setcolor'}
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
           ~any(strcmp(varargin{4}(:)', {'bgcol', 'hicol', 'locol'}))
            return;
        end
        clbt = varargin{4}(:)';

        % update color button
        clbtn = struct( ...
            'bgcol', 'Background color', ...
            'hicol', 'Full-intensity color', ...
            'locol', 'Null-intensity color');
        newcol = colorpicker(cc.(clbt), {clbtn.(clbt)});

        % update color
        if numel(newcol) == 3 && ...
           ~isequal(newcol, cc.(clbt))

            % set in config
            ne_gcfg.fcfg.Render.(clbt) = newcol;
            tags.(sprintf('BT_render_%s', clbt)).CData = ...
                repmat(reshape(uint8(newcol), [1, 1, 3]), 16, 22);
        end

        % update view
        if ne_gcfg.fcfg.page == 4
            ne_render_setview;
        end

    % dis-/enable filename
    case {'show'}
        tags.ED_vismontage_filename.Enable = 'off';
        vm.VMFig.RadioGroupSetOne('VisMOut', 1);
    case {'write'}
        tags.ED_vismontage_filename.Enable = 'on';
        vm.VMFig.RadioGroupSetOne('VisMOut', 2);
end

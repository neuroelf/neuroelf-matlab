function ne_setoption(varargin)
% ne_setoption  - set option for NeuroElf GUI
%
% FORMAT:       ne_setoption(SRC, EVT, optname [, optval])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       optname     option (name) to be altered (names see in examples)
%       optval      option value (if not given, either flip or request)
%
% No output fields.
%
% Examples:
%
%       ne_setoption(0, 0, 'ctableadd', true); % add table to existing VOIs
%       ne_setoption(0, 0, 'ctablelupcrd', 'cog'); % set to center-of-gravity
%       ne_setoption(0, 0, 'ctablescsizes', true); % show local maxima
%       ne_setoption(0, 0, 'ctablesort', 'x'); % sort cluster table by X

% Version:  v1.1
% Build:    16061500
% Date:     Jun-15 2016, 12:19 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011 - 2016, Jochen Weber
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

% disallow invalid calls
if nargin < 3 || ~ischar(varargin{3}) || isempty(varargin{3})
    return;
end
opt = lower(varargin{3}(:)');
if nargin > 3
    if ischar(varargin{4})
        optval = varargin{4}(:)';
    else
        optval = varargin{4};
    end
    if nargin > 4
        if ischar(varargin{5})
            optval2 = varargin{5}(:)';
        else
            optval2 = varargin{5};
        end
    else
        optval2 = [];
    end
else
    optval = [];
    optval2 = [];
end

% global variable(s)
global ne_gcfg;
ch = ne_gcfg.h;
ci = ne_gcfg.c.ini;
cpage = ne_gcfg.fcfg.page;
ts = ch.MainFig.TagStruct;

% option
switch (opt)

    % cluster-table "add"
    case 'ctableadd'

        % test optval
        if ~islogical(optval) || numel(optval) ~= 1
            return;
        end

        % make setting
        if optval
            strnew = 'off';
            stradd = 'on';
        else
            strnew = 'on';
            stradd = 'off';
        end
        ci.Statistics.ClusterTableAdd = optval;
        ts.UIM_NeuroElf_CTableNewVOI.Checked = strnew;
        ts.UIM_NeuroElf_CTableAddToVOI.Checked = stradd;

    % cluster-table lookup coordinate
    case 'ctablelupcrd'

        % test optval
        if ~ischar(optval) || ~any(strcmp(optval, {'center', 'cog', 'peak'}))
            return;
        end

        % make setting
        strcenter = 'off';
        strcog = 'off';
        strpeak = 'off';
        switch (optval)
            case 'peak'
                strpeak = 'on';
            case 'center'
                strcenter = 'on';
            case 'cog'
                strcog = 'on';
        end
        ci.Statistics.LookupCoord = optval;
        ts.UIM_NeuroElf_CTablePeak.Checked = strpeak;
        ts.UIM_NeuroElf_CTableCenter.Checked = strcenter;
        ts.UIM_NeuroElf_CTableCOG.Checked = strcog;

    % flip cluster-table sub-cluster sizes print flag
    case 'ctablescsizes'

        % swap flag
        if ~isempty(optval) && islogical(optval) && numel(optval) == 1
            ne_gcfg.fcfg.localmaxsz = optval;
        else
            ne_gcfg.fcfg.localmaxsz = ~ne_gcfg.fcfg.localmaxsz;
        end
        ci.Statistics.LocalMaxSizes = ne_gcfg.fcfg.localmaxsz;
        if ne_gcfg.fcfg.localmaxsz
            ts.UIM_NeuroElf_CTableSCSizes.Checked = 'on';
        else
            ts.UIM_NeuroElf_CTableSCSizes.Checked = 'off';
        end

    % set table sort option
    case 'ctablesort'

        % bad input
        if isempty(optval) || ~ischar(optval) || ...
           ~any(strcmp(optval, {'maxstat', 'maxstats', 'size', 'x', 'y', 'z'}))
            return;
        end

        % make setting
        ne_gcfg.fcfg.clsort = optval;
        ci.Statistics.Sorting = optval;
        ts.UIM_NeuroElf_CTableSortMax.Checked = 'off';
        ts.UIM_NeuroElf_CTableSortMaxS.Checked = 'off';
        ts.UIM_NeuroElf_CTableSortSize.Checked = 'off';
        ts.UIM_NeuroElf_CTableSortX.Checked = 'off';
        ts.UIM_NeuroElf_CTableSortY.Checked = 'off';
        ts.UIM_NeuroElf_CTableSortZ.Checked = 'off';
        switch (optval)
            case 'maxstat'
                ts.UIM_NeuroElf_CTableSortMax.Checked = 'on';
            case 'maxstats'
                ts.UIM_NeuroElf_CTableSortMaxS.Checked = 'on';
            case 'size'
                ts.UIM_NeuroElf_CTableSortSize.Checked = 'on';
            case 'x'
                ts.UIM_NeuroElf_CTableSortX.Checked = 'on';
            case 'y'
                ts.UIM_NeuroElf_CTableSortY.Checked = 'on';
            case 'z'
                ts.UIM_NeuroElf_CTableSortZ.Checked = 'on';
        end

    % drawing
    case 'drawon'

        % only for certain optvals
        if ~ischar(optval) || ...
            isempty(optval) || ...
           ~isfield(ci.Drawing, ['On' optval])
            return;
        end
        optval = ['On' optval];

        % field content
        if ~islogical(optval2) || ...
            numel(optval2) ~= 1
            optval2 = ~ci.Drawing.(optval);
        end

        % update field
        ci.Drawing.(optval) = optval2;

        % update menu item
        if optval2
            ts.(['UIM_NeuroElf_D' optval]).Checked = 'on';
        else
            ts.(['UIM_NeuroElf_D' optval]).Checked = 'off';
        end

    % echo calls
    case {'echocalls'}

        % no input
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = ~ci.MainFig.EchoCalls;
        end

        % set option
        ci.MainFig.EchoCalls = optval;
        ne_gcfg.c.echo = optval;
        if optval
            lstate = 'on';
        else
            lstate = 'off';
        end
        ts.UIM_NeuroElf_EchoCalls.Checked = lstate;

    % extended map names
    case {'extmapnames'}

        % old value
        oldval = ci.Statistics.ExtendedMapNames;

        % no input
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = ~ci.Statistics.ExtendedMapNames;
        end

        % set option
        ci.Statistics.ExtendedMapNames = optval;
        ne_gcfg.c.extmapnames = optval;
        if optval
            lstate = 'on';
        else
            lstate = 'off';
        end
        ts.UIM_NeuroElf_ExtMapNames.Checked = lstate;

        % update
        if ~isequal(optval, oldval) && ...
            isxff(ne_gcfg.fcfg.StatsVar, true)
            ne_gcfg.fcfg.StatsVar.Browse;
        end

    % beta extraction on selection
    case {'extonselect'}

        % test optval
        if ~ischar(optval) || ...
           ~any(strcmp(optval, {'manual', 'multi', 'single'}))
            return;
        end

        % make setting
        ci.Statistics.ExtractOnSelect = optval;
        ts.UIM_NeuroElf_CTableExtManual.Checked = 'off';
        ts.UIM_NeuroElf_CTableExtMulti.Checked = 'off';
        ts.UIM_NeuroElf_CTableExtSingle.Checked = 'off';
        switch (optval)
            case {'manual'}
                ts.UIM_NeuroElf_CTableExtManual.Checked = 'on';
            case {'multi'}
                ts.UIM_NeuroElf_CTableExtMulti.Checked = 'on';
            case {'single'}
                ts.UIM_NeuroElf_CTableExtSingle.Checked = 'on';
        end

    % subject-IDs with extraction
    case {'extsepchars'}

        % no input
        if ~isa(optval, 'double') || ...
            isempty(optval) || ...
            numel(optval) > 2 || ...
            any(isinf(optval) | isnan(optval) | optval < 9 | optval > 63)
            return;
        end
        optval = optval(:)';

        % set option
        ci.Statistics.ExtractSepChars = optval;
        spcstate = 'off';
        tabstate = 'off';
        comstate = 'off';
        smcstate = 'off';
        if isequal(optval, 32) || ...
            isequal(optval, [32, 32])
            spcstate = 'on';
        elseif isequal(optval, 9)
            tabstate = 'on';
        elseif isequal(optval, 44)
            comstate = 'on';
        elseif isequal(optval, 59)
            smcstate = 'on';
        end
        ts.UIM_NeuroElf_CTableExtSepSpc.Checked = spcstate;
        ts.UIM_NeuroElf_CTableExtSepTab.Checked = tabstate;
        ts.UIM_NeuroElf_CTableExtSepCom.Checked = comstate;
        ts.UIM_NeuroElf_CTableExtSepSmc.Checked = smcstate;

    % extract from transio-based GLMs
    case {'exttransio'}

        % no input
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = ~ci.Statistics.ExtractTransIO;
        end

        % set option
        ci.Statistics.ExtractTransIO = optval;
        if optval
            lstate = 'on';
        else
            lstate = 'off';
        end
        ts.UIM_NeuroElf_CTableExtTransIO.Checked = lstate;
        
    % subject-IDs with extraction
    case {'extwithsids'}

        % no input
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = ~ci.Statistics.ExtractWithSubIDs;
        end

        % set option
        ci.Statistics.ExtractWithSubIDs = optval;
        if optval
            lstate = 'on';
        else
            lstate = 'off';
        end
        ts.UIM_NeuroElf_CTableExtWithSID.Checked = lstate;

    % GLM conditions to extract
    case {'glmxconds'}

        % make sure object is a GLM
        if nargin < 4 || ...
            numel(optval) ~= 1 || ...
           ~isxff(optval, 'glm')
            optval = ne_gcfg.fcfg.StatsVar;
        end
        if ~isxff(optval, 'glm')
            return;
        end
        glmpreds = optval.SubjectPredictors;

        % get list to present
        if nargin < 5 || ...
           ~isa(optval2, 'double') || ...
            any(isinf(optval2(:)) | isnan(optval2(:)) | optval2(:) < 1)
            optrtv = optval.RunTimeVars;
            if ~isfield(optrtv, 'ExtractConds') || ...
               ~isa(optrtv.ExtractConds, 'double')
                optval2 = 1:numel(glmpreds);
            else
                optval2 = optrtv.ExtractConds(:)';
            end
        end
        optval2 = unique(min(numel(glmpreds), max(1, round(optval2(:)))))';

        % present list
        [optval2, selok] = listdlg( ...
            'ListString',    glmpreds, ...
            'SelectionMode', 'multiple', ...
            'InitialValue',  optval2, ...
            'Name',          'NeuroElf - GLM beta extraction selection', ...
            'PromptString',  'Please select the conditions/predictors to extract from the GLM:');
        if ~isequal(selok, 1)
            return;
        end

        % set new list
        optval.RunTimeVars.ExtractConds = optval2(:);

    % gray-LUT
    case {'graylut'}

        % invalid input
        if ~ischar(optval) || ...
            isempty(optval) || ...
           ~any(strcmpi(optval(:)', {'bw', 'bwb', 'bwbwb', 'color', 'wb', 'wbw'}))
            return;
        end

        % set value
        optval = lower(optval);
        switch (optval)

            % default
            case {'bw'}
                ne_gcfg.fcfg.graylut = [];

            % one band black-white-black
            case {'bwb'}
                ne_gcfg.fcfg.graylut = flexinterpn([0, 0, 0; 255, 255, 255; 0, 0, 0], ...
                    [Inf, Inf; 1, 1; 2/255, 1; 3, 3]);

            % one band black-white-black
            case {'bwbwb'}
                ne_gcfg.fcfg.graylut = flexinterpn([0; 255; 0; 255; 0] * ones(1, 3), ...
                    [Inf, Inf; 1, 1; 4/255, 1; 5, 3]);

            % color (with black as 0)
            case {'color'}
                ne_gcfg.fcfg.graylut = [0, 0, 0; ...
                    double(hsvconv([(0:1/277:11/12)', ones(254, 2)], 1)); ...
                    255, 255, 255];

            % inverse
            case {'wb'}
                ne_gcfg.fcfg.graylut = flexinterpn([255, 255, 255; 0, 0, 0], ...
                    [Inf, Inf; 1, 1; 1/255, 1; 2, 3]);

            % one band black-white-black
            case {'wbw'}
                ne_gcfg.fcfg.graylut = flexinterpn([255; 0; 255] * ones(1, 3), ...
                    [Inf, Inf; 1, 1; 2/255, 1; 3, 3]);
        end

        % update
        if cpage < 3
            ne_setslicepos;
        elseif cpage == 4
            ne_render_setview;
        end
        drawnow;

    % instant seed correlations
    case {'instseedcorr'}

        % reverse setting (without arguments)
        if nargin < 4
            optval = ~ci.Statistics.InstantSeedCorr;
        end

        % test optval
        if (~islogical(optval) && ~isa(optval, 'double')) || numel(optval) ~= 1
            return;
        end
        if ~islogical(optval)
            optval = (optval ~= 0);
        end

        % make setting
        ci.Statistics.InstantSeedCorr = optval;
        ne_gcfg.fcfg.instscorr = optval;
        if optval
            chflag = 'on';
            if isxff(ne_gcfg.fcfg.SliceVar, {'hdr', 'head', 'vtc'}) && ...
                ne_gcfg.fcfg.SliceVar.NrOfVolumes >= 20
                ne_gcfg.fcfg.instscvar = ne_gcfg.fcfg.SliceVar;
            end
            ne_setslicepos;
            ci.Statistics.InstantSeedCorr = optval;
            ne_gcfg.fcfg.instscorr = optval;
        else
            chflag = 'off';
        end
        ts.UIM_NeuroElf_InstantSeed.Checked = chflag;

    % join max-dist for 2 maps
    case {'joinmd2'}

        % reverse setting (without arguments)
        if nargin < 4
            optval = ~ci.Statistics.JoinMapsMaxDist;
        end

        % test optval
        if (~islogical(optval) && ...
            ~isa(optval, 'double')) || ...
            numel(optval) ~= 1
            return;
        end
        if ~islogical(optval)
            optval = (optval ~= 0);
        end

        % make setting
        ci.Statistics.JoinMapsMaxDist = optval;
        ne_gcfg.fcfg.joinmd2 = optval;
        if optval
            chflag = 'on';
        else
            chflag = 'off';
        end
        ts.UIM_NeuroElf_StatsJoinMaxDist.Checked = chflag;

    % underlay join
    case {'joinulay'}
        if numel(optval) == 1 && ...
            isa(optval, 'double') && ...
           ~isinf(optval) && ...
           ~isnan(optval) && ...
            optval >= 0 && ...
            optval <= 6
            ne_gcfg.fcfg.joinulay = optval;

            % update
            ne_setcvar;

            % menu checked
            ts.UIM_NeuroElf_ULayBlendOLay.Checked = 'off';
            ts.UIM_NeuroElf_ULayBlendOLayF.Checked = 'off';
            ts.UIM_NeuroElf_ULayBlendOLayW.Checked = 'off';
            ts.UIM_NeuroElf_ULayBlendMix.Checked = 'off';
            ts.UIM_NeuroElf_ULayBlendULayW.Checked = 'off';
            ts.UIM_NeuroElf_ULayBlendULayF.Checked = 'off';
            ts.UIM_NeuroElf_ULayBlendULay.Checked = 'off';
            switch (optval)
                case {0}
                    ts.UIM_NeuroElf_ULayBlendULay.Checked = 'on';
                case {1}
                    ts.UIM_NeuroElf_ULayBlendULayF.Checked = 'on';
                case {2}
                    ts.UIM_NeuroElf_ULayBlendULayW.Checked = 'on';
                case {3}
                    ts.UIM_NeuroElf_ULayBlendMix.Checked = 'on';
                case {4}
                    ts.UIM_NeuroElf_ULayBlendOLayW.Checked = 'on';
                case {5}
                    ts.UIM_NeuroElf_ULayBlendOLayF.Checked = 'on';
                case {6}
                    ts.UIM_NeuroElf_ULayBlendOLay.Checked = 'on';
            end
        end

    % MKDA option
    case {'mkda'}

        % test optval
        if ~ischar(optval) || ...
           ~any(strcmp(optval, {'lookuponcluster', 'lookuponcursor'}))
            return;
        end
        optval = sprintf('LookupOn%s%s', upper(optval(9)), optval(10:end));
        if nargin < 5 || ...
           ~islogical(optval2) || ...
            numel(optval(2)) ~= 1
            optval2 = ~ci.MKDA.(optval);
        end

        % make setting
        ci.MKDA.(optval) = optval2;
        if optval2
            chkd = 'on';
        else
            chkd = 'off';
        end
        switch lower(optval)
            case {'oncluster'}
                ts.UIM_NeuroElf_PLPLUCluster.Checked = chkd;
            case {'oncursor'}
                ts.UIM_NeuroElf_PLPLUCursor.Checked = chkd;
        end

    % new VMR resolution
    case {'newvmrres'}

        % test optval
        if ~isa(optval, 'double') || ...
            numel(optval) ~= 1 || ...
            isinf(optval) || ...
            isnan(optval) || ...
            optval < 0.1 || ...
            optval > 1
            if ischar(optval) && ...
                strcmpi(optval, 'res05')
                optval = 0.5 + 0.5 .* strcmpi(ts.UIM_NeuroElf_VMRRes05.Checked, 'on');
            else
                optval = inputdlg({'Create new VMRs with a resolution of:'}, ...
                    'NeuroElf - user input', 1, {'  1'});
                if isempty(optval) || ...
                   ~iscell(optval) || ...
                    numel(optval) ~= 1 || ...
                   ~ischar(optval{1}) || ...
                    isempty(ddeblank(optval{1}))
                    return;
                end
                try
                    optval = str2double(optval{1});
                    if numel(optval) ~= 1 || ...
                        isinf(optval) || ...
                        isnan(optval) || ...
                        optval < 0.1 || ...
                        optval > 1
                        return;
                    end
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                    return;
                end
            end
        end

        % set optval
        ne_gcfg.c.ini.MainFig.NewVMRRes = optval;
        if optval == 0.5
            m = 'on';
        else
            m = 'off';
        end
        ts.UIM_NeuroElf_VMRRes05.Checked = m;
        
    % orientation
    case {'orientation'}

        % test optval
        if ~ischar(optval) || ...
            isempty(optval) || ...
           ~any(lower(optval(1)) == 'nr')
            return;
        end
        optval = lower(optval);

        % for the MainFig
        if ~ischar(optval2) || ...
           ~isfield(ne_gcfg.cc, optval2)

            % set the orientation
            ci.MainFig.Orientation = upper(optval);
            ne_gcfg.fcfg.orient = optval;

            % update flag in menu
            switch (optval)
                case {'n'}
                    ts.UIM_NeuroElf_OrientNeuro.Checked = 'on';
                    ts.UIM_NeuroElf_OrientRadio.Checked = 'off';
                case {'r'}
                    ts.UIM_NeuroElf_OrientNeuro.Checked = 'off';
                    ts.UIM_NeuroElf_OrientRadio.Checked = 'on';
            end

            % update slice pos
            if nargin < 6 || ...
               ~islogical(varargin{6}) || ...
                numel(varargin{6}) ~= 1 || ...
                varargin{6}
                ne_setslicepos;
            end

        % for satellite figure
        else

            % set the orientation
            ne_gcfg.cc.(optval2).Config.orient = optval;

            % update slice pos
            if nargin < 6 || ...
               ~islogical(varargin{6}) || ...
                numel(varargin{6}) ~= 1 || ...
                varargin{6}
                ne_setsatslicepos(0, 0, optval2);
            end
        end

    % renderer
    case {'renderer'}

        % test optval
        if ~ischar(optval) || ...
           ~any(strcmpi(optval, {'OpenGL', 'zbuffer'}))
            return;
        end

        % set the renderer
        ci.MainFig.Renderer = optval;
        ne_gcfg.fcfg.renderer = optval;
        ch.MainFig.Renderer = optval;

        % update flag in menu
        switch (optval)
            case {'OpenGL'}
                ts.UIM_NeuroElf_RendOpenGL.Checked = 'on';
                ts.UIM_NeuroElf_Rendzbuffer.Checked = 'off';
            case {'zbuffer'}
                ts.UIM_NeuroElf_RendOpenGL.Checked = 'off';
                ts.UIM_NeuroElf_Rendzbuffer.Checked = 'on';
        end

        % draw again
        drawnow;

    % set scaling window
    case {'scalingwindow'}

        % only valid for good SliceVar
        svar = ne_gcfg.fcfg.SliceVar;
        if numel(svar) ~= 1 || ...
           ~isxff(svar, {'dmr', 'fmr', 'hdr', 'head', 'vmr', 'vtc'})
            return;
        end

        % get current scaling window
        rtv = svar.RunTimeVars;
        if isxff(svar, 'vmr') && ...
            rtv.ShowV16
            slim = rtv.ScalingWindowLim16;
            swin = rtv.ScalingWindow16;
        else
            slim = rtv.ScalingWindowLim;
            swin = rtv.ScalingWindow;
        end

        % use inputdlg
        rval = inputdlg( ...
            {'Selectable minimum value:', 'Selectable maximum value:', ...
             'Selected minimum value:', 'Selected maximum value:'}, ...
            'NeuroElf - user input', 1, ...
            {sprintf('%5g', slim(1)), sprintf('%5g', slim(2)), ...
             sprintf('%5g', swin(1)), sprintf('%5g', swin(2))});
        if ~iscell(rval) || ...
            numel(rval) ~= 4
            return;
        end

        % make sure input is good
        try
            slim = [str2double(rval{1}), str2double(rval{2})];
            if any(isinf(slim) | isnan(slim))
                error('GOING_AUTOMODE');
            end
            if slim(2) <= slim(1)
                slim(2) = slim(1) + 1;
            end
            swin = limitrangec([str2double(rval{3}), str2double(rval{4})], ...
                slim(1), slim(2), slim(1));
            if swin(2) <= swin(1)
                swin(2) = slim(2);
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            svar.SetScalingWindow([-Inf, Inf], true);
            ne_setcvar;
            return;
        end

        % update values
        if isxff(svar, 'vmr') && rtv.ShowV16
            svar.RunTimeVars.ScalingWindowLim16 = slim;
        else
            svar.RunTimeVars.ScalingWindowLim = slim;
        end
        svar.SetScalingWindow(slim, true);
        svar.RunTimeVars.ScalingWindow = swin;
        ne_setcvar;

    % show thresh bars
    case 'showthreshbars'

        % get optval from menu item of necessary
        if ~islogical(optval) || numel(optval) ~= 1
            optval = ~strcmpi(ts.UIM_NeuroElf_ShowThreshBars.Checked, 'on');
        end

        % make setting
        if optval
            show = 'on';
        else
            show = 'off';
        end
        ci.Statistics.ShowThreshBars = optval;
        ts.UIM_NeuroElf_ShowThreshBars.Checked = show;

        % update screen for pages 1 and 2 (unless suppressed)
        set(ch.CorStatsText, 'Visible', 'off');
        set(ch.ZoomStatsText, 'Visible', 'off');
        if cpage == 1 && ci.Statistics.ShowThreshText
            set(ch.CorStatsText, 'Visible', show);
        elseif cpage == 2 && ci.Statistics.ShowThreshText
            set(ch.ZoomStatsText, 'Visible', show);
        end
        if cpage < 3 && isempty(optval2)
            ne_setslicepos;
            drawnow;
        elseif cpage == 3 && isempty(optval2)
            ne_setsurfpos(0, 0, 1);
            drawnow;
        end

        % update other screens as well
        ccs = fieldnames(ne_gcfg.cc);
        for ccname = ccs(:)'
            satcfg = ne_gcfg.cc.(ccname{1});
            if isstruct(satcfg) && isfield(satcfg, 'Config') && ...
                isstruct(satcfg.Config) && isfield(satcfg.Config, 'sattype')
                if strcmpi(satcfg.Config.sattype, 'slice')
                    ne_showsatpage(0, 0, ccname{1}, satcfg.Config.page);
                elseif strcmpi(satcfg.Config.sattype, 'surf')
                    ne_setcsrfstatbars(0, 0, ccname{1});
                end
            end
        end

    % show thresh bar stats limits
    case 'showthreshtext'

        % get optval from menu item of necessary
        if ~islogical(optval) || numel(optval) ~= 1
            optval = ~strcmpi(ts.UIM_NeuroElf_ShowThreshText.Checked, 'on');
        end

        % make setting
        if optval
            show = 'on';
        else
            show = 'off';
        end
        ci.Statistics.ShowThreshText = optval;
        ts.UIM_NeuroElf_ShowThreshText.Checked = show;

        % update screen for pages 1 and 2 (unless suppressed)
        set(ch.CorStatsText, 'Visible', 'off');
        set(ch.ZoomStatsText, 'Visible', 'off');
        if cpage == 1 && ci.Statistics.ShowThreshText
            set(ch.CorStatsText, 'Visible', show);
        elseif cpage == 2 && ci.Statistics.ShowThreshText
            set(ch.ZoomStatsText, 'Visible', show);
        end
        if cpage < 3 && isempty(optval2)
            ne_setslicepos;
            drawnow;
        elseif cpage == 3 && isempty(optval2)
            ne_setsurfpos(0, 0, 1);
            drawnow;
        end

        % update other screens as well
        ccs = fieldnames(ne_gcfg.cc);
        for ccname = ccs(:)'
            satcfg = ne_gcfg.cc.(ccname{1});
            if isstruct(satcfg) && isfield(satcfg, 'Config') && ...
                isstruct(satcfg.Config) && isfield(satcfg.Config, 'sattype')
                if strcmpi(satcfg.Config.sattype, 'slice')
                    ne_showsatpage(0, 0, ccname{1}, satcfg.Config.page);
                elseif strcmpi(satcfg.Config.sattype, 'surf')
                    ne_setcsrfstatbars(0, 0, ccname{1});
                end
            end
        end

    % show V16 data on slicevar
    case {'showv16'}

        % only valid for VMR with correct data
        svar = ne_gcfg.fcfg.SliceVar;
        if ~isxff(svar, 'vmr') || ...
           ~isequal(size(svar.VMRData), size(svar.VMRData16))
            return;
        end

        % get value
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = (ts.BT_NeuroElf_showv16.Value ~= 0);
        end

        % set option value
        svar.RunTimeVars.ShowV16 = optval;
        ne_setcvar;
        if cpage ~= 3
            drawnow;
        end

    % set SPM normalization file
    case {'spmsn'}

        % only valid for SliceVar and StatsVar
        if ~ischar(optval) || ...
           ~any(strcmp(optval, {'SliceVar', 'StatsVar'}))
            return;
        end

        % not a valid xff
        if ~isxff(ne_gcfg.fcfg.(optval), true)
            return;
        end

        % request file
        if isempty(optval2)
            [spmsnf, spmsnp] = uigetfile( ...
                {'*.mat', 'SPM spatial normalization file (*_sn.mat)'}, ...
                'Please select an SPM spatial normalization file...');
            if isequal(spmsnf, 0) || ...
                isequal(spmsnp, 0) || ...
                isempty(spmsnf)
                return;
            end

            % try to load file
            try
                optval2 = load([spmsnp '/' spmsnf]);
                if ~isfield(optval2, 'VG') || ...
                   ~isfield(optval2, 'VF') || ...
                   ~isfield(optval2, 'Tr') || ...
                   ~isfield(optval2, 'Affine')
                    warning( ...
                        'neuroelf:BadFileContent', ...
                        'Invalid SPM *_sn.mat file loaded.' ...
                    );
                    return;
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                return;
            end
        end

        % set content
        ne_gcfg.fcfg.(optval).RunTimeVars.SPMsn = optval2;

        % update
        if cpage < 3
            ne_setslicepos;
            drawnow;
        end

    % surface background color
    case {'surfbgcol'}

        % request color
        ocol = ci.Surface.BackgroundColor;
        if ~isa(optval, 'double') || ...
           ~isequal(size(optval), [1, 3]) || ...
            any(isinf(optval) | isnan(optval) | optval < 0 | optval > 255)
            ncol = colorpicker(ocol, {'Surface view background color'});
        else
            ncol = round(optval);
        end

        % update
        if ~isequal(ncol, ocol)
            srfbcl = (1 / 255) .* ncol;
            ci.Surface.BackgroundColor = ncol;
            ne_gcfg.fcfg.SurfBackColor = srfbcl;
            set(ch.Surface, 'Color', srfbcl);
            ne_setsurfpos(0, 0, 1);
        end

    % surface reconstruction 4 triangles per surface
    case {'surfreco4tps'}

        % option not given
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = (ci.Surface.RecoTriPerVoxelFace ~= 4);
        end

        % update
        if optval
            ci.Surface.RecoTriPerVoxelFace = 4;
            ts.UIM_NeuroElf_SOReco4TPS.Checked = 'on';
        else
            ci.Surface.RecoTriPerVoxelFace = 2;
            ts.UIM_NeuroElf_SOReco4TPS.Checked = 'off';
        end

    % surface reconstruction colors
    case {'surfrecocol'}

        % request colors
        ocol = ci.Surface.RecoColors;
        if ~isa(optval, 'double') || ...
           ~isequal(size(optval), [2, 3]) || ...
            any(isinf(optval(:)) | isnan(optval(:)) | optval(:) < 0 | optval(:) > 255)
            ncol = colorpicker(ocol, {'Surface convex color'; 'Surface concave color'});
        else
            ncol = round(optval);
        end

        % update
        if ~isequal(ocol, ncol)
            ci.Surface.RecoColors = ncol;
        end

    % surface reconstruct single surface option
    case {'surfrecoonesurf'}

        % option not given
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = ~ci.Surface.RecoOneSurfaceOnly;
        end

        % update
        ci.Surface.RecoOneSurfaceOnly = optval;
        if optval
            ts.UIM_NeuroElf_SORecoOneSurf.Checked = 'on';
        else
            ts.UIM_NeuroElf_SORecoOneSurf.Checked = 'off';
        end

    % surface reuse XSM mapping
    case {'surfreusexsm'}

        % option not given
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = ~ci.Surface.ReuseSphereMapping;
        end

        % update
        ci.Surface.ReuseSphereMapping = optval;
        if optval
            ts.UIM_NeuroElf_SOReuseXSMMapping.Checked = 'on';
        else
            ts.UIM_NeuroElf_SOReuseXSMMapping.Checked = 'off';
        end

    % use FDR thresholds on VMP p-val dropdown
    case {'vmpusefdr'}

        % test optval
        if ~ischar(optval) || ...
            isempty(optval) || ...
           ~any(strcmpi(optval, {'indpos', 'nonpar', 'raw'}))
            return;
        end
        optval = lower(optval(:)');

        % make setting
        strind = 'off';
        strnon = 'off';
        strraw = 'off';
        switch (optval)
            case {'indpos'}
                strind = 'on';
            case {'nonpar'}
                strnon = 'on';
            case {'raw'}
                strraw = 'on';
        end
        ci.Statistics.VMPUseFDR = optval;
        ts.UIM_NeuroElf_VMPOPTRawThresh.Checked = strraw;
        ts.UIM_NeuroElf_VMPOPTFDRIndPos.Checked = strind;
        ts.UIM_NeuroElf_VMPOPTFDRNonPar.Checked = strnon;
end

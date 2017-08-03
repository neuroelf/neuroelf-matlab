% FUNCTION ne_mkda_setoption: set option for MKDA in NeuroElf
function ne_mkda_setoption(varargin)

% Version:  v0.9c
% Build:    12031512
% Date:     Dec-02 2011, 1:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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

% global variable(s)
global ne_gcfg;

% disallow access without MKDA being loaded
if ~isstruct(ne_gcfg.fcfg.MKDA)
    return;
end

% shorthands
ci = ne_gcfg.c.ini;
ts = ne_gcfg.h.MKDA.MKDAFig.TagStruct;
mcfg = ci.MKDA;

% disallow invalid calls
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~isfield(mcfg, varargin{3}(:)');
    return;
end
opt = varargin{3}(:)';
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

% option
switch (opt)

    % border color
    case {'BorderColor'}

        % requires valid PLP object
        ch = ne_gcfg.h.MKDA.h;
        plps = ch.PLPs;
        plpud = plps.UserData;
        plpid = plps.Value;
        try
            plp = plpud{plpid, 3};
            if numel(plp) ~= 1 || ...
               ~isxff(plp, 'plp')
                return;
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % test optval
        if numel(optval) ~= 3 || ...
            any(isinf(optval(:)) | isnan(optval(:)) | optval(:) < 0 | optval(:) > 255)

            % get if necessary
            bcolor = plp.RunTimeVars.Config.BorderColor;
            if isempty(bcolor)
                bcolor = [0, 0, 0];
            end
            [optval, optcncl] = colorpicker(bcolor, {'PLP points label color'});
            if optcncl
                optval = [];
            end
        else
            optval = round(optval(:)');
        end

        % set option
        plp.RunTimeVars.Config.BorderColor = optval;

    % apply mask to results
    case {'ApplyMask'}

        % get optval
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = strcmp(ts.UIM_NeuroElf_MKDA_OPT_mskr.Checked, 'off');
        end

        % make setting
        if optval
            checked = 'on';
        else
            checked = 'off';
        end
        ts.UIM_NeuroElf_MKDA_OPT_mskr.Checked = checked;

    % contrast computation
    case {'ContrastComp'}

        % test optval
        if ~ischar(optval) || ...
           ~any(strcmp(optval, {'diff', 'excl', 'wexcl'}))
            return;
        end

        % make setting
        strdiff = 'off';
        strexcl = 'off';
        strwexcl = 'off';
        ts.UIM_NeuroElf_MKDA_OCNT_swex.Enable = 'off';
        switch (optval)
            case {'diff'}
                strdiff = 'on';
            case {'excl'}
                strexcl = 'on';
            case {'wexcl'}
                strwexcl = 'on';
            ts.UIM_NeuroElf_MKDA_OCNT_swex.Enable = 'on';
        end
        ts.UIM_NeuroElf_MKDA_OCNT_diff.Checked = strdiff;
        ts.UIM_NeuroElf_MKDA_OCNT_excl.Checked = strexcl;
        ts.UIM_NeuroElf_MKDA_OCNT_wex.Checked = strwexcl;

    % contrast computation exclusion weight
    case {'ContrastCompExclWeight'}

        % get optval if necessary
        if isempty(optval) || ...
            numel(optval) ~= 1 || ...
            isinf(optval) || ...
            isnan(optval) || ...
            optval <= 0 || ...
            optval >= 1
            optval = inputdlg({'Relative exclusion threshold:   '}, ...
                'NeuroElf - user input', 1, {' 0.5'});
            if ~iscell(optval) || ...
                numel(optval) ~= 1 || ...
               ~ischar(optval{1}) || ...
                isempty(ddeblank(optval{1}))
                return;
            end
            try
                optval = str2double(ddeblank(optval{1}));
                if ~isa(optval, 'double') || ...
                    numel(optval) ~= 1 || ...
                    isinf(optval) || ...
                    isnan(optval) || ...
                    optval <= 0 || ...
                    optval >= 1
                    return;
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                return;
            end
        end

        % update text
        ts.UIM_NeuroElf_MKDA_OCNT_swex.Label = sprintf( ...
            'Set relative exclusion weight... (T = %g)', optval);

    % group map computation
    case {'GroupMapComp'}

        % bad input
        if isempty(optval) || ...
           ~ischar(optval) || ...
           ~any(strcmp(optval, {'ost', 'sum', 'wsum'}))
            return;
        end

        % make setting
        strsum = 'off';
        strwsum = 'off';
        switch (optval)
            case {'sum'}
                strsum = 'on';
            case {'wsum'}
                strwsum = 'on';
        end
        ts.UIM_NeuroElf_MKDA_OGRP_sum.Checked = strsum;
        ts.UIM_NeuroElf_MKDA_OGRP_wsum.Checked = strwsum;

    % join blobs computation
    case {'JoinBlobComp'}

        % bad input
        if isempty(optval) || ...
           ~ischar(optval) || ...
           ~any(strcmp(optval, {'max', 'rsum'}))
            return;
        end

        % make setting
        strmax = 'off';
        strrsum = 'off';
        switch (optval)
            case {'max'}
                strmax = 'on';
            case {'rsum'}
                strrsum = 'on';
        end
        ts.UIM_NeuroElf_MKDA_OJBM_max.Checked = strmax;
        ts.UIM_NeuroElf_MKDA_OJBM_rsum.Checked = strrsum;

    % keep individual maps
    case {'KeepIndivMaps'}

        % get optval
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = strcmp(ts.UIM_NeuroElf_MKDA_OPT_imap.Checked, 'off');
        end

        % make setting
        if optval
            checked = 'on';
        else
            checked = 'off';
        end
        ts.UIM_NeuroElf_MKDA_OPT_imap.Checked = checked;

    % label color
    case {'LabelColor'}

        % requires valid PLP object
        ch = ne_gcfg.h.MKDA.h;
        plps = ch.PLPs;
        plpud = plps.UserData;
        plpid = plps.Value;
        try
            plp = plpud{plpid, 3};
            if numel(plp) ~= 1 || ...
               ~isxff(plp, 'plp')
                return;
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % test optval
        if numel(optval) ~= 3 || ...
            any(isinf(optval(:)) | isnan(optval(:)) | optval(:) < 0 | optval(:) > 255)

            % get if necessary
            optval = colorpicker( ...
                plp.RunTimeVars.Config.LabelColor, {'PLP points label color'});
        else
            optval = round(optval(:)');
        end

        % set option
        plp.RunTimeVars.Config.LabelColor = optval;

    % label column name
    case {'LabelColumn'}

        % requires valid PLP object
        ch = ne_gcfg.h.MKDA.h;
        plps = ch.PLPs;
        plpud = plps.UserData;
        plpid = plps.Value;
        try
            plp = plpud{plpid, 3};
            if numel(plp) ~= 1 || ...
               ~isxff(plp, 'plp')
                return;
            end
            cnames = plp.ColumnNames(:);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % test optval
        if ~ischar(optval) || ...
           (~isempty(optval) && ...
            ~any(strcmpi(cnames, optval(:)')))

            % ask for value
            optval = listdlg( ...
                'ListString',    [{'<none>'}; cnames], ...
                'SelectionMode', 'single', ...
                'ListSize',      [min(640, max(320, 10 * size(char(cnames), 2))), 300], ...
                'InitialValue',  1, ...
                'Name',          'NeuroElf - user input', ...
                'PromptString',  'Please select the column for labeling points...');
            if numel(optval) ~= 1 || ...
               ~isa(optval, 'double') || ...
                optval < 1 || ...
                optval > (numel(cnames) + 1)
                return;
            end
            if optval == 1
                optval = '';
            else
                optval = cnames{optval - 1};
            end
        elseif isempty(optval)
            optval = '';
        end

        % set option
        if ~isempty(optval)
            optval = cnames{findfirst(strcmpi(cnames, optval(:)'))};
        end
        plp.RunTimeVars.Config.LabelColumn = optval;
        if isempty(optval)
            ch.LabelColumnMenu.Label = 'Label points by column (none)';
        else
            ch.LabelColumnMenu.Label = sprintf('Label points by column (%s)', optval);
        end

    % points-per-study weighting function
    case {'PPSWeighting'}

        % test optval
        if ~ischar(optval) || ...
           ~any(strcmp(optval, ...
                {'confidence', 'logpoints', 'none', 'points', 'sqrtpoints'}))
            return;
        end

        % make setting
        strconf = 'off';
        strlog  = 'off';
        strnone = 'off';
        strnpts  = 'off';
        strsqrt = 'off';
        switch (optval)
            case {'confidence'}
                strconf = 'on';
            case {'logpoints'}
                strlog = 'on';
            case {'none'}
                strnone = 'on';
            case {'points'}
                strnpts = 'on';
            case {'sqrtpoints'}
                strsqrt = 'on';
        end
        ts.UIM_NeuroElf_MKDA_OWPS_conf.Checked = strconf;
        ts.UIM_NeuroElf_MKDA_OWPS_log.Checked  = strlog;
        ts.UIM_NeuroElf_MKDA_OWPS_none.Checked = strnone;
        ts.UIM_NeuroElf_MKDA_OWPS_npts.Checked = strnpts;
        ts.UIM_NeuroElf_MKDA_OWPS_sqrt.Checked = strsqrt;

    % spatial null (full/near)
    case {'SpatialNull'}

        % bad input
        if isempty(optval) || ...
           ~ischar(optval) || ...
           ~any(strcmp(optval, {'full', 'near'}))
            return;
        end

        % make setting
        strfull = 'off';
        strnear = 'off';
        switch (optval)
            case {'full'}
                strfull = 'on';
            case {'near'}
                strnear = 'on';
        end
        ts.UIM_NeuroElf_MKDA_OSPN_full.Checked = strfull;
        ts.UIM_NeuroElf_MKDA_OSPN_near.Checked = strnear;

    % generate VMP with copies of summary maps
    case {'SummaryVMP'}

        % get optval
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = strcmp(ts.UIM_NeuroElf_MKDA_OPT_svmp.Checked, 'off');
        end

        % make setting
        if optval
            checked = 'on';
        else
            checked = 'off';
        end
        ts.UIM_NeuroElf_MKDA_OPT_svmp.Checked = checked;

    % unique unit points
    case {'UniqueUnitPoints'}

        % get optval
        if ~islogical(optval) || ...
            numel(optval) ~= 1
            optval = strcmp(ts.UIM_NeuroElf_MKDA_OPT_uniq.Checked, 'off');
        end

        % make setting
        if optval
            checked = 'on';
        else
            checked = 'off';
        end
        ts.UIM_NeuroElf_MKDA_OPT_uniq.Checked = checked;

end

% update option in ini
ci.MKDA.(opt) = optval;

% and PLP
if numel(ne_gcfg.fcfg.plp) == 1 && ...
    isxff(ne_gcfg.fcfg.plp, 'plp')
    ne_gcfg.fcfg.plp.RunTimeVars.Config = ci.MKDA;
end

% FUNCTION ne_rm_select: select a mediation model component
function ne_rm_select(varargin)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-10 2011, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
cf = ne_gcfg.h.RM.RMFig;
ch = ne_gcfg.h.RM.h;
ci = ne_gcfg.c.ini;

% disallow invalid calls
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3})
    return;
end
arg = lower(varargin{3}(:)');

% get GLM
glm = cf.UserData.lastglm;
if ~isxff(glm, 'glm')
    ne_rm_closeui;
    return;
end

% what to select
switch (arg)

    % type of X
    case {'xtype'}

        % set value of list to 1 (to avoid invalid setting)
        ch.XList.Value = 1;

        % depending on Value
        if ch.XType.Value == 1
            ch.XList.String = glm.RunTimeVars.Contrasts(:, 1);
            tval = 1;
        else
            ch.XList.String = glm.RunTimeVars.CovariatesNames(:);
            tval = 2;
        end

        % update ini file
        ci.Mediation.XType = tval;

        % call XList setting as well
        ne_rm_select(0, 0, 'xlist');

    % type of Y
    case {'ytype'}
        ch.YList.Value = 1;
        if ch.YType.Value == 1
            ch.YList.String = glm.RunTimeVars.Contrasts(:, 1);
            tval = 1;
        else
            ch.YList.String = glm.RunTimeVars.CovariatesNames(:);
            tval = 2;
        end
        ci.Mediation.YType = tval;
        ne_rm_select(0, 0, 'ylist');

    % entry within X list
    case {'xlist'}

        % get current xtype
        tval = ch.XType.Value;

        % selection
        sval = ch.XList.Value;

        % depending on value
        switch (tval)

            % contrasts
            case {1}

                % make sure the same is not set in M/C
                mval = ch.MCons.Value;
                if any(mval == sval)
                    mval(mval == sval) = [];
                    ch.MCons.Value = mval;
                end
                mval = ch.CCons.Value;
                if any(mval == sval)
                    mval(mval == sval) = [];
                    ch.CCons.Value = mval;
                end

                % ensure goodness of M list
                ne_rm_select(0, 0, 'mcons');

            % covariates
            case {2}
                mval = ch.MCovs.Value;
                if any(mval == sval)
                    mval(mval == sval) = [];
                    ch.MCovs.Value = mval;
                end
                mval = ch.CCovs.Value;
                if any(mval == sval)
                    mval(mval == sval) = [];
                    ch.CCovs.Value = mval;
                end
                ne_rm_select(0, 0, 'mcovs');
        end

        % give shout if Y is the same
        if tval == ch.YType.Value && ...
            ch.XList.Value == ch.YList.Value
            uiwait(warndlg('X and Y cannot be the same term!', ...
                'NeuroElf - warning', 'modal'));
        end

    % entry within Y list
    case {'ylist'}
        tval = ch.YType.Value;
        sval = ch.YList.Value;
        switch (tval)
            case {1}
                mval = ch.MCons.Value;
                if any(mval == sval)
                    mval(mval == sval) = [];
                    ch.MCons.Value = mval;
                end
                mval = ch.CCons.Value;
                if any(mval == sval)
                    mval(mval == sval) = [];
                    ch.CCons.Value = mval;
                end
                ne_rm_select(0, 0, 'mcons');
            case {2}
                mval = ch.MCovs.Value;
                if any(mval == sval)
                    mval(mval == sval) = [];
                    ch.MCovs.Value = mval;
                end
                mval = ch.CCovs.Value;
                if any(mval == sval)
                    mval(mval == sval) = [];
                    ch.CCovs.Value = mval;
                end
                ne_rm_select(0, 0, 'mcovs');
        end
        if tval == ch.XType.Value && ...
            ch.XList.Value == ch.YList.Value
            uiwait(warndlg('X and Y cannot be the same term!', ...
                'NeuroElf - warning', 'modal'));
        end

    % M contrast selection
    case {'mcons'}

        % adapted? edited?
        a = false;
        e = false;

        % ensure that the same is not selected in the X/Y list
        sval = ch.MCons.Value;
        if ch.XType.Value == 1 && ...
            any(sval == ch.XList.Value)
            sval(sval == ch.XList.Value) = [];
            ch.MCons.Value = sval;
            e = true;
        end
        if ch.YType.Value == 1 && ...
            any(sval == ch.YList.Value)
            sval(sval == ch.YList.Value) = [];
            ch.MCons.Value = sval;
            e = true;
        end
        if ~isempty(intersect(sval, ch.CCons.Value)) || ...
           (isempty(sval) && ...
            isempty(ch.MCovs.Value))
            a = true;
        end

        % warning
        if a
            if e
                uiwait(warndlg('The mediator term (M) was altered and needs attention', ...
                    'NeuroElf - warning', 'modal'));
            else
                uiwait(warndlg('The mediator term (M) needs (re-)selection', ...
                    'NeuroElf - warning', 'modal'));
            end
        end

    % M covariate selection
    case {'mcovs'}
        a = false;
        e = false;
        sval = ch.MCovs.Value;
        if ch.XType.Value == 2 && ...
            any(sval == ch.XList.Value)
            sval(sval == ch.XList.Value) = [];
            ch.MCovs.Value = sval;
            e = true;
        end
        if ch.YType.Value == 2 && ...
            any(sval == ch.YList.Value)
            sval(sval == ch.YList.Value) = [];
            ch.MCovs.Value = sval;
            e = true;
        end
        if ~isempty(intersect(sval, ch.CCovs.Value)) || ...
           (isempty(sval) && ...
            isempty(ch.MCons.Value))
            a = true;
        end
        if a
            if e
                uiwait(warndlg('The mediator term (M) was altered and needs attention', ...
                   'NeuroElf - warning', 'modal'));
            else
                uiwait(warndlg('The mediator term (M) needs (re-)selection', ...
                   'NeuroElf - warning', 'modal'));
            end
        end

    % C contrast selection
    case {'ccons'}
        a = false;
        e = false;
        sval = ch.CCons.Value;
        if ch.XType.Value == 1 && ...
            any(sval == ch.XList.Value)
            sval(sval == ch.XList.Value) = [];
            ch.CCons.Value = sval;
            e = true;
        end
        if ch.YType.Value == 1 && ...
            any(sval == ch.YList.Value)
            sval(sval == ch.YList.Value) = [];
            ch.CCons.Value = sval;
            e = true;
        end
        if ~isempty(intersect(sval, ch.MCons.Value)) || ...
           (isempty(sval) && ...
            isempty(ch.CCovs.Value))
            a = true;
        end
        if a
            if e
                uiwait(warndlg('The covariate term (C) was altered and might need attention', ...
                    'NeuroElf - warning', 'modal'));
            else
                uiwait(warndlg('The covariate term (C) needs alteration', ...
                    'NeuroElf - warning', 'modal'));
            end
        end

    % C covariate selection
    case {'ccovs'}
        a = false;
        e = false;
        sval = ch.CCovs.Value;
        if ch.XType.Value == 2 && ...
            any(sval == ch.XList.Value)
            sval(sval == ch.XList.Value) = [];
            ch.CCovs.Value = sval;
            e = true;
        end
        if ch.YType.Value == 2 && ...
            any(sval == ch.YList.Value)
            sval(sval == ch.YList.Value) = [];
            ch.CCovs.Value = sval;
            e = true;
        end
        if ~isempty(intersect(sval, ch.MCovs.Value)) || ...
           (isempty(sval) && ...
            isempty(ch.CCons.Value))
            a = true;
        end
        if a
            if e
                uiwait(warndlg('The covariate term (C) was altered and might need attention', ...
                    'NeuroElf - warning', 'modal'));
            else
                uiwait(warndlg('The covariate term (C) needs alteration', ...
                    'NeuroElf - warning', 'modal'));
            end
        end

    % number of bootstrap samples
    case {'numsmp'}

        % make sure the value is numeric and between 100 and 100000
        try
            numsmp = str2double(ch.BootNumSmp.String);
            if ~isa(numsmp, 'double') || ...
                numel(numsmp) ~= 1 || ...
                isinf(numsmp) || ...
                isnan(numsmp) || ...
                numsmp <= 0
                numsmp = 1000;
            else
                numsmp = fix(max(500, min(1000000, real(numsmp))));
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            numsmp = 1000;
        end
        ch.BootNumSmp.String = sprintf('%d', numsmp);
        ci.Mediation.BootNumSmp = numsmp;

    % bootstrapping
    case {'boot'}

        % set button state
        cf.RadioGroupSetOne('abmeth', 1);

        % update ini
        ci.Mediation.Strategy = 'boot';

        % set groups
        cf.SetGroupEnabled('robust', 'off');
        cf.SetGroupEnabled('numsmp', 'on');
        cf.SetGroupEnabled('booted', 'on');

    % MCMAM
    case {'mcmam'}
        cf.RadioGroupSetOne('abmeth', 2);
        ci.Mediation.Strategy = 'mcmam';
        cf.SetGroupEnabled('robust', 'on');
        cf.SetGroupEnabled('numsmp', 'on');
        cf.SetGroupEnabled('booted', 'off');

    % Sobel-Test
    case {'sobel'}
        cf.RadioGroupSetOne('abmeth', 3);
        ci.Mediation.Strategy = 'sobel';
        cf.SetGroupEnabled('robust', 'on');
        cf.SetGroupEnabled('numsmp', 'off');
        cf.SetGroupEnabled('booted', 'off');

    % robust checkbox
    case {'robust'}

        % update ini
        ci.Mediation.Robust = (ch.Robust.Value > 0);

    % z-trans checkbox
    case {'ztrans'}
        ci.Mediation.ZTrans = (ch.Ztrans.Value > 0);

    % bootstrap sample reuse checkbox
    case {'breuse'}
        ci.Mediation.BootReuse = (ch.BootReuse.Value > 0);

    % bootstrap SE computation: percentile radiobutton
    case {'perc'}
        cf.RadioGroupSetOne('bstrap', 1);
        ci.Mediation.BootSE = 'perc';

    % bootstrap SE computation: percentile radiobutton
    case {'bca'}
        cf.RadioGroupSetOne('bstrap', 2);
        ci.Mediation.BootSE = 'bca';

    % bootstrap SE computation: percentile radiobutton
    case {'var'}
        cf.RadioGroupSetOne('bstrap', 3);
        ci.Mediation.BootSE = 'var';
end

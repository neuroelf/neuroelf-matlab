% FUNCTION ne_setstatthrpval: set thresholds from p < dropdown
function varargout = ne_setstatthrpval(varargin)

% Version:  v1.1
% Build:    16052921
% Date:     May-29 2016, 9:57 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% get configuration and handles
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
cinist = ne_gcfg.c.ini.Statistics;

% and current stats var and index
stvar = cc.StatsVar;
stvix = cc.StatsVarIdx;

% get selected p < -value
if nargin < 3 || ~isa(varargin{3}, 'double') || numel(varargin{3}) ~= 1 || ...
    isinf(varargin{3}) || isnan(varargin{3}) || varargin{3} <= 0 || varargin{3} > 0.2
    pval = str2double(ch.Stats.PThresh.String{ch.Stats.PThresh.Value});
else
    pval = varargin{3};
end

% for two-tailed stats
if ch.Stats.NegTail.Value ~= 0 && ch.Stats.PosTail.Value ~= 0

    % certain stats need to apply factor 0.5
    ttld = 0.5;

% but for (up to) one tail, we take the face value
else
    ttld = 1;
end

% special VMP-based cases
if numel(stvar) == 1 && isxff(stvar, 'vmp') && numel(stvix) == 1 && stvix >= 1 && ...
    stvix <= numel(stvar.Map) && stvar.Map(stvix).EnableClusterCheck == 0
    stmap = stvar.Map(stvix);

    % for VMP and FDR (without clustersize check!)
    if ~isempty(stmap.FDRThresholds) && ...
        any(stmap.FDRThresholds(:, 1) == (pval * ttld)) && ...
        any(strcmp(cinist.VMPUseFDR, {'indpos', 'nonpar'}))

        % get lower threshold
        fdrcol = 2 + double(strcmp(cinist.VMPUseFDR, 'nonpar'));
        lthresh = stmap.FDRThresholds(stmap.FDRThresholds(:, 1) == (pval * ttld), fdrcol);

        % get estimate of upper threshold
        vals = sort(double(abs(stmap.VMPData(abs(stmap.VMPData) >= lthresh))));
        if numel(vals) < 2
            uthresh = lthresh + 0.5;
        else
            uthresh = vals(ceil(0.95 * numel(vals)));
        end

        % set thresholds
        stvar.Map(stvix).LowerThreshold = lthresh;
        stvar.Map(stvix).UpperThreshold = uthresh;

        % and the text controls
        ch.Stats.LThresh.String = sprintf('%.4f', lthresh);
        ch.Stats.UThresh.String = sprintf('%.4f', uthresh);

        % update display then return
        ne_setslicepos;
        return;

    % for SVC-correction AND p<0.05
    elseif pval == 0.05 && isfield(stmap, 'RunTimeVars') && ...
        isstruct(stmap.RunTimeVars) && numel(stmap.RunTimeVars) == 1 && ...
        isfield(stmap.RunTimeVars, 'SVCResels') && ...
        isa(stmap.RunTimeVars.SVCResels, 'double') && ...
        numel(stmap.RunTimeVars.SVCResels) == 1 && stmap.RunTimeVars.SVCResels >= 1 && ...
       ~isinf(stmap.RunTimeVars.SVCResels) && ~isnan(stmap.RunTimeVars.SVCResels)

        % compute threshold
        pval = pval / stmap.RunTimeVars.SVCResels;
    end
end

% compute upper threshold (with prange factor)
ttud = cc.prange * ttld;

% computation depends on type of map
lttxt = ' ';
uttxt = ' ';
sdtxt = sprintf('%%.%df', cinist.ShowThreshDecimals);
switch (cc.StatsVarPar{1})

    % for t-maps
    case 't'

        % use (custom) tinv
        lthr = -sdist('tinv', ttld * pval, cc.StatsVarPar{2});
        uthr = -sdist('tinv', ttud * pval, cc.StatsVarPar{2});
        lttxt = sprintf(['t_{[%d]}>=' sdtxt], cc.StatsVarPar{2}, lthr);
        uttxt = sprintf(['t_{[%d]}<' sdtxt], cc.StatsVarPar{2}, uthr);

    % for F-maps
    case 'F'

        % use (custom) finv function
        lthr = sdist('finv', 1 - pval, cc.StatsVarPar{2}, cc.StatsVarPar{3});
        uthr = sdist('finv', 1 - cc.prange * pval, cc.StatsVarPar{2}, cc.StatsVarPar{3});
        lttxt = sprintf(['F_{[%d,%d]}>=' sdtxt], cc.StatsVarPar{2}, cc.StatsVarPar{3}, lthr);
        uttxt = sprintf(['F_{[%d,%d]}<' sdtxt], cc.StatsVarPar{2}, cc.StatsVarPar{3}, uthr);

    % for r-maps
    case 'r'

        % pass p-value through (custom) tinv and then set r-threshold
        lthr = correlinvtstat(-sdist('tinv', ttld * pval, cc.StatsVarPar{2}), ...
            cc.StatsVarPar{2} + 2);
        uthr = correlinvtstat(-sdist('tinv', ttud * pval, cc.StatsVarPar{2}), ...
            cc.StatsVarPar{2} + 2);
        lttxt = sprintf(['r_{df=%d}>=' sdtxt], cc.StatsVarPar{2}, lthr);
        uttxt = sprintf(['r_{df=%d}<' sdtxt], cc.StatsVarPar{2}, uthr);

    % for MKDA maps
    case 'm'

        % try to locate correct value
        if numel(stvar) == 1 && isxff(stvar, 'vmp') && numel(stvix) == 1 && ...
            stvar.Map(stvix).Type == 9
            strtv = stvar.Map(stvix).RunTimeVars;
            if pval == 0.05 && isfield(strtv, 'FWEMaxDist') && ...
                stvar.Map(stvix).EnableClusterCheck == 0
                lthr = strtv.FWEMaxDist(round(1 + 0.95 * numel(strtv.FWEMaxDist)));
                uthr = double(max(stvar.Map(stvix).VMPData(:)));
            else
                if isfield(strtv, 'AlphaSimMKDA') && iscell(strtv.AlphaSimMKDA) && ...
                    any(abs(strtv.AlphaSimMKDA{2}(:, 1) - pval) < eps)
                    lthr = strtv.AlphaSimMKDA{2}(abs(strtv.AlphaSimMKDA{2}(:, 1) - pval) < eps, 2);
                    uthr = (2 + 2 * (pval ^ 0.2)) * lthr;
                else
                    return;
                end
            end
        else
            return;
        end

    % for AvgVTCs
    case 'a'

        % rescale
        lthr = -sdist('tinv', ttld * pval, cc.StatsVarPar{2});
        uthr = -sdist('tinv', 0.01 * ttld * pval, cc.StatsVarPar{2});
        stvar.RunTimeVars.ConditionThresholds(stvix, 2, 1) = lthr;
        stvar.RunTimeVars.ConditionThresholds(stvix, 2, 2) = uthr;

        % update and return
        ne_setslicepos;
        return;

    % for unknown stats
    otherwise

        % simply leave
        return;
end

% and the text controls
ch.Stats.LThresh.String = sprintf('%.4f', lthr);
ch.Stats.UThresh.String = sprintf('%.4f', uthr);
set([ch.CorStatsText(1), ch.ZoomStatsText(1)], 'String', lttxt);
set([ch.CorStatsText(2), ch.ZoomStatsText(2)], 'String', uttxt);

% if one VMP map is selected
if numel(stvar) == 1 && isxff(stvar, {'hdr', 'head', 'vmp'}) && numel(cc.StatsVarIdx) == 1

    % get current threshold (for comparison)
    stvmt = stvar.Map(stvix).Type;
    oldt = [stvar.Map(stvix).LowerThreshold, stvar.Map(stvix).UpperThreshold];

    % then set thresholds
    stvar.Map(stvix).LowerThreshold = lthr;
    stvar.Map(stvix).UpperThreshold = uthr;

    % and get new values for comparison
    newt = [stvar.Map(stvix).LowerThreshold, stvar.Map(stvix).UpperThreshold];

    % if different and clustering requested
    if any(oldt ~= newt) && ...
        stvar.Map(stvix).EnableClusterCheck > 0 && ...
        isxff(stvar, 'vmp')

        % also adapt cluster size
        ntails = 2 - double(stvar.Map(stvix).Type == 4);
        utails = min(ntails, max(1, stvar.Map(stvix).ShowPositiveNegativeFlag - 1));
        if stvmt == 9
            fld = 'AlphaSimMKDA';
        else
            fld = sprintf('AlphaSim%d%d', utails, ntails);
        end
        if isfield(stvar.Map, 'RunTimeVars') && ...
            isstruct(stvar.Map(stvix).RunTimeVars) && ...
            isfield(stvar.Map(stvix).RunTimeVars, fld) && ...
            iscell(stvar.Map(stvix).RunTimeVars.(fld)) && ...
           ~isempty(stvar.Map(stvix).RunTimeVars.(fld)) && ...
            any(abs(stvar.Map(stvix).RunTimeVars.(fld){2}(:, 1) - pval) < eps)
            if stvmt == 9
                stvar.Map(stvix).ClusterSize = stvar.Map(stvix).RunTimeVars.(fld){2}( ...
                    abs(stvar.Map(stvix).RunTimeVars.(fld){2}(:, 1) - pval) < eps, 3);
            else
                stvar.Map(stvix).ClusterSize = stvar.Map(stvix).RunTimeVars.(fld){2}( ...
                    abs(stvar.Map(stvix).RunTimeVars.(fld){2}(:, 1) - pval) < eps, 2);
            end
            set(ch.Stats.kThresh, 'String', ...
                sprintf('%d', stvar.Map(stvix).ClusterSize));
        end

        % re-cluster with @VMP::ClusterTable
        stvar.Map(stvix).VMPDataCT = [];
        stvar.ClusterTable(stvix, []);
    end
end

% and then update window
ne_setslicepos;

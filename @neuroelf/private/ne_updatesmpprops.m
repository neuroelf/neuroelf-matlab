% FUNCTION ne_updatesmpprops: update SMP map properties
function varargout = ne_updatesmpprops(varargin)

% Version:  v1.1
% Build:    16051916
% Date:     May-19 2016, 4:36 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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

% get object/index
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h.SrfStats;
stvar = cc.SurfStatsVar;
stvix = cc.SurfStatsVarIdx;
if ~isxff(stvar, {'fsmf', 'glm', 'mtc', 'smp'}) || numel(stvix) ~= 1
    return;
end
strtv = stvar.RunTimeVars;
try
    smp = stvar.Map(stvix);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% update from dropdown?
if nargin > 2 && ischar(varargin{3}) && strcmpi(varargin{3}(:)', 'pval')

    % get selected p < -value
    pval = str2double(ch.PThresh.String{ch.PThresh.Value});

    % for two-tailed stats
    if ch.NegTail.Value ~= 0 && ch.PosTail.Value ~= 0

        % certain stats need to apply factor 0.5
        ttld = 0.5;

    % but for (up to) one tail, we take the face value
    else
        ttld = 1;
    end

    % compute upper threshold (with prange factor)
    ttud = cc.prange * ttld;

    % computation depends on type of map
    switch (smp.Type)

        % for t-maps
        case 1

            % use (custom) tinv
            lthr = -sdist('tinv', ttld * pval, smp.DF1);
            uthr = -sdist('tinv', ttud * pval, smp.DF1);

        % for F-maps
        case 4

            % use (custom) finv function
            lthr = sdist('finv', 1 - pval, smp.DF1, smp.DF2);
            uthr = sdist('finv', 1 - cc.prange * pval, smp.DF1, smp.DF2);

        % for r-maps
        case 2

            % pass p-value through (custom) tinv and then set r-threshold
            lthr = correlinvtstat(-sdist('tinv', ttld * pval, smp.DF1), smp.DF1 + 2);
            uthr = correlinvtstat(-sdist('tinv', ttud * pval, smp.DF1), smp.DF1 + 2);

        % for average MTC maps
        case 30

            % pass p-value through sdist for norminv
            lthr = [];
            uthr = [];
            zthrs = -sdist('norminv', [1, 0.01] .* (ttld * pval), 0, 1);

            % update RunTimeVars
            if isfield(strtv, 'AvgMTC') && islogical(strtv.AvgMTC) && strtv.AvgMTC
                stvar.RunTimeVars.ConditionThresholds(stvix, 2, :) = zthrs;
            end

        % for unknown stats
        otherwise

            % simply leave it as is
            lthr = [];
            uthr = [];
    end

    % and the text controls
    if ~isempty(lthr) && ~isempty(uthr)
        ch.LThresh.String = sprintf('%.4f', lthr);
        ch.UThresh.String = sprintf('%.4f', uthr);
    end
end

% try to update items
try
    stvar.Map(stvix).LowerThreshold = limitrangec(str2double(ch.LThresh.String), -1000, 1000, 1);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
try
    stvar.Map(stvix).UpperThreshold = limitrangec(str2double(ch.UThresh.String), -1000, 1000, 1);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
if stvar.Map(stvix).UpperThreshold <= stvar.Map(stvix).LowerThreshold
    stvar.Map(stvix).UpperThreshold = stvar.Map(stvix).LowerThreshold + 1;
end
try
    stvar.Map(stvix).ClusterSize = ...
        limitrangec(str2double(strrep(lower(get(ch.kThresh, 'String')), ...
        'm', '')), 0, 100000, stvar.Map(stvix).ClusterSize);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
stvar.Map(stvix).ShowPositiveNegativeFlag = ...
    double(ch.PosTail.Value > 0) + 2 * double(ch.NegTail.Value > 0);
stvar.Map(stvix).EnableClusterCheck = double(ch.UsekThr.Value > 0);
ch.LThresh.String = sprintf('%.4f', stvar.Map(stvix).LowerThreshold);
ch.UThresh.String = sprintf('%.4f', stvar.Map(stvix).UpperThreshold);
set(ch.kThresh, 'String', sprintf('%.1fmm', stvar.Map(stvix).ClusterSize));
if stvar.Map(stvix).EnableClusterCheck > 0
    stvar.Map(stvix).SMPDataCT = [];
end

% update map coloring
ne_setcsrfstatmap;

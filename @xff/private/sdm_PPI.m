function xo = sdm_PPI(xo, ppitc, ppicond, ppiname)
% SDM::PPI  - add PPI regressors
%
% FORMAT:       [sdm = ] sdm.PPI(ppitc, ppicond [, ppiname]);
%
% Input fields:
%
%       ppitc       time-course from VOI/phys data (in TR resolution)
%       ppicond     list of condition(s) to interact ppitc with
%       ppiname     name of region from which ppitc was derived ('ppi')
%
% Output fields:
%
%       sdm         SDM object (design matrix)
%
% Using: ddeblank, splittocell, ztrans.

% Version:  v1.1
% Build:    16021209
% Date:     Feb-12 2016, 9:51 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library
global ne_methods;
ddeblank    = ne_methods.ddeblank;
splittocell = ne_methods.splittocell;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'sdm') || ...
   ~isnumeric(ppitc) || isempty(ppitc) || numel(ppitc) ~= size(ppitc, 1) || ...
    any(isinf(ppitc) | isnan(ppitc)) || ~iscell(ppicond) || isempty(ppicond)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if numel(ppitc) ~= size(bc.SDMMatrix, 1)
    error('neuroelf:xff:badArgument', 'Invalid ppitc.');
end
ppitc = ne_methods.ztrans(ppitc);
if nargin < 4 || ~ischar(ppiname) || isempty(ppiname)
    ppiname = 'ppi';
else
    ppiname = ppiname(:)';
end

% get a bit more of content
sdml = bc.PredictorColors;
sdmn = bc.PredictorNames;
sdm = bc.SDMMatrix;

% insertion position
fc = bc.FirstConfoundPredictor;
if fc > (size(sdm, 2) + 1)
    fc = size(sdm, 2) + 1;
elseif fc <= 1
    error('neuroelf:xff:badFileContent', 'Invalid SDM configuration.');
else
    fc = round(fc);
end

% regress variance out of PPI-tc that we don't want to look at
sdmx = sdm(:, sum(abs(sdm)) > 0);
ppitc = ppitc - sdmx * ((pinv(sdmx' * sdmx) * sdmx') * ppitc);

% initialize marker
mark = cell2struct(cell(numel(ppicond), 1, 2), {'pos', 'neg'}, 3);

% initialize vector with conditions to sum
sumcond = false(1, numel(sdmn));

% iterate over condition list
for cc = numel(mark):-1:1

    % set marker
    mark(cc).pos = false(1, numel(sdmn));
    mark(cc).neg = false(1, numel(sdmn));

    % get name
    cname = ppicond{cc};

    % contrast?
    if any(cname == '>')

        % split
        parts = splittocell(cname, '>');
        parts{1} = ddeblank(parts{1});
        parts{2} = ddeblank(parts{2});
    else
        parts = {cname};
    end

    % split parts even further
    if numel(parts{1}) > 2 && parts{1}(1) == '(' && parts{1}(end) == ')'
        parts{1} = parts{1}(2:end-1);
    end
    parts{1} = splittocell(parts{1}, '+');
    for cpc = 1:numel(parts{1})
        parts{1}{cpc} = ddeblank(parts{1}{cpc});
        condi = strcmpi(parts{1}{cpc}, sdmn);
        if ~any(condi)
            error('neuroelf:xff:badArgument', ...
                'PPI condition particle not found: %s.', parts{1}{cpc});
        end
        mark(cc).pos = mark(cc).pos | condi(:)';
        sumcond = sumcond | condi(:)';
    end
    if numel(parts) > 1
        if numel(parts{2}) > 2 && parts{2}(1) == '(' && parts{2}(end) == ')'
            parts{2} = parts{2}(2:end-1);
        end
        parts{2} = splittocell(parts{2}, '+');
        for cpc = 1:numel(parts{2})
            parts{2}{cpc} = ddeblank(parts{2}{cpc});
            condi = strcmpi(parts{2}{cpc}, sdmn);
            if ~any(condi)
                error('neuroelf:xff:badArgument', ...
                    'PPI condition particle not found: %s.', parts{2}{cpc});
            end
            mark(cc).neg = mark(cc).neg | condi(:)';
            sumcond = sumcond | condi(:)';
        end
    end
end

% patch SDM
sumcond = find(sumcond);
if numel(sumcond) <= numel(mark)
    error('neuroelf:xff:badArgument', 'PPI conditions/contrasts not independent.');
end
sumridx = 1:max(sumcond);
sumridx(sumcond) = 1:numel(sumcond);
sdmsum = sdm(:, sumcond);
sdmlsum = sdml(sumcond, :);
sdm(:, sumcond(1)) = sum(sdmsum, 2);
sdmn{sumcond(1)} = sprintf('%s + ', sdmn{sumcond});
sdmn{sumcond(1)}(end-2:end) = [];
sdml(sumcond(1), :) = round(mean(sdml(sumcond, :)));

% increase size
ppitc(numel(ppitc), 1 + numel(mark)) = 0;

% ppi names
ppin = cell(1, size(ppitc, 2));
ppin{1} = ppiname;
ppil = floor(255.99 .* rand(1, 3));
ppil(1 + numel(mark), 3) = 0;

% for each condition/contrast
for cc = 1:numel(mark)

    % get pos/neg columns in smaller matrix
    mpos = sumridx(mark(cc).pos);
    mneg = sumridx(mark(cc).neg);

    % build new regressor
    sdm(:, sumcond(cc+1)) = sum(sdmsum(:, mpos), 2) - sum(sdmsum(:, mneg), 2);

    % set new name
    sdmn{sumcond(cc+1)} = ppicond{cc};
    ppin{cc + 1} = sprintf('%s-x-%s', ppiname, ppicond{cc});

    % and color
    sdml(sumcond(cc+1), :) = round(min(255, max(0, ...
        mean(sdmlsum(mpos, :), 1) - mean(sdmlsum(mneg, :), 1))));

    % and interact
    ppitc(:, cc + 1) = ppitc(:, 1) .* sdm(:, sumcond(cc+1));
    ppil(cc + 1, :) = round(0.5 * (sdml(sumcond(cc+1), :) + ppil(1, :)));
end

% remove conditions if necessary
if (numel(mark) + 1) < numel(sumcond)
    sdm(:, sumcond(cc+2:end)) = [];
    sdmn(sumcond(cc+2:end)) = [];
    sdml(sumcond(cc+2:end), :) = [];
end

% add interacted timecourses
sdm = [sdm(:, 1:(fc-1)), ppitc, sdm(:, fc:end)];
sdmn = [sdmn(:, 1:(fc-1)), ppin, sdmn(:, fc:end)];
sdml = [sdml(1:(fc-1), :); ppil; sdml(fc:end, :)];
fc = fc + size(ppitc, 2);

% put into output
bc.PredictorColors = sdml;
bc.PredictorNames = sdmn;
bc.FirstConfoundPredictor = fc;
bc.SDMMatrix = sdm;
bc.RTCMatrix = sdm(:, 1:fc-1);
xo.C = bc;

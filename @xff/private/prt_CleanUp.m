function xo = prt_CleanUp(xo, varargin)
% PRT::CleanUp  - clean up protocols
%
% FORMAT:       [prt] = prt.CleanUp([fullclean]);
%
% Input fields:
%
%       fullclean   if given and true, perform full cleanup
%
% Output fields:
%
%       prt         altered PRT
%
% Using: gluetostring.

% Version:  v1.1
% Build:    16021018
% Date:     Feb-10 2016, 6:05 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin > 1 && (isa(varargin{1}, 'double') || islogical(varargin{1})) && ...
   ~isempty(varargin{1}) && varargin{1}(1)
    fullclean = true;
else
    fullclean = false;
end

% get settings
ncon = numel(bc.Cond);
if ~isempty(bc.ResolutionOfTime) && lower(bc.ResolutionOfTime(1)) == 'm'
    tresms = 1;
    tresvol = 0;
else
    tresms = 0;
    tresvol = 1;
end

% loop over number of conditions
lastend = 0;
for cc = ncon:-1:1

    % get condition, update last end
    c = bc.Cond(cc);
    c.NrOfOnOffsets = size(c.OnOffsets, 1);
    if size(c.Weights, 1) ~= size(c.OnOffsets, 1)
        c.Weights = zeros(size(c.OnOffsets, 1), 0);
    end

    % skip empty conditions
    if isempty(c.OnOffsets)
        if fullclean
            bc.Cond(cc) = [];
        end
        continue;
    end

    % sort onsets
    [oos, oosi] = sort(c.OnOffsets, 1);
    c.OnOffsets = oos;
    if ~isempty(c.Weights)
        c.Weights = c.Weights(oosi, :);
    end

    % update last index and set changed flag
    lastend = max(lastend, max(c.OnOffsets(:, 2)));

    % loop over onsets, starting at the end
    for oc = c.NrOfOnOffsets:-1:2

        % check intervals
        if (c.OnOffsets(oc - 1, 2) + 1) >= c.OnOffsets(oc, 1)

            % give warning ?
            if (c.OnOffsets(oc - 1, 2) + tresvol) > c.OnOffsets(oc, 1)
                warning('neuroelf:xff:protocolWarning', ...
                    'Overlapping stimulus durations in same condition.');
            end

            % concatenate epochs ?
            if all(c.Weights(oc, :) == c.Weights(oc - 1, :))
                c.OnOffsets(oc, 1) = c.OnOffsets(oc - 1, 1);
                c.OnOffsets(oc - 1, :) = [];
                c.Weights(oc - 1, :) = [];

            % otherwise, set beginning of second to end of first
            else
                c.OnOffsets(oc, 1) = c.OnOffsets(oc - 1, 2);
            end
        end
    end

    % put back
    c.NrOfOnOffsets = size(c.OnOffsets, 1);
    bc.Cond(cc) = c;
end

% update number of conditions
if fullclean
    ncon = numel(bc.Cond);
    bc.NrOfConditions = ncon;
end

% that's it ?
if ~fullclean || ncon < 1
    xo.C = bc;
    return;
end

% get old conditions
ocond = bc.Cond;
oc1 = ocond(1);

% otherwise empty protocol conditions
nbc = bc;
nbc.Cond(:) = [];

% get combined list and matrix of on/offsets
iact = uint8(0);
iact(1:ncon, 1:lastend) = iact(1);
for cc = 1:ncon
    c = ocond(cc);
    onum = c.NrOfOnOffsets;
    for oc = 1:onum
        iact(cc,(c.OnOffsets(oc,1)+tresms):(c.OnOffsets(oc,2))) = uint8(1);
    end
end
numc = sum(iact);

% resting condition ?
rest = (numc == 0);
if any(rest)

    % get rest on/offsets
    ron  = 1 + find(diff(rest) > 0);
    roff = 1 + find(diff(rest) < 0);
    if isempty(roff)
        roff = length(rest);
    end
    if isempty(ron) || ron(1) > roff(1)
         ron = [1, ron];
    end
    if roff(end) < ron(end)
        roff(end + 1) = length(rest);
    end
    numrest = length(ron);
    if numrest ~= length(roff)
        error('neuroelf:xff:internalError', 'Mishap when figuring out rest phases.');
    end
    ronoff = zeros(numrest, 2);
    for rc = 1:numrest
        ronoff(cc, :) = [ron, roff - tresvol];
    end

    % add rest to new file (first condition)
    fprintf(' -> adding ''Rest'' condition with %d onsets\n', numrest);
    oc1.ConditionName = {'Rest'};
    oc1.NrOfOnOffsets = size(ronoff, 1);
    oc1.OnOffsets = ronoff;
    oc1.Weights = zeros(size(ronoff, 1), 0);
    oc1.Color = [32, 32, 32];
    nbc.Cond(end + 1) = oc1;
else

    % add rest with no onsets
    disp(' -> adding ''Rest'' condition with 0 onsets');
    oc1.ConditionName = {'Rest'};
    oc1.NrOfOnOffsets = 0;
    oc1.OnOffsets = zeros(0, 2);
    oc1.Weights = [];
    oc1.Color = [32, 32, 32];
    nbc.Cond(end + 1) = oc1;
end
xo.C = nbc;

% find overlap between one .. ncon conditions
for cc = 1:max(numc)

    % scan conditions
    for scc = 1:ncon

        % find where numc == cc and condition(scc) is on
        ciact = diff([false, (numc == cc & (iact(scc, :) > 0)), false]);
        ciaon = find(ciact > 0) - tresms;
        ciaoff = find(ciact < 0) - 1;

        % check length
        if length(ciaon) ~= length(ciaoff)
            xo.C = bc;
            error('neuroelf:xff:internalError', 'Error finding non-interaction periods.');
        end

        % do nothing on emtpy array
        if isempty(ciaon)
            continue;
        end

        % no interaction
        if cc == 1

            % put into new file
            fprintf(' -> keeping condition ''%s'' with %d onset(s)\n', ...
                deblank(ocond(scc).ConditionName{1}), length(ciaon));
            prt_AddCond(xo, deblank(ocond(scc).ConditionName{1}), ...
                [ciaon(:), ciaoff(:)], ocond(scc).Color, ocond(scc).Weights);
            iact(cc, (numc == cc & (iact(scc, :) > 0))) = false;
            continue;
        end
        nbc = xo.C;

        % repeat until interactions of level cc with condition(scc) OK
        while ~isempty(ciaon)

            % get "on" conditions on first onset
            oncon = find(iact(:, ciaon(1) + tresms) > 0);
            if length(oncon) ~= cc
                xo.C = bc;
                error('neuroelf:xff:internalError', 'Error looking up interaction conditions.');
            end

            % find phases where sum == cc and sum(oncon) == cc !
            siact = diff([false, (numc == cc & sum(iact(oncon(:)', :)) == cc), false]);
            siaon = find(siact > 0) - tresms;
            siaoff = find(siact < 0) - 1;

            % check length
            if length(siaon) ~= length(siaoff) || isempty(siaon)
                xo.C = bc;
                error('neuroelf:xff:internalError', 'Error finding interaction periods.');
            end

            % building new name
            newname = cell(1, numel(oncon));
            newcol = zeros(1, 3);
            for sscc = oncon(:)'
                newname{sscc} = deblank(ocond(sscc).ConditionName{1});
                newcol = newcol + ocond(sscc).Color;
            end
            newname = sprintf('Interaction (%s)', ne_methods.gluetostring(newname, ' & '));
            newcol = round(newcol / cc);

            % put into new file
            fprintf(' -> building interaction ''%s'' with %d onset(s)\n', newname, length(siaon));
            oc1.Name = newname;
            oc1.NrOfOnOffsets = numel(siaon);
            oc1.OnOffsets = [siaon(:), siaoff(:)];
            oc1.Weights = zeros(numel(siaon), 0);
            oc1.Color = newcol;
            nbc.Cond(end + 1) = oc1;
            xo.C = nbc;

            % remove from array
            for soc = 1:length(siaon)
                iact(oncon(:)', siaon(soc)+tresms:siaoff(soc)) = false;
            end

            % do find again to get loop right
            numc = sum(iact);
            ciact = diff([false, (numc == cc & (iact(scc, :) > 0)), false]);
            ciaon = find(ciact > 0) - tresms;
            ciaoff = find(ciact < 0) - 1;

            % check length
            if length(ciaon) ~= length(ciaoff)
                xo.C = bc;
                error('neuroelf:xff:internalError', 'Error finding non-interaction periods.');
            end
        end
    end
end

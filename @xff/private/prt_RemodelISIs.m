function xo = prt_RemodelISIs(xo, remod)
% PRT::RemodelISIs  - re-model ISI lengths between off- and onsets
%
% FORMAT:       [prt =] prt.RemodelISIs(remod)
%
% Input fields:
%
%       remod       Rx2 or Rx3 cell array with remodelling parameters
%                   with content {condition_specifier, ISIvalue, relative}
%                   where the relative flag determines whether ISIvalue is
%                   set as the new value (false) or is added to the
%                   existing ISI (true, default)
%                   the condition_specifier can be a specific name or a
%                   regular expression pattern to match several conditions
%
% Output fields:
%
%       prt         altered PRT
%
% Examples:
%
%   prt.RemodelISIs({'.*_(now|later)^', -1000, true; 'resp', 3000, false});
%
%   this will relatively shift the onsets following stimuli in conditions
%   ending in "now" or "later" by 1000ms closer to the preceding stimulus
%   and will set the ISI following stimuli in the resp condition to 3000ms.

% Version:  v1.1
% Build:    16021017
% Date:     Feb-10 2016, 5:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2013, 2014, 2016, Jochen Weber
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt') || ...
   ~iscell(remod) || ~any(size(remod, 2) == [2, 3]) || ndims(remod) ~= 2
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if isempty(remod)
    return;
end
if size(remod, 2) == 2
    remod(:, 3) = {true};
end
bc = xo.C;

% only works for millisecond protocols
if ~strcmpi(bc.ResolutionOfTime, 'msec')
    error('neuroelf:xff:badArgument', ...
        'This method only works for millisecond resolution protocols.');
end

% get number of conditions
ncond = numel(bc.Cond);

% generate list of onsets in order of time
ons = cat(1, bc.Cond.OnOffsets);
ons(:, 3) = 0;
oc = 0;
for cc = 1:ncond
    noc = size(bc.Cond(cc).OnOffsets, 1);
    ons(oc+1:oc+noc, 3) = cc;
    oc = oc + noc;
end
cnames = cat(1, bc.Cond.ConditionName);

% sort by time
[sons, soni] = sort(ons(:, 1));
ons = ons(soni, :);
lons = size(ons, 1) - 1;

% for each condition
for cc = 1:ncond

    % determine if condition fits to one of the patterns
    mc = 0;
    for rc = 1:size(remod, 1)
        if (any(remod{rc, 1} == '*' | remod{rc, 1} == '+' | remod{rc, 1} == '?' | remod{rc, 1} == '.') && ...
            ~isempty(regexpi(cnames{cc}, remod{rc, 1}))) || strcmpi(cnames{cc}, remod{rc, 1})
            mc = rc;
            break;
        end
    end

    % no match found
    if mc == 0
        continue;
    end

    % relative
    if remod{rc, 3}

        % iterate over onsets
        for oc = 1:lons

            % if condition matches
            if ons(oc, 3) == cc

                % shift subsequent event times accordingly
                ons(oc+1:end, 1:2) = ons(oc+1:end, 1:2) + remod{rc, 2};
            end
        end

    % absolute
    else

        % iterate over onsets
        for oc = 1:lons

            % if condition matches
            if ons(oc, 3) == cc

                % compute time difference
                tdiff = remod{rc, 2} - (ons(oc+1, 1) - ons(oc, 2));

                % then shift subsequent event times accordingly
                ons(oc+1:end, 1:2) = ons(oc+1:end, 1:2) + tdiff;
            end
        end
    end
end

% set on/offsets again
for cc = 1:ncond
    bc.Cond(cc).OnOffsets = ons(ons(:, 3) == cc, 1:2);
end

% set new content
xo.C = bc;

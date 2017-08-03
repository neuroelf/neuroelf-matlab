function xo = prt_Collapse(xo, cfrom, cto, ncol)
% PRT::Collapse  - collapse a set of conditions into one
%
% FORMAT:       [prt] = prt.Collapse(cfrom, newname [, newcolor])
%
% Input fields:
%
%       cfrom       from-condition specification
%                   can be either a 1xN double array or a regexp pattern
%       newname     name for collapsed condition
%       newcolor    1x3 optional color (will be mixed otherwise)
%
% Output fields:
%
%       prt         altered PRT
%
% Examples:
%
%   prt.Collapse([1, 3, 5, 7], 'Odd conditions', [255, 0, 0]);
%   prt.Collapse([2, 3, 4, 5], 'Even conditions', [0, 255, 0]);
%
%  - or -
%
%   prt.Collapse('.*odd.*',  'Odd conditions');
%   prt.Collapse('.*even.*', 'Even conditions');
%
% Using: makelabel, multimatch.

% Version:  v1.1
% Build:    16021017
% Date:     Feb-10 2016, 5:43 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library
global ne_methods;

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt') || ...
    isempty(cfrom) || isempty(cto) || ~ischar(cto)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% get number of conditions
ncon = numel(bc.Cond);

% calc new color ?
if nargin > 3 && isa(ncol, 'double') && numel(ncol) == 3 && ...
    ~any(isinf(ncol) | isnan(ncol) | ncol < 0 | ncol > 255)

    % accept
    ncol = fix(ncol(:)');
    usenc = true;
else
    ncol = [0, 0, 0];
    usenc = false;
end

% resolve cell with names to numbers
if iscell(cfrom)
    try

        % resolve conditions to numbers
        cfrom = unique(ne_methods.multimatch(cfrom(:), prt_ConditionNames(xo)));
        if cfrom(1) < 1
            error('neuroelf:xff:badArgument', 'Source (from) condition not found.');
        end
    catch xfferror
        rethrow(xfferror);
    end
end

% if parametric weights are used
if bc.ParametricWeights
    mw = 0;
    for cc = 1:numel(bc.Cond)
        mw = max(mw, size(bc.Cond(cc).Weights, 2));
    end
    for cc = 1:numel(bc.Cond)
        if size(bc.Cond(cc).Weights, 2) < mw
            bc.Cond(cc).Weights = nan(size(bc.Cond(cc).Weights, 1), mw);
        end
    end
end

% what kind of input
if isa(cfrom, 'double')

    % check min, max, etc...
    cfrom = cfrom(:)';
    if any(isinf(cfrom) | isnan(cfrom) | cfrom < 1 | cfrom > ncon | fix(cfrom) ~= cfrom)
        error('neuroelf:xff:badArgument', 'Bad from specification given.');
    end

    % get conditions
    selected = bc.Cond(cfrom);

    % combine onsets
    oos = zeros(0, 2);
    wgh = zeros(0, size(selected(1).Weights, 2));
    for cc = 1:numel(cfrom)
        oos = [oos; selected(cc).OnOffsets];
        if size(selected(cc).Weights, 2) == size(wgh, 2) && ...
            size(selected(cc).Weights, 1) == size(selected(cc).OnOffsets, 1)
            wgh = [wgh; selected(cc).Weights];
        elseif ~isempty(wgh)
            wgh(end+1:end+size(selected(cc).OnOffsets, 1), :) = 0;
        end
        if ~usenc
            ncol = ncol + selected(cc).Color;
        end
    end
    if ~usenc
        ncol =ncol / numel(cfrom);
    end

    % sort onsets
    [ooso{1:2}] = sort(oos(:, 1));
    oos = oos(ooso{2}, :);
    if ~isempty(wgh)
        wgh = wgh(ooso{2}, :);
    else
        wgh = zeros(numel(ooso{2}), 0);
    end

    % get only first one
    selected = selected(1);
    selected.ConditionName = {cto};
    selected.NrOfOnOffsets = size(oos, 1);
    selected.OnOffsets = oos;
    selected.Weights = wgh;
    selected.Color = fix(ncol);

    % put back into PRT
    bc.Cond(cfrom(1)) = selected;

    % remove collapsed
    bc.Cond(cfrom(2:end)) = [];

    % set new number of conditions
    bc.NrOfConditions = numel(bc.Cond);

% character input
elseif ischar(cfrom) && ~isempty(cfrom) && size(cfrom, 2) == numel(cfrom)

    % create empty condition array
    newconds = bc.Cond;
    newconds(:) = [];

    % build list of matches
    condlist = struct;
    condnums = struct;
    for cc = 1:numel(bc.Cond)

        % get condition name
        condname = bc.Cond(cc).ConditionName{1};

        % match against pattern
        [cdmatch{1:3}] = regexpi(condname, cfrom);
        cdmatch = cdmatch{3};

        % if no match, simply store in newconds!
        if isempty(cdmatch)
            newconds(end + 1) = bc.Cond(cc);
            continue;
        end

        % otherwise get name right!
        newcondname = regexprep(condname, cfrom, cto, 'ignorecase');

        % then create a valid label
        condtag = ne_methods.makelabel(newcondname);

        % and check whether it's already in the list
        if isfield(condlist, condtag)

            % get condition number to combine with
            cmbcnum = condlist.(condtag);

            % combine onsets and weights
            newons = [newconds(cmbcnum).OnOffsets; bc.Cond(cc).OnOffsets];
            if numel(bc.Cond(cc).Weights) == size(bc.Cond(cc).OnOffsets, 1) && ...
                size(bc.Cond(cc).Weights, 2) == size(newconds(cmbcnum).Weights, 2)
                newwgh = [newconds(cmbcnum).Weights; bc.Cond(cc).Weights];
            else
                newwgh = [newconds(cmbcnum).Weights; ...
                    zeros(size(bc.Cond(cc).OnOffsets, 1), size(newconds(cmbcnum).Weights, 2))];
            end

            % sort onsets
            ooso = cell(1, 2);
            [ooso{1:2}] = sort(newons(:, 1));
            newons = newons(ooso{2}, :);
            if ~isempty(newwgh)
                newwgh = newwgh(ooso{2}, :);
            else
                newwgh = zeros(numel(ooso{2}), 0);
            end

            % write onsets into newconds
            newconds(cmbcnum).OnOffsets = newons;
            newconds(cmbcnum).Weights = newwgh;

            % recalculate color
            if ~usenc
                newconds(cmbcnum).Color = (condnums.(condtag) * newconds(cmbcnum).Color + ...
                    bc.Cond(cc).Color) / (condnums.(condtag) + 1);
            end

        % otherwise add it as well
        else
            newconds(end + 1) = bc.Cond(cc);
            if size(newconds(end).Weights, 1) ~= size(newconds(end).OnOffsets, 1)
                newconds(end).Weights = zeros(size(newconds(end).OnOffsets, 1), 0);
            end

            % and exchange the name
            newconds(end).ConditionName{1} = newcondname;

            % and set the correct condition number!
            condlist.(condtag) = numel(newconds);
            condnums.(condtag) = 1;

            % set new color ?
            if usenc
                newconds(end).Color = ncol;
            end
        end
    end

    % fix colors to integral numbers
    if ~usenc
        for cc = 1:numel(newconds)
            newconds(cc).Color = fix(newconds(cc).Color);
        end
    end

    % set newconds into content
    bc.Cond = newconds;
    bc.NrOfConditions = numel(newconds);

% reject the rest for the moment
else
    error('neuroelf:xff:notYetImplemented', 'Collapsing over names not yet implemented.');
end

% reduce empty/non-useful weights arrays again
if bc.ParametricWeights
    for cc = 1:numel(bc.Cond)
        if ~isempty(bc.Cond(cc).Weights) && ~isempty(bc.Cond(cc).OnOffsets)
            nw = size(bc.Cond(cc).Weights, 1);
            for wc = size(bc.Cond(cc).Weights, 2):-1:1
                if isequal(bc.Cond(cc).Weights(ones(1, nw), wc), bc.Cond(cc).Weights(:, wc))
                    bc.Cond(cc).Weights(:, wc) = [];
                end
            end
        end
    end
end

% set new content
xo.C = bc;

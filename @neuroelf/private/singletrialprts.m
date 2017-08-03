function stlist = singletrialprts(prts, subjids, sngtskip, allc)
% singletrialprts  - convert list of PRTs to single-trial PRTs
%
% FORMAT:       stlist = singletrialprts(prts, subjids [, sngtskip [, allc]])
%
% Input fields:
%
%       prts        Px1 cell array of PRT objects
%       subjids     Px1 cell array with matching subject IDs
%       sngtskip    1xC list of trial names to skip during conversion
%       allc        boolean flag, PRTs have all possible conditions (false)
%
% Output fields:
%
%       stlist      Tx1 cell array with names of single trial conditions
%
% Note: given that the objects are altered (in memory), there is no
%       need to return them; please note that if an error occurs
%       the state of the PRT objects is not determined (they should be
%       discarded subsequently!)

% Version:  v1.1
% Build:    16040717
% Date:     Apr-07 2016, 5:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
if nargin < 2 || ...
   ~iscell(prts) || ...
    numel(prts) ~= size(prts, 1) || ...
   ~iscell(subjids) || ...
    numel(prts) ~= numel(subjids)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~iscell(sngtskip)
    sngtskip = {};
end
if nargin < 4 || ...
   ~islogical(allc) || ...
    numel(allc) ~= 1
    allc = false;
end
subjids = subjids(:);
sngtskip = sngtskip(:);
numstudy = numel(prts);
for sc = 1:numstudy
    if ~isxff(prts{sc}, 'prt')
        error( ...
            'neuroelf:BadArgument', ...
            'prts{%d} is not a valid PRT object.', ...
            sc ...
        );
    end
    if ~ischar(subjids{sc}) || ...
        isempty(subjids{sc})
        error( ...
            'neuroelf:BadArgument', ...
            'subjids{%d} is invalid.', ...
            sc ...
        );
    end
    subjids{sc} = subjids{sc}(:)';
end
tpcs = unique(subjids);
tpci = cell(numel(tpcs), 1);
tpcc = cell(0, 2);

% iterate over studies and convert protocols
prtpo = prts;
prtpn = prts;
for sc = 1:numstudy

    % get corresponding subject id index
    sidx = find(strcmp(tpcs, subjids{sc}));
    tpcis = tpci{sidx};

    % ensure format
    if ~iscell(tpcis)
        tpcis = cell(0, 2);
    end

    % get conditions
    conds = prts{sc}.Cond;

    % make sure each condition has a value in sidx call
    skipidx = ones(1, numel(conds));
    for pcc = 1:numel(conds)
        cname = conds(pcc).ConditionName{1};
        if isempty(tpcc) || ...
           ~any(strcmpi(cname, tpcc(:, 1)))
            tpcc(end+1, :) = {cname, conds(pcc).Color};
        end
        if ~isempty(cname)
            cnfound = find(strcmpi(tpcis(:, 1), cname));
            if isempty(cnfound)
                tpcis(end + 1, :) = {cname, 0};
            else
                skipidx(pcc) = tpcis{cnfound, 2} + 1;
            end
        end
    end

    % convert protocol (but keep original order)
    prtpo{sc} = prts{sc}.ConditionNames;
    prts{sc}.ConvertToSingleTrial(struct( ...
        'digits', 3, 'sidx', skipidx, 'skip', {sngtskip}));
    prtpn{sc} = prts{sc}.ConditionNames;

    % update list of skip indices
    for pcc = 1:numel(conds)
        cname = conds(pcc).ConditionName{1};
        if ~isempty(cname)
            cnfound = find(strcmpi(tpcis(:, 1), cname));
            tpcis{cnfound, 2} = tpcis{cnfound, 2} + ...
                size(conds(pcc).OnOffsets, 1);
        end
    end
    tpci{sidx} = tpcis;
end

% create full list of possible condition names
stlisto = unique(cat(1, prtpn{:}));
stlistr = regexprep(stlisto, '_T\d+$', '');
stlist = cell(numel(stlisto), 1);
stcols = zeros(numel(stlist), 3);
stli = 1;
for sc = 1:numel(prtpo)
    stplist = prtpo{sc};
    for pcc = 1:numel(stplist)
        if isempty(stplist{pcc})
            continue;
        end
        stpmatch = strcmp(stlistr, stplist{pcc});
        stpadd = sort(stlisto(stpmatch));
        stlist(stli:stli+numel(stpadd)-1) = stpadd;
        stcols(stli:stli+numel(stpadd)-1, :) = ...
            ones(numel(stpadd), 1) * prts{sc}.Cond(pcc).Color;
        stli = stli + numel(stpadd);
        stplist(strcmp(stplist, stplist{pcc})) = {''};
        stlistr(stpmatch) = {''};
    end
    if stli > numel(stlist)
        break;
    end
end

% all conditions
if allc

    % for each protocol, create the full set
    for sc = 1:numstudy

        % get a match
        cmatch = multimatch(stlist, prts{sc}.ConditionNames);
        pconds = prts{sc}.Cond;
        prts{sc}.Cond(:) = [];

        % iterate over list
        for pcc = 1:numel(stlist)

            % condition in the original protocol
            if cmatch(pcc) > 0

                % set condition
                prts{sc}.Cond(pcc) = pconds(cmatch(pcc));

            % otherwise
            else

                % add a new condition
                prts{sc}.AddCond(stlist{pcc}, zeros(0, 2), stcols(pcc, :));
            end
        end
    end

% just the conditions we need
else

    % for each protocol, bring conditions into this order
    for sc = 1:numstudy
        conds = prts{sc}.Cond;
        condn = cat(1, conds.ConditionName);
        [condo, condoi] = sort(multimatch(condn, stlist));
        if any(condo < 1) || ...
            numel(condo) ~= numel(unique(condo))
            clearxffobjects(prts);
            error( ...
                'xff:BadArgument', ...
                'Error resolving condition names for single trial SDMs.' ...
            );
        end
        prts{sc}.Cond = conds(condoi);
    end
end

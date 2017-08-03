function xo = vtc_CombineAvgTCs(xo, cfunc)
% VTC::CombineAvgTCs  - combine two (or more) conditions of an AvgVTC
%
% FORMAT:       [vtc = ] vtc.CombineAvgTCs(cfunc)
%
% Input fields:
%
%       cfunc       combination function, e.g. 'Cond 1 > Cond 2'
%
% Output fields:
%
%       vtc         VTC with added (combined) time course
%
% Note: functions can be written as
%       C1 + C2 (weighted adding of functions)
%       C1 - C2 (or C1 > C2, difference of conditions)
%       C1 & C1 (conjunction, based on statistic and weighted mean)
%
% Using: ddeblank, lsqueeze, multimatch, splittocellc, tstfrommsw.

% Version:  v1.1
% Build:    16021321
% Date:     Feb-13 2016, 9:28 PM EST
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

% neuroelf library
global ne_methods;
ddeblank     = ne_methods.ddeblank;
multimatch   = ne_methods.multimatch;
splittocellc = ne_methods.splittocellc;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ~ischar(cfunc) || isempty(cfunc)
    error('neuroelf:xff:badArgument', 'Invalid object in call.');
end
bc = xo.C;
if ~isfield(bc.RunTimeVars, 'AvgVTC') || ~islogical(bc.RunTimeVars.AvgVTC) || ...
    numel(bc.RunTimeVars.AvgVTC) ~= 1 || ~bc.RunTimeVars.AvgVTC
    error('neuroelf:xff:badArgument', 'Only valid for Average-VTCs.');
end
cnames = lower(ddeblank(bc.RunTimeVars.ConditionNames));
nconsets = bc.RunTimeVars.NrOfConditionOnsets;
ntcpc = bc.RunTimeVars.NrOfTCsPerCondition;
nvptc = bc.RunTimeVars.NrOfVolumesPerTC;
nvpc = nvptc * ntcpc;

% parse conjunctions
cfunc = cfunc(:)';
cfuno = cfunc;
if ~any(cfunc == '&')
    cfunc = {cfunc};
else
    cfunc = ddeblank(splittocellc(cfunc, '&', true, false));
    cfunc(cellfun('isempty', cfunc)) = [];
    if isempty(cfunc)
        error('neuroelf:xff:badArgument', 'Invalid contrast function specification.');
    end
end

% for each conjunction particle
np = 0;
for cc = 1:numel(cfunc)

    % no > or -
    if ~any(cfunc{cc} == '-' | cfunc{cc} == '>')

        % no negative particle
        cfunc{cc} = {ddeblank(splittocellc(cfunc{cc}, '+', true))};
    else
        cfunc{cc} = ddeblank(splittocellc(cfunc{cc}, '->', true, true));
        cfunc{cc}(cellfun('isempty', cfunc{cc})) = [];
        if numel(cfunc{cc}) ~= 2
            error('neuroelf:xff:badArgument', 'Invalid contrast function specification.');
        end
        cfunc{cc}{1} = ddeblank(splittocellc(cfunc{cc}{1}, '+', true));
        cfunc{cc}{2} = ddeblank(splittocellc(cfunc{cc}{2}, '+', true));
    end

    % check that condition particles exist
    for cpc = 1:numel(cfunc{cc})
        if numel(cfunc{cc}{cpc}) ~= numel(unique(cfunc{cc}{cpc}))
            error('neuroelf:xff:badArgument', 'Double particles not allowed.');
        end
        cfunc{cc}{cpc} = multimatch(lower(cfunc{cc}{cpc}), cnames);
        if any(cfunc{cc}{cpc} == 0)
            error('neuroelf:xff:badArgument', 'Unknown particle in contrast function.');
        end
        np = np + numel(cfunc{cc}{cpc});
    end
    if cpc == 2 && ~isempty(intersect(cfunc{cc}{1}(:), cfunc{cc}{2}(:)))
        error('neuroelf:xff:badArgument', 'Particles may not appear on positive and negative side.');
    end
end
if np < 2
    error('neuroelf:xff:badArgument', 'At least two particles total required.');
end

% get mask
msk = ne_methods.lsqueeze(all(bc.VTCData(ceil(0.5 * nvptc):nvpc:end, :, :, :) ~= 0, 1));
smsk = sum(msk);

% for each conjunction particle
for cc = numel(cfunc):-1:1

    % deal with positive particle
    pp = cfunc{cc}{1};

    % only one
    if numel(pp) == 1

        % get data
        npo = nconsets(pp);
        pmdata = bc.VTCData((pp-1)*nvpc+1:(pp-1)*nvpc+nvptc, msk);
        pedata = bc.VTCData((pp-1)*nvpc+nvptc+1:(pp-1)*nvpc+2*nvptc, msk);
        if ntcpc > 2
            pwdata = npo .* bc.VTCData((pp-1)*nvpc+2*nvptc+1:pp*nvpc, msk);
        else
            pwdata = npo .* ...
                (bc.RunTimeVars.TCOnsetWeights((pp-1)*nvpc+1:(pp-1)*nvpc+nvptc) * ones(1, smsk));
        end

        % get auxiliary data
        ccolor = bc.RunTimeVars.ConditionColors(pp, :);
        consets = bc.RunTimeVars.ConditionOnsets(:, pp);
        cthresh = bc.RunTimeVars.ConditionThresholds(pp, :, :);
        pplcolor = bc.RunTimeVars.Map(pp).RGBLowerThreshPos;
        ppucolor = bc.RunTimeVars.Map(pp).RGBUpperThreshPos;
        pnlcolor = bc.RunTimeVars.Map(pp).RGBLowerThreshNeg;
        pnucolor = bc.RunTimeVars.Map(pp).RGBUpperThreshNeg;
        spnflag = bc.RunTimeVars.Map(pp).ShowPositiveNegativeFlag;
        ptcweights = npo .* bc.RunTimeVars.TCOnsetWeights((pp-1)*nvpc+1:pp*nvpc);

    % average particles
    else
        error('neuroelf:xff:notYetImplemented', 'Weighting not yet implemented.');
    end

    % deal with negative particle
    if numel(cfunc{cc}) > 1
        pp = cfunc{cc}{2};

        % only one
        if numel(pp) == 1

            % get data
            nno = nconsets(pp);
            nmdata = bc.VTCData((pp-1)*nvpc+1:(pp-1)*nvpc+nvptc, msk);
            nedata = bc.VTCData((pp-1)*nvpc+nvptc+1:(pp-1)*nvpc+2*nvptc, msk);
            if ntcpc > 2
                nwdata = nno .* bc.VTCData((pp-1)*nvpc+2*nvptc+1:pp*nvpc, msk);
            else
                nwdata = nno .* ...
                    (bc.RunTimeVars.TCOnsetWeights((pp-1)*nvpc+1:(pp-1)*nvpc+nvptc) * ones(1, smsk));
            end

            % update colors
            pnlcolor = bc.RunTimeVars.Map(pp).RGBLowerThreshPos;
            pnucolor = bc.RunTimeVars.Map(pp).RGBUpperThreshPos;
            nconsets = bc.RunTimeVars.ConditionOnsets(:, pp);
            ntcweights = nno .* bc.RunTimeVars.TCOnsetWeights((pp-1)*nvpc+1:pp*nvpc);

        % average particles
        else
            error('neuroelf:xff:notYetImplemented', 'Weighting not yet implemented.');
        end

        % combine
        npo = npo + nno;
        mdata = pmdata - nmdata;
        [ttdata, utdata, edata] = ...
            ne_methods.tstfrommsw(pmdata, nmdata, pedata, nedata, pwdata, nwdata);
        edata(isinf(edata) | edata <= 0) = max(edata(~isinf(edata) & edata > 0));
        wdata = pwdata + nwdata;
        tcweights = ptcweights + ntcweights;

        % for each onset list
        for oc = 1:numel(consets)

            % combine list
            clist = [consets{oc}; -nconsets{oc}];

            % sort index
            [slist, sindex] = sort(abs(clist(:, 1)));

            % sort with sign
            consets{oc} = clist(sindex, :);
        end

    % otherwise simply copy
    else
        mdata = pmdata;
        edata = pedata;
        wdata = pwdata;
    end

    % re-weigh
    wdata = (1 / npo) .* wdata;
    tcweights = (1 / npo) .* tcweights;

    % more than one particle
    if numel(cfunc) > 1

        % first particle
        if cc == numel(cfunc)
            cmdata = mdata;
            cedata = edata;
            cwdata = wdata;
            cnpo = npo;
            ctcweights = tcweights;

        % conjunction
        else
        end

        % last particle
        if cc == 1
            mdata = cmdata;
            edata = cedata;
            wdata = (1 / numel(cfunc)) .* cwdata;
            npo = cnpo;
            tcweights = ctcweights;
        end
    end
end

% extend data
tv = size(bc.VTCData, 1);
bc.VTCData(tv+1:tv+nvpc, :, :, :) = 0;
if ntcpc > 2
    bc.VTCData(tv+1:tv+nvpc, msk) = cat(1, mdata, edata, wdata);
else
    bc.VTCData(tv+1:tv+nvpc, msk) = cat(1, mdata, edata);
end
bc.RunTimeVars.Map(end+1) = bc.RunTimeVars.Map(1);
bc.RunTimeVars.Map(end).Name = cfuno;
bc.RunTimeVars.Map(end).RGBLowerThreshPos = pplcolor;
bc.RunTimeVars.Map(end).RGBUpperThreshPos = ppucolor;
bc.RunTimeVars.Map(end).RGBLowerThreshNeg = pnlcolor;
bc.RunTimeVars.Map(end).RGBUpperThreshNeg = pnucolor;
bc.RunTimeVars.Map(end).DF1 = npo - 1;
bc.RunTimeVars.Map(end).ShowPositiveNegativeFlag = spnflag;
bc.RunTimeVars.Map(end).OverlayColors = zeros(0, 3);
bc.RunTimeVars.NrOfConditions = bc.RunTimeVars.NrOfConditions + 1;
bc.RunTimeVars.NrOfConditionOnsets(1, end+1) = npo;
bc.RunTimeVars.ConditionColors(end+1, :) = ccolor;
bc.RunTimeVars.ConditionNames{end+1} = cfuno;
bc.RunTimeVars.ConditionOnsets(:, end+1) = consets;
bc.RunTimeVars.ConditionThresholds(end+1, :, :) = cthresh;
bc.RunTimeVars.TCOnsetWeights(tv+1:tv+nvpc, 1) = tcweights;

% set back
xo.C = bc;

% update color scheme
try
    aft_SetColors(xo, numel(bc.RunTimeVars.Map), 'xauto');
catch xfferror
    neuroelf_lasterr(xfferror);
end

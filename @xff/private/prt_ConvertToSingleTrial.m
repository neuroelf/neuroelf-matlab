function xo = prt_ConvertToSingleTrial(xo, opts)
% PRT::ConvertToSingleTrial  - convert to a single-trial version
%
% FORMAT:       [prt =] prt.ConvertToSingleTrial([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .keepempty flag, keep empty conditions (default: false)
%        .keepparam flag, keep parametric value after split (default: false)
%        .sidx      start-index (either Cx1 double or 1x1 doubles in struct)
%        .skip      either cell array of names or 1xN double vector
%                   with conditions to skip conversion (default: [])
%
% Output fields:
%
%       prt         altered PRT
%
% Examples:
%
%   prt.ConvertToSingleTrial;
%

% Version:  v1.1
% Build:    16040710
% Date:     Apr-07 2016, 10:21 AM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
ncon = numel(bc.Cond);
try
    bcond = {bc.Cond(:).ConditionName};
    for cc = 1:ncon
        bcond{cc} = bcond{cc}{1}(:)';
        if ~ischar(bcond{cc}) || isempty(bcond{cc})
            error('BADCONDNAME');
        end
    end
catch xfferror
    error('neuroelf:xff:badObject', ...
        'Invalid ConditionName field for at least one condition: %s.', xfferror.message);
end
if ncon == 0
    return;
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'keepempty') || ~islogical(opts.keepempty) || numel(opts.keepempty) ~= 1
    opts.keepempty = false;
end
if ~isfield(opts, 'keepparam') || ~islogical(opts.keepparam) || numel(opts.keepparam) ~= 1
    opts.keepparam = false;
end
if ~isfield(opts, 'sidx') || ((~isa(opts.sidx, 'double') || numel(opts.sidx) ~= ncon || ...
    any(isinf(opts.sidx(:)) | isnan(opts.sidx(:)) | opts.sidx(:) < 1 | opts.sidx(:) ~= fix(opts.sidx(:)))) && ...
    (~isstruct(opts.sidx) || numel(opts.sidx) ~= 1))
    opts.sidx = ones(1, ncon);
end
if isstruct(opts.sidx)
    sidx = opts.sidx;
    opts.sidx = ones(1, ncon);
    for fc = 1:ncon
        if isfield(sidx, bcond{fc})
            fcidx = sidx.(bcond{fc});
            if isa(fcidx, 'double') && numel(fcidx) == 1 && ~isinf(fcidx) && ...
               ~isnan(fcidx) && fcidx >= 1 && fcidx == fix(fcidx)
                opts.sidx(fc) = fcidx;
            end
        end
    end
end
sidx = opts.sidx - 1;
if ~isfield(opts, 'skip') || isempty(opts.skip) || ...
   (~isa(opts.skip, 'double') && ~ischar(opts.skip) && ~iscell(opts.skip))
    opts.skip = [];
elseif isa(opts.skip, 'double')
    opts.skip = opts.skip(:);
elseif ischar(opts.skip)
    opts.skip = {opts.skip};
end
if iscell(opts.skip)
    for sc = numel(opts.skip):-1:1
        if ~ischar(opts.skip{sc}) || isempty(opts.skip{sc})
            opts.skip(sc) = [];
        else
            opts.skip{sc} = find(~cellfun('isempty', regexpi(bcond, opts.skip{sc})));
            opts.skip{sc} = opts.skip{sc}(:);
        end
    end
    opts.skip = unique(cat(1, opts.skip{:}));
end
if isempty(opts.skip)
    opts.skip = [];
else
    opts.skip(isinf(opts.skip) | isnan(opts.skip) | opts.skip < 1 | opts.skip > ncon) = [];
    opts.skip = unique(round(opts.skip));
end

% work is all but skip
work = 1:ncon;
work(opts.skip) = [];

% get number of single trial conditions
too = 0;
for cc = 1:numel(work)
    if isempty(bc.Cond(work(cc)).OnOffsets) && opts.keepempty
        too = too + 1;
    else
        too = too + size(bc.Cond(work(cc)).OnOffsets, 1);
    end
end

% create new conditions
nc = bc.Cond(1);
nc.NrOfOnOffsets = 1;
nc.Weights = zeros(1, 0);
nc = nc(ones(1, too));

% now iterate over the conditions we need to split
too = 0;
for cc = 1:numel(work)

    % get condition and counting offset
    sc = bc.Cond(work(cc));
    sco = sidx(work(cc));

    % for empty conditions and keepempty
    if isempty(sc.OnOffsets) && opts.keepempty

        % set number of onsets to 0
        too = too + 1;
        nc(too).ConditionName = sc.ConditionName;
        nc(too).NrOfOnOffsets = 0;
        nc(too).OnOffsets = zeros(0, 2);
        nc(too).Weights = [];
        nc(too).Color = sc.Color;
        continue;
    end

    % condition contains single-trial numbers
    if isfield(sc, 'RunTimeVars') && isstruct(sc.RunTimeVars) && numel(sc.RunTimeVars) == 1 && ...
        isfield(sc.RunTimeVars, 'SingleTrialNums') && isa(sc.RunTimeVars.SingleTrialNums, 'double') && ...
        numel(sc.RunTimeVars.SingleTrialNums) == size(sc.OnOffsets, 1)
        stnums = sc.RunTimeVars.SingleTrialNums;
        if isfield(sc.RunTimeVars, 'SingleTrialDigits') && isa(sc.RunTimeVars.SingleTrialDigits, 'double') && ...
            numel(sc.RunTimeVars.SingleTrialDigits) == 1 && any(1:8 == sc.RunTimeVars.SingleTrialDigits)
            stdigits = sc.RunTimeVars.SingleTrialDigits;
        else
            stdigits = ceil(log10(max(stnums) + 0.5));
        end
        stdigits = sprintf('%%s_T%%0%dd', stdigits);
    else
        stnums = [];
    end

    % otherwise
    for scc = 1:size(sc.OnOffsets, 1)

        % increase counter
        too = too + 1;

        % then create new name
        if isempty(stnums)
            nc(too).ConditionName = {sprintf('%s_T%03d', sc.ConditionName{1}, scc + sco)};
        else
            nc(too).ConditionName = {sprintf(stdigits, sc.ConditionName{1}, stnums(scc))};
        end

        % get on/offset of trial
        nc(too).OnOffsets = sc.OnOffsets(scc, :);

        % and possibly weight
        if size(sc.Weights, 1) == size(sc.OnOffsets, 1) && opts.keepparam
            nc(too).Weights = sc.Weights(scc, :);
        end

        % and color
        nc(too).Color = sc.Color;
    end
end

% set new content
if ~isempty(opts.skip)
    kc = bc.Cond(opts.skip);
    bc.Cond = [kc(:)', nc(:)'];
else
    bc.Cond = nc(:)';
end
bc.NrOfConditions = numel(bc.Cond);
xo.C = bc;

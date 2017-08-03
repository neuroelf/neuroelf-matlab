function xo = prt_ImportStimChannels(xo, schan, freq, cspec, opts)
% PRT::ImportStimChannels  - import EEG/MEG experiment stimulus data
%
% FORMAT:       [prt] = prt.ImportStimChannels(schan, freq, cspec [, opts])
%
% Input fields:
%
%       schann      signal channel data
%       freq        frequency of channel data
%       cspec       1xC condition specification with fields
%        .Name      - condition name
%        .Color     - corresponding color
%        .Pattern   - 1xS, 2xS, or 4xS list(s) of 0, 1, and NaNs;
%                     e.g. [0, 1, 0, NaN, 0] means channels 1, 3, and 5
%                     must be off, channel 2 must be on, and the state of
%                     channel 4 does not matter
%            NOTE   * for 1xS list, the channels must remain in this
%                     configuration for the entire stimulation period
%                   * for 2xS lists, the first list is used to detect
%                     the begin, the second list to detect the end of
%                     the stimulation periods
%                   * for 4xS lists, the third pattern means that all
%                     and for the forth pattern any channel combinations
%                     between begin and end must be present for a match
%        .Duration  - if given and > 0, does not look for an
%                     "end of pattern" but uses this for all events (usec)
%       opts        optional 1x1 struct with fields
%        .ExpName   - experiment name
%        .NrOfChan  - if schan is a combined channel (binary!), gives the
%                     number (e.g. 6: values in schan from 0 to 63)
%        .Offset    - value in samples, added to all on/offsets
%        .Smooth    - false|true, smooth durations of events (if close)
%        .Threshold - pulse detection threshold (default: auto)
%        .Tolerance - when to assume that two quickly following channel
%                     state combination changes belong together
%                     (default: auto)
%
% Output field:
%
%       prt         PRT object with added conditions
%
% Using: emptystruct.

% Version:  v1.1
% Build:    16021017
% Date:     Feb-10 2016, 5:35 PM EST
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
if nargin < 4 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt') || ...
   ~isnumeric(schan) || isempty(schan) || ~isa(freq, 'double') || ...
    numel(freq) ~= 1 || isinf(freq) || isnan(freq) || freq <= 0 || freq > 1e6 || ...
   ~isstruct(cspec) || isempty(cspec) || numel(fieldnames(cspec)) ~= 4 || ...
   ~all(strcmp(fieldnames(cspec), {'Name'; 'Color'; 'Pattern'; 'Duration'}))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
ssz = size(schan);
if ssz(1) < ssz(2)
    schan = schan';
    ssz = size(schan);
end
NrOfSamples = ssz(1);
NrOfChannels = ssz(2);
NrOfAddConds = numel(cspec);
if nargin < 5 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if isfield(opts, 'ExpName') && ischar(opts.ExpName) && ~isempty(opts.ExpName)
    bc.Experiment = opts.ExpName(:)';
end
if isfield(opts, 'NrOfChan') && isa(opts.NrOfChan, 'double') && numel(opts.NrOfChan) == 1 && ...
   ~isinf(opts.NrOfChan) && ~isnan(opts.NrOfChan) && opts.NrOfChan > 1 && ...
    opts.NrOfChan <= 16 && opts.NrOfChan == fix(opts.NrOfChan)
    NrOfChannels = opts.NrOfChan;
    opts.Threshold = 1;
end
if ~isfield(opts, 'Offset') || ~isa(opts.Offset, 'double') || ...
    numel(opts.Offset) ~= 1 || isinf(opts.Offset) || isnan(opts.Offset)
    opts.Offset = 0;
end
if ~isfield(opts, 'Smooth') || ~islogical(opts.Smooth) || numel(opts.Smooth) ~= 1
    opts.Smooth = false;
end
if ~isfield(opts, 'Threshold') || ~isa(opts.Threshold, 'double') || ...
    numel(opts.Threshold) ~= 1 || isinf(opts.Threshold)
    opts.Threshold = NaN;
end
if ~isfield(opts, 'Tolerance') || ~isa(opts.Tolerance, 'double') || ...
    numel(opts.Tolerance) ~= 1 || isinf(opts.Tolerance)
    opts.Tolerance = NaN;
end

% check condition specs
for cc = 1:NrOfAddConds
    cs = cspec(cc);
    if ~ischar(cs.Name) || isempty(cs.Name) || ~isa(cs.Color, 'double') || numel(cs.Color) ~= 3 || ...
        any(isinf(cs.Color) | isnan(cs.Color) | cs.Color < 0 | cs.Color > 256) || ~isa(cs.Pattern, 'double') || ...
       ~any(numel(cs.Pattern) == [NrOfChannels, 2 * NrOfChannels, 4 * NrOfChannels]) || ...
       ~all(cs.Pattern(:) == 0 | cs.Pattern(:) == 1 | isnan(cs.Pattern(:))) || ...
       ~isa(cs.Duration, 'double') || ~any(numel(cs.Duration) == [0, 1])
        error('neuroelf:xff:badArgument', 'Invalid cspec (%d).', cc);
    end
    cs.Name = cs.Name(:)';
    cs.Color = cs.Color(:)';
    if any(cs.Color ~= fix(cs.Color))
        cs.Color = fix(max(cs.Color * 255, 255));
    end
    if numel(cs.Pattern) > NrOfChannels
        if size(cs.Pattern, 2) ~= NrOfChannels && any(size(cs.Pattern, 2) == [1, 2, 4])
            cs.Pattern = cs.Pattern';
        else
            cs.Pattern = reshape(cs.Pattern, [numel(cs.Pattern) / NrOfChannels, NrOfChannels]);
        end
    else
        cs.Pattern = cs.Pattern(:)';
    end
    cs.Duration = fix(cs.Duration);
    cspec(cc) = cs;
end

% what threshold
if isnan(opts.Threshold)
    scmed = median(schan);
    if any(scmed ~= 0)
        schan = schan - ones(NrOfSamples, 1) * scmed;
    end
    opts.Threshold = max(schan) / 2;
    schan = single(schan >= (ones(NrOfSamples, 1) * opts.Threshold));
else
    if size(schan, 2) == NrOfChannels
        schan = single(schan >= opts.Threshold);
    end
end

% find places in stimulus channels that are of interest
iint = find(any([schan(1, :); diff(schan)] ~= 0, 2));
ichan = schan(iint, :);

% expand combined channel ?
if size(schan, 2) == 1 && NrOfChannels > 1

    % pre-check ichan
    badic = find(ichan > (2 ^ NrOfChannels));
    if ~isempty(badic)
        warning('neuroelf:xff:outOfBounds', ...
            'Too high combined channel value, e.g. at sample %d.', iint(badic(1)));
        ichan(badic) = mod(ichan(badic), 2 ^ NrOfChannels);
    end

    % create combined matrix
    cchan = single(0);
    cchan(size(ichan, 1), NrOfChannels) = 0;

    % fill with values
    for cc = NrOfChannels:-1:1
        cchan(:, cc) = ichan >= (2 ^ (cc - 1));
        ichan = ichan - cchan(:, cc) * (2 ^ (cc - 1));
    end

    % post check ichan
    badic = find(ichan ~= 0);
    if ~isempty(badic)
        warning('neuroelf:xff:outOfBounds', ...
            'Invalid combined channel value, e.g. at sample %d.', iint(badic(1)));
    end

    % re-write ichan
    ichan = cchan;
end

% find entries that are below the tolerance
dint = [NrOfSamples; diff(iint)];
if isnan(opts.Tolerance)
    sdint = sort(dint);
    ddint = find([1; diff(sdint)] > 3);
    if isempty(ddint)
        error('neuroelf:xff:automationError', 'Error finding timing tolerance level.');
    end
    opts.Tolerance = floor((sdint(ddint(1) - 1) + sdint(ddint(1))) / 2);
end
intol = find(dint < opts.Tolerance);

% replace data in samples before that by those
for ic = numel(intol):-1:1
    itc = intol(ic);
    if itc > 1
        ichan(itc - 1, :) = ichan(itc, :);
    end
end

% remove from iint and ichan
iint(intol) = [];
ichan(intol, :) = [];
numi = size(ichan, 1);

% create and loop over conditions to add
aconds = ne_methods.emptystruct({'ConditionName', 'NrOfOnOffsets', 'OnOffsets', 'Weights', 'Color'});
aconds(NrOfAddConds).OnOffsets = [];
for cc = 1:NrOfAddConds

    % get spec
    thisspec = cspec(cc);

    % set basic properties
    aconds(cc).ConditionName = {thisspec.Name};
    aconds(cc).Color = thisspec.Color;

    % find matching "start" onsets
    thispatt = thisspec.Pattern(1, :);
    tpi = ~isnan(thispatt);
    bmatch = find(all(ones(numi, 1) * thispatt(tpi) == ichan(:, tpi), 2));

    % in not by duration
    if thisspec.Duration == 0

        % any non matching next entries mean end
        if size(thisspec.Pattern, 1) < 2

            % get end time points and build OnOffsets
            if bmatch(end) < numi
                ematch = iint(bmatch + 1);
            else
                ematch = iint(bmatch(1:end-1) + 1);
            end

        % end is defined by channel combination
        else

            % look for end matches first
            thispatt = thisspec.Pattern(2, :);
            tpi = ~isnan(thispatt);
            ematch = find(all(ones(numi, 1) * thispatt(tpi) == ichan(:, tpi), 2));

            % match begins with ends
            sst = intersect(bmatch, ematch);
            if ~isempty(sst)
                bmatch = setdiff(bmatch, sst);
                ematch = setdiff(ematch, sst);
            end
            if isempty(bmatch) || isempty(ematch)
                error('neuroelf:xff:internalError', 'Error matching begins and ends, empty array.');
            end

            % linearize bmatch and find corresponding end
            bmatch = bmatch(:);
            bmatch(:, 2) = NaN;
            for bmc = 1:size(bmatch, 1)
                emc = find(ematch > bmatch(bmc, 1));
                if isempty(emc) && bmc < size(bmatch, 1)
                    error('neuroelf:xff:internalError', ...
                        'No END match found for period starting at sample %d.', iint(bmatch(bmc, 1)));
                end
                if ~isempty(emc)
                    bmatch(bmc, 2) = ematch(emc(1));
                else
                    bmatch(bmc, 2) = numi + 1;
                end
            end
            ematch = bmatch(:, 2);
            bmatch = bmatch(:, 1);

            % possibly remove interfering entries
            if size(thisspec.Pattern, 1) > 3

                mallpatt = thisspec.Pattern(3, :);
                manypatt = thisspec.Pattern(4, :);
                alli = ~isnan(mallpatt);
                anyi = ~isnan(manypatt);
                % iterate over matches
                for mc = numel(bmatch):-1:1

                    % we believe it matches all but not any (!)
                    mall = true;
                    many = false;

                    % get from and to indices
                    if bmatch(mc) >= numi
                        continue;
                    end
                    fromi = bmatch(mc) + 1;
                    if ematch(mc) > numi
                        toi = numi;
                    else
                        toi = ematch(mc) - 1;
                    end

                    % iterate over pattern combinations
                    for pc = fromi:toi

                        % first check "all" pattern
                        if any(mallpatt(alli) ~= ichan(pc, alli))
                            mall = false;
                            break;
                        end

                        % then check "any" pattern
                        if all(manypatt(anyi) == ichan(pc, alli))
                            many = true;
                        end
                    end

                    % if not mall or not many, delete
                    if ~mall || ~many
                        bmatch(mc) = [];
                        ematch(mc) = [];
                    end
                end
            end

            % check last ematch against numi
            if ~isempty(ematch) && ematch(end) > numi
                ematch(end) = [];
            end
        end
    end

    % get time points
    btimes = 1000000 * (iint(bmatch) + opts.Offset) / freq;

    % duration is expressly given
    if thisspec.Duration > 0

        % build Offsets with duration
        etimes = btimes + thisspec.Duration;

    % otherwise
    else

        % determine from ematch
        if numel(ematch) == numel(bmatch)
            etimes = 1000000 * (iint(ematch) + opts.Offset) / freq;
        else
            etimes = 1000000 * (iint(ematch) + opts.Offset) / freq;
            etimes(end + 1) = 1000000 * NrOfSamples / freq;
        end

        % smooth ?
        if opts.Smooth && numel(btimes) > 1
            dtimes = etimes - btimes;
            if std(dtimes) < (2000000 * opts.Tolerance / freq)
                etimes = btimes + mean(dtimes);
            end
        end

    end

    % set (NrOf) OnOffsets + standard weights
    aconds(cc).OnOffsets = round([btimes(:), etimes(:)]);
    aconds(cc).Weights = zeros(numel(btimes), 0);
    aconds(cc).NrOfOnOffsets = numel(btimes);
end

% put into PRT
bc.Cond(end + 1:end + NrOfAddConds) = aconds;
bc.NrOfConditions = numel(bc.Cond);
xo.C = bc;

% clean up PRT
prt_CleanUp(xo);

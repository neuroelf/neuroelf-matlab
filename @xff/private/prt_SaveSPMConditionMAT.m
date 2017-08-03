function mcont = prt_SaveSPMConditionMAT(xo, filename, opts)
% PRT::SaveSPMConditionMAT  - save onsets into SPM-compatible MAT file
%
% FORMAT:       [mcont] = prt.SaveSPMConditionMAT(filename [, opts]);
%
% Input fields:
%
%       filename    output filename (use 'struct' if not to be saved)
%       opts        optional settings
%        .dropconds cell array with conditions to drop from MAT file
%        .erthresh  event-related threshold (default: 200ms)
%        .toffset   time-offset, added to all onset/offsets (in PRT unit)
%
% Output fields:
%
%       mcont       contents of SPM-compatible mat file
%
% Using: multimatch.

% Version:  v1.1
% Build:    16021016
% Date:     Feb-10 2016, 4:45 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2012, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt') || ...
   ~ischar(filename) || isempty(filename)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
filename = filename(:)';
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dropconds') || ~iscell(opts.dropconds)
    opts.dropconds = {};
end
if ~isfield(opts, 'erthresh') || ~isa(opts.erthresh, 'double') || numel(opts.erthresh) ~= 1 || ...
    isinf(opts.erthresh) || isnan(opts.erthresh) || opts.erthresh < 0
    if lower(bc.ResolutionOfTime(1)) == 'm'
        opts.erthresh = 200;
    else
        opts.erthresh = 0;
    end
end
if ~isfield(opts, 'toffset') || ~isa(opts.toffset, 'double') || ...
    numel(opts.toffset) ~= 1 || isinf(opts.toffset) || isnan(opts.toffset)
    opts.toffset = 0;
end

% time in MS -> convert to second
if lower(bc.ResolutionOfTime(1)) == 'm'
    tfac = 0.001;
    opts.erthresh = tfac * opts.erthresh;
else
    tfac = 1;
    opts.toffset = opts.toffset - 1;
end

% drop conditions as requested
conds = bc.Cond;
names = {conds.ConditionName};
for cc = 1:numel(names)
    try
        names{cc} = names{cc}{1};
    catch xfferror
        rethrow(xfferror);
    end
end
if ~isempty(opts.dropconds)
    try
        dconds = ne_methods.multimatch(opts.dropconds(:), names(:));
        dconds = dconds(dconds > 0);
        conds(dconds) = [];
        names(dconds) = [];
    catch xfferror
        rethrow(xfferror);
    end
end

% create output arrays
nconds = numel(conds);
onsets = cell(1, nconds);
durations = cell(1, nconds);
pmod = repmat(struct('name', [], 'param', [] , 'poly', []), [1, nconds]);

% fill arrays
wc = 0;
for cc = 1:nconds

    % get condition details
    c = conds(cc);

    % get on/offsets (with global offset)
    oo = tfac .* (c.OnOffsets + opts.toffset);

    % compute durations
    dur = oo(:, 2) - oo(:, 1);

    % set ER durations
    dur(dur <= opts.erthresh) = 0;

    % cut oo
    oo(:, 2) = [];

    % add to arrays
    onsets{cc} = oo;
    durations{cc} = dur;

    % parametric modulation?
    if ~isempty(c.Weights)

        % get weights
        w = c.Weights;
        nw = size(w, 2);

        % setup structure
        pmod(cc).name = cell(1, nw);
        pmod(cc).param = cell(1, nw);
        pmod(cc).poly = repmat({1}, [1, nw]);

        % for each weight
        for wc = 1:size(w, 2)

            % set name
            pmod(cc).name{wc} = sprintf('%s-pmod%d', names{cc}, wc);

            % set weights
            pmod(cc).param{wc} = w(:, wc);
        end
    end
end

% return as fields
mcont = struct('names', {names}, 'onsets', {onsets}, 'durations', {durations}, 'pmod', pmod);

% save
if ~strcmpi(filename, 'struct')

    % try to save
    try
        if wc > 0
            save(filename, 'names', 'onsets', 'durations', 'pmod', '-v6');
        else
            save(filename, 'names', 'onsets', 'durations', '-v6');
        end
    catch xfferror
        rethrow(xfferror);
    end
end

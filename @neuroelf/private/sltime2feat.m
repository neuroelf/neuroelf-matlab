function [f, l, oos, unames] = sltime2feat(t, oos, varargin)
%SLTIME2FUNC  Extract features from time courses for Searchlight analyses.
%   [F, L] = SLTIME2FEAT(TC, OOS) extracts features F (samples-by-features
%   matrix) along the labels, by applying the temporal information in OOS.
%
%   TC must be a Time-by-SourceFeatures matrix, and if several runs are
%   given, the onsets must either be corrected for run-length, or the
%   length of each run must be given either by supplying PRT contents where
%   the RunTimeVars.NrOfTimePoints is set, or by supplying the .runlens
%   field in the options (see below).
%
%   OOS can be either the concantenated content of PRT protocol objects,
%   supporting several runs (in which case run-length must be given as an
%   additional field in PRT(RUN).RunTimeVars.NrOfTimePoints, and to convert
%   millisecond protocols to volumes, the TR must be given in additional
%   field PRT(RUN).RunTimeVars.TR). Alternatively, OOS can be an array of
%   onsets, in the form of a matrix wherein each row is coded as
%   [RUN, COND, ONSET, OFFSET]. This is also the format any PRT(s) input
%   will be converted into (and returned as the third, optional output).
%
%   Numeric OOS must be given in unit of samples (TRs), with the first
%   timepoint of TC being indexed by 1 (MATLAB-based indexing). If OPTS is
%   not given, the function extracts the average of the 2nd and 3rd time
%   point following the onset (assuming a TR of 2 seconds, and rounding the
%   onsets to full integers). See below for more options.
%
%   [F, L] = SLTIME2FEAT(TC, OOS, NUIS, ICNT) uses the matrices NUIS and
%   ICNT (which must have a size of NxTime) to remove nuisance from the
%   data in TC whereas ICNT must be the inverse of COV(NUIS) * NUIS' to
%   compute least-squares estimates of fit, and NUIS is then the matrix
%   of regressors to be removed (if certain effects are to be retained,
%   simply set those columns to all 0, which in effect estimates the
%   de-correlated betas from TC but leaves those effects in the data).
%
%   [F, L] = SLTIME2FEAT(..., OPTS) allows to set additional options by
%   supplying a 1x1 struct with settings as the last input. Recognized
%   options are
%
%    .baseline  for avg, stack, and wavg computation, pick of baseline
%               volume(s) relative to onset; if multiple values are given,
%               average those, default: 0
%    .comp      which computation to perform on source features, one of
%               {'avg'}, 'estim', 'hrfcov', 'stack', 'wavg'
%    .hrf       HRF function to estimate covariance or weighted average
%    .imeth     interpolation method (if needed, default: lanczos3)
%    .round     round onsets (default: true if comp is '(w)avg' or 'stack')
%               if round is false, the function internally makes calls to
%               flexinterpn (with pre-made kernels to reduce overhead)
%    .runlens   an Rx1 set of numbers indicating the length of each
%               functional run; importantly its sum must match the total
%               number of volumes (by default: consider the data to be
%               from a single run)
%    .timepts   time points to draw on for averaging (default: [2, 3])
%
%   [F, L, OOS] = SLTIME2FEAT(...) also returns a OOS array suitable for
%   subsequent calls.
%
%   [F, L, OOS, UNAMES] = SLTIME2FEAT(TC, OOS, ...) also returns the unique
%   names of the conditions labeled in the second column of OOS if (and
%   only if) the OOS input to SLTIME2FEAT is given as an array of PRT
%   objects. Otherwise an error will be produced by the call (unassigned
%   output UNAMES).

% Version:  v1.1
% Build:    16040114
% Date:     Apr-01 2016, 2:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% persistent interpolation kernels
persistent t2fik;
if numel(t2fik) ~= 1
    t2fik = struct;
end

% first arguments check
if nargin < 2 || (~isa(t, 'double') && ~isa(t, 'single')) || isempty(t) || ndims(t) > 2 || size(t, 1) < 4
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
numt = size(t, 1);
numf = size(t, 2);
if nargin > 2 && isstruct(varargin{end}) && numel(varargin{end}) == 1
    opts = varargin{end};
else
    opts = struct;
end
if ~isfield(opts, 'comp') || ~ischar(opts.comp) || ~any(strcmpi(opts.comp(:)', ...
   {'a', 'avg', 'e', 'estim', 'h', 'hrfcov', 's', 'stack', 'w', 'wavg'}))
    comp = 'a';
else
    comp = lower(opts.comp(1));
end
if isfield(opts, 'round') && islogical(opts.round) && numel(opts.round) == 1
    doround = opts.round;
else
    doround = any(comp == 'asw');
end
if ~doround
    if ~isfield(opts, 'imeth') || ~ischar(opts.imeth) || isempty(opts.imeth) || ~isrealvarname(opts.imeth(:)')
        imeth = 'lanczos3';
    else
        imeth = opts.imeth(:)';
    end
    if ~isfield(t2fik, imeth)
        [ndata, t2fik.(imeth)] = flexinterpn_method(zeros(10, 1), [Inf; 2.5; 5; 10], imeth);
    end
end
if ~isfield(opts, 'timepts') || ~isa(opts.timepts, 'double') || isempty(opts.timepts) || ...
    any(isinf(opts.timepts(:)) | isnan(opts.timepts(:)) | opts.timepts(:) < 0)
    timepts = [2; 3];
else
    timepts = unique(opts.timepts(:));
    if numel(timepts) > 8
        error('neuroelf:general:tooManyValues', 'Too many time points to average.');
    end
end
if nargin > 3 && isa(varargin{1}, 'double') && isa(varargin{2}, 'double') && ...
    ndims(varargin{1}) == 2 && ndims(varargin{2}) == 2 && size(varargin{1}, 1) == numt && ...
    size(varargin{2}, 1) == size(varargin{1}, 2) && size(varargin{2}, 2) == numt
    nsdm = varargin{1};
    isdm = varargin{2};
else
    nsdm = [];
end

% PRT content
if isstruct(oos) && isfield(oos, 'Cond')
    
    % total number of conditions (match labels later)
    try

        % store in prts (we will repurpose oos later)
        prts = oos;
        prtv = {prts.RunTimeVars};
        rlen = 0;

        % PRTs still in MS?
        for pc = 1:numel(prts)
            if lower(prts(pc).ResolutionOfTime(1)) ~= 'v'
                volfac = 1 / prtv{pc}.TR;
                voldiff = 1 + rlen;
            else
                volfac = 1;
                voldiff = rlen;
            end
            if volfac ~= 1 || voldiff ~= 0
                for cc = 1:numel(prts(pc).Cond)
                    prts(pc).Cond(cc).OnOffsets = ...
                        voldiff + volfac .* prts(pc).Cond(cc).OnOffsets;
                end
            end

            % add to rlen if needed
            if pc < numel(prts)
                rlen = rlen + prtv{pc}.NrOfTimePoints;
            end
        end

        % get all conditions
        tconds = lsqueezecells({prts.Cond});

        % conditions/run
        rconds = cellfun('prodofsize', tconds);

        % concatenate
        tconds = cat(1, tconds{:});

        % get the number of onsets per condition
        tnoos = cellfun('size', {tconds.OnOffsets}, 1);

        % create required matrix and names for resolving to labels
        oos = zeros(sum(tnoos), 4);
        cnames = lower(cat(1, tconds.ConditionName));
        [unames, unamei, unameir] = unique(cnames);

        % iterate while keeping track of runs and conditions
        ci = 1;
        rv = 1;
        rc = 1;
        for cc = 1:numel(tconds)
            oos(ci:ci+tnoos(cc)-1, :) = ...
                [ones(tnoos(cc), 1) * [rv, unameir(cc)], tconds(cc).OnOffsets(:, 1:2)];

            % keep track of run (value)
            rc = rc + 1;
            if rc > rconds(rv)
                rc = 1;
                rv = rv + 1;
            end
            ci = ci + tnoos(cc);
        end

        % any invalid onsets
        if any(oos(:, 3) < 1)
            error('Invalid onset(s).');
        end

    catch ne_eo;
        error('neuroelf:xff:badObjectContent', 'Invalid PRT contents.');
    end

% check for onsets
elseif ~isa(oos, 'double') || isempty(oos) || ndims(oos) > 2 || size(oos, 2) ~= 4 || ...
    any(any(isinf(oos), 2) | any(isnan(oos), 2) | any(oos < 1, 2) | any(oos(:, 1:2) ~= round(oos(:, 1:2)), 2))
    error('neuroelf:general:badArgument', 'Bad on/offsets argument.');
end

% remove nuisance
if ~isempty(nsdm)
    tb = isdm * t;
    t = t - nsdm * tb;
end

% round
if doround
    oos = round(oos);
end
numoos = size(oos, 1);

% output features
numof = numf;
if comp(1) == 's'
    numof = numof * numel(timepts);
end

% create outputs
f = NaN(numoos, numof);
l = zeros(numoos, 1);

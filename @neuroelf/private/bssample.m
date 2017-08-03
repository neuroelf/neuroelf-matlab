function bs = bssample(n, opts)
% bssample  - bootstrap sample generator
%
% FORMAT:       bs = bssample(n [, opts])
%
% Input fields:
%
%       n           number of samples to draw
%       opts        optional settings
%        .maxsmp    maximum number of same index in sample
%        .minsmp    minimum number of samples to cover
%        .numparm   number of parameters to sample (next dimension)
%        .numsmp    number of re-samples to create (next dimension)
%        .perm      permutation forced (maxsmp = minsmp = n)
%        .uniquem   NxP model predictors, for which drawn samples must
%                   be valid and unique
%        .uniquemr  flag, if model is a combination (one parameter) and
%                   each regressor must be uniquely identified
%
% Output fields:
%
%       bs          NxPxS indices to sample data with

% Version:  v0.9b
% Build:    10082816
% Date:     Aug-24 2010, 12:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
if nargin < 1 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n < 1
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
n = floor(real(n));

% no/invalid options
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1

    % handle default case quickly!
    bs = ceil(n .* rand(n, 1));
    return;
end

% options
if ~isfield(opts, 'maxsmp') || ...
   ~isa(opts.maxsmp, 'double') || ...
    numel(opts.maxsmp) ~= 1 || ...
    isinf(opts.maxsmp) || ...
    isnan(opts.maxsmp) || ...
    opts.maxsmp < 1
    opts.maxsmp = n;
else
    opts.maxsmp = ceil(min(n, opts.maxsmp));
end
if ~isfield(opts, 'minsmp') || ...
   ~isa(opts.minsmp, 'double') || ...
    numel(opts.minsmp) ~= 1 || ...
    isinf(opts.minsmp) || ...
    isnan(opts.minsmp) || ...
    opts.minsmp < 1
    opts.minsmp = 1;
else
    opts.minsmp = floor(min(n, opts.minsmp));
end
if ~isfield(opts, 'numparm') || ...
   ~isa(opts.numparm, 'double') || ...
    numel(opts.numparm) ~= 1 || ...
    isinf(opts.numparm) || ...
    isnan(opts.numparm) || ...
    opts.numparm < 1
    opts.numparm = 1;
else
    opts.numparm = floor(opts.numparm);
end
if ~isfield(opts, 'numsmp') || ...
   ~isa(opts.numsmp, 'double') || ...
    numel(opts.numsmp) ~= 1 || ...
    isinf(opts.numsmp) || ...
    isnan(opts.numsmp) || ...
    opts.numsmp < 1
    opts.numsmp = [];
else
    opts.numsmp = floor(opts.numsmp);
    if opts.numsmp == 1
        opts.numsmp = [];
    end
end
if ~isfield(opts, 'perm') || ...
   ~islogical(opts.perm) || ...
    numel(opts.perm) ~= opts.numparm
    opts.perm = false(1, opts.numparm);
else
    opts.perm = opts.perm(:);
end
if ~isfield(opts, 'uniquem') || ...
   ~isa(opts.uniquem, 'double') || ...
    ndims(opts.uniquem) > 2 || ...
    size(opts.uniquem, 1) ~= n || ...
   (opts.numparm > 1 && ...
    any(size(opts.uniquem, 2) ~= [1, opts.numparm])) || ...
    any(isinf(opts.uniquem(:)) | isnan(opts.uniquem(:)))
    opts.uniquem = [];
end
unim = opts.uniquem;

% sample size
np = opts.numparm;
ns = opts.numsmp;
pp = opts.perm;
ps = [n, 1, ns];
ss = [n, np, ns];
nu = size(unim, 2);

% only one sample
if ss(1) == 1

    % all indices must be 1, depending on number of parameters
    bs = ones(ss);
    return;
end

% no need to heed maxsmp/minsmp/uniquem
maxs = opts.maxsmp;
mins = opts.minsmp;
if maxs == n && ...
    mins == 1 && ...
   ~any(pp) && ...
    nu == 0

    % sample away!
    bs = ceil(n .* rand(ss));
    return;
end

% forced perm for all?
if maxs == 1 || ...
    mins == n
    pp(:) = true;
end

% create indices
bs = zeros(ss);

% iterate over parameters
for pc = 1:np

    % for permutation parameters, take a different turn
    if pp(pc)

        % use sort to create a unique mapping
        [t, bs(:, pc, :)] = sort(rand(ps), 1, 'ascend');

        % needs unique sampling
        if nu > 0

            % model has numparm columns
            if nu == np

                % generate a random variable
                rr = randn(n, 1);

                % regress permutations against a random variable
                [t, rb] = mmregress(reshape(subsref(unim(:, pc), struct( ...
                    'type', '()', 'subs', {{squeeze(bs(:, pc, :))}})), ps), rr);

                % remove invalid items
                rb(isinf(rb) | isnan(rb)) = 0;

            % model has more columns
            else

                % preset rb
                rb = ones(1, 1, ns);

                % generate random variables
                if opts.uniquemr
                    rr = repmat(randn(n, 1), 1, nu);
                else
                    rr = randn(n, nu);
                end

                % iterate over model columns
                for mc = 1:nu

                    % regress permutations against a random variable
                    [t, rbp] = mmregress(reshape(subsref(unim(:, mc), struct( ...
                        'type', '()', 'subs', {{squeeze(bs(:, pc, :))}})), ps), rr(:, mc));

                    % remove invalid items
                    rbp(isinf(rbp) | isnan(rbp)) = 0;

                    % multiply to get unique rb
                    rb = rb .* rbp + rbp;
                end
            end

            % sort values
            rb = sort(rb(:));

            % find indices that are not unique
            gi = (rb ~= 0 & diff([rb; rb(1)]) ~= 0);
            t = find(~gi);

            % we need to find replacement for those
            if ~isempty(t)

                % try to estimate number of unique combinations from data
                nt = numel(t);
                ng = ns - nt;
                xd = (nt + 1) / ng;
                unic = floor(0.6 * ng / xd);

                % generously estimate number of additionally required draws
                addc = 50 + ceil(1.01 * nt * ((1 + (nt * unic) / (ng * (unic - ng))) .^ 2));

                % way too many? then this is not working
                if addc > (2 * ns)
                    warning( ...
                        'neuroelf:BadArgument', ...
                       ['Too hard to find unique permutations for parameter %d ' ...
                        '(with %d estimated unique samples.)'], ...
                        pc, unic ...
                    );
                    continue;
                end

                % create additional samples
                [t, bsa] = sort(rand([n, 1, addc]), 1, 'ascend');

                % concatenate samples
                bsa = cat(3, bs(:, pc, :), bsa);
                nsa = ns + addc;

                % needs unique sampling
                if nu > 0

                    % model has numparm columns
                    if nu == np

                        % regress permutations against a random variable
                        [t, rb] = mmregress(reshape(subsref(unim(:, pc), ...
                            struct('type', '()', 'subs', ...
                            {{squeeze(bsa)}})), [n, 1, nsa]), rr);

                        % remove invalid items
                        rb(isinf(rb) | isnan(rb)) = 0;

                    % model has more columns
                    else

                        % preset rb
                        rb = ones(1, 1, addc);

                        % iterate over model columns
                        for mc = 1:nu

                            % regress permutations against a random variable
                            [t, rbp] = mmregress(reshape(subsref(unim(:, mc), ...
                                struct('type', '()', 'subs', ...
                                {{squeeze(bsa)}})), [n, 1, nsa]), rr(:, mc));

                            % remove invalid items
                            rbp(isinf(rbp) | isnan(rbp)) = 0;

                            % multiply to get unique rb
                            rb = rb .* rbp + rbp;
                        end
                    end

                    % sort *all* values
                    [rb, rbi] = sort(rb(:));

                    % find indices that are unique
                    gi = (rb ~= 0 & diff([rb; rb(1)]) ~= 0);

                    % get the indices of the good ones
                    rbig = sort(rbi(gi));

                    % test whether we have enough
                    if numel(rbig) >= ns

                        % get the first ns items
                        bs(:, pc, :) = bsa(:, 1, rbig(1:ns));
                    else

                        % give a warning
                        warning( ...
                            'neuroelf:InternalError', ...
                           ['%d samples not unique for parameter %d ' ...
                            '(%d samples requested, %d unique samples estimated, ' ...
                            '%d samples found on first pass, %d additional generated)'], ...
                            ns - numel(rbig), pc, ns, unic, ng, addc ...
                        );

                        % add some of the non-unique ones (*sigh*)
                        rbit = sort(rbi(~gi));

                        % combine the samples
                        bs(:, pc, :) = cat(3, bsa(:, 1, rbig), ...
                            bsa(:, 1, rbit(1:ns - numel(rbig))));
                    end
                end
            end
        end

    % standard bootstrapping sampling (with replacement)
    else

        % no further conditions to be heeded
        if maxs == n && ...
            mins == 1

            % fill indices
            bs(:, pc, :) = ceil(n .* rand([n, 1, ns]));

            % and move on
            continue;
        end

        % create boolean flag list
        bb = true(1, ns);

        % repeat sampling until conditions satisfactory
        while any(bb)

            % fill indices
            bs(:, pc, bb) = ceil(n .* rand([n, 1, sum(bb)]));

            % get indices (to set state)
            bbi = find(bb);

            % compute histogram (with indices as boundaries)
            bsh = histcount(bs(:, pc, bb), 1, n, 1, 1);

            % eliminate those who fail the tests
            if maxs < n && ...
                mins > 1
                bb(bbi(all(bsh <= maxs, 1) & ...
                    (sum(bsh > 0, 1) >= mins))) = false;
            elseif maxs < n
                bb(bbi(all(bsh <= maxs, 1))) = false;
            else
                bb(bbi(sum(bsh > 0, 1) >= mins)) = false;
            end
        end
    end
end

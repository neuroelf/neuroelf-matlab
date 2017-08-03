function [words, solved] = solvewordbrain(puzzle, sizes, opts, varargin)
%SOLVEWORDBRAIN  Solve WordBrain puzzle.
%   WORDS = SOLVEWORDBRAIN(PUZZLE, SIZES) solves the puzzle in PUZZLE with
%   the given word lengths in sizes (1xW double).
%
%   WORDS = SOLVEWORDBRAIN(PUZZLE, SIZES, OPTS) allows to set options in a
%   1x1 struct with fields
%
%     .knowns       only use words from the known list of words (false)
%     .matches      known words (in order of sizes)
%     .numwords     default as many as sizes
%     .simplify     simplify (top level only, default: false)
%     .smart        try a number of techniques (default: true)
%     .timeout      number of seconds to spend at most (default: Inf)
%
%   [WORDS, SOLVED] = SOLVEWORDBRAIN(PUZZLE, SIZES, OPTS) will also return
%   an image with colored letters showing the paths graphically, if (and
%   only if) a solution can be found.

% Version:  v1.1
% Build:    16060708
% Date:     Jun-07 2016, 8:08 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% persistent dictionary
persistent wbdict;
if isempty(wbdict) || ~isstruct(wbdict) || isempty(fieldnames(wbdict))
    try

        % prepare dictionary
        wbdict = struct;

        % read a dictionary (wordlist)
        wbd = unique(lower(splittocell(asciiread('/usr/share/dict/words'), char(10))));
        wbd = wbd(:);

        % also read the "short" list (known words from puzzles)
        shlist = splittocellc(lower(asciiread( ...
            [neuroelf_path('contrib') '/wordbrain_list.txt'])), ...
            char(10), true, true);
        shlist = shlist(:);
        shlist(~cellfun('isempty', regexp(shlist, '(\?|\d)'))) = [];

        % only keep 2- to 10-letter words
        wbd(cellfun('prodofsize', wbd) < 2 | cellfun('prodofsize', wbd) > 10) = [];
        shlist(cellfun('prodofsize', shlist) < 2) = [];

        % remove illegal character words
        wbd(~cellfun('isempty', regexpi(wbd, '\-'))) = [];

        % and remove numbers (used as markers)
        wbd = regexprep(wbd, '[0-9]+$', '');

        % patch dictionary with missing entries
        wbd = [wbd; shlist];

        % and order again
        wbd = unique(wbd);
        wln = cellfun('prodofsize', wbd);

        % for each length
        for lc = 2:10

            % add to dictionary
            wbdict.(['f' char(95 + lc)]) = gendict(char(wbd(wln == lc)));
        end

        % at this point the dictionary is a struct such that for the
        % word list car, carpark, cart, etc.
        % wbdict.c.a.r.X, wbdict.c.a.r.p.a.r.k.X, wbdict.c.a.r.t.X, etc.
        % such that .X is a field meaning a word terminated at that letter

        % and also create shorter dictionary, ensuring it has all letters
        shscore = struct;
        for l1 = 1:numel(shlist)
            if ~isfield(shscore, shlist{l1})
                shscore.(shlist{l1}) = 1;
            else
                shscore.(shlist{l1}) = shscore.(shlist{l1}) + 1;
            end
        end
        shlist = unique(shlist(:));
        shlen = cellfun('prodofsize', shlist);
        for lc = 2:max(shlen)
            wbdict.(['s' char(95 + lc)]) = gendict(char(shlist(shlen == lc)));
        end

        % store as short list
        wbdict.score = shscore;
        
    catch ne_eo;
        wbdict(:) = [];
        error('neuroelf:dependency:dictFileError', 'Error processing dictionary file.');
    end
end

% argument check
onarg = nargin;
if nargin < 2 || isempty(puzzle) || (~ischar(puzzle) && ~iscell(puzzle)) || ...
   ~isa(sizes, 'double') || isempty(sizes) || any(isinf(sizes(:)) | isnan(sizes(:)) | sizes(:) < 2)
    if nargin == 1 && ischar(puzzle) && strcmpi(puzzle(:)', 'ask')
        puzzle = input('Puzzle: ', 's');
        psz = numel(puzzle);
        if ~any(psz == ((2:10) .^ 2))
            error('neuroelf:input:badInput', 'Invalid puzzle string.');
        end
        sizes = input('Sizes: ', 's');
        if any(sizes ~= ' ' & (sizes < '2' | sizes > '9'))
            error('neuroelf:input:badInput', 'Invalid sizes.');
        end
        sizes = u8str2double(sizes);
        if sum(sizes) ~= psz
            error('neuroelf:input:badInput', 'Invalid sizes (sum).');
        end
    elseif nargin > 0 && ~isempty(puzzle) && any('0123456789' == puzzle(end))
        if onarg > 2
            varargin = [{opts}, varargin];
        end
        if onarg > 1
            opts = sizes;
        end
        onarg = onarg + 1;
        puzzle = lower(puzzle(:)');
        sizes = double(regexprep(puzzle, '^.*[a-z](\d+)$', '$1')) - 48;
        if ~any(sum(sizes) == ((2:10) .^ 2))
            error('neuroelf:general:badArgument', 'Bad or missing argument');
        end
        puzzle = puzzle(puzzle >= 'a' & puzzle <= 'z');
        if sum(sizes) ~= numel(puzzle)
            error('neuroelf:general:badArgument', 'Bad or missing argument');
        end
    else
        error('neuroelf:general:badArgument', 'Bad or missing argument');
    end
end
if onarg < 3 || ~isstruct(opts) || numel(opts) ~= 1
    if onarg > 2 && iscell(opts)
        opts = struct('knowns', false, 'matches', {opts(:)'}, ...
            'simplify', false, 'smart', false, 'timeout', 60);
        if onarg > 3 && ischar(varargin{1}) && strcmpi(varargin{1}(:)', 'knowns')
            opts.knowns = true;
        end
    elseif onarg > 2 && ischar(opts) && strcmpi(opts(:)', 'knowns')
        opts = struct('knowns', true, 'matches', {{}}, 'simplify', false, ...
            'smart', false, 'timeout', 120);
    else
        opts = struct;
    end
end
if ~isfield(opts, 'knowns') || ~islogical(opts.knowns) || numel(opts.knowns) ~= 1
    opts.knowns = false;
end
if ~isfield(opts, 'matches') || ~iscell(opts.matches) || isempty(opts.matches)
    opts.matches = {};
else
    opts.matches = opts.matches(:);
end
if ~isfield(opts, 'numwords') || ~isa(opts.numwords, 'double') || numel(opts.numwords) ~= 1 || ...
    isinf(opts.numwords) || isnan(opts.numwords) || opts.numwords < 1 || opts.numwords > numel(sizes)
    opts.numwords = numel(sizes);
else
    opts.numwords = floor(opts.numwords);
end
if numel(opts.matches) < numel(sizes)
    opts.matches(end+1:numel(sizes)) = {''};
end
if ~isfield(opts, 'simplify') || ~islogical(opts.simplify) || numel(opts.simplify) ~= 1
    opts.simplify = false;
end
if ~isfield(opts, 'smart') || ~islogical(opts.smart) || numel(opts.smart) ~= 1
    opts.smart = true;
end
if ~isfield(opts, 'timeout') || ~isa(opts.timeout, 'double') || numel(opts.timeout) ~= 1 || ...
    isinf(opts.timeout) || isnan(opts.timeout) || opts.timeout <= 0
    opts.timeout = Inf;
else
    opts.timeout = now + opts.timeout / 86400;
end

% parse puzzle if necessary
if iscell(puzzle)
    psz = numel(puzzle{1});
    if numel(puzzle) ~= psz || any(cellfun('prodofsize', puzzle) ~= psz)
        error('neuroelf:general:badArgument', 'Bad puzzle cell argument.');
    end
    for pc = 1:psz
        puzzle{pc} = puzzle{pc}(:)';
    end
    puzzle = cat(1, puzzle{:});

% check size
else
    if size(puzzle, 1) ~= size(puzzle, 2)
        if any(numel(puzzle) == ((1:10) .^ 2))
            puzzle = reshape(puzzle, round(sqrt(numel(puzzle))) .* [1, 1])';
        else
            error('neuroelf:general:badArgument', 'Bad puzzle char argument.');
        end
    end
    psz = size(puzzle, 1);
end

% make sure it's all lower-case letters
puzzle = lower(puzzle);

% all valid characters
if any(puzzle(:) ~= ' ' & puzzle(:) ~= '-' & (puzzle(:) < 'a' | puzzle(:) > 'z'))
    error('neuroelf:general:badArgument', 'Invalid letter(s) in puzzle.');
end

% check total sizes (of remaining letters)
sizes = round(sizes(:));
if sum(sizes) ~= sum(puzzle(:) >= 'a')
    error('neuroelf:general:badArgument', 'Cannot solve puzzle with given sizes.');
end
nsz = numel(sizes);
solved = [];

% smart mode
if opts.smart
    
    % try to solve with known words first
    oknowns = opts.knowns;
    opts.knowns = true;
    opts.numwords = nsz;
    opts.simplify = false;
    opts.smart = false;
    opts.timeout = 15;
    try
        [words, solved] = solvewordbrain(puzzle, sizes, opts);
    catch ne_eo;
        rethrow(ne_eo);
    end
    wf = fieldnames(words);
    swf = regexprep(wf, '_.*$', '');
    
    % solved (and simplified)
    if numel(wf) == nsz && isa(words.(wf{1}), 'double') && ~isempty(words.(wf{1}))
        return;
        
    % still solved? (test further)
    elseif numel(unique(swf)) == 1
        twf = wf;
        uwords = words;
        twords = cell(nsz, 1);
        for wc = 1:nsz
            if wc < nsz
                swords = {};
                for swc = 1:numel(twf)
                    if ~isstruct(uwords.(twf{swc}))
                        swords = {};
                        break;
                    end
                    swords = unique(cat(1, swords, regexprep(fieldnames(uwords.(twf{swc})), '_.*$', '')));
                end
                if numel(swords) ~= 1
                    break;
                end
            end
            twords(wc) = swf(1);
            if wc < nsz
                uwords = uwords.(twf{1});
                twf = fieldnames(uwords);
                swf = regexprep(twf, '_.*$', '');
            end
        end

        % all filled in
        if ~any(cellfun('isempty', twords))
            try
                [words, solved] = solvewordbrain(puzzle, sizes, twords(:)');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            return;
        end
    end
    
    % we didn't find a smart solution, simply report what we can get easily
    opts.knowns = oknowns;
    opts.matches = {};
    opts.numwords = 1;
    opts.simplify = true;
    opts.timeout = 2;
    try
        if nsz == 1
            words = solvewordbrain(puzzle, sizes, opts);
        else
            uwords = struct;
            [psizes, px, py] = pairs(sizes);
            psizes = permute(psizes, [1, 3, 2]);
            psizes = [psizes; psizes(:, [2, 1])];
            pxy = [px; py];
            py = [py; px];
            px = pxy;
            [psizes, pu] = unique(psizes, 'rows');
            px = px(pu);
            py = py(pu);
            for sc = 1:numel(px)
                newsizes = [psizes(sc, :)'; sizes(setdiff(1:nsz, [px(sc), py(sc)]), 1)];
                words = solvewordbrain(puzzle, newsizes, opts);
                wf = fieldnames(words);
                for wc = 1:numel(wf)
                    if isstruct(words.(wf{wc}))
                        swords = words.(wf{wc});
                        swf = fieldnames(swords);
                        for swc = 1:numel(swf)
                            uwords.(wf{wc}).(swf{swc}) = swords.(swf{swc});
                        end
                    elseif ~isfield(uwords, wf{wc})
                        uwords.(wf{wc}) = [];
                    end
                end
            end
            wf = sort(fieldnames(uwords));
            words = struct;
            for wc = 1:numel(wf)
                words.(wf{wc}) = uwords.(wf{wc});
            end
        end

        % sort
        wf = fieldnames(words);
        wfs = zeros(numel(wf), 1);
        for wc = 1:numel(wf)
            if isfield(wbdict.score, wf{wc})
                wfs(wc) = wbdict.score.(wf{wc});
            end
        end
        [wfs, wfsi] = sort(wfs);
        wf = wf(wfsi);
        uwords = struct;
        for wc = 1:numel(wf)
            uwords.(wf{wc}) = words.(wf{wc});
        end
        words = uwords;
        return;
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% try to solve with given sizes
try
    words = wbsolve(psz, puzzle, nsz, sizes, wbdict, opts);
catch ne_eo;
    rethrow(ne_eo);
end

% nothing returned
if isempty(fieldnames(words))
    return;
end

% try to see if a unique solution came out
uwords = struct;
twords = words;
scf = false;
for sc = 1:nsz
    uw = fieldnames(twords);
    uwt = unique(regexprep(uw, '_.*$', ''));
    if numel(uwt) == 1 && numel(uw) > 1
        uwx = cell(numel(uw), 1);
        for wc = 1:numel(uw)
            uwv = strrep(strrep(strrep(uw{wc}, [uwt{1} '_'], ''), 'r', ''), '_', ',');
            uwv = reshape(u8str2double(uwv), 2, numel(uwt{1}));
            uwx{wc} = sort(psz * uwv(1, :) + uwv(2, :));
        end
        if ~any(any(diff(cat(1, uwx{:})) ~= 0))
            uw = uw{1};
        end
    end
    if numel(uw) == 1
        uws = regexprep(uw{1}, '_.*$', '');
        uwv = strrep(strrep(strrep(uw{1}, [uws '_'], ''), 'r', ''), '_', ',');
        twords = twords.(uw{1});
        if ~isstruct(twords) && sc < nsz
            break;
        end
        uwv = reshape(u8str2double(uwv), 2, numel(uws));
        uwords.(uws) = uwv;
        if sc == nsz
            scf = true;
        end
    else
        break;
    end
end
if scf
    if nargout > 1
        solved = gensolved(psz, puzzle, nsz, sizes, words);
    end
    words = uwords;
    return;
end

% sort fields
wf = sort(fieldnames(words));
uwords = struct;
for fc = 1:numel(wf)
    uwords.(wf{fc}) = words.(wf{fc});
end
words = uwords;

% simplify
if opts.simplify
    wf = fieldnames(words);
    swf = regexprep(wf, '_.*$', '');
    uwords = struct;
    for fc = 1:numel(wf)
        wsc = words.(wf{fc});
        if isstruct(wsc)
            wfs = fieldnames(wsc);
            for sfc = 1:numel(wfs)
                uwords.(swf{fc}).(wfs{sfc}) = wsc.(wfs{sfc});
            end
        elseif ~isfield(uwords, swf{fc})
            uwords.(swf{fc}) = [];
        end
    end
    if nargout > 1
        solved = gensolved(psz, puzzle, nsz, sizes, words);
    end
    words = uwords;
elseif nargout > 1
    solved = gensolved(psz, puzzle, nsz, sizes, words);
end



% internal function
function w = wbsolve(psz, p, nsz, s, d, opts)

% generate output
w = struct;

% inputs
aw = opts.numwords > 0;
opts.numwords = opts.numwords - 1;
sw = opts.matches{1};
opts.matches(1) = [];

% word candidates with letters to remove
if opts.knowns
    [wc, wrm] = wbwords(psz, p, s(1), d.(['s' char(95 + s(1))]), sw);
else
    [wc, wrm] = wbwords(psz, p, s(1), d.(['f' char(95 + s(1))]), sw);
end

% sort by known status (always do high scoring first, in case of timeout)
if ~isinf(opts.timeout) && numel(wc) > 1
    wcs = zeros(numel(wc), 1);
    for cc = 1:numel(wc)
        if isfield(d.score, wc{cc})
            wcs(cc) = d.score.(wc{cc});
        end
    end
    if any(wcs > 0)
        [wcs, wcsi] = sort(wcs, 'descend');
        wc = wc(wcsi);
        wrm = wrm(wcsi, :);
    end
end

% solve entire puzzle
if aw && nsz > 1 && ~isempty(wc) && (isempty(sw) || any(~cellfun('isempty', regexpi(wc, sw))))

    % only candidate?
    if ~isempty(sw)
        wch = ~cellfun('isempty', regexpi(wc, sw));
        wc = wc(wch);
        wrm = wrm(wch, :);
    end

    % for each candidate
    for cc = 1:numel(wc)

        % solve updated puzzle if requested and time available
        if opts.timeout > now && opts.numwords > 0
            up = uppuzzle(psz, p, wrm(cc, :));
            swc = wbsolve(psz, up, nsz - 1, s(2:end), d, opts);

        % otherwise stop (empty return)
        else
            swc = [];
        end

        % create path string
        y = 1 + mod(wrm(cc, :) - 1, psz);
        x = 1 + round((wrm(cc, :) - y) / psz);
        xy = [x(:)'; y(:)'];
        rms = sprintf('r%d_%d_', xy(:));
        rms(end) = [];

        % solved?
        if opts.numwords == 0 || (isstruct(swc) && ~isempty(fieldnames(swc)))

            % store information
            w.([wc{cc} '_' rms]) = swc;
        elseif opts.timeout <= now
            w.([wc{cc} '_' rms]) = 'timeout';
        end
    end

% only look for single words
elseif (~aw || nsz == 1) && ~isempty(wc)

    % only candidate?
    if ~isempty(sw)
        wch = ~cellfun('isempty', regexpi(wc, sw));
        wc = wc(wch);
        wrm = wrm(wch, :);
    end

    % store information
    for cc = 1:numel(wc)
        y = 1 + mod(wrm(cc, :) - 1, psz);
        x = 1 + round((wrm(cc, :) - y) / psz);
        xy = [x(:)'; y(:)'];
        rms = sprintf('r%d_%d_', xy(:));
        rms(end) = [];
        w.([wc{cc} '_' rms]) = [];
    end
end



% word candidates from position
function [wc, wrm] = wbwords(psz, p, sz, d, sw)

% initialize arrays
psz2 = psz * psz;
wc = repmat(' ', psz2, sz);
wrm = zeros(psz2, sz);
wf = 0;
swl = '';
if numel(sw) ~= sz
    sw = '';
else
    swl = sw(1);
    sw(1) = [];
end

% begin looking for words
for p1 = 1:psz
    for p2 = 1:psz

        % continue if nothing at this position or letter not in dictionary
        pl = p(p1, p2);
        if pl < 'a' || (~isempty(swl) && pl ~= swl) || ~isfield(d, pl)
            continue;
        end

        % get all letter combinations from this position
        [lcs, swrm] = candidates(psz, p, sz, d.(pl), p1, p2, sw);

        % add to list
        if ~isempty(lcs)
            nlcs = size(lcs, 1);
            wc(wf+1:wf+nlcs, :) = lcs;
            wrm(wf+1:wf+nlcs, :) = swrm;
            wf = wf + nlcs;
        end
    end
end

if wf < psz2
    wc(wf+1:end, :) = [];
    wrm(wf+1:end, :) = [];
end
if wf > 0
    wc = cellstr(wc);

    % simplify here if all letters the same
    if numel(wc) > 1 && numel(unique(wc)) == 1
        wrmt = sort(wrm, 2);
        if ~any(any(diff(wrmt)))
            wc = wc(1);
            wrm = wrm(1, :);
        end
    end
else
    wc = {};
end



% letter candidates
function [lcs, lcp] = candidates(psz, p, sz, d, p1, p2, sw)

% last letter
if sz == 1
    lcs = p(p1, p2);
    lcp = psz * (p2 - 1) + p1;
    return;
end

% generate list and fill
psz2 = psz * psz;
lcs = repmat(' ', psz2, sz);
lcp = zeros(psz2, sz);
ti = 0;
cp = p(p1, p2);
p(p1, p2) = ' ';
if ~isempty(sw)
    swl = sw(1);
    sw(1) = [];
else
    swl = '';
end

% sub lists
if isempty(swl)
    if p1 > 1
        if p2 > 1
            pl = p(p1 - 1, p2 - 1);
            if pl >= 'a' && isfield(d, pl)
                [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 - 1, p2 - 1, sw);
                if ~isempty(slcs)
                    nsl = size(slcs, 1);
                    lcs(ti+1:ti+nsl, 2:end) = slcs;
                    lcp(ti+1:ti+nsl, 2:end) = slcp;
                    ti = ti + nsl;
                end
            end
        end
        pl = p(p1 - 1, p2);
        if pl >= 'a' && isfield(d, pl)
            [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 - 1, p2, sw);
            if ~isempty(slcs)
                nsl = size(slcs, 1);
                lcs(ti+1:ti+nsl, 2:end) = slcs;
                lcp(ti+1:ti+nsl, 2:end) = slcp;
                ti = ti + nsl;
            end
        end
        if p2 < psz
            pl = p(p1 - 1, p2 + 1);
            if pl >= 'a' && isfield(d, pl)
                [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 - 1, p2 + 1, sw);
                if ~isempty(slcs)
                    nsl = size(slcs, 1);
                    lcs(ti+1:ti+nsl, 2:end) = slcs;
                    lcp(ti+1:ti+nsl, 2:end) = slcp;
                    ti = ti + nsl;
                end
            end
        end
    end
    if p2 > 1
        pl = p(p1, p2 - 1);
        if pl >= 'a' && isfield(d, pl)
            [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1, p2 - 1, sw);
            if ~isempty(slcs)
                nsl = size(slcs, 1);
                lcs(ti+1:ti+nsl, 2:end) = slcs;
                lcp(ti+1:ti+nsl, 2:end) = slcp;
                ti = ti + nsl;
            end
        end
    end
    if p2 < psz
        pl = p(p1, p2 + 1);
        if pl >= 'a' && isfield(d, pl)
            [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1, p2 + 1, sw);
            if ~isempty(slcs)
                nsl = size(slcs, 1);
                lcs(ti+1:ti+nsl, 2:end) = slcs;
                lcp(ti+1:ti+nsl, 2:end) = slcp;
                ti = ti + nsl;
            end
        end
    end
    if p1 < psz
        if p2 > 1
            pl = p(p1 + 1, p2 - 1);
            if pl >= 'a' && isfield(d, pl)
                [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 + 1, p2 - 1, sw);
                if ~isempty(slcs)
                    nsl = size(slcs, 1);
                    lcs(ti+1:ti+nsl, 2:end) = slcs;
                    lcp(ti+1:ti+nsl, 2:end) = slcp;
                    ti = ti + nsl;
                end
            end
        end
        pl = p(p1 + 1, p2);
        if pl >= 'a' && isfield(d, pl)
            [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 + 1, p2, sw);
            if ~isempty(slcs)
                nsl = size(slcs, 1);
                lcs(ti+1:ti+nsl, 2:end) = slcs;
                lcp(ti+1:ti+nsl, 2:end) = slcp;
                ti = ti + nsl;
            end
        end
        if p2 < psz
            pl = p(p1 + 1, p2 + 1);
            if pl >= 'a' && isfield(d, pl)
                [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 + 1, p2 + 1, sw);
                if ~isempty(slcs)
                    nsl = size(slcs, 1);
                    lcs(ti+1:ti+nsl, 2:end) = slcs;
                    lcp(ti+1:ti+nsl, 2:end) = slcp;
                    ti = ti + nsl;
                end
            end
        end
    end
else
    if p1 > 1
        if p2 > 1
            pl = p(p1 - 1, p2 - 1);
            if pl == swl && isfield(d, pl)
                [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 - 1, p2 - 1, sw);
                if ~isempty(slcs)
                    nsl = size(slcs, 1);
                    lcs(ti+1:ti+nsl, 2:end) = slcs;
                    lcp(ti+1:ti+nsl, 2:end) = slcp;
                    ti = ti + nsl;
                end
            end
        end
        pl = p(p1 - 1, p2);
        if pl == swl && isfield(d, pl)
            [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 - 1, p2, sw);
            if ~isempty(slcs)
                nsl = size(slcs, 1);
                lcs(ti+1:ti+nsl, 2:end) = slcs;
                lcp(ti+1:ti+nsl, 2:end) = slcp;
                ti = ti + nsl;
            end
        end
        if p2 < psz
            pl = p(p1 - 1, p2 + 1);
            if pl == swl && isfield(d, pl)
                [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 - 1, p2 + 1, sw);
                if ~isempty(slcs)
                    nsl = size(slcs, 1);
                    lcs(ti+1:ti+nsl, 2:end) = slcs;
                    lcp(ti+1:ti+nsl, 2:end) = slcp;
                    ti = ti + nsl;
                end
            end
        end
    end
    if p2 > 1
        pl = p(p1, p2 - 1);
        if pl == swl && isfield(d, pl)
            [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1, p2 - 1, sw);
            if ~isempty(slcs)
                nsl = size(slcs, 1);
                lcs(ti+1:ti+nsl, 2:end) = slcs;
                lcp(ti+1:ti+nsl, 2:end) = slcp;
                ti = ti + nsl;
            end
        end
    end
    if p2 < psz
        pl = p(p1, p2 + 1);
        if pl == swl && isfield(d, pl)
            [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1, p2 + 1, sw);
            if ~isempty(slcs)
                nsl = size(slcs, 1);
                lcs(ti+1:ti+nsl, 2:end) = slcs;
                lcp(ti+1:ti+nsl, 2:end) = slcp;
                ti = ti + nsl;
            end
        end
    end
    if p1 < psz
        if p2 > 1
            pl = p(p1 + 1, p2 - 1);
            if pl == swl && isfield(d, pl)
                [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 + 1, p2 - 1, sw);
                if ~isempty(slcs)
                    nsl = size(slcs, 1);
                    lcs(ti+1:ti+nsl, 2:end) = slcs;
                    lcp(ti+1:ti+nsl, 2:end) = slcp;
                    ti = ti + nsl;
                end
            end
        end
        pl = p(p1 + 1, p2);
        if pl == swl && isfield(d, pl)
            [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 + 1, p2, sw);
            if ~isempty(slcs)
                nsl = size(slcs, 1);
                lcs(ti+1:ti+nsl, 2:end) = slcs;
                lcp(ti+1:ti+nsl, 2:end) = slcp;
                ti = ti + nsl;
            end
        end
        if p2 < psz
            pl = p(p1 + 1, p2 + 1);
            if pl == swl && isfield(d, pl)
                [slcs, slcp] = candidates(psz, p, sz - 1, d.(pl), p1 + 1, p2 + 1, sw);
                if ~isempty(slcs)
                    nsl = size(slcs, 1);
                    lcs(ti+1:ti+nsl, 2:end) = slcs;
                    lcp(ti+1:ti+nsl, 2:end) = slcp;
                    ti = ti + nsl;
                end
            end
        end
    end
end

% remove all other elements
if ti < psz2
    lcs(ti+1:end, :) = [];
    lcp(ti+1:end, :) = [];
end

% set first letter and removal element
lcs(:, 1) = cp;
lcp(:, 1) = psz * (p2 - 1) + p1;



% update puzzle
function up = uppuzzle(psz, up, r)

% remove letters (replace by ' ')
up(r) = ' ';
for cc = 1:psz
    col = up(psz:-1:1, cc);
    if any(col == ' ')
        col(col == ' ') = [];
        col(end+1:psz) = ' ';
        up(psz:-1:1, cc) = col;
    end
end


% generate sub-dictionary
function d = gendict(wlist)

% nothing else to do
if isempty(wlist)
    d = [];
    return;

% single word left
elseif size(wlist, 1) == 1
    d.(wlist(1)) = gendict(wlist(1, 2:end));
    return;
end

% work
lb = 1 + [0; find(diff(wlist(:, 1)))];
le = [lb(2:end) - 1; size(wlist, 1)];
d = struct;
for lc = 1:numel(lb)
    d.(wlist(lb(lc), 1)) = gendict(wlist(lb(lc):le(lc), 2:end));
end


% generate solution
function solved = gensolved(psz, puzzle, nsz, sizes, words)

% corner shape
corner = repmat(round(255 .* ...
    [.0, .0, .0, .1, .4; ...
     .0, .0, .3, .6, 1 ; ...
     .0, .3, .9, 1 , 1 ; ...
     .1, .6, 1 , 1 , 1 ; ...
     .4, 1 , 1 , 1 , 1 ]), [1, 1, 3]);
ncrn = size(corner, 1);
lsize = 100;
lsmls = lsize - 4;

% colors
if nsz <= 5
    hues = 0:(1/nsz):(1 + 0.001 - 1/nsz);
    cols = hsvconv([hues(:), 0.5 * ones(nsz, 1), ones(nsz, 1)], 1);
else
    hues = 0.2:(1.5/nsz):(1.7 + 0.001 - 1.5/nsz);
    hues = hues - floor(hues);
    sats = min(1, 1.5 - 0.1 .* (1:nsz)');
    vals = min(1, 0.5 + 0.1 .* (1:nsz)');
    cols = hsvconv([hues(:), sats, vals], 1);
end
cols = repmat(reshape(cols', [1, 1, 3, nsz]), lsize, lsize);
colc = 1;

% create image of puzzle
dpuzzle = double(puzzle);
upuzzle = unique(dpuzzle);
upuzzle(end+1:end+9) = (49:57)';
lpuzzle = cell(1, 255);
for uc = 1:numel(upuzzle)
    lnum = upuzzle(uc);
    letter = image_font(upper(char(lnum)), 'Calibri', 72);
    if lnum < 64
        letter = image_resize(letter, 0.4);
        letter(end+1:lsize-12, :, :) = 255;
        letter(:, end+1:lsize-12, :) = 255;
    end
    slet = size(letter);
    letbox = 255 .* ones(lsmls, lsmls, 3);
    blet = 1 + floor(0.5 * (lsmls - slet(1:2)));
    letbox(blet(1):blet(1)+slet(1)-1, blet(2):blet(2)+slet(2)-1, :) = letter;
    letbox(1:ncrn, 1:ncrn, :) = corner;
    letbox(1:ncrn, end+1-ncrn:end, :) = corner(:, end:-1:1, :);
    letbox(end+1-ncrn:end, 1:ncrn, :) = corner(end:-1:1, :, :);
    letbox(end+1-ncrn:end, end+1-ncrn:end, :) = corner(end:-1:1, end:-1:1, :);
    letbox = cat(2, zeros(lsize, 2, 3), ...
        cat(1, zeros(2, lsmls, 3), letbox, zeros(2, lsmls, 3)), zeros(lsize, 2, 3));
    lpuzzle{upuzzle(uc)} = uint8(letbox);
end

% figure out how to pack words below the image
nrows = 1;
rowlc = zeros(6, 1);
for wc = 1:nsz
    if (rowlc(nrows) + sizes(wc) + 1) > (2.5 * psz)
        rowlc(nrows) = rowlc(nrows) - 1;
        nrows = nrows + 1;
    end
    rowlc(nrows) = rowlc(nrows) + sizes(wc) + 1;
end
rowlc(nrows) = rowlc(nrows) - 1;

% create puzzle image
solved = uint8(0);
solved(psz * lsize + ceil((nrows + 1) * (lsize * 0.5)), psz * lsize, 3) = 0;
solvedwr = ceil(1 + psz * lsize + 0.4 * lsize);
for cc = 1:psz
    for rc = 1:psz
        solved((rc-1)*lsize+1:rc*lsize, (cc-1)*lsize+1:cc*lsize, :) = ...
            lpuzzle{dpuzzle(rc, cc)};
    end
end

% to color it, begin by creating a char version with row and column indices
row = char((1:psz)' * ones(1, psz) + 48);

% now begin coloring until no more words are found
wf = fieldnames(words);
tw = '';
for wc = 1:numel(wf)
    if colc == nsz || isstruct(words.(wf{wc}))
        tw = wf{wc};
        break;
    end
end
rowc = 1;
rowis = cell(1, rowlc(rowc));
wci = 1;
while ~isempty(tw)

    % extract path
    p = u8str2double(strrep(strrep(regexprep(tw, '^[^_]+', ''), '_r', ','), '_', ','));
    p = reshape(p, 2, round(0.5 * numel(p)));

    % color solved image
    for cc = 1:size(p, 2)

        % get the column and row of position at time of solving
        lc = p(1, cc);
        lr = p(2, cc);

        % update row to original row
        lr = double(row(lr, lc)) - 48;

        % color letter
        limage = solved((lr-1)*lsize+1:lr*lsize, (lc-1)*lsize+1:lc*lsize, :);
        solved((lr-1)*lsize+1:lr*lsize, (lc-1)*lsize+1:lc*lsize, :) = ...
            min(min(cols(:, :, :, colc), lpuzzle{48+cc}), limage);

        % add to row
        rowis{wci} = solved((lr-1)*lsize+1:lr*lsize, (lc-1)*lsize+1:lc*lsize, :);
        wci = wci + 1;
    end

    % add space
    if wci < rowlc(rowc)
        rowis{wci} = uint8(zeros(lsize, ceil(0.8 * lsize), 3));
        wci = wci + 1;

    % place image
    else
        rowletters = image_resize(cat(2, rowis{:}), 0.375);
        rowlsz = size(rowletters);
        solvedwc = floor(1 + 0.5 * (psz * lsize - rowlsz(2)));
        solved(solvedwr:solvedwr+rowlsz(1)-1, solvedwc:solvedwc+rowlsz(2)-1, :) = rowletters;
        solvedwr = solvedwr + ceil(0.45 * lsize);
        rowc = rowc + 1;
        rowis = cell(1, rowlc(rowc));
        wci = 1;
    end

    % update the puzzle (and row representation)
    rl = psz .* (p(1, :) - 1) + p(2, :);
    puzzle = uppuzzle(psz, puzzle, rl);
    row = uppuzzle(psz, row, rl);

    % increase color counter
    colc = colc + 1;

    % get new words struct
    if colc <= nsz
        words = words.(tw);
    else
        words = struct;
    end

    % re-get next word
    wf = fieldnames(words);
    tw = '';
    for wc = 1:numel(wf)
        if colc == nsz || isstruct(words.(wf{wc}))
            tw = wf{wc};
            break;
        end
    end
end

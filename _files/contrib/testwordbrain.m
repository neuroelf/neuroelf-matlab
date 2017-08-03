% neuroelf library and solver functions
n = neuroelf;
ffirst = n.findfirst;
solver = n.solvewordbrain;

% load puzzles (not words)
puzzles = n.splittocellc(lower(n.asciiread( ...
    [neuroelf_path('contrib') '/wordbrain_list.txt'])), char(10), true, true);
puzzles = puzzles(:);

% remove empty rows and rows with letters + digits
puzzles(cellfun('isempty', puzzles)) = [];
puzzles(~cellfun('isempty', regexpi(puzzles, '^[a-z]+.*\d$'))) = [];

% progress bar
pbar = xprogress;
xprogress(pbar, 'settitle', 'Testing WordBrain solver...');
xprogress(pbar, 0, 'Testing puzzles...', 'visible', 0, numel(puzzles));

% start from the beginning
mt = 0;
pc = 1;
pn = 0;
ps = 0;
pa = 0;
pf = 0;
pe = 0;
while pc < numel(puzzles)

    % number of words in puzzle
    nw = numel(regexprep(puzzles{pc}, '.*[a-z]+(\d+)$', '$1'));
    pwords = puzzles(pc+1:pc+nw);
    pn = pn + 1;
    
    % try to solve puzzles with known words
    try
        psucc = 'failure';
        xprogress(pbar, pc);
        nt = now;
        w = solver(puzzles{pc}, struct('knowns', true, 'timeout', 30));
        wf = fieldnames(w);
        if numel(wf) == nw && isa(w.(wf{1}), 'double') && ~isempty(w.(wf{1})) && ...
            all(strcmp(wf, pwords))
            psucc = 'complete';
            nnt = now;
            if (nnt - nt) >= (2 / 86400)
                fprintf('%s: time to solve %.2fs\n', puzzles{pc}, 86400 * (nnt - nt));
            end
            ps = ps + 1;
        else
            twords = pwords;
            tw = w;
            twf = wf;
            twfi = strcmp(regexprep(twf, '_.*$', ''), twords{1});
            while ~isempty(twords) && any(twfi)
                twfi = ffirst(twfi);
                tw = w.(wf{twfi});
                twords(1) = [];
                if ~isempty(twords) && isstruct(tw)
                    twf = fieldnames(tw);
                end
            end
            if isempty(twords)
                psucc = 'complete';
                nnt = now;
                if (nnt - nt) >= (10 / 86400)
                    fprintf('%s: time to solve %.2fs\n', puzzles{pc}, 86400 * (nnt - nt));
                end
                ps = ps + 1;
            else
                w = solver(puzzles{pc}, puzzles(pc+1:pc+nw)');
                wf = fieldnames(w);
                if ~isempty(wf)
                    nnt = now;
                    if (nnt - nt) >= (2/86400)
                        fprintf('%s: time to almost-solve %.2fs\n', puzzles{pc}, 86400 * (nnt - nt));
                    end
                    psucc = 'almost';
                    pa = pa + 1;
                end
            end
        end
        nt = now - nt;
        if nt > mt
            mt = nt;
            mtp = puzzles{pc};
        end
    catch ne_eo;
        fprintf('%s: %s\n', puzzles{pc}, ne_eo.message);
        pe = pe + 1;
        pc = pc + nw + 1;
        continue;
    end
    if psucc(1) == 'f'
        fprintf('%s: %s\n', puzzles{pc}, psucc);
        pf = pf + 1;
    end

    % advance pointer
    pc = pc + nw + 1;
end

% stats
closebar(pbar);
fprintf('Success: %.2f%%\nAlmost:  %.2f%%\nFailure: %.2f%%\nErrors:  %.2f%%\n', ...
    (100 / pn) .* [ps, pa, pf, pe]);

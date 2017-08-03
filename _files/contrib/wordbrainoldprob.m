function [a, oldwords] = wordbrainoldprob(smk)
if nargin < 1 || ~isa(smk, 'double') || numel(smk) ~= 1
    smk = 25;
end

% neuroelf library and solver functions
n = neuroelf;

% load puzzle words (not puzzles or headers)
puzzles = n.splittocellc(lower(n.asciiread( ...
    [neuroelf_path('contrib') '/wordbrain_list.txt'])), char(10), true, true);
puzzles = puzzles(:);
puzzles(cellfun('isempty', puzzles)) = [];
puzzles(~cellfun('isempty', regexpi(puzzles, '\d'))) = [];

% now compute for each word whether it appeared already
oldword = zeros(numel(puzzles), 1);
oldwords = struct;
for wc = 1:numel(puzzles)
    if ~isfield(oldwords, puzzles{wc})
        oldwords.(puzzles{wc}) = wc;
    else
        oldword(wc) = 1;
    end
end

% and also get the lengths of all/new words
allwlen = cellfun('prodofsize', puzzles);
newwords = fieldnames(oldwords);
newwlen = cellfun('prodofsize', newwords);
newwpos = struct2cell(oldwords);
newwpos = cat(1, newwpos{:});
newwlps = NaN .* zeros(size(oldword));
newwlps(newwpos) = newwlen;
newwlpn = double(isnan(newwlps));

% smooth a bit
oldwordsm = n.flexinterpn(oldword, [inf;1;1;numel(oldword)], n.smoothkern(smk), 1);
allwlensm = 0.125 .* n.flexinterpn(allwlen, [inf;1;1;numel(allwlen)], n.smoothkern(smk), 1);
newwlpssm = 0.125 .* n.flexinterpn(newwlps, [inf;1;1;numel(newwlps)], n.smoothkern(smk), 1);
newwlpnsm = n.flexinterpn(newwlpn, [inf;1;1;numel(newwlpn)], n.smoothkern(smk), 1);
newwlpssm(newwlpnsm > 0.99) = NaN;

% plot
figure;
a = axes;
p = plot([oldwordsm, allwlensm, newwlpssm]);
set(p, 'LineWidth', 4);
lg = legend(a, 'P(oldword)', 'Len(all words) / 8', 'Len(new words) / 8', 'Location', 'SouthEast');
set(lg, 'FontSize', 16)

% add size breaks (2x2 -> 3x3 -> 4x4 -> 5x5 -> 6x6)
brs = [20.5, 180.5, 553.5, 1114.5, 2047.5, 2962.5];
hold(a, 'on');
for bc = 1:numel(brs)
    l(bc) = line(brs(bc) .* [1; 1], [0; 1]);
end
set(l, 'Color', [0,0,0]);

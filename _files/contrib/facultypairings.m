function [numids, authors, singleids, pairids, fig, ax] = facultypairings(filename, year)
%FACULTYPAIRINGS Create a matrix of PubMed IDs shared by pairs of faculty.
%   [NUMIDS, AUTHORS, SINGLEIDS, PAIRIDS] = FACULTY(FILENAME) calculates
%   the collaboration matrix of a list of faculty in file FILENAME. The
%   first row must contain the name of the institution.
%
%   The function searches PubMed, and publications not available in this
%   archive will not be part of the returned results!
%
%   [NUMIDS, AUTHORS, SINGLEIDS, PAIRIDS, FIG, AX] = FACULTY(FILENAME) also
%   returns the figure and axes handles.

% (c) 2019 Dmitriy

% color threshold
cthr = 10;

% test for second input
if nargin > 1 && ~isempty(year)
    if numel(year) == 1
        year = sprintf('%d[year]', year);
    else
        year = sprintf('%d[year] OR ', year);
        year = ['(' year(1:end-4) ')'];
    end
    year = [' AND ' year];
else
    year = '';
end

% load and process list of faculty
fid = fopen(filename);
iname = fgetl(fid);
t = '';
authors = {};
while ischar(t)
    t = fgetl(fid);
    if ischar(t)
        authors{end+1} = t; %#ok<AGROW>
    end
end
authors = sort(authors);
numauthors = numel(authors);
fclose(fid);

% prepare other outputs
singleids = cell(numauthors, 1);
pairids = Inf .* ones(numauthors + 1);
numids = cell(numauthors, numauthors);

% set up a figure
fig = figure;
ax = axes;
numplus = numauthors + 1;
img = surf(1:numplus, 1:numplus, pairids);
view([0, 90]);
strAbbr = cellfun(@(x) x(1:end-2), authors, 'uniformoutput', false);
set(gca, ...
    'Position', [.15 0.01 .84 0.84], ...
    'ytick', (1:numel(authors))+0.5, 'yticklabel', strAbbr, ...
    'xtick', (1:numel(authors))+0.5, 'xticklabel', strAbbr, ...
    'xticklabelrotation', 90, 'xaxislocation', 'top', ...
    'tickdir', 'out', 'ticklength', [0 0],'fontsize',12, 'ydir','reverse','clim',[0 cthr+1]);
set(gcf, ...
    'position', [40, 40, 800, 800], ...
    'colormap', [0, 0, 0; spring(cthr); 0, 1, 1], 'color', 'w');
axis tight;
axis equal;

% Loop through faculty
warning off
for a = 1:numel(authors)
    isgood = false;
    while ~isgood
        try % Sometimes URL reading fails; keep trying if it does
            % Search pubmed for Name1[au] AND Name2[au] AND Columbia
            %term = [authors{a} '[au] AND ' authors{b} '[au] AND ' iname];
            fprintf('Fetching results for %s...\n', authors{a});
            term = [authors{a} '[au] AND ' iname year];
            result = urlread(['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=200&term=' term]);
            [st, en] = regexp(result,'<Id>[0-9]+</Id>');
            singleids{a} = NaN .* zeros(numel(st), 1);
            for j = 1:numel(st)
                singleids{a}(j) = str2double(result(st(j)+4:en(j)-5));
            end
            pause(0.3);
            isgood = true;
        catch
            pause(5);
            drawnow;
        end
    end
end
for a = 1:numel(authors)
    for b = 1:a-1
        st = intersect(singleids{a}, singleids{b});
        pairids(a,b) = min(10, numel(st));
        pairids(b,a) = pairids(a,b);
        if ~isempty(st)
            numids{a,b} = st;
            numids{b,a} = numids{a,b};
        end
    end
end
mm = pairids;
mm(isinf(mm)) = cthr+1;
set(img,'zdata', mm);

% shorten array
pairids(end, :) = [];
pairids(:, end) = [];

% reset warning
warning on

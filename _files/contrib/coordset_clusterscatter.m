% code for meeting with Jocelyn, Ajay

% column codes
studycol = 3;
spacecol = 20;
xcol = 25;
ycol = 26;
zcol = 27;

% read in text-version of Excel file (
eac = acsvread( ...
    'EmotionalAwarenessCoordinates3.txt', ...
    char(9), ...
    struct('asmatrix', 1, 'convert', true, 'noquotes', 1));

% get coordinates in one Px3 matrix
eacc = cat(2, ...
    cat(1, eac{2:end, xcol}), ...
    cat(1, eac{2:end, ycol}), ...
    cat(1, eac{2:end, zcol}));

% convert TAL coordinates to MNI (using tal2icbm, spelling mistake...)
eacc(strcmpi(eac(2:end, spacecol), 'talairach'), :) = ...
    tal2icbm(eacc(strcmpi(eac(2:end, spacecol), 'talairach'), :));

% compute distances with new function
pd = clusterdist(eacc, eac(2:end, studycol));

% number of unique studies
ns = numel(unique(eac(2:end, studycol)));

% scatter accordingly
maxdist = 10 * ns * max(pd);
scsize = max(5, maxdist - (12 * ns) .* pd);
maxsize = max(scsize);
sccol = max(0.1, scsize ./ maxsize);
sccol = sccol * ([0, 0, 1]) + (1 - sccol) * ones(1, 3);
scatter3(eacc(:,1), eacc(:,2), eacc(:,3), scsize, sccol);

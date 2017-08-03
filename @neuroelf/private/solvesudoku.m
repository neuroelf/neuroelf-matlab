function [s, si] = solvesudoku(s, bfd)
% solvesudoku  - solve a sudoku puzzle
%
% FORMAT:       [s, sm] = solvesudoku(s [, bfd])
%
% Input fields:
%
%       s           either a 9x9 double array with 1...9 for known entries
%                   or a 81x9 logical array with true for potential entries
%       bfd         brute-force depth (number of guesses, default: 4)
%
% Output fields:
%
%       s           output as input with solution
%       sm          bit field (for pencil markings)
%
% Note: as the brute-force algorithm picks the *first* occurrence of a cell
%       that has the fewest remaining candidates, transposing or reordering
%       the puzzle *can* make a difference in run-time (but not solution)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% persistent memory
persistent ss_pm;
if isempty(ss_pm)

    % indices into rows, columns and 3x3 packets
    ss_pm.ri = {1:9:81, 2:9:81, 3:9:81, 4:9:81, 5:9:81, 6:9:81, 7:9:81, 8:9:81, 9:9:81};
    ss_pm.ci = {1:9, 10:18, 19:27, 28:36, 37:45, 46:54, 55:63, 64:72, 73:81};
    ss_pm.bi = {[1:3, 10:12, 19:21], [4:6, 13:15, 22:24], [7:9, 16:18, 25:27], ...
        [28:30, 37:39, 46:48], [31:33, 40:42, 49:51], [34:36, 43:45, 52:54], ...
        [55:57, 64:66, 73:75], [58:60, 67:69, 76:78], [61:63, 70:72, 79:81]};

    % combined array with all 27 "packets"
    ss_pm.ai = [ss_pm.ri, ss_pm.ci, ss_pm.bi];

    % reverse logic indices
    ss_pm.ii = cell(1, 81);
    for c = 1:27
        ci = ss_pm.ai{c};
        for cc = 1:9
            ii = true(1, 9);
            ii(cc) = false;
            ss_pm.ii{ci(cc)} = cat(2, ss_pm.ii{ci(cc)}, ci(ii));
        end
    end
    for c = 1:81
        ss_pm.ii{c} = unique(ss_pm.ii{c});
    end

    % combinations of three of up to 9 remaining cells
    ss_pm.t9 = cell(1, 84);
    c = 1;
    for c1 = 1:7
        for c2 = (c1 + 1):8
            for c3 = (c2 + 1):9
                ss_pm.t9{c} = [c1, c2, c3];
                c = c + 1;
            end
        end
    end

    % combinations of three of up to 8 remaining cells
    ss_pm.t8 = cell(1, 56);
    c = 1;
    for c1 = 1:6
        for c2 = (c1 + 1):7
            for c3 = (c2 + 1):8
                ss_pm.t8{c} = [c1, c2, c3];
                c = c + 1;
            end
        end
    end

    % combinations of three of up to 7 remaining cells
    ss_pm.t7 = cell(1, 35);
    c = 1;
    for c1 = 1:5
        for c2 = (c1 + 1):6
            for c3 = (c2 + 1):7
                ss_pm.t7{c} = [c1, c2, c3];
                c = c + 1;
            end
        end
    end

    % combinations of three of up to 6 remaining cells
    ss_pm.t6 = {[1, 2, 3], [1, 2, 4], [1, 2, 5], [1, 2, 6], [1, 3, 4], ...
                [1, 3, 5], [1, 3, 6], [1, 4, 5], [1, 4, 6], [1, 5, 6], ...
                [2, 3, 4], [2, 3, 5], [2, 3, 6], [2, 4, 5], [2, 4, 6], ...
                [2, 5, 6], [3, 4, 5], [3, 4, 6], [3, 5, 6], [4, 5, 6]};

    % combinations of three of up to 5 remaining cells
    ss_pm.t5 = {[1, 2, 3], [1, 2, 4], [1, 2, 5], [1, 3, 4], [1, 3, 5], ...
                [1, 4, 5], [2, 3, 4], [2, 3, 5], [2, 4, 5], [3, 4, 5]};
end

% argument check
if nargin < 1 || ...
   ((~islogical(s) || ...
     ~isequal(size(s), [81, 9])) && ...
    (~isa(s, 'double') || ...
     ~isequal(size(s), [9, 9])))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument supplied.' ...
    );
end
if nargin < 2 || ...
   ~isa(bfd, 'double') || ...
    numel(bfd) ~= 1 || ...
    isinf(bfd) || ...
    isnan(bfd)
    bfd = 4;
else
    bfd = min(64, max(0, round(bfd)));
end

% convert to logical array
if ~islogical(s)

    % keep track of input format
    sconv = true;

    % set invalid entries to 0
    s(isinf(s) | isnan(s)) = 0;

    % round
    s = round(s);

    % find indices between 1 and 9
    fs = find(s(:) > 0 & s(:) < 10);

    % create logical representation
    si = true(81, 9);

    % iterate over known numbers
    for c = 1:numel(fs)

        % get position
        fc = fs(c);

        % get number
        sc = s(fc);

        % set all but number in cell to false
        si(fc, :) = false;
        si(fc, sc) = true;

        % set all affected cells for that number to false
        si(ss_pm.ii{fc}, sc) = false;
    end

% already logical
else
    % keep track of input format
    sconv = false;
    si = s;
end

% initiate 1x9 logical indexing
un = true(1, 9);

% compute current score of determination
sf = sum(si(:));
sf2 = sf;

% go on while not solved
while sf > 81

    % display?
    %disp(sudokufrombit(si));

    % iterate over all cells
    for c = 1:81

        % wherever only one solution remains
        if sum(si(c, :)) == 1

            % remove from affected cells
            si(ss_pm.ii{c}, si(c, :)) = false;
        end
    end

    % if sum already lower
    sf2 = sum(si(:));
    if sf2 < sf

        % re-set sum
        sf = sf2;

        % and continue
        continue;
    end

    % check all 27 cells for places where only one is possible
    for c = 1:27

        % get indices of packet
        ti = ss_pm.ai{c};

        % find uniquely specified numbers
        ui = find(sum(si(ti, :)) == 1);

        % iterate over these
        for uc = 1:numel(ui)

            % get number
            ux = ui(uc);

            % prepare indexing
            un(:) = true;
            un(ux) = false;

            % set all other numbers to false
            ci = ti(si(ti, ux));
            si(ci, un) = false;

            % and also reflect in affected cells
            si(ss_pm.ii{ci}, ux) = false;
        end
    end

    % if sum already lower
    sf2 = sum(si(:));
    if sf2 < sf

        % re-set sum
        sf = sf2;

        % and continue
        continue;
    end

    % iterate over 27 packets again to check for pairs
    for c = 1:27

        % get indices of packet
        ti = ss_pm.ai{c};

        % iterate over first 8 indices
        for tc1 = 1:8

            % if sum of these is 2
            if sum(si(ti(tc1), :)) == 2

                % iterate over remaining indices
                for tc2 = (tc1 + 1):9

                    % if the set indices are the same
                    if all(si(ti(tc2), :) == si(ti(tc1), :))

                        % prepare logical indexing
                        un(:) = true;
                        un(tc1) = false;
                        un(tc2) = false;

                        % set to false in other boxes
                        si(ti(un), si(ti(tc1), :)) = false;
                    end
                end
            end
        end

        % and check for "hidden" pairs as well
        ui = si(ti, :);

        % iterate over numbers
        for tc1 = 1:8

            % if only in two boxes
            if sum(ui(:, tc1)) == 2

                % iterate over remaining numbers
                for tc2 = (tc1 + 1):9

                    % same two numbers
                    if all(ui(:, tc2) == ui(:, tc1))

                        % prepare logical indexing
                        un(:) = true;
                        un(tc1) = false;
                        un(tc2) = false;

                        % remove from remaining cells
                        si(ti(ui(:, tc1)), un) = false;
                        ui = si(ti, :);
                    end
                end
            end
        end
    end

    % if sum already lower
    sf2 = sum(si(:));
    if sf2 < sf

        % re-set sum
        sf = sf2;

        % and continue
        continue;
    end

    % iterate over 27 packets again to check for triplets
    for c = 1:27

        % get indices of packet and packet
        ti = ss_pm.ai{c};
        ui = si(ti, :);

        % check whether anything "good" is to be expected from this
        if sum(sum(ui, 2) == 1) > 4
            continue;
        end

        % get vector with cells to consider and combinations of these
        tc = find(sum(ui, 2) > 1);
        t0 = ss_pm.(sprintf('t%d', numel(tc)));

        % iterate over all combinations of cells
        for cc = 1:numel(t0)

            % get all candidates in triplet
            pci = any(ui(tc(t0{cc}), :), 1);

            % check whether fields share the same candidates
            if sum(pci) == 3

                % remove those from other cells
                un(:) = true;
                un(tc(t0{cc})) = false;
                si(ti(un), pci) = false;
                ui = si(ti, :);
            end
        end

        % also check for hidden triplets
        tc = find(sum(ui, 1) > 1);
        t0 = ss_pm.(sprintf('t%d', numel(tc)));

        % iterate over all combinations of numbers
        for cc = 1:numel(t0)

            % get all cells with those three numbers in triplet
            pci = any(ui(:, tc(t0{cc})), 2);

            % check whether these three appear only in three cells
            if sum(pci) == 3

                % remove other numbers from those cells
                un(:) = true;
                un(tc(t0{cc})) = false;
                si(ti(pci), un) = false;
                ui = si(ti, :);
            end
        end
    end

    % if sum lower
    sf2 = sum(si(:));
    if sf2 < sf

        % re-set sum
        sf = sf2;

        % and continue
        continue;
    end

    % make a soundness check
    if ~all(any(si, 2))
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid Sudoku grid.' ...
        );
    end

    % break
    break;
end

% make sure brute-force is disabled if solved
if sf2 < 82
    bfd = 0;
end

% brute force
if bfd > 0

    % find the first cell with two (or three) candidates
    bfc = findfirst(sum(si, 2) == 2);
    if isempty(bfc)
        bfc = findfirst(sum(si, 2) == 3);
        if isempty(bfc)
            bfc = findfirst(sum(si, 2) == 4);
        end
    end

    % get candidates
    cn = find(si(bfc, :));

    % make copy
    st = si;

    % iterate over candidates
    for cc = 1:numel(cn)

        % set to candidate
        st(bfc, :) = false;
        st(bfc, cn(cc)) = true;

        % try to solve
        try
            stn = solvesudoku(st, bfd - 1);
            if all(sum(stn, 2) == 1)
                si = stn;
                break;
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
end

% sanity check
if ~all(any(si, 2))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid Sudoku grid.' ...
    );
end

% convert back?
if sconv
    s = sudokufrombit(si);

% otherwise make sure to return bit-field
else
    s = si;
end

% function to convert bitfield back
function s = sudokufrombit(si)
    s = zeros(9, 9);
    for c = 1:81
        try
            s(c) = find(si(c, :));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
% end of function sudokufrombit

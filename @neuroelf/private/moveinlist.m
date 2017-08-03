function [tidx, moveto] = moveinlist(lsize, sel, moveto, mskip)
% moveinlist  - compute new indices after a subselection is moved
%
% FORMAT:       [tidx, moveto] = moveinlist(lsize, sel, moveto [, mskip])
%
% Input fields:
%
%       lsize       size of list (1x1 double, integer value >= 1)
%       sel         selection (indices to be moved)
%       moveto      relative move
%       mskip       boolean flag, skip over unselected at the end (false)
%
% Output fields:
%
%       tidx        target indices, so that newlist = oldlist(tidx)
%       moveto      moved-to indices (e.g. to update a listbox selection)

% Version:  v0.9b
% Build:    11041215
% Date:     Apr-12 2011, 3:05 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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

% input check
if nargin < 3 || ...
   ~isa(lsize, 'double') || ...
    numel(lsize) ~= 1 || ...
    isinf(lsize) || ...
    isnan(lsize) || ...
    lsize < 1 || ...
    lsize ~= fix(lsize) || ...
   ~isa(sel, 'double') || ...
   ~any(numel(sel) == size(sel)) || ...
    any(isinf(sel) | isnan(sel) | sel < 1 | sel > lsize | sel ~= fix(sel)) || ...
    numel(sel) ~= numel(unique(sel)) || ...
   ~isa(moveto, 'double') || ...
   ~any(numel(moveto) == [1, numel(sel)]) || ...
   ~any(numel(moveto) == size(moveto)) || ...
    any(isinf(moveto) | isnan(moveto) | moveto <= -lsize | moveto == 0 | moveto > lsize) || ...
    numel(moveto) ~= numel(unique(moveto)) || ...
   (numel(moveto) > 1 && ...
    any(moveto < 0))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input.' ...
    );
end
if nargin < 4 || ...
   ~islogical(mskip) || ...
    numel(mskip) ~= 1
    mskip = false;
end

% create lists for summing
sel = sel(:);
moveto = moveto(:);

% one moveto argument
if numel(moveto) == 1

    % sort selection
    sel = sort(sel);

    % target is relative to input
    moveto = sel + moveto;

    % skipping allowed
    if mskip

        % values below 1
        if any(moveto < 1)

            % push within limits
            moveto(moveto < 1) = moveto(moveto < 1) - (min(moveto) - 1);

            % for multiple entries
            if numel(moveto) > 1

                % ensure that no double indices occur
                dmoveto = diff(moveto);
                if any(dmoveto == 0)
                    moveto = moveto(1) + [0; cumsum(max(1, dmoveto))];
                end
            end

        % or greater lsize
        elseif any(moveto > lsize)

            % push within limits
            moveto(moveto > lsize) = moveto(moveto > lsize) - (max(moveto) - lsize);

            % for multiple entries
            if numel(moveto) > 1

                % ensure that no double indices occur
                dmoveto = diff((lsize + 1) - moveto(end:-1:1));
                if any(dmoveto == 0)
                    dmove0 = find(dmoveto == 0);
                    if dmove0(end) == numel(dmoveto)
                        dmove0(end) = [];
                    end
                    dmoveto(dmove0) = 1;
                    dmoveto(dmove0+1) = min(1, dmoveto(dmove0+1) - 1);
                    dmoveto = moveto(end) - [0; cumsum(dmoveto)];
                    moveto = dmoveto(end:-1:1);
                end
            end
        end

    % skipping not allowed
    else

        % values below 1
        if any(moveto < 1)

            % push all back the same amount
            moveto = moveto - (min(moveto) - 1);

        % or greater lsize
        elseif any(moveto > lsize)

            % push all back the same amount
            moveto = moveto - (max(moveto) - lsize);
        end
    end
end

% create lists
sidx = 1:lsize;
tidx = sidx;

% get indices for remaining items (from and to)
rindex = setdiff(sidx, sel');
oindex = setdiff(sidx, moveto');

% compile and convert new array, set value
tidx(oindex) = sidx(rindex);
tidx(moveto) = sidx(sel);

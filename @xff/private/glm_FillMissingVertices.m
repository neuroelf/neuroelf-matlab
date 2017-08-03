function xo = glm_FillMissingVertices(xo, srf, method)
% GLM::FillMissingVertices  - fill missing vertices in beta maps
%
% FORMAT:       [glm] = glm.FillMissingVertices(srf [, method])
%
% Input fields:
%
%       srf         surface object (needed for neighborhood information)
%       method      either of 'mean', {'nearest'}, 'dist-weighted'
%              where
%                   mean - interpolate missing by mean of valid neighbors
%                   dist - interpolate by an 1/distance-weighted mean
%                   nearest - set missing value to nearest valid neighbor
%
% Output fields:
%
%       glm         altered GLM structure

% Version:  v1.1
% Build:    16020314
% Date:     Feb-03 2016, 2:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || numel(srf) ~= 1 || ...
   ~xffisobject(xo, true, 'glm') || ~xffisobject(srf, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
srfc = srf.C;
if bc.ProjectType ~= 2 || bc.ProjectTypeRFX == 0 || bc.NrOfStudies < 2 || ...
    bc.NrOfVertices ~= srfc.NrOfVertices
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~ischar(method) || isempty(method) || ...
   ~any(strcmpi(method(:)', {'all', 'dist', 'dist-weighted', 'mean', 'mean-weighted'}))
    method = 0;
else
    switch (lower(method(:)'))
        %case {'all'}
        %    method = 3;
        case {'dist', 'dist-weighted'}
            method = 2;
        case {'mean', 'mean-weighted'}
            method = 1;
        otherwise
            method = 0;
    end
end

% get numbers, coordinates, and neighbors
numsubs = numel(bc.GLMData.Subject);
nummaps = size(bc.GLMData.Subject(1).BetaMaps, 2);
pcoords = srfc.VertexCoordinate;
neilist = srfc.Neighbors(:, 2);
myeps = eps;

% iterate over subjects
for sc = 1:numsubs

    % get BetaMaps
    filltry = 0;
    bmaps = bc.GLMData.Subject(sc).BetaMaps(:, :);

    % get indices to fill
    fillidx = find(all(bmaps == 0, 2));

    % while indices to fill
    while ~isempty(fillidx)

        % get neighborhood information on indices
        fillnei = neilist(fillidx);

        % for mean, don't use distance information
        if method == 1

            % iterate over indices
            for ic = 1:numel(fillidx)

                % get neighbors and their betas
                inei = fillnei{ic};
                bnei = bmaps(inei, :);

                % fill into missing
                try
                    knei = find(any(bnei ~= 0, 2));
                    if ~isempty(knei)
                        bmaps(fillidx(ic), :) = mean(bnei(knei, :), 1);
                    end
                catch xfferror
                    neuroelf_lasterr(xfferror);
                end
            end

        % either nearest or distance weighted
        else

            % iterate over indices
            for ic = 1:numel(fillidx)

                % get neighbors and their betas
                inei = fillnei{ic};
                bnei = bmaps(inei, :);

                % remove invalid entries
                rnei = find(all(bnei == 0, 2));
                if ~isempty(rnei)
                    inei(rnei) = [];
                    if isempty(inei)
                        continue;
                    end
                    bnei(rnei, :) = [];
                end

                % get coordinate positions
                neipos = pcoords(inei, :);
                mispos = repmat(pcoords(fillidx(ic), :), [size(neipos, 1), 1]);

                % calculate distance
                pdist = sqrt(sum((neipos - mispos) .^ 2, 2));

                % nearest neighbor
                if method == 0

                    % minimum position
                    [srcpos{1:2}] = min(pdist);

                    % replace beta values
                    bmaps(fillidx(ic), :) = bnei(srcpos{2}, :);

                % distance weighted
                else

                    % set betas to weighted mean
                    pdist = 1 ./ (pdist + myeps);
                    bmaps(fillidx(ic), :) = ...
                        sum(double(bnei) .* repmat(pdist, [1, nummaps]), 1) ./ sum(pdist);
                end
            end
        end

        % increase fill try
        filltry = filltry + 1;
        if filltry > 3
            warning('neuroelf:xff:badFileContent', ...
                'Too many beta values missing for subject %d.', sc);
            break;
        end

        % re-get fillidx
        fillidx = find(all(bmaps == 0, 2));
    end

    % set BetaMaps back
    if filltry > 0
        bc.GLMData.Subject(sc).BetaMaps = bmaps;
    end
end

% put back
xo.C = bc;

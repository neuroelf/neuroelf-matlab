function xo = srf_MarkPOI(xo, poi, npoi, colcode, mix)
% SRF::MarkPOI  - mark a patch-of-interest
%
% FORMAT:       [srf] = srf.MarkPOI(poi, npoi, colcode [, mix])
%
% Input fields:
%
%       poi         POI xff object
%       npoi        either a name (does nothing if not found) or number
%       colcode     1x3 RGB color code
%       mix         if given and between 0 .. 1, mix colcode with current
%
% Output fields:
%
%       srf         altered SRF object

% Version:  v1.1
% Build:    16021119
% Date:     Feb-11 2016, 7:47 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% check arguments
if nargin < 4 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
    numel(poi) ~= 1 || ~xffisobject(poi, true, 'poi') || ...
    isempty(npoi) || (~ischar(npoi) && ~isa(npoi, 'double')) || ~isa(colcode, 'double') || ...
   (numel(colcode) ~= 1 && numel(colcode) ~= 3) || ...
    any(isinf(colcode) | isnan(colcode) | colcode < 0 | colcode > 255)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
poic = poi.C;
if all(colcode <= 1) && any(colcode ~= fix(colcode))
    colcode = round(colcode(:)' * 255);
else
    colcode = round(colcode(:)');
end
if nargin < 5 || ~isa(mix, 'double') || numel(mix) ~= 1 || ...
    isinf(mix) || isnan(mix) || mix <= 0 || mix > 1
    mix = 1;
end

% check npoi
pois = poic.POI;
numpois = numel(pois);
if ischar(npoi)
    npoi = npoi(:)';
    for pc = 1:numpois
        if strcmpi(pois(pc).Name, npoi)
            npoi = pc;
            break;
        end
    end
end

% char POI not found
if ~isa(npoi, 'double')
    return;
end

% check double POI
if numel(npoi) ~= 1 || isinf(npoi) || isnan(npoi) || npoi ~= fix(npoi) || npoi < 1 || npoi > numpois
    return;
end

% get selected poi's vertices
v = pois(npoi).Vertices(:)';

% check validity
nvert = size(bc.VertexColor, 1);
if any(v > nvert)
    warning('neuroelf:xff:badFileContents', 'Highest vertex in POI surpasses SRF vertex index.');
    return;
end

% if mixing...
if mix < 1

    % get original colors
    cconvex = [0, round(bc.ConvexRGBA(1:3) * 255)];
    cconcave = [0, round(bc.ConcaveRGBA(1:3) * 255)];
    ocol = bc.VertexColor(v, :);

    % set correct convex/concave color code
    idxcc = isnan(ocol(:, 1));
    idxcx = idxcc & ocol(:, 1) == 0;
    idxcv = idxcc & ~idxcx;
    ocol(idxcx, :) = repmat(cconvex, [sum(idxcx), 1]);
    ocol(idxcv, :) = repmat(cconcave, [sum(idxcv), 1]);

    % mix colors
    ocol(:, 1) = NaN;
    ocol(:, 2:4) = round((1 - mix) * ocol(:, 2:4) + mix * repmat(colcode, [numel(v), 1]));

% otherwise
else

    % create color array
    ocol = [(NaN * ones(numel(v), 1)); repmat(colcode, [numel(v), 1])];
end

% set colors
bc.VertexColor(v, :) = ocol;
xo.C = bc;

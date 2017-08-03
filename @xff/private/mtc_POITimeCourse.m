function [poitc, poin, wr] = mtc_POITimeCourse(xo, poi, opts)
% MTC::POITimeCourse  - extract POI time course data
%
% FORMAT:       [poitc, poin] = mtc.POITimeCourse(poi [, opts])
%
% Input fields:
%
%       poi         POI object or vertex indices
%       opts        optional settings
%        .poisel    poi selection (either name or number)
%        .subsel    subject ID of current object (default: from filename)
%        .subpois   if poi is a POI object, one of 'sub_', '_sub', {'poi'}
%        .weight    if set to 2, get a cell array of TxC arrays
%
% Output fields:
%
%       poitc       TxV time course of poi(s) (or 1xV cell with TxC double)
%       poin        1xV poinames (without subject ID)
%
% Using: findfirst, limitrangec.

% Version:  v1.1
% Build:    16040518
% Date:     Apr-05 2016, 6:14 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% neuroelf library
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mtc') || ...
   (~all(xffisobject(poi(:), true, 'poi')) && ~isa(poi, 'double')) || isempty(poi)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin > 2 && isa(opts, 'double') && numel(opts) == 1 && ~isnan(opts)
    opts = struct('poisel', [], 'weight', opts);
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
bc = xo.C;
if isa(poi, 'double')
    poiidx = poi(:);
    poi = xff('new:poi');
    poic = poi.C;
    poic.POI(1).Vertices = poiidx;
    poic.POI.NrOfVertices = numel(poiidx);
    poic.POI.Name = 'temp';
    poi.C = poic;
    fromdouble = true;
    opts.poisel = 'temp';
    opts.subsel = '';
    opts.subpois = 'poi';
else
    fromdouble = false;
    poic = poi.C;
end
pn = {poic.POI.Name};
[sfp, sf] = fileparts(xo.F);
if ~isfield(opts, 'subsel') || ~ischar(opts.subsel) || isempty(opts.subsel)
    opts.subsel = '';
end
if ~isfield(opts, 'subpois') || ~ischar(opts.subpois) || ...
   ~any(strcmpi(opts.subpois(:)', {'sub_', '_sub', 'poi'}))
    opts.subpois = 'p';
else
    opts.subpois = lower(opts.subpois(1));
    if opts.subpois ~= 'p'
        if isempty(opts.subsel)
            if isempty(sf) || ~any(sf == '_') || sf(1) == '_'
                error('neuroelf:xff:badArgument', 'Cannot resolve subject ID for POI access.');
            end
            opts.subsel = regexprep(sf, '^([^_]+)_.*$', '$1');
        end
    end
end
if ~isfield(opts, 'poisel') || isempty(opts.poisel) || ...
   (~isa(opts.poisel, 'double') && ~ischar(opts.poisel) && ~iscell(opts.poisel))
    if opts.subpois == 'p'
        opts.poisel = 1:numel(poic.POI);
    end
elseif isa(opts.poisel, 'double')
    opts.poisel = opts.poisel(:)';
    if any(isinf(opts.poisel) | isnan(opts.poisel) | ...
        opts.poisel < 1 | opts.poisel > numel(poic.POI) | opts.poisel ~= fix(opts.poisel)) || ...
        numel(opts.poisel) ~= numel(unique(opts.poisel))
        error('neuroelf:xff:badArgument', 'Invalid numeric POI selection.');
    end
else
    if ischar(opts.poisel)
        opts.poisel = {opts.poisel(:)'};
    else
        opts.poisel = opts.poisel(:);
    end
    try
        for pc = numel(opts.poisel):-1:1
            switch (opts.subpois)
                case {'p'}
                    opts.poisel{pc} = ne_methods.findfirst(strcmpi(opts.poisel{pc}(:)', pn));
                case {'_'}
                    opts.poisel{pc} = ...
                        ne_methods.findfirst(strcmpi([opts.poisel{pc}(:)' '_' opts.subsel], pn));
                case {'s'}
                    opts.poisel{pc} = ...
                        ne_methods.findfirst(strcmpi([opts.subsel '_' opts.poisel{pc}(:)'], pn));
            end
        end
        opts.poisel = cat(1, opts.poisel{:});
    catch xfferror
        neuroelf_lasterr(xfferror);
        error('neuroelf:xff:badArgument', 'Invalid POI name selection: %s.', xfferror.message);
    end
end
numpois = numel(opts.poisel);
if ~isfield(opts, 'weight') || ...
   (~isa(opts.weight, 'double') && ~iscell(opts.weight)) || ...
   (isa(opts.weight, 'double') && (numel(opts.weight) ~= 1 || ...
     ~isa(opts.weight, 'double') || isinf(opts.weight) || isnan(opts.weight) || ...
     ~any(opts.weight == (0:3)))) || ...
   (iscell(opts.weight) && ~any(numel(opts.weight) == [1, numpois]))
    opts.weight = 1;
end
if isa(opts.weight, 'double')
    opts.weight = {opts.weight};
end
if numel(opts.weight) ~= numpois
    opts.weight = opts.weight(1, ones(1, numpois));
end
poig = cell(1, numpois);
pois = cell(1, numpois);
cellout = false;
for pc = 1:numpois
    pois{pc} = unique(round(poic.POI(opts.poisel(pc)).Vertices));
    if numel(opts.weight{pc}) == 1 && opts.weight{pc} == 2
        cellout = true;
    end
end
if nargout > 1
    poin = cell(1, numpois);
    for pc = 1:numpois
        switch (opts.subpois)
            case 'p'
                poin{pc} = pn{opts.poisel(pc)};
            case '_'
                poin{pc} = regexprep(pn{opts.poisel(pc)}, ['_' opts.subsel '^'], '', 'ignorecase');
            case 's'
                poin{pc} = regexprep(pn{opts.poisel(pc)}, ['$' opts.subsel '_'], '', 'ignorecase');
        end
    end
    if nargout > 2
        wr = cell(1, numpois);
    end
end

% clear temporary object
if fromdouble
    delete(poi);
end

% iterate
nvert = size(bc.MTCData, 2);
for pc = 1:numpois
    poig{pc} = ~(isnan(pois{pc}) | pois{pc} < 1 | pois{pc} > nvert);
    if nargout > 2
        wr{pc} = pois{pc} .* double(poig{pc});
        wr{pc}(isnan(wr{pc})) = 0;
    end
end

% prepare output
poitc = cell(1, numpois);
numtp = size(bc.MTCData, 1);
for pc = 1:numpois
    poitc{pc} = zeros(numtp, size(pois{pc}, 1));
end

% access
for pc = 1:numpois
    poitc{pc}(:, poig{pc}) = bc.MTCData(:, pois{pc}(poig{pc}));
    if numel(opts.weight{pc}) == size(poitc{pc}, 2)
        poitc{pc} = diag(opts.weight{pc}(:)) * poitc{pc};
        poig{pc} = opts.weight{pc}(:) .* double(poig{pc}(:));
    end
end

% average within POIs
if ~cellout
    for pc = 1:numpois
        poitc{pc} = sum(poitc{pc}, 2) ./ sum(poig{pc});
    end
    poitc = cat(2, poitc{:});
end

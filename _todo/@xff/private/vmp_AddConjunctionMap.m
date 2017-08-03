function xo = vmp_AddConjunctionMap(xo, maps)
% VMP::AddConjunctionMap  - add a conjunction map to a VMP
%
% FORMAT:       [vmp] = vmp.AddConjunctionMap(maps)
%
% Input fields:
%
%       maps        1xN vector specifying the maps to conjugate
%
% Output fields:
%
%       vmp         VMP with added map
%
% Using: correlinvtstat, correltstat, sdist.

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:12 PM EST
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

% neuroelf library
global ne_methods;
correlinvtstat = ne_methods.correlinvtstat;
correltstat    = ne_methods.correltstat;
sdist          = ne_methods.sdist;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp') || ...
   ~isa(maps, 'double') || numel(maps) < 2 || any(isinf(maps(:)) | isnan(maps(:))) || ...
    any(maps(:) < 1 | maps(:) ~= fix(maps(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
maps = maps(:)';
if any(maps > numel(bc.Map))
    error('neuroelf:xff:badArgument', 'Given map number(s) out of bounds.');
end

% create empty map
tmap = zeros(size(bc.Map(1).VMPData));
gmap = nan(size(bc.Map(1).VMPData));

% iterate over conjugate maps
for cc = maps

    % get stats
    smap = double(bc.Map(cc).VMPData);
    smap(isinf(smap) | isnan(smap)) = 0;
    sgmap = sign(smap);

    % compute p-values
    switch (bc.Map(cc).Type)

        % t-Map
        case 1
            pmap = sdist('tcdf', -abs(smap), bc.Map(cc).DF1);

        % r-Map
        case 2
            pmap = sdist('tcdf', -correltstat(abs(smap), bc.Map(cc).DF1 + 2), bc.Map(cc).DF1);

        % F-Map
        case 4
            pmap = sdist('fcdf', smap, bc.Map(cc).DF1, bc.Map(cc).DF2);
            sgmap(:) = nan;

    end

    % conjunction
    tmap = max(tmap, pmap);
    gmap(isnan(gmap) & ~isnan(sgmap)) = sgmap(isnan(gmap) & ~isnan(sgmap));
    cgmap = ~isnan(gmap) & ~isnan(sgmap);
    gmap(cgmap) = gmap(cgmap) .* (gmap(cgmap) == sgmap(cgmap));
end

% set nan in sign map to 1
gmap(isnan(gmap)) = 1;

% convert back to stats of first map
map1 = bc.Map(maps(1));
switch (map1.Type)

    % t-Map
    case 1
        tmap = gmap .* (-sdist('tinv', tmap, map1.DF1));

    % r-Map
    case 2
        tmap = gmap .* correlinvtstat(-sdist('tinv', tmap, map1.DF1), map1.DF1 + 2);

    % F-Map
    case 4
        tmap = sdist('finv', tmap, map1.DF1, map1.DF2);

end

% add conjugated map to VMP
bc.Map(end + 1) = map1;
bc.Map(end).Name = sprintf('%s conjugated with%s', bc.Map(end).Name, sprintf(' %d', maps(2:end)));
bc.Map(end).VMPData = single(tmap);
xo.C = bc;

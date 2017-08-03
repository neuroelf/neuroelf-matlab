function xo = smp_SplitCCRetinoMap(xo, mapno, factor)
% SMP::SplitCCRetinoMap  - add a CC retinotopic map with split info
%
% FORMAT:       smp.SplitCCRetinoMap(mapno [,factor]);
%
% Input fields:
%
%       mapno       number of map(s) to split
%       factor      either 2 or 4 (number of visual fields, default: 4)
%
% No output fields.

% Version:  v1.1
% Build:    16021210
% Date:     Feb-12 2016, 10:10 AM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'smp') || ...
   ~isa(mapno, 'double') || isempty(mapno) || numel(mapno) ~= length(mapno) || ...
    any(isinf(mapno) | isnan(mapno) | mapno < 1 | mapno ~= fix(mapno))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if any(mapno > numel(bc.Map))
    error('neuroelf:xff:badArgument', 'Selected map number(s) out of bounds.');
end
mapno = unique(mapno(:)');
if nargin < 3 || ~isa(factor, 'double') || numel(factor) ~= 1 || ~any(factor == [2, 4])
    factor = 4;
end

% iterate over maps to add
added = false;
for mc = mapno

    % get and check map type
    map = bc.Map(mc);
    if map.Type ~= 3
        warning('neuroelf:xff:invalidOption', 'CC-splitting only valid for CC maps.');
        continue;
    end

    % get data and split lag/r values
    mdat = map.SMPData;
    mlag = fix(mdat) / 1000;
    if any(mlag ~= fix(mlag))
        error('neuroelf:xff:invalidFile', 'Bad lag values in map.');
    end
    mr = mdat - 1000 * mlag;
    if any(mr < 0 | mr > 1)
        error('neuroelf:xff:invalidFile', 'Bad correlation value in map.');
    end

    % threshold r
    delr = (mr < map.LowerThreshold);
    mr = mr / factor;

    % get max number of lags
    maxl = map.NrOfLags - 1;

    % split into two halves anyway
    spi = find(mlag >= (maxl / 2));
    mr(spi) = mr(spi) + 0.5;
    mlag(spi) = maxl - mlag(spi);
    maxl = max(mlag);

    % only split again for factor 4
    if factor == 4
        spi = find(mlag >= (maxl / 2));
        mr(spi) = mr(spi) + 0.25;
        mlag(spi) = maxl - mlag(spi);
    end
    mr(delr) = 0;

    % build new map
    map.SMPData = mr + 1000 * mlag;
    map.MaxLag = max(mlag);
    map.NrOfLags = map.MaxLag + 1;
    map.CCOverlay = 0;
    map.LowerThreshold = map.LowerThreshold / factor;
    map.UpperThreshold = (factor + map.UpperThreshold - 1) / factor;
    map.Name = sprintf('Split (%d): %s', factor, map.Name);
    bc.Map(end + 1) = map;
    added = true;
end

% added maps?
if added
    bc.NrOfMaps = numel(bc.Map);
end

% set back
xo.C = bc;

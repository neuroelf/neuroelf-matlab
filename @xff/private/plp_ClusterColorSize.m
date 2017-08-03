function xo = plp_ClusterColorSize(xo, opts)
% PLP::ClusterColorSize  - set color and size values by clustering
%
% FORMAT:       [plp = ] plp.ClusterColorSize([opts])
%
% Input fields:
%
%       opts        1x1 options struct
%        .clcolumn  clustering column (name or number, default: 'Study')
%        .colhigh   color for high clustering (default: [64, 255, 64]);
%        .collow    color for low clustering (default: [0, 0, 0])
%
% Output fields:
%
%       plp         altered object
%
% Using: clusterdist, findfirst.

% Version:  v1.1
% Build:    16021210
% Date:     Feb-12 2016, 10:24 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
findfirst = ne_methods.findfirst;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'plp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
cns = bc.ColumnNames;
if nargin < 2 || numel(opts) ~= 1 || ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'clcolumn') || isempty(opts.clcolumn) || ...
   (~ischar(opts.clcolumn) && (~isa(opts.clcolumn, 'double') || ...
     numel(opts.clcolumn) ~= 1 || isinf(opts.clcolumn) || isnan(opts.clcolumn) || ...
     opts.clcolumn < 1 || opts.clcolumn > size(bc.Points, 2) || opts.clcolumn ~= fix(opts.clcolumn)))
    opts.clcolumn = 'Study';
end
if ischar(opts.clcolumn)
    opts.clcolumn = findfirst(strcmpi(opts.clcolumn(:)', cns));
    if isempty(opts.clcolumn)
        error('neuroelf:xff:badArgument', 'Missing named column in PLP object.');
    end
end
if ~isfield(opts, 'colhigh') || ~isa(opts.colhigh, 'double') || numel(opts.colhigh) ~= 3 || ...
    any(isinf(opts.colhigh) | isnan(opts.colhigh))
    opts.colhigh = [0.25, 1, 0.25];
else
    opts.colhigh = min(255, max(0, opts.colhigh(:)'));
    if any(opts.colhigh > 1)
        opts.colhigh = (1 / 255) .* opts.colhigh;
    end
end
if ~isfield(opts, 'collow') || ~isa(opts.collow, 'double') || numel(opts.collow) ~= 3 || ...
    any(isinf(opts.collow) | isnan(opts.collow))
    opts.collow = [0, 0, 0];
else
    opts.collow = min(255, max(0, opts.collow(:)'));
    if any(opts.collow > 1)
        opts.collow = (1 / 255) .* opts.collow;
    end
end

% get other columns we need
xc = findfirst(strcmpi(cns, 'x'));
yc = findfirst(strcmpi(cns, 'y'));
zc = findfirst(strcmpi(cns, 'z'));
if isempty(xc) || isempty(yc) || isempty(zc)
    error('neuroelf:xff:badObject', 'PLP object is missing X/Y/Z columns.');
end

% computation
pcrd = bc.Points(:, [xc, yc, zc]);
pidx = ~any(isinf(pcrd) | isnan(pcrd), 2) | ...
    bc.Points(:, opts.clcolumn) < 1 | bc.Points(:, opts.clcolumn) > numel(bc.Labels);
pd = ne_methods.clusterdist(pcrd(pidx, :), bc.Labels(round(bc.Points(pidx, opts.clcolumn))));
ns = numel(unique(bc.Points(pidx, opts.clcolumn)));
maxdist = max(4 * ns, 0.5 * ns * max(pd));
scsize = 4 .* max(2, maxdist - (0.75 * ns) .* pd);
maxsize = max(scsize);
sccol = max(0.1, scsize ./ maxsize);
sccol = sccol * opts.colhigh + (1 - sccol) * opts.collow;
sccol = max(0, min(255, round(255 .* sccol)));

% set or add
szcol = findfirst(strcmpi(cns, 'size'));
clcol = findfirst(strcmpi(cns, 'color'));
if isempty(szcol)
    bc.Points(:, end+1) = 0;
    bc.ColumnNames{end+1} = 'Size';
    szcol = size(bc.Points, 2);
    bc.NrOfColumns = szcol;
end
if isempty(clcol)
    bc.Points(:, end+1) = 0;
    bc.ColumnNames{end+1} = 'Color';
    clcol = size(bc.Points, 2);
    bc.NrOfColumns = clcol;
end

% set back
bc.Colors = sccol;
bc.Points(~pidx, [clcol, szcol]) = 1;
bc.Points(pidx, clcol) = 1:numel(scsize);
bc.Points(pidx, szcol) = scsize;
xo.C = bc;

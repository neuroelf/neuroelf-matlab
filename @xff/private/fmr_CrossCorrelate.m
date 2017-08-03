function ccmap = fmr_CrossCorrelate(xo, xo2, opts)
% FMR::CrossCorrelate  - create CC map of two FMRs
%
% FORMAT:       ccmap = fmr.CrossCorrelate(fmr2 [, opts])
%
% Input fields:
%
%       fmr2        second FMR (must match in dims and layout)
%       opts        options settings
%        .reverse   reverse time courses of second VTC
%
% Output fields:
%
%       ccmap       cross-correlation r-VMP
%
% Note: the toolbox internal cov_nd function is used which gives
%       slightly different r values than corrcoef.
%
% Using: correlinvtstat, cov_nd, sdist.

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:47 AM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr') || ...
    numel(xo2) ~= 1 || ~xffisobject(xo2, true, 'fmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
fmr1f = xo.F;
fmr2f = xo2.F;
if isempty(fmr1f)
    fmr1f = '<Unknown 1>';
end
if isempty(fmr2f)
    fmr2f = '<Unknown 2>';
end
fmr1 = xo.C;
fmr2 = xo2.C;
if fmr1.ResolutionX ~= fmr2.ResolutionX || fmr1.ResolutionY ~= fmr2.ResolutionY || ...
    fmr1.NrOfSlices ~= fmr2.NrOfSlices || fmr1.NrOfVolumes ~= fmr2.NrOfVolumes
    error('neuroelf:xff:badArgument', 'Dimension/Layout mismatch.');
end
if fmr1.NrOfVolumes < 3
    error('neuroelf:xff:badArgument', 'FMRs must have at least 3 volumes each.');
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'reverse') || isempty(opts.reverse) || ~islogical(opts.reverse)
    opts.reverse = false;
else
    opts.reverse = opts.reverse(1);
end
if opts.reverse
    revstr = ' TC-reversed';
else
    revstr = '';
end

% create map
df1 = fmr1.NrOfVolumes - 2;
ccmap = xff('new:map');
map = ccmap.C;
map.Type = 1;
map.NrOfSlices = fmr1.NrOfSlices;
map.CombinedTypeSlices = 10000 * map.Type + map.NrOfSlices;
map.DimX = fmr1.ResolutionX;
map.DimY = fmr1.ResolutionY;
map.ClusterSize = 1;
map.LowerThreshold = ne_methods.correlinvtstat( ...
    -ne_methods.sdist('tinv', 0.005, df1), fmr1.NrOfVolumes);
map.UpperThreshold = ne_methods.correlinvtstat( ...
    -ne_methods.sdist('tinv', 0.0001, df1), fmr1.NrOfVolumes);
map.DF1 = df1;
map.DF2 = 0;
map.NameOfSDMFile = sprintf('<CC %s <-> %s%s>', fmr1f, fmr2f, revstr);
map.Map = struct('Number', 1, 'Data', single(zeros(fmr1.ResolutionY, fmr1.ResolutionX)));
map.Map = map.Map(1, ones(1, fmr1.NrOfSlices));

% iterate over last spatial dim
cov_nd = ne_methods.cov_nd;
for z = 1:fmr1.NrOfSlices

    % get components for cov_nd
    if fmr1.DataStorageFormat < 2
        r1 = double(fmr1.Slice(z).STCData(:, :, :));
    else
        r1 = double(fmr1.Slice.STCData(:, :, :, z));
    end
    if fmr2.DataStorageFormat < 2
        r2 = double(fmr2.Slice(z).STCData(:, :, :));
    else
        r2 = double(fmr2.Slice.STCData(:, :, :, z));
    end
    if opts.reverse
        r2 = r2(:, :, end:-1:1);
    end

    % compute r value
    [cc, cr] = cov_nd(r1, r2);
    cr(isinf(cr) | isnan(cr)) = 0;
    cr = sign(cr) .* (1 - abs(cr));
    cr(cr == 1) = 0;
    map.Map(z).Number = z;
    map.Map(z).Data = single(cr);
end

% set data to VMP
ccmap.C = map;

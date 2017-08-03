function sd = smoothned(dc, dv, k, kd, sc)
% smoothned  - smooth non-equidistant data with a gaussian kernel
%
% FORMAT:       sd = smoothned(dc, dv, k [, kd [, sc]])
%
% Input fields:
%
%       dc      data coordinates
%       dv      1-d data
%       k       gaussian kernel size or 1-d kernel
%       kd      kernel value distance
%       sc      sample coordinates (default: at dc)
%
% Output fields:
%
%       sd      smoothed data

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
if nargin < 3 || ...
   ~isa(dc, 'double') || ...
   ~isnumeric(dv) || ...
   ~isequal(size(dc), size(dv)) || ...
    numel(dc) ~= max(size(dc)) || ...
    isempty(dv) || ...
    any(isinf(dc) | isnan(dc)) || ...
   ~isa(k, 'double') || ...
    isempty(k) || ...
    numel(k) ~= max(size(k)) || ...
   (numel(k) == 1 && ...
    any(isinf(k) | isnan(k) | k < 0)) || ...
   (numel(k) > 1 && ...
    any(isinf(k) | isnan(k) | k < 0 | k > 1))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument.' ...
    );
end
if ~isa(dv, 'double')
    dv = double(dv);
end
dv = dv(:);
if nargin < 4 || ...
    numel(kd) ~= 1 || ...
   ~isa(kd, 'double') || ...
    isinf(kd) || ...
    isnan(kd) || ...
    kd < 1 || ...
    kd ~= round(kd)
    kd = 1024;
end
if numel(k) == 1
    f = k / sqrt(8 * log(2));
    md = ceil(4 * f);
    k = exp(- (-md:1/kd:md)' .^ 2 ./ (2 * f .^ 2));
elseif ((numel(k) - 1) / kd) ~= round((numel(k) - 1) / kd)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid kernel spacing supplied.' ...
    );
end
kw = 0.5 * (numel(k) - 1) / kd;
k = [k(1); k(:); k(end)];
if nargin < 5 || ...
   ~isa(sc, 'double') || ...
    numel(sc) ~= max(size(sc))
    sc = dc;
elseif ndims(sc) > 2
    sc = sc(:);
end

% generate output data first
sd = zsz(sc);
wd = sd;

% remove bad samples
bv = (isinf(dv) | isnan(dv));
dc(bv) = [];
dv(bv) = [];

% sort coordinates
nc = numel(dc);
[dc, dci] = sort(dc(:));
dv = dv(dci);
[sc, sci] = sort(sc(:));

% set first sample to take into account as first coordinate
fs = 1;

% for each sample coordinate
for c = 1:numel(sc)

    % get minimum coordinate to weight
    bg = sc(c) - kw;
    ed = sc(c) + kw;

    % find first coordinate that matches this criterion
    fs = findfirst(dc >= bg, fs);

    % iterate over those coordinates
    for s = fs:nc

        % if coordinate out of range
        if dc(s) > ed

            % leave loop
            break;
        end

        % otherwise get weight
        diffc = sc(c) - dc(s);
        diffs = 2 + kd * (diffc + kw);
        diffr = diffs - floor(diffs);
        sw = (1 - diffr) * k(floor(diffs)) + diffr * k(ceil(diffs));

        % add to sample and weight
        si = sci(c);
        sd(si) = sd(si) + sw * dv(s);
        wd(si) = wd(si) + sw;
    end
end

% re-weight to average
sd = sd ./ wd;
sd(isinf(sd)) = NaN;

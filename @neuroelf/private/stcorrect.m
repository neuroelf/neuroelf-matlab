function [d, w] = stcorrect(d, shifts, opts)
% stcorrect  - correct slice time differences in 4D data
%
% FORMAT:       d = stcorrect(d, shifts [, opts])
%
% Input fields:
%
%       d           4D data (or transio object)
%       shifts      slice-wise shift values (must match size(d, 3))
%       opts        optional settings
%        .dtype     output datatype, either of {'single'} or 'double'
%        .interp    interpolation method (default: lanczos8)
%        .snoise    detect spatial noise (within slice, default: true)
%        .tfcutoff  temporal filter cutoff in seconds (default: 60)
%        .tr        TR in secones (default: 2)
%
% Output fields:
%
%       d           corrected data
%       w           weighting matrix (pseudo certainty, single!)
%
% Note: if all shifts are 0, only outlier removal is performed.

% Version:  v0.9b
% Build:    11022719
% Date:     Feb-26 2011, 4:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% argument check
if nargin < 2 || ...
   ~isnumeric(d) || ...
    numel(size(d)) ~= 4 || ...
    isempty(d) || ...
   ~isa(shifts, 'double') || ...
    numel(shifts) ~= size(d, 3) || ...
    any(isinf(shifts(:)) | isnan(shifts(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
shifts = shifts(:)';
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dtype') || ...
   ~ischar(opts.dtype) || ...
    isempty(opts.dtype) || ...
    lower(opts.dtype(1)) ~= 'd'
    d = single(d);
else
    d = double(d);
end
if ~isfield(opts, 'interp') || ...
   ~ischar(opts.interp) || ...
    isempty(regexpi(opts.interp(:)', '^(cubic|linear|lanczos\d)$'))
    opts.interp = 'lanczos8';
else
    opts.interp = lower(opts.interp(:)');
end
if ~isfield(opts, 'snoise') || ...
   ~islogical(opts.snoise) || ...
    numel(opts.snoise) ~= 1
    opts.snoise = true;
end
if ~isfield(opts, 'tfcutoff') || ...
   ~isa(opts.tfcutoff, 'double') || ...
    numel(opts.tfcutoff) ~= 1 || ...
    isinf(opts.tfcutoff) || ...
    isnan(opts.tfcutoff)
    opts.tfcutoff = 60;
else
    opts.tfcutoff = max(15, min(240, opts.tfcutoff));
end
if ~isfield(opts, 'tr') || ...
   ~isa(opts.tr, 'double') || ...
    numel(opts.tr) ~= 1 || ...
    isinf(opts.tr) || ...
    isnan(opts.tr) || ...
    opts.tr < 0.05 || ...
    opts.tr > 30
    opts.tr = 2;
end

% prepare weight output if required
if nargout > 1
    w = single(0);
    w(numel(d)) = 0;
    w = reshape(w, size(d));
end

% for spatial functionality mean image is required
md = double(meannoinfnan(d, 4));
gm = sum(md(:)) / numel(md);
bthr = 0.5 .* gm;

% from which a histogram is computed
mdmm = minmaxmean(md(md ~= 0));
mdh = histcount(md(:), mdmm(1), mdmm(2));

% generate filtering matrices
nt = size(d, 4);
fthr = 0.75 .* nt;
ntinv = 1 ./ nt;
nti = [1, 1:nt, nt];
[d0, fm] = tempfilter(zeros(nt, 1), struct('tempdct', opts.tfcutoff / opts.tr));
fm(:, end+1) = 1;
[d0, fmd] = tempfilter(zeros(nt - 1, 1), ...
    struct('tempdct', ((nt - 1) / nt) * opts.tfcutoff / opts.tr));
fmd(:, end+1) = 1;
fmixx = pinv(fm' * fm) * fm';
fmdixx = pinv(fmd' * fmd) * fmd';

% iterate over slices
nx = size(d, 1);
ny = size(d, 2);
nxy = nx * ny;
nxyinv = 1 ./ nxy;
for sc = 1:numel(shifts)

    % get slice's time course
    stc = squeeze(double(d(:, :, sc, :)));
    szs = size(stc);

    % set slice weights = 1
    sw = single(true(size(stc)));

    % background voxels
    bvox = find(lsqueeze((md(:, :, sc) < bthr) & (sum(stc < bthr, 3) >= fthr)));
    fvox = find(lsqueeze((md(:, :, sc) > gm) & (sum(stc > gm, 3) >= fthr)));

    % detect outliers -> spatial noise
    if opts.snoise

        % compute get transposed version
        sd = reshape(stc, nxy, nt)';

        % compute derivative
        sdd = diff(sd);

        % temporally filter both
        sd = sd - fm * (fmixx * sd);
        sdd = sdd - fmd * (fmdixx * sdd);

        % remove over-noisy voxels from estimate
        sdh = ntinv .* sum(sd .* sd, 1);
        sd = limitrangec(sd, 0, 2 .* sqrt(median(sdh(sdh > 0))), 0);
        sddh = ntinv .* sum(sdd .* sdd, 1);
        sdd = limitrangec(sdd, 0, 2 .* sqrt(median(sddh(sddh > 0))), 0);

        % estimate smoothness as sum of squares / number voxel
        sdsm = nxyinv .* sum(sd .* sd, 2);
        sddsm = nxyinv .* sum(sdd .* sdd, 2);
    end

    % extend time course by one at the beginning and the end
    stc = stc(:, :, nti);
    stcx = stc(:, :, 2:7);
    for ic = 1:4
        stcx = flexinterpn_method(stcx, [Inf * ones(1, 3); ...
            1, 1, 0.75; ones(1, 3); size(stcx)], 0, 'lanczos5');
    end
    stc(:, :, 1) = stcx(:, :, 1);
    stcx = stc(:, :, end-6:end-1);
    szs(3) = szs(3) + 1;
    for ic = 1:4
        stcx = flexinterpn_method(stcx, [Inf * ones(1, 3); ...
             1, 1, 1.25; ones(1, 3); size(stcx) + [0, 0, 1]], 0, 'lanczos5');
    end
    stc(:, :, end) = stcx(:, :, end);

    % don't do anything for a shift of 0!
    if opts.tshift(sc) == 0
        continue;
    end

    % interpolate new data
    scslc = flexinterpn_method(stcd(:, :, :, opts.order(sc)), ...
        [Inf * ones(1, 3); 1, 1, 2 - opts.tshift(sc); ones(1, 3); szs(1:3)], ...
        0, opts.interp);

    % put back into data
    if bc.FileVersion > 4 && ...
        numel(bc.Slice) == 1
        bc.Slice.STCData(:, :, :, opts.order(sc)) = max(0, scslc);
    else
        if isui16
            bc.Slice(opts.order(sc)).STCData = uint16(round(max(0, scslc)));
        else
            bc.Slice(opts.order(sc)).STCData = single(max(0, scslc));
        end
    end
end

% set into content
xffsetcont(hfile.L, bc);

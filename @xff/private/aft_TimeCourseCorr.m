function [stcc, st, gso] = aft_TimeCourseCorr(xo, opts)
% AFT::TimeCourseCorr  - compute slice-based cross-correlation
%
% FORMAT:       [stcc, st, gso] = obj.TimeCourseCorr([opts])
%
% Input fields:
%
%       opts        optional fields
%        .fgthresh  foreground threshold (default: 2)
%        .remmean   remove global mean (default: true)
%        .tfilter   temporally high-pass filter (default: 30, in TRs)
%
% Output fields:
%
%       stcc        slice timecourse correlation matrix
%       st          slice timecourses
%       gso         guessed slice acquisition order
%
% TYPES: FMR, HDR, HEAD, VTC
%
% Using: calcbetas, limitrangec, maxpos, meannoinfnan, minmaxmean,
%        tempfilter, ztrans.

% Version:  v1.1
% Build:    16012511
% Date:     Jan-25 2016, 11:15 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012 - 2014, 2016, Jochen Weber
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

% using many library functions
global ne_methods;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fmr', 'hdr', 'head', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'fgthresh') || ~isa(opts.fgthresh, 'double') || numel(opts.fgthresh) ~= 1 || ...
    isinf(opts.fgthresh) || isnan(opts.fgthresh) || opts.fgthresh <= 0
    opts.fgthresh = 2;
end
if ~isfield(opts, 'remmean') || ~islogical(opts.remmean) || numel(opts.remmean) ~= 1
    opts.remmean = true;
end
if ~isfield(opts, 'tfilter') || ~isa(opts.tfilter, 'double') || numel(opts.tfilter) ~= 1
    opts.tfilter = 30;
else
    opts.filter = ne_methods.limitrangec(opts.tfilter, 0, 1e6, 0);
end

% get data and temporal/slice dim for object
ot = lower(xo.S.Extensions{1});
bc = xo.C;
switch (ot)

    % for FMRs
    case 'fmr'
        tdim = 3;
        sdim = 4;

    % for Analyze/NIftI
    case 'hdr'
        tdim = 4;
        sdim = 3;
        tc = bc.VoxelData;

    % for AFNI/HEAD
    case 'head'
        tdim = 4;
        sdim = 3;

    % for BVQX VTCs
    case 'vtc'
        tdim = 1;
        sdim = 3;
        tc = bc.VTCData;
end

% check that the data can be used
if size(tc, tdim) < 5 || ndims(tc) ~= 4
    error('neuroelf:xff:badArgument', ...
        'Less than 5 functional volumes detected or data not 4D.');
end
stc = size(tc);

% resolve and make double precision
if istransio(tc) || ~isa(tc, 'double')
    tc = double(tc);
end

% slice-based access
sa = {':', ':', ':', ':'};
ta = {':', ':', ':', ':'};

% get mask
mm = mean(tc, tdim);
mmsk = (mm > (opts.fgthresh .* mean(mm(:))));

% refine mask by removing voxels with too large range
minm = min(tc, [], tdim);
maxm = max(tc, [], tdim);
difm = maxm - minm;
difm(difm == 0 | isinf(difm)) = NaN;
difmm = ne_methods.minmaxmean(difm(~isnan(difm)), 1);
mmsk = mmsk & (difm > 0) & (difm < (difmm(3) + 2 * sqrt(difmm(6))));

% refine mask by removing voxels with too high variance
stdtc = stc;
stdtc(tdim) = 1;
stdtc = zeros(stdtc);
for c = 1:stc(tdim)
    ta{tdim} = c;
    tcp = tc(ta{:});
    stdtc = stdtc + (mm - tcp) .^ 2;
end
stdtc = sqrt(stdtc ./ (stc(tdim) - 1));
mmsk = mmsk & (stdtc < mean(stdtc(mmsk)));

% remove mean foreground
if opts.remmean
    mtc = zeros(stc(tdim), 1);
    for c = 1:stc(tdim)
        ta{tdim} = c;
        tcp = tc(ta{:});
        mtc(c) = ne_methods.meannoinfnan(tcp(mmsk(:)));
    end
    mtc = [ne_methods.ztrans(mtc), ones(stc(tdim), 1)];
    b = ne_methods.calcbetas(mtc, tc, tdim);
    for c = 1:stc(tdim)
        ta{tdim} = c;
        tc(ta{:}) = tc(ta{:}) - mtc(c) .* reshape(b(:, :, :, 1), size(tcp));
    end
end

% create time courses
nslc = stc(sdim);
st = zeros(stc(tdim), nslc);

% iterate over slices
for sc = 1:nslc

    % set access
    sa{sdim} = sc;

    % get mask and time course data
    smsk = mmsk(sa{:});
    sltc = tc(sa{:});

    % iterate over time points
    for c = 1:stc(tdim)
        ta{tdim} = c;
        sltp = sltc(ta{:});
        st(c, sc) = mean(mean(sltp(smsk(:))));
    end
end

% filter time courses
if opts.tfilter
    st = ne_methods.tempfilter(st, struct('tdim', 1, 'temp', true, 'tempdct', opts.tfilter));
end

% compute cross correlation
stcc = corrcoef(st);
stcc(isnan(stcc)) = 0;

% guess slice order
if nargout > 2

    % create output
    gso = zeros(nslc, 1);

    % compute for how many "lags" we check
    lc = ceil(sqrt(nslc));
    lcc = zeros(lc, 1);

    % look at the "middle" slice and look for typical signatures
    mss = stcc(:, round(0.5 * nslc));

    % compute correlations for typical values
    for sc = 1:lc
        mscc = corrcoef(mss(1:end-sc), mss(1+sc:end));
        lcc(sc) = mscc(2);
    end

    % guess that the lag is the max of these
    lag = ne_methods.maxpos(lcc);

    % if correlation < 0.5, stop
    if lcc(lag) < 0.5
        return;
    end

    % check with this lag for all slices
    slcc = zeros(nslc, 1);
    for sc = 1:nslc
        mscc = corrcoef(stcc(1:end-lag, sc), stcc(1+lag:end, sc));
        slcc(sc) = mscc(2);
    end
    if mean(slcc) < 0.5
        return;
    end

    % figure out order with given lag
end

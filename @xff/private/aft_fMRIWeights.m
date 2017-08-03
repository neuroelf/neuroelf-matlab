function xo = aft_fMRIWeights(xo, opts)
% AFT::fMRIWeights  - estimate goodness-of-data weights for fMRI data
%
% FORMAT:       [obj = ] obj.fMRIWeights([opts])
%
% Input fields:
%
%       opts        optional fields
%        .diff      compute 1st derivative and use bisquare weight (true)
%        .slicetc   compute slice-based time course (false)
%
% Output fields:
%
%       obj         object with .RunTimeVars.fMRIWeights set
%
% TYPES: FMR, HDR, HEAD, VTC
%
% Using: flexinterpn_method, limitrangec, lsqueeze.

% Version:  v1.0
% Build:    16012413
% Date:     Jan-24 2016, 1:57 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012 - 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fmr', 'hdr', 'head', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'diff') || ~islogical(opts.diff) || numel(opts.diff) ~= 1
    opts.diff = true;
end
if ~isfield(opts, 'slicetc') || ~islogical(opts.slicetc) || numel(opts.slicetc) ~= 1
    opts.slicetc = false;
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
        sdim = 4;
        tc = bc.VTCData;
end
xydim = setdiff(1:4, [sdim, tdim]);

% check that the data can be used
if size(tc, tdim) < 5 || ndims(tc) ~= 4
    error('neuroelf:xff:badArgument', ...
        'Less than 5 functional volumes detected or data not 4D.');
end
stc = size(tc);

% nothing to do
if ~opts.diff && ~opts.slicetc

    % set weights to 1
    bc.fMRIWeights = single(1);
    bc.fMRIWeights(1:stc(1), 1:stc(2), 1:stc(3), 1:stc(4)) = 1;
    xo.C = bc;

    % and return early
    return;
end

% resolve and make double precision
if istransio(tc) || ~isa(tc, 'double')
    tc = double(tc);
end

% remove mean
trm = ones(1, 4);
trm(tdim) = stc(tdim);
mtc = mean(tc, tdim);
tc = tc - repmat(mtc, trm);

% create initial weights
tcw = zeros(size(tc));

% diff-based
if opts.diff
end

% slice timecourse-based
if opts.slicetc

    % slice-based access
    sa = {':', ':', ':', ':'};

    % iterate over slices
    for sc = 1:stc(sdim)

        % compute slice-based timecourse
        sa{sdim} = sc;
        sltc = tc(sa{:});
        msltc = ne_methods.lsqueeze(mean(mean(sltc, xydim(1)), xydim(2)));
        ssltc = zeros(size(msltc));

        % compute per-volume noise estimate
        ta = {':', ':', ':', ':'};
        for c = 1:stc(tdim)
            ta{tdim} = c;
            ssltc(c) = std(ne_methods.lsqueeze(sltc(ta{:})));
        end

        % compute discrete diff of slice time course
        sdiff = abs(ne_methods.flexinterpn_method(abs(diff( ...
            msltc)), [inf; 0.5; 1; stc(tdim)], 'cubic'));

        % detect outliers on diff
        msdiff = mean(sdiff);
        ssdiff = std(sdiff);
        sdiff(sdiff < msdiff) = msdiff;
        sdiff = max(0, 1 + (msdiff - sdiff) ./ (4.685 * ssdiff));

        % detect outliers on standard deviation
        mssltc = mean(ssltc);
        sssltc = std(ssltc);
        ssltc(ssltc < mssltc) = mssltc;
        sdiff = min(sdiff, max(0, 1 + (mssltc - ssltc) ./ (4.685 * sssltc)));

        % add to weights
        ta = sa;
        for c = 1:stc(tdim)
            ta{tdim} = c;
            tcw(ta{:}) = tcw(ta{:}) + sdiff(c);
        end
    end
end

% scale weights
tcw = (1 / (double(opts.diff) + double(opts.slicetc))) .* tcw;

% store weights
bc.RunTimeVars.fMRIWeights = single(ne_methods.limitrangec(tcw, 0, 1, 0));
xo.C = bc;

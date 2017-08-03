function [tc, ftc, ci] = guesssliceorder
% guesssliceorder  - algorithm to guess slice order
%
% FORMAT:       [tc, ftc, ci] = guesssliceorder
%
% No input fields. (file selector)
%
% Output fields:
%
%        tc         slice time course correlation matrices (3D)
%        ftc        fourier transform over correlation matrices (3D)
%        ci         confidence index (1D vector)
%
% Note: the order of planes/CI values is:
%       - unchanged
%       - ascending (non-interleaved)
%       - descending (non-interleaved)
%       - ascending interleaved (odd first)
%       - ascending interleaved (even first)
%       - descending interleaved (odd first)
%       - descending interleaved (even first)
%       - ascending interleaved (sqrt(N) ordering)
%       - descending interleaved (sqrt(N) ordering)

% Version:  v0.9c
% Build:    13020213
% Date:     Feb-02 2013, 1:14 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, 2013, Jochen Weber
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

% try to load HDR/NII file
try
    h = xff('*.(hdr|nii)');
    if isempty(h)
        return;
    end
    if ~isxff(h, 'hdr')
        error( ...
            'neuroelf:BadFileSelected', ...
            'Bad file selected.' ...
        );
    end
catch ne_eo;
    rethrow(ne_eo);
end

% load data
h.LoadVoxelData;

% compute uncorrected features
[u_tc, u_ftc, u_ci] = hdr_tc_ci(h);

% test all reasonable combinations
[asc_tc, asc_ftc, asc_ci] = hdr_tc_ci(h, 'asc');
[desc_tc, desc_ftc, desc_ci] = hdr_tc_ci(h, 'des');
[ai1_tc, ai1_ftc, ai1_ci] = hdr_tc_ci(h, 'aint1');
[ai2_tc, ai2_ftc, ai2_ci] = hdr_tc_ci(h, 'aint2');
[di1_tc, di1_ftc, di1_ci] = hdr_tc_ci(h, 'dint1');
[di2_tc, di2_ftc, di2_ci] = hdr_tc_ci(h, 'dint2');
[asq_tc, asq_ftc, asq_ci] = hdr_tc_ci(h, 'asqr');
[dsq_tc, dsq_ftc, dsq_ci] = hdr_tc_ci(h, 'dsqr');

% clear original object
h.ClearObject;

% create output figure
if nargout < 1

% or create output variables
else
    tc = cat(3, u_tc, asc_tc, desc_tc, ai1_tc, ai2_tc, di1_tc, di2_tc, asq_tc, dsq_tc);
    ftc = cat(3, u_ftc, asc_ftc, desc_ftc, ai1_ftc, ai2_ftc, di1_ftc, di2_ftc, asq_ftc, dsq_ftc);
    ci = cat(1, u_ci, asc_ci, desc_ci, ai1_ci, ai2_ci, di1_ci, di2_ci, asq_ci, dsq_ci);
end

% sub functions



function [tc, ftc, ci] = hdr_tc_ci(h, so)

% slice time correction?
if nargin > 1
    h = h.CopyObject;
    h.SliceTiming(struct('order', so));
end

% get time course cross-correlation matrix
tc = h.TimeCourseCorr;

% remove temporary object
if nargin > 1
    h.ClearObject;
end

% compute fft and normalize
ftc = abs(fftshift(fftn(tc)));
ftc = ftc ./ sum(ftc(:));
nx = size(ftc, 1);

% compute composite index
[x, y] = ndgrid(1:nx, 1:nx);
xy = (x - (nx + 1)) .^ 2 + (y - (nx + 1)) .^ 2;
ci = sum(ftc(:) .* xy(:));

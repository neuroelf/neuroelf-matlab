function xo3 = sdm_Concatenate(xo, xo2, combine)
% SDM::Concatenate  - concatenate SDMs
%
% FORMAT:       csdm = sdm1.Concatenate(sdm2 [, combine]);
%
% Input fields:
%
%       sdm2        second SDM
%       combine     if given and evaluates to true only concatenate in Y
%                   otherwise blocked
%
% Output fields:
%
%       csdm        combined/concatenated RTC

% Version:  v1.1
% Build:    16021210
% Date:     Feb-12 2016, 10:13 AM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'sdm') || ...
    numel(xo2) ~= 1 || ~xffisobject(xo2, true, 'sdm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc1 = xo.C;
bc2 = xo2.C;
if nargin > 2 && numel(combine) == 1 && (islogical(combine) || isnumeric(combine))
    if combine
        combine = true;
    else
        combine = false;
    end
else
    combine = false;
end

% get information
pcol1 = bc1.PredictorColors;
pcol2 = bc2.PredictorColors;
pnam1 = bc1.PredictorNames(:)';
pnam2 = bc2.PredictorNames(:)';
sdm1 = bc1.SDMMatrix;
sdm2 = bc2.SDMMatrix;

% split confounds and no confounds
if bc1.IncludesConstant || (bc1.FirstConfoundPredictor > 0 && ...
    bc1.FirstConfoundPredictor <= size(sdm1, 2))
    fcp1 = bc1.FirstConfoundPredictor;
else
    fcp1 = size(sdm1, 2) + 1;
end
if bc2.IncludesConstant || (bc2.FirstConfoundPredictor > 0 && ...
    bc2.FirstConfoundPredictor <= size(sdm2, 2))
    fcp2 = bc2.FirstConfoundPredictor;
else
    fcp2 = size(sdm2, 2) + 1;
end
if bc1.IncludesConstant > 0 && bc2.IncludesConstant == 0
    pcol2(end+1, :) = [255, 255, 255];
    pnam2(end+1) = {'Constant'};
    sdm2(:, end+1) = 1;
elseif bc2.IncludesConstant > 0 && bc1.IncludesConstant == 0
    pcol1(end+1, :) = [255, 255, 255];
    pnam1(end+1) = {'Constant'};
    sdm1(:, end+1) = 1;
end
ccol1 = pcol1(fcp1:end, :);
pcol1(fcp1:end, :) = [];
cnam1 = pnam1(fcp1:end);
pnam1(fcp1:end) = [];
csdm1 = sdm1(:, fcp1:end);
sdm1(:, fcp1:end) = [];
ccol2 = pcol2(fcp2:end, :);
pcol2(fcp2:end, :) = [];
cnam2 = pnam2(fcp2:end);
pnam2(fcp2:end) = [];
csdm2 = sdm2(:, fcp2:end);
sdm2(:, fcp2:end) = [];

% combine or block
if combine

    % combine confounds
    ccol1 = [ccol1; ccol2];
    cnam1 = [cnam1, cnam2];
    csdm1 = [[csdm1; zeros(size(csdm2, 1), size(csdm1, 2))], ...
             [zeros(size(csdm1, 1), size(csdm2, 2)); csdm2]];

    % find matching predictor pairs
    nidx = numel(pnam1) + 1;
    tidx2 = zeros(1, numel(pnam2));
    for ixc = 1:numel(tidx2)
        ixf = find(strcmpi(pnam2{ixc}, pnam1));
        if ~isempty(ixf)
            tidx2(ixc) = ixf(1);
        else
            tidx2(ixc) = nidx;
            nidx = nidx + 1;
        end
    end
    nidx = nidx - 1;

    % combine
    sdm = zeros(size(sdm1, 1) + size(sdm2, 1), nidx);
    sdm(1:size(sdm1, 1), 1:size(sdm1, 2)) = sdm1;
    for ixc = 1:numel(tidx2)
        sdm((size(sdm1, 1)+1):end, tidx2(ixc)) = sdm2(:, ixc);
    end
    pcol1(end+1:end+sum(tidx2 > numel(pnam1)), :) = pcol2(tidx2 > numel(pnam1), :);
    pnam1(end+1:end+sum(tidx2 > numel(pnam1))) = pnam2(tidx2 > numel(pnam1));

    % add confounds back to combined SDM
    pcol1 = [pcol1; ccol1];
    pnam1 = [pnam1, cnam1];
    sdm = [sdm, csdm1];
    fcp = fcp1;
else

    % block
    sdm = ...
       [[sdm1; zeros(size(sdm2, 1), size(sdm1, 2))], ...
        [zeros(size(sdm1, 1), size(sdm2, 2)); sdm2], ...
        [csdm1; zeros(size(sdm2, 1), size(csdm1, 2))], ...
        [zeros(size(sdm1, 1), size(csdm2, 2)); csdm2]];
    pcol1 = [pcol1; pcol2; ccol1; ccol2];
    pnam1 = [pnam1, pnam2, cnam1, cnam2];

    % get names lists
    onam1 = pnam1;

    % put names of list2 into list1
    for p2c = 2:length(pnam1)

        % name exists?
        nnc = 1;
        while any(strcmp(pnam1{p2c}, pnam1(1:p2c-1)))

            % build new name
            nnc = nnc + 1;
            pnam1{p2c} = sprintf('%s - S%d', onam1{p2c}, nnc);
        end
    end
    fcp = fcp1 + fcp2 - 1;
end

% put into new RTC
xo3 = xff('new:sdm');
bc3 = xo3.C;
bc3.NrOfPredictors = size(sdm, 2);
bc3.NrOfDataPoints = size(sdm, 1);
if bc1.IncludesConstant || bc2.IncludesConstant
    bc3.IncludesConstant = 1;
else
    bc3.IncludesConstant = 0;
end
bc3.FirstConfoundPredictor = fcp;
bc3.PredictorColors = pcol1;
bc3.PredictorNames = pnam1;
bc3.SDMMatrix = sdm;
bc3.RTCMatrix = sdm(:, 1:fcp-1);
xo3.C = bc3;

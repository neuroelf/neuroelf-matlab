function glm = mtc_CreateGLM(xo, sdm, glmopts)
% MTC::CreateGLM  - create a single study GLM
%
% FORMAT:       glm = mtc.CreateGLM(sdm [, glmopts])
%
% Input fields:
%
%       sdm         SDM object with same number of time points
%       glmopts     optional fields for GLM calculation
%        .ssm       SSM object
%        .threshval threshold value to compute a GLM (default: 100)
%        .threshvol threshold vols (mean value, default: [fix(N/2)+1,N])
%        .tctrans   time course transformation before estimation
%                   'none' or '' no transformation
%                   'z' z-transformation
%                   'psc' %-transformation
%        .tsm       TSM object (as alternative to SSM)
%        .tsmsrf    TSM sphere/srf object used for mapping (for triangles)
%
% Output fields:
%
%       glm         GLM object/structure
%
% Using: psctrans, varc, ztrans.

% to-do:
%        .arcorr    auto-correlation correction (default: 0)

% Version:  v1.1
% Build:    16020917
% Date:     Feb-09 2016, 5:06 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mtc') || ...
    numel(sdm) ~= 1 || ~xffisobject(sdm, true, 'sdm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get MTC object for some more checks
bc = xo.C;
szmtc = size(bc.MTCData);
numtp = szmtc(1);
szmtc(1) = [];

% further checks on SDM object
bcsdm = sdm.C;
sdmfile = sdm.F;
if isempty(sdmfile)
    sdmfile = 'Interactive';
end
if bcsdm.NrOfDataPoints ~= numtp
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% options
if nargin < 3 || ~isstruct(glmopts) || numel(glmopts) ~= 1
    glmopts = struct;
end
xsm = [];
xsf = '<none>';
if ~isfield(glmopts, 'ssm') || numel(glmopts.ssm) ~= 1 || ~xffisobject(glmopts.ssm, true, 'ssm')
    glmopts.ssm = [];
else
    xsm = glmopts.ssm;
    xsmc = xsm.C;
    xsf = xsm.F;
    if xsmc.NrOfSourceVertices ~= szmtc
        error('neuroelf:xff:badArgument', 'SSM and MTC don''t match.');
    end
end
if ~isfield(glmopts, 'tctrans') || ~ischar(glmopts.tctrans) || ...
    numel(glmopts.tctrans) ~= size(glmopts.tctrans, 2) || ...
   ~any(strcmpi(glmopts.tctrans, {'psc', 'z'}))
    glmopts.tctrans = 'none';
else
    glmopts.tctrans = lower(glmopts.tctrans);
end
switch (glmopts.tctrans)
    case 'none'
        transval = 0;
    case 'psc'
        transval = 3;
    case 'z'
        transval = 1;
end
if ~isfield(glmopts, 'threshval') || ~isa(glmopts.threshval, 'double') || ...
    numel(glmopts.threshval) ~= 1 || isinf(glmopts.threshval) || isnan(glmopts.threshval)
    glmopts.threshval = 100;
else
    glmopts.threshval = real(glmopts.threshval);
end
if ~isfield(glmopts, 'threshvol') || ~isa(glmopts.threshvol, 'double') || ...
    any(isinf(glmopts.threshvol(:)) | isnan(glmopts.threshvol(:)) | ...
    glmopts.thresvol(:) < 1 | glmopts.threshvol(:) > numtp)
    glmopts.threshvol = [fix(numtp / 2) + 1, numtp];
else
    glmopts.threshvol = fix(real(glmopts.threshvol(:)'));
end
if ~isfield(glmopts, 'tsm') || numel(glmopts.tsm) ~= 1 || ~xffisobject(glmopts.tsm, true, 'tsm') || ...
   ~isfield(glmopts, 'tsmsrf') || numel(glmopts.tsmsrf) ~= 1 || ~xffisobject(glmopts.tsmsrf, true, 'srf')
    glmopts.tsm = [];
else
    xsm = glmopts.tsm;
    xsmc = xsm.C;
    xsf = xsm.F;
    if xsmc.NrOfSourceVertices ~= szmtc
        error('neuroelf:xff:badArgument', 'TSM and MTC don''t match.');
    end
    srfc = tsmsrf.C;
    if xsmc.NrOfSourceTriangles ~= srfc.NrOfTriangles
        error('neuroelf:xff:badArgument', 'TSM and TSM-SRF don''t match.');
    end
    xsmc.TriangleVertex = srfc.TriangleVertex(xsmc.SourceOfTarget, :);
end

% resample timecourse
if ~isempty(xsm)

    % for SSM
    if xffisobject(xsm, true, 'ssm')

        % simply look up sources of target vertices
        bc.MTCData = bc.MTCData(:, xsmc.SourceOfTarget(:)');

    % for TSM
    else

        % create target matrix
        mtcd = single(zeros(numtp, 1));
        mtcd(end, xsmc.NrOfTargetVertices) = 0;

        % get necessary arguments
        tv = xsmc.TriangleVertex;
        tl = xsmc.TriangleEdgeLengths';

        % iterate over time points
        for tc = 1:numtp

            % get data of three triangle vertices
            trid = reshape(bc.MTCData(tc, tv), size(tv))';
            mtcd(tc, :) = trid(1, :) + tl(1, :) .* (trid(2, :) - trid(1, :)) + ...
                tl(2, :) .* (trid(3, :) - trid(1, :));
        end

        % put into bc
        bc.MTCData = mtcd;
    end

    % update bc
    bc.NrOfVertices = size(bc.MTCData, 2);
end

% define mask
if ~isempty(glmopts.threshvol)
    tvol = glmopts.threshvol;
    tval = glmopts.threshval;
    maski = (mean(bc.MTCData(tvol, :)) > tval);
else
    maski = true([1, szmtc]);
end

% get iXX
[sdmret{1:2}] = sdm_CalcBetas(sdm, rand(numtp, 1));
iXX = sdmret{2};

% check constant confound
sdmc = sdm.C;
psdm = bcsdm.SDMMatrix;
numpred = size(iXX, 2);
if numpred > size(psdm, 2)
    psdm(:, end + 1) = 1;
    sdmc.PredictorNames{end + 1} = 'Constant';
    sdmc.PredictorColors(end + 1, :) = 255;
end

% build GLM structure
glm = xff('new:glm');
glmc = glm.C;

% build complete GLM structure
glmc.ProjectType = 2;
glmc.ProjectTypeRFX = 0;
glmc.FileVersion = 4;
glmc.NrOfSubjects = 1;
glmc.NrOfSubjectPredictors = numpred;
glmc.NrOfTimePoints = size(bc.MTCData, 1);
glmc.NrOfPredictors = numpred;
glmc.NrOfConfounds = 1;
glmc.NrOfStudies = 1;
glmc.NrOfStudiesWithConfounds = 1;
glmc.NrOfConfoundsPerStudy = glmc.NrOfConfounds;
glmc.SeparatePredictors = 0;
glmc.TransformationType = transval;
glmc.SerialCorrelation = 0;
glmc.MeanARPre = 0;
glmc.MeanARPost = 0;
glmc.NrOfVertices = szmtc;
glmc.CortexBasedStatistics = 0;
glmc.NrOfVoxelsForBonfCorrection = sum(maski(:));
glmc.CortexBasedStatisticsMaskFile = '';
glmc.Study = struct;
glmc.Study.NrOfTimePoints = numtp;
glmc.Study.NameOfAnalyzedFile = xo.F;
glmc.Study.NameOfSSMFile = xsf;
glmc.Study.NameOfSDMFile = sdmfile;
glmc.Predictor = struct;
for pc = 1:numpred
    glmc.Predictor(pc).Name1 = sprintf('Predictor: %d', pc);
    glmc.Predictor(pc).Name2 = sdmc.PredictorNames{pc};
    glmc.Predictor(pc).RGB = [255, 255, 255; zeros(3, 3)];
end
glmc.Predictor(end).Name1 = 'Constant confound';
glmc.DesignMatrix = psdm;
glmc.GLMData.MultipleRegressionR = single(zeros(szmtc, 1));
glmc.GLMData.MCorrSS = single(zeros(szmtc, 1));
glmc.GLMData.BetaMaps = single(zeros(szmtc, numpred));
glmc.GLMData.XY = single(zeros(szmtc, numpred));
glmc.GLMData.TimeCourseMean = single(zeros(szmtc, 1));

% only work within mask!
maskx = find(maski);

% drop field
keep = 1:size(bc.MTCData, 1);
drop = [];
if isfield(bc.RunTimeVars, 'Discard') && isa(bc.RunTimeVars.Discard, 'double')
    drop = bc.RunTimeVars.Discard(:)';
    if ~any(isinf(drop) | isnan(drop) | drop < 1)
        drop = unique(round(max(1, min(numel(keep), drop))));
    end
end
if ~isempty(drop)
    keep(drop) = [];
end

% get masked data
switch transval
    case 0
        data = single(bc.MTCData(:, maskx));
    case 1
        data = ne_methods.ztrans(bc.MTCData(:, maskx), 1, keep);
    case 3
        data = ne_methods.psctrans(bc.MTCData(:, maskx), 1, keep);
end

% calculate betas and PTC
[beta, isdm, ptc] = sdm_CalcBetas(sdm, data, struct('drop', drop));
glmc.iXX = isdm;

% calculate MC R and total SS and set fields
vartc = ne_methods.varc(data);
sdmc.SDMMatrix(:, end+1:size(beta, 2)) = 1;
xy = (sdmc.SDMMatrix' * data)';

% create values
glmc.GLMData.MultipleRegressionR(maskx) = std(ptc) ./ sqrt(vartc);
glmc.GLMData.MCorrSS(maskx) = (numtp - 1) .* vartc;
for betac = 1:size(beta, 2)
    glmc.GLMData.BetaMaps((betac - 1) .* szmtc + maskx) = beta(:, betac);
    glmc.GLMData.XY((betac - 1) .* szmtc + maskx) = xy(:, betac);
end
glmc.GLMData.TimeCourseMean(maskx) = (1 / numtp) .* sum(data);

% put glm data into object
glm.C = glmc;

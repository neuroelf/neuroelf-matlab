function glm = vtc_CreateGLM(xo, sdm, glmopts)
% VTC::CreateGLM  - create a single study GLM
%
% FORMAT:       glm = vtc.CreateGLM(sdm [, glmopts])
%
% Input fields:
%
%       sdm         SDM object with same number of time points
%       glmopts     optional fields for GLM calculation
%        .threshval threshold value to compute a GLM (default: 100)
%        .threshvol threshold vols (mean value, default: [fix(N/2)+1,N])
%        .tctrans   time course transformation before estimation
%                   'none' or '' no transformation (default)
%                   'z' z-transformation
%                   'psc' %-transformation
%
% Output fields:
%
%       glm         GLM object/structure
%
% Note: instead of an SDM object, also a PRT object can be given, in
%       which case all options of PRT::CreateSDM are supported as well
%
% Using: psctrans, varc, ztrans.

% to-do:
%        .arcorr    auto-correlation correction (default: 0)

% Version:  v1.1
% Build:    16021321
% Date:     Feb-13 2016, 9:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
varc = ne_methods.varc;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ...
    numel(sdm) ~= 1 || ~xffisobject(sdm, true, {'prt', 'sdm'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isstruct(glmopts) || numel(glmopts) ~= 1
    glmopts = struct;
end

% get VTC object for some more checks
bc = xo.C;
szvtc = size(bc.VTCData);
numtp = szvtc(1);
szvtc(1) = [];
nlvtc = prod(szvtc);

% not from SDM
if ~xffisobject(sdm, true, 'sdm')
    fromsdm = false;
    try
        glmopts.nvol = size(bc.VTCData, 1);
        glmopts.prtr = bc.TR;
        prt = sdm;
        sdm = prt_CreateSDM(prt, glmopts);
        if ~xffisobject(sdm, true, 'sdm')
            error('neuroelf:xff:internalError', 'Couldn''t create SDM from PRT.');
        end
    catch xfferror
        rethrow(xfferror);
    end
else
    fromsdm = true;
end

% further checks on SDM object
bcsdm = sdm.C;
sdmfile = sdm.F;
if isempty(sdmfile)
    sdmfile = 'Interactive';
end
if bcsdm.NrOfDataPoints ~= numtp
    if ~fromsdm
        delete(sdm);
    end
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if isfield(glmopts, 'robust') && islogical(glmopts.robust) && numel(glmopts.robust) == 1
    robust = glmopts.robust;
else
    robust = false;
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
        psctrans = ne_methods.psctrans;
    case 'z'
        transval = 1;
        ztrans = ne_methods.ztrans;
end
if ~isfield(glmopts, 'threshval') || ~isa(glmopts.threshval, 'double') || ...
    numel(glmopts.threshval) ~= 1 || isinf(glmopts.threshval) || isnan(glmopts.threshval)
    glmopts.threshval = 100;
else
    glmopts.threshval = real(glmopts.threshval);
end
if ~isfield(glmopts, 'threshvol') || ~isa(glmopts.threshvol, 'double') || ...
    any(isinf(glmopts.threshvol(:)) | isnan(glmopts.threshvol(:)) | ...
    glmopts.threshvol(:) < 1 | glmopts.threshvol(:) > numtp)
    glmopts.threshvol = [fix(numtp / 2) + 1, numtp];
else
    glmopts.threshvol = fix(real(glmopts.threshvol(:)'));
end

% define mask
maski = false([1, szvtc]);
if ~isempty(glmopts.threshvol)
    tvol = glmopts.threshvol;
    tval = glmopts.threshval;
    for slc = 1:szvtc(3)
        sldat = bc.VTCData(:, :, :, slc);
        maski(1, :, :, slc) = mean(sldat(tvol, :, :), 1) > tval;
    end
else
    for slc = 1:szvtc(3)
        sldat = bc.VTCData(:, :, :, slc);
        maski(1, :, :, slc) = sum(abs(sldat), 1) > 0;
    end
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
glmc.ProjectType = 1;
glmc.ProjectTypeRFX = 0;
glmc.FileVersion = 4;
glmc.NrOfSubjects = 1;
glmc.NrOfSubjectPredictors = numpred;
glmc.NrOfTimePoints = size(bc.VTCData, 1);
glmc.NrOfPredictors = numpred;
glmc.NrOfConfounds = 1;
glmc.NrOfStudies = 1;
glmc.NrOfStudiesWithConfounds = 1;
glmc.NrOfConfoundsPerStudy = glmc.NrOfConfounds;
glmc.SeparatePredictors = 0;
glmc.TransformationType = transval;
glmc.Resolution = bc.Resolution;
glmc.SerialCorrelation = 0;
glmc.MeanARPre = 0;
glmc.MeanARPost = 0;
glmc.XStart = bc.XStart;
glmc.XEnd = bc.XEnd;
glmc.YStart = bc.YStart;
glmc.YEnd = bc.YEnd;
glmc.ZStart = bc.ZStart;
glmc.ZEnd = bc.ZEnd;
glmc.CortexBasedStatistics = 0;
glmc.NrOfVoxelsForBonfCorrection = sum(maski(:));
glmc.CortexBasedStatisticsMaskFile = '';
glmc.Study = struct;
glmc.Study.NrOfTimePoints = numtp;
glmc.Study.NameOfAnalyzedFile = xo.F;
glmc.Study.NameOfSDMFile = sdmfile;
glmc.Predictor = struct;
for pc = 1:numpred
    glmc.Predictor(pc).Name1 = sprintf('Predictor: %d', pc);
    glmc.Predictor(pc).Name2 = sdmc.PredictorNames{pc};
    glmc.Predictor(pc).RGB = [255, 255, 255; zeros(3, 3)];
end
glmc.Predictor(end).Name1 = 'Constant confound';
glmc.DesignMatrix = psdm;
glmc.iXX = iXX;
glmc.GLMData.MultipleRegressionR = single(zeros(szvtc));
glmc.GLMData.MCorrSS = single(zeros(szvtc));
glmc.GLMData.BetaMaps = single(zeros([szvtc, numpred]));
glmc.GLMData.XY = single(zeros([szvtc, numpred]));
glmc.GLMData.TimeCourseMean = single(zeros(szvtc));

% only work within mask!
maskxa = find(maski);

% drop field
if isfield(bc.RunTimeVars, 'Discard')
    drop = bc.RunTimeVars.Discard;
else
    drop = [];
end

% iterate over voxels as necessary
numvox = floor((2 ^ 24) / (numpred * size(bc.VTCData, 1)));
for vc = 1:numvox:numel(maskxa)

    % get sub-mask
    maskx = maskxa(vc:min(vc + numvox - 1, numel(maskxa)));

    % get masked data
    switch transval
        case 0
            data = single(bc.VTCData(:, maskx));
        case 1
            data = ztrans(bc.VTCData(:, maskx));
        case 3
            data = psctrans(bc.VTCData(:, maskx));
    end

    % calculate betas and PTC
    [beta, isdm, ptc] = ...
        sdm_CalcBetas(sdm, data, struct('drop', drop, 'pbar', [], 'robust', robust));

    % calculate MC R and total SS and set fields
    vartc = varc(data);
    sdmc.SDMMatrix(:, end+1:size(beta, 2)) = 1;
    xy = (sdmc.SDMMatrix' * data)';

    % create values
    glmc.GLMData.MultipleRegressionR(maskx) = std(ptc) ./ sqrt(vartc);
    glmc.GLMData.MCorrSS(maskx) = (numtp - 1) .* vartc;
    for betac = 1:size(beta, 2)
        glmc.GLMData.BetaMaps((betac - 1) .* nlvtc + maskx) = beta(:, betac);
        glmc.GLMData.XY((betac - 1) .* nlvtc + maskx) = xy(:, betac);
    end
    glmc.GLMData.TimeCourseMean(maskx) = (1 / numtp) .* sum(data);
end

% store additional options
if robust
    glmc.RunTimeVars.Robust = true;
end

% add further fields
if isfield(bc.RunTimeVars, 'TrfPlus')
    glmc.AutoSave = true;
    glmc.RunTimeVars.TrfPlus = bc.RunTimeVars.TrfPlus;
end

% put glm data in global storage
glm.C = glmc;

% remove temporary SDM
if ~fromsdm
    delete(sdm);
end

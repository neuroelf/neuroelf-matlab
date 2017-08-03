function xo3 = glm_JoinFFX(xo, xo2)
% GLM::JoinFFX  - join two fixed effects GLMs
%
% FORMAT:       combined = glm1.JoinFFX(glm2);
%
% Input fields:
%
%       glm2        second FFX GLM to be added to glm1
%
% Joins the Fixed Effects GLM results of glm1 and glm2 to one
% single structure.
%
% Using: catstruct, findfirst, multimatch.

% Version:  v1.1
% Build:    17021613
% Date:     Feb-16 2017, 1:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, 2017, Jochen Weber
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
catstruct  = ne_methods.catstruct;
findfirst  = ne_methods.findfirst;
multimatch = ne_methods.multimatch;

% argument check
if nargin < 2 || numel(xo) ~= 1 || numel(xo2) ~= 1 || ...
   ~xffisobject(xo, true, 'glm') || ~xffisobject(xo2, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc1 = xo.C;
bc2 = xo2.C;
if bc1.ProjectTypeRFX ~= 0 || bc2.ProjectTypeRFX ~= 0 || ~any(bc1.SeparatePredictors == [1, 2])
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% check type and sizes of GLMs
if ~isempty(bc1.GLMData.RFXGlobalMap) || ~isempty(bc2.GLMData.RFXGlobalMap) || ...
    bc1.ProjectType ~= bc2.ProjectType || bc1.FileVersion ~= bc2.FileVersion || ...
    bc1.SeparatePredictors ~= bc2.SeparatePredictors || ...
    bc1.TransformationType ~= bc2.TransformationType || ...
    bc1.Resolution ~= bc2.Resolution || ...
    bc1.XStart ~= bc2.XStart || bc1.XEnd ~= bc2.XEnd || ...
    bc1.YStart ~= bc2.YStart || bc1.YEnd ~= bc2.YEnd || ...
    bc1.ZStart ~= bc2.ZStart || bc1.ZEnd ~= bc2.ZEnd || ...
   ~isequal(bc1.RunTimeVars.TrfPlus, bc2.RunTimeVars.TrfPlus)
    error('neuroelf:xff:badObject', ...
        'Invalid object(s) given in call. Crucial fields mismatch.');
end

% disallow identical subjects (for SPSB GLMs)
SPS = bc1.SeparatePredictors;
if SPS == 2 && any(multimatch(glm_Subjects(xo), glm_Subjects(xo2)) > 0)
    error('neuroelf:xff:badArgument', 'GLMs must contain different subjects to combine.');
end

% allow to use transio data (only if both are transio)
bmtio = false;
xytio = false;
if istransio(bc1.GLMData.BetaMaps) || istransio(bc2.GLMData.BetaMaps)
    if ~istransio(bc1.GLMData.BetaMaps)
        try
            bc2.GLMData.BetaMaps = resolve(bc2.GLMData.BetaMaps);
        catch xfferror
            error('neuroelf:xff:outOfMemory', ...
                'Error resolving transio for BetaMaps, both GLMs must have stored data to work.');
        end
    elseif ~istransio(bc2.GLMData.BetaMaps)
        try
            bc1.GLMData.BetaMaps = resolve(bc1.GLMData.BetaMaps);
        catch xfferror
            error('neuroelf:xff:outOfMemory', ...
                'Error resolving transio for BetaMaps, both GLMs must have stored data to work.');
        end
    else
        bmtio = true;
        tbms1 = struct(bc1.GLMData.BetaMaps);
        tbmf1 = tbms1.FileName;
        nmaps = tbms1.DataDims(end);
        if tbms1.LittleND
            tbme1 = 'ieee-le';
        else
            tbme1 = 'ieee-be';
        end
        if ~iscell(tbmf1)
            tbmf1 = repmat({tbmf1}, nmaps, 1);
            smap = prod(tbms1.DataDims(1:end-1)) * tbms1.TypeSize;
            tbms1.IOOffset = tbms1.IOOffset + (0:smap:(nmaps-1)*smap);
            bc1.GLMData.BetaMaps = transio(tbmf1, tbme1, tbms1.DataType, ...
                tbms1.IOOffset, tbms1.DataDims);
        end
        tbms2 = struct(bc2.GLMData.BetaMaps);
        tbmf2 = tbms2.FileName;
        nmaps = tbms2.DataDims(end);
        if ~iscell(tbmf2)
            tbmf2 = repmat({tbmf2}, nmaps, 1);
            if tbms2.LittleND
                tbme2 = 'ieee-le';
            else
                tbme2 = 'ieee-be';
            end
            smap = prod(tbms2.DataDims(1:end-1)) * tbms2.TypeSize;
            tbms2.IOOffset = tbms2.IOOffset + (0:smap:(nmaps-1)*smap);
            bc2.GLMData.BetaMaps = transio(tbmf2, tbme2, tbms2.DataType, ...
                tbms2.IOOffset, tbms2.DataDims);
        end
        if tbms1.LittleND ~= tbms2.LittleND || ...
           ~strcmp(tbms1.DataType, tbms2.DataType) || ...
           ~isequal(tbms1.DataDims(1:end-1), tbms2.DataDims(1:end-1))
            error('neuroelf:xff:badArgument', ...
                'Invalid transio combination (little endian, datatype, size).');
        end
    end
end
if istransio(bc1.GLMData.XY) || istransio(bc2.GLMData.XY)
    if ~istransio(bc1.GLMData.XY)
        try
            bc2.GLMData.XY = resolve(bc2.GLMData.XY);
        catch xfferror
            error('neuroelf:xff:outOfMemory', ...
                'Error resolving transio for XY, both GLMs must have stored data to work.');
        end
    elseif ~istransio(bc2.GLMData.XY)
        try
            bc1.GLMData.XY = resolve(bc1.GLMData.XY);
        catch xfferror
            error('neuroelf:xff:outOfMemory', ...
                'Error resolving transio for XY, both GLMs must have stored data to work.');
        end
    else
        xytio = true;
        txys1 = struct(bc1.GLMData.XY);
        txyf1 = txys1.FileName;
        nmaps = txys1.DataDims(end);
        if txys1.LittleND
            txye1 = 'ieee-le';
        else
            txye1 = 'ieee-be';
        end
        if ~iscell(txyf1)
            txyf1 = repmat({txyf1}, nmaps, 1);
            smap = prod(txys1.DataDims(1:end-1)) * txys1.TypeSize;
            txys1.IOOffset = txys1.IOOffset + (0:smap:(nmaps-1)*smap);
            bc1.GLMData.XY = transio(txyf1, txye1, txys1.DataType, ...
                txys1.IOOffset, txys1.DataDims);
        end
        txys2 = struct(bc2.GLMData.XY);
        txyf2 = txys2.FileName;
        nmaps = txys2.DataDims(end);
        if txys2.LittleND
            txye2 = 'ieee-le';
        else
            txye2 = 'ieee-be';
        end
        if ~iscell(txyf2)
            txyf2 = repmat({txyf2}, nmaps, 1);
            smap = prod(txys2.DataDims(1:end-1)) * txys2.TypeSize;
            txys2.IOOffset = txys2.IOOffset + (0:smap:(nmaps-1)*smap);
            bc2.GLMData.XY = transio(txyf2, txye2, txys2.DataType, ...
                txys2.IOOffset, txys2.DataDims);
        end
        if txys1.LittleND ~= txys2.LittleND || ...
           ~strcmp(txys1.DataType, txys2.DataType) || ...
           ~isequal(txys1.DataDims(1:end-1), txys2.DataDims(1:end-1))
            error('neuroelf:xff:badArgument', ...
                'Invalid transio combination (little endian, datatype, size).');
        end
    end
end

% build new file
ct = sqrt(eps);
DM1 = bc1.DesignMatrix;
DM1(abs(DM1) < ct) = 0;
if (sum(DM1(:) == 0) / numel(DM1)) > 0.875
    DM1 = sparse(DM1);
end
DM2 = bc2.DesignMatrix;
DM2(abs(DM2) < ct) = 0;
if (sum(DM2(:) == 0) / numel(DM2)) > 0.875
    DM2 = sparse(DM2);
end
NP1 = bc1.NrOfPredictors;
NP2 = bc2.NrOfPredictors;
NPS = NP1 + NP2;
NS1 = bc1.NrOfStudies;
NS2 = bc2.NrOfStudies;
if SPS == 1
    if NS1 == 1
        NC1 = bc1.NrOfConfounds;
    else
        pnames = {bc1.Predictor.Name2};
        s1names = ~cellfun('isempty', regexpi(pnames(:), '^study\s+1\:'));
        NC1 = NP1 - findfirst(~s1names(1:end-1) & s1names(2:end));
    end
    if NS2 == 1
        NC2 = bc2.NrOfConfounds;
    else
        pnames = {bc2.Predictor.Name2};
        s1names = ~cellfun('isempty', regexpi(pnames(:), '^study\s+1\:'));
        NC2 = NP2 - findfirst(~s1names(1:end-1) & s1names(2:end));
    end
else
    pnames = {bc1.Predictor.Name2};
    NC1 = sum(cellfun('isempty', regexpi(pnames(:), '^subject\s+')));
    pnames = {bc2.Predictor.Name2};
    NC2 = sum(cellfun('isempty', regexpi(pnames(:), '^subject\s+')));
end
NSS = NS1 + NS2;
NT1 = bc1.NrOfTimePoints;
NT2 = bc2.NrOfTimePoints;
NTP = NT1 + NT2;
VA1 = (bc1.GLMData.MCorrSS > 0);
VA2 = (bc2.GLMData.MCorrSS > 0);
VAC = (VA1 & VA2);
xo3 = aft_CopyObject(xo);
xo3.F = '';
bc3 = xo3.C;
bc3.NrOfTimePoints = NTP;
bc3.NrOfPredictors = NPS;
bc3.NrOfConfounds = NC1 + NC2;
bc3.NrOfStudies = NSS;
bc3.NrOfStudiesWithConfounds = bc1.NrOfStudiesWithConfounds + bc2.NrOfStudiesWithConfounds;
bc3.NrOfConfoundsPerStudy = [bc1.NrOfConfoundsPerStudy, bc2.NrOfConfoundsPerStudy];
bc3.NrOfVoxelsForBonfCorrection = sum(VAC(:));
bc3.Study(1:end) = [];
bc3.Predictor = catstruct( ...
    bc1.Predictor(1:(NP1-NC1)), bc2.Predictor(1:(NP2-NC2)), ...
    bc1.Predictor((end+1-NC1):end), bc2.Predictor((end+1-NC2):end));

% concatenate design matrices
try
    bc3.DesignMatrix = [ ...
        DM1(:, 1:end-NC1), ...
        zeros(size(DM1, 1), size(DM2, 2) - NC2), ...
        DM1(:, end+1-NC1:end), ...
        zeros(size(DM1, 1), NC2); ...
        ...
        zeros(size(DM2, 1), size(DM1, 2) - NC1), ...
        DM2(:, 1:end-NC2), ...
        zeros(size(DM2, 1), NC1), ...
        DM2(:, end+1-NC2:end)];
    bc3.iXX = inv(bc3.DesignMatrix' * bc3.DesignMatrix);
catch xfferror
    error('neuroelf:xff:badObjects', ...
        'Error pseudo-inverting concatenated design matrix: %s.', xfferror.message);
end

% copy studys
bc3.Study = catstruct(bc1.Study(:), bc2.Study(:))';

% copy predictors
tpc = NP1 + 1 - NC1;
for pc = 1:(NP2 - NC2)
    bc3.Predictor(tpc) = bc2.Predictor(pc);
    bc3.Predictor(tpc).Name1 = sprintf('Predictor: %d', tpc);
    if SPS == 1
        prname = bc3.Predictor(tpc).Name2;
        [pm_m{1:3}] = regexpi(prname, '^study\s+(\d+)\:\s+(.*)$');
        if isempty(pm_m{1}) || isempty(pm_m{3}) || ~all(size(pm_m{3}{1}) == 2)
            error('neuroelf:xff:internalError', 'Error extracting study predictor name.');
        end
        pm_m = pm_m{3};
        bc3.Predictor(tpc).Name2 = sprintf('Study %d: %s', ...
            str2double(prname(pm_m{1}(1, 1):pm_m{1}(1, 2))) + NS1, ...
            prname(pm_m{1}(2, 1):pm_m{1}(2, 2)));
    end
    tpc = tpc + 1;
end
for pc = (NP1 + 1 - NC1):NP1
    bc3.Predictor(tpc) = bc1.Predictor(pc);
    bc3.Predictor(tpc).Name1 = sprintf('Predictor: %d', tpc);
    tpc = tpc + 1;
end
for pc = (NP2 + 1 - NC2):NP2
    bc3.Predictor(tpc) = bc2.Predictor(pc);
    bc3.Predictor(tpc).Name1 = sprintf('Predictor: %d', tpc);
    prname = bc3.Predictor(tpc).Name2;
    [pm_m{1:3}] = regexpi(prname, '^study\s+(\d+)\:\s+(.*)$');
    if isempty(pm_m{1}) || isempty(pm_m{3}) || ~all(size(pm_m{3}{1}) == 2)
        error('neuroelf:xff:internalError', 'Error extracting study predictor name.');
    end
    pm_m = pm_m{3};
    bc3.Predictor(tpc).Name2 = sprintf('Study %d: %s', ...
        str2double(prname(pm_m{1}(1, 1):pm_m{1}(1, 2))) + NS1, ...
        prname(pm_m{1}(2, 1):pm_m{1}(2, 2)));
    tpc = tpc + 1;
end

% create GLMData fields
bc3.GLMData.MCorrSS = bc1.GLMData.MCorrSS + bc2.GLMData.MCorrSS;
bc3.GLMData.MultipleRegressionR = ...
   ((bc1.GLMData.MultipleRegressionR .* bc1.GLMData.MCorrSS) + ...
    (bc2.GLMData.MultipleRegressionR .* bc2.GLMData.MCorrSS)) ./ bc3.GLMData.MCorrSS;
bc3.GLMData.MultipleRegressionR(~VAC) = 0;
maporder = [1:1:(NP1 - NC1), -1:-1:-(NP2 - NC2), (NP1 + 1 - NC1):1:NP1, -(NP2 + 1 - NC2):-1:-NP2];
maps1 = find(maporder > 0);
maps2 = find(maporder < 0);
betasz = size(bc3.GLMData.BetaMaps);
betasz(end) = [];
numvox = prod(betasz);
if bmtio
    tbmf = [tbmf1(:); tbmf2(:)];
    tbmo = zeros(1, numel(tbmf));
    tbmf(maps1) = tbmf1(maporder(maps1));
    tbmo(maps1) = tbms1.IOOffset(maporder(maps1));
    tbmf(maps2) = tbmf2(-maporder(maps2));
    tbmo(maps2) = tbms2.IOOffset(-maporder(maps2));
    tbmsz = tbms1.DataDims;
    tbmsz(end) = tbmsz(end) + tbms2.DataDims(end);
    bc3.GLMData.BetaMaps = transio(tbmf, tbme1, tbms1.DataType, tbmo, tbmsz);
else
    bc1.GLMData.BetaMaps = reshape(bc1.GLMData.BetaMaps, numvox, NP1);
    bc2.GLMData.BetaMaps = reshape(bc2.GLMData.BetaMaps, numvox, NP2);
    bc3.GLMData.BetaMaps = single(0);
    bc3.GLMData.BetaMaps(numvox, NPS) = 0;
    bc3.GLMData.BetaMaps(:, maps1) = bc1.GLMData.BetaMaps(:, maporder(maps1));
    bc3.GLMData.BetaMaps(:, maps2) = bc2.GLMData.BetaMaps(:, -maporder(maps2));
    bc3.GLMData.BetaMaps(~VAC(:), :) = 0;
    bc3.GLMData.BetaMaps = reshape(bc3.GLMData.BetaMaps, [betasz, NPS]);
end
if xytio
    txyf = [txyf1(:); txyf2(:)];
    txyo = zeros(1, numel(txyf));
    txyf(maps1) = txyf1(maporder(maps1));
    txyo(maps1) = txys1.IOOffset(maporder(maps1));
    txyf(maps2) = txyf2(-maporder(maps2));
    txyo(maps2) = txys2.IOOffset(-maporder(maps2));
    txysz = txys1.DataDims;
    txysz(end) = txysz(end) + txys2.DataDims(end);
    bc3.GLMData.XY = transio(txyf, txye1, txys1.DataType, txyo, txysz);
else
    bc1.GLMData.XY = reshape(bc1.GLMData.XY, numvox, NP1);
    bc2.GLMData.XY = reshape(bc2.GLMData.XY, numvox, NP2);
    bc3.GLMData.XY = single(0);
    bc3.GLMData.XY(numvox, NPS) = 0;
    bc3.GLMData.XY(:, maps1) = bc1.GLMData.XY(:, maporder(maps1));
    bc3.GLMData.XY(:, maps2) = bc2.GLMData.XY(:, -maporder(maps2));
    bc3.GLMData.XY(~VAC(:), :) = 0;
    bc3.GLMData.XY = reshape(bc3.GLMData.XY, [betasz, NPS]);
end
bc3.GLMData.TimeCourseMean = (1 / NTP) .* ...
    (NT1 .* bc1.GLMData.TimeCourseMean + NT2 .* bc2.GLMData.TimeCourseMean);
bc3.GLMData.TimeCourseMean(~VAC) = 0;

% clean up fields
if issparse(bc3.DesignMatrix)
    bc3.DesignMatrix = full(bc3.DesignMatrix);
end
if issparse(bc3.iXX)
    bc3.iXX = full(bc3.iXX);
end
bc3.GLMData.MultipleRegressionR(abs(bc3.GLMData.MultipleRegressionR) < ct) = 0;
bc3.GLMData.MCorrSS(abs(bc3.GLMData.MCorrSS) < ct) = 0;
bc3.GLMData.TimeCourseMean(abs(bc3.GLMData.TimeCourseMean) < ct) = 0;

% deal with RunTimeVars
if isfield(bc1.RunTimeVars, 'Groups')
    bc3.RunTimeVars = rmfield(bc3.RunTimeVars, 'Groups');
end
if isfield(bc1.RunTimeVars, 'MotionParameters') && ...
    iscell(bc1.RunTimeVars.MotionParameters) && ...
    numel(bc1.RunTimeVars.MotionParameters) == NS1 && ...
    isfield(bc2.RunTimeVars, 'MotionParameters') && ...
    iscell(bc2.RunTimeVars.MotionParameters) && ...
    numel(bc2.RunTimeVars.MotionParameters) == NS2
    bc3.RunTimeVars.MotionParameters = [ ...
        bc1.RunTimeVars.MotionParameters(:); bc2.RunTimeVars.MotionParameters(:)];
    bc3.RunTimeVars.AutoSave = true;
elseif isfield(bc1.RunTimeVars, 'MotionParameters')
    bc3.RunTimeVars = rmfield(bc3.RunTimeVars, 'MotionParameters');
end

% put back
xo3.C = bc3;

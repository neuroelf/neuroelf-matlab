function map = glm_FFX_tMap(xo, c, mapopts)
% GLM::FFX_tMap  - calculate a t contrast map
%
% FORMAT:       map = glm.FFX_tMap([c, mapopts])
%
% Input fields:
%
%       c           PxC contrast vector (default: full model and main eff)
%       mapopts     structure with optional fields
%        .arlag     use BV-style ARLag maps (default: true)
%        .arlagv2w  use V/W matrix operation for SE term (default: false)
%        .interp    mesh-based interpolation (default: true)
%        .names     1xC cell array with names for contrasts
%        .srf       surface file, required for interpolation
%
% Output fields:
%
%       map         VMP/SMP object with C maps
%
% Using: glmtstat, lsqueeze, multimatch, ne_sqrtm, newnatresvmp.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2012, 2014, 2016, Jochen Weber
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
glmtstat   = ne_methods.glmtstat;
lsqueeze   = ne_methods.lsqueeze;
multimatch = ne_methods.multimatch;
ne_sqrtm   = ne_methods.ne_sqrtm;
sdist      = ne_methods.sdist;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
rtv = bc.RunTimeVars;
glmfile = xo.F;
glmid = xo.L;
if isempty(glmfile)
    glmfile = glmid;
end
if bc.ProjectTypeRFX > 0
    error('neuroelf:xff:badArgument', 'FFX contrasts available only for FFX-GLMs.');
end
if bc.SeparatePredictors == 2
    ffxsubs = glm_Subjects(xo);
end
ffxpred = {bc.Predictor.Name2};
ffxpred = ffxpred(:);
ffxspred = glm_SubjectPredictors(xo);
numtp = bc.NrOfTimePoints;
numpred = bc.NrOfPredictors;
nval = numtp - numpred;
szmap = size(bc.GLMData.MCorrSS);
if nargin < 2 || ~isa(c, 'double') || isempty(c) || size(c, 1) > numpred || ...
    any(isinf(c(:)) | isnan(c(:)))
    if numpred > 2
        c = [[eye(numpred - 1);zeros(1, numpred - 1)], [ones(numpred - 1, 1);0]];
    else
        c = [1; 0];
    end
elseif numel(c) == size(c, 2)
    if numel(c) > numpred
        c((numpred + 1):end) = [];
    end
    c = real(c');
else
    c = real(c);
end
nummaps = size(c, 2);
if size(c, 1) < numpred
    c = [c; zeros(numpred - size(c, 1), size(c, 2))];
end
if nargin < 3 || ~isstruct(mapopts) || numel(mapopts) ~= 1
    mapopts = struct;
end
if ~isfield(mapopts, 'arlag') || ~islogical(mapopts.arlag) || numel(mapopts.arlag) ~= 1
    mapopts.arlag = ~isempty(bc.GLMData.ARLag);
end
if ~isfield(mapopts, 'arlagv2w') || ~islogical(mapopts.arlagv2w) || numel(mapopts.arlagv2w) ~= 1
    mapopts.arlagv2w = false;
end
ipo = true;
if isfield(mapopts, 'interp') && ~isempty(mapopts.interp) && ...
   (isnumeric(mapopts.interp) || islogical(mapopts.interp))
    ipo = ipo && mapopts.interp(1);
end
if ~isfield(mapopts, 'names') || ~iscell(mapopts.names) || numel(mapopts.names) ~= size(c, 2)
    mapopts.names = cell(1, size(c, 2));
    for cc = 1:size(c, 2)
        if sum(c(:, cc) ~= 0) == 1 && sum(c(:, cc)) > 0
            mapopts.names(cc) = ffxpred(c(:, cc) ~= 0);
        else
            mapopts.names{cc} = sprintf('Contrast:%s', sprintf(' %d', c(:, cc)'));
        end
    end
else
    for cc = 1:size(c, 2)
        if ~ischar(mapopts.names{cc}) || numel(mapopts.names{cc}) ~= size(mapopts.names{cc}, 2) || ...
            isempty(mapopts.names{cc})
            if sum(c(:, cc) ~= 0) == 1 && sum(c(:, cc)) > 0
                mapopts.names(cc) = ffxpred(c(:, cc) ~= 0);
            else
                mapopts.names{cc} = sprintf('Contrast:%s', sprintf(' %d', c(:, cc)'));
            end
        end
    end
end

% get study lengths
stp = cat(1, bc.Study.NrOfTimePoints);

% for serial correlation correction
if mapopts.arlag && numel(bc.Study) > 1

    % and compute study boundaries
    r1r = cumsum(stp(1:end-1, 1));
    r2r = sort([r1r(:) - 1; r1r(:)]);

    % load helper file
    if bc.SerialCorrelation == 2
        %load
    end

% for a single study
else

    % no boundaries!
    r1r = [];
end

% set additional data
artv = struct('SourceGLM', glmfile, 'SourceGLMID', glmid, 'Contrast', [], ...
    'Covariates', [], 'Groups', {{}}, 'MeanRem', false, ...
    'RFXGLM', (bc.SeparatePredictors == 2), 'Robust', false, 'SubPreds', {ffxspred});

% generate map according to type
iXX = double(bc.iXX);
switch (bc.ProjectType)

    % currently unsupported types
    case 0
        error('neuroelf:xff:notYetSupported', 'FMR/MAP generation not yet supported.');

    % VTCs
    case 1

        % create VMP
        bbox = aft_BoundingBox(xo);
        map = ne_methods.newnatresvmp(bbox.BBox, bc.Resolution, ones(1, size(c, 2)));
        mapc = map.C;
        mapc.Resolution = bc.Resolution;
        mapc.XStart = bc.XStart;
        mapc.XEnd = bc.XEnd;
        mapc.YStart = bc.YStart;
        mapc.YEnd = bc.YEnd;
        mapc.ZStart = bc.ZStart;
        mapc.ZEnd = bc.ZEnd;

        % copy TrfPlus
        mapc.RunTimeVars.TrfPlus = rtv.TrfPlus;

        % add subject-specific normalization
        if bc.NrOfSubjects == 1 && isfield(rtv, 'SubjectSPMsn') && ...
            isstruct(rtv.SubjectSPMsn) && numel(rtv.SubjectSPMsn) == 1 && ...
            numel(fieldnames(rtv.SubjectSPMsn)) == 1
            ffxtsid = fieldnames(rtv.SubjectSPMsn);
            mapc.RunTimeVars.SPMsn = rtv.SubjectSPMsn.(ffxtsid{1});
        end
        if bc.NrOfSubjects == 1 && isfield(rtv, 'SubjectTrfPlus') && ...
            isstruct(rtv.SubjectTrfPlus) && numel(rtv.SubjectTrfPlus) == 1 && ...
            numel(fieldnames(rtv.SubjectTrfPlus)) == 1
            ffxtsid = fieldnames(rtv.SubjectTrfPlus);
            mapc.RunTimeVars.TrfPlus = rtv.TrfPlus * rtv.SubjectTrfPlus.(ffxtsid{1});
        end

        % calculation of SE map
        semap = sqrt((1 - (double(bc.GLMData.MultipleRegressionR) .^ 2)) .* ...
            double(bc.GLMData.MCorrSS) / nval);
        semap(semap == 0) = Inf;

        % serial correlation correction?
        if mapopts.arlag
            semap = sqrt(nval / (nval - (bc.SerialCorrelation * bc.NrOfStudies))) .* semap;
        end

        % generate tmaps
        for cc = 1:size(c, 2)

            % set additional info
            artv.Contrast = c(:, cc);

            % subject selection
            if bc.SeparatePredictors == 2
                subsel = sort(multimatch(unique(regexprep(ffxpred(c(:, cc) ~= 0), ...
                    '^Subject\s+(\w[^\:]+)\:.*$', '$1')), ffxsubs));
                if ~isempty(subsel)
                    if subsel(1) == 0
                        subsel(1) = [];
                    end
                    if ~isempty(subsel)
                        artv.SubSel= subsel(:);

                        % adapt contrast
                        artv.Contrast = c(~cellfun('isempty', regexp(ffxpred, ...
                            ['^Subject\s+' ffxsubs{subsel(1)}])), cc);
                    end
                end
            end

            % map options
            mapc.Map(cc).DF1 = nval - (bc.SerialCorrelation * bc.NrOfStudies);
            mapc.Map(cc).DF2 = 0;
            mapc.Map(cc).BonferroniValue = bc.NrOfVoxelsForBonfCorrection;
            mapc.Map(cc).UseRGBColor = 0;
            mapc.Map(cc).RunTimeVars = artv;
            mapc.Map(cc).LowerThreshold = -sdist('tinv', 1e-3, mapc.Map(cc).DF1);
            mapc.Map(cc).UpperThreshold = -sdist('tinv', 1e-6, mapc.Map(cc).DF1);

            % complicated scheme
            if mapopts.arlag

                % lags distribution
                a1lag = 0.05 * round(20 .* double(bc.GLMData.ARLag(:, :, :, 1)));
                a1min = 21 + round(20 * min(lsqueeze(a1lag)));
                a1max = 21 + round(20 * max(lsqueeze(a1lag)));
                if bc.SerialCorrelation > 1
                    a2lag = 0.05 * round(20 .* double(bc.GLMData.ARLag(:, :, :, 2)));
                    a2min = 21 + round(20 * min(lsqueeze(a2lag)));
                    a2max = 21 + round(20 * max(lsqueeze(a2lag)));
                    spi = (1:size(bc.DesignMatrix, 1))';
                    spin = numel(spi);
                    spi1 = ones(spin, 1);
                else
                    a2lag = [];
                    a2min = 21;
                    a2max = 21;
                end

                % reshape remainder so we can work
                cci = find(c(:, cc) ~= 0);
                bmaps = reshape(bc.GLMData.BetaMaps(:, :, :, cci), numel(a1lag), numel(cci));
                semap = semap(:);

                % initialize VMPData
                mapc.Map(cc).VMPData = zeros(size(a1lag));

                % iterate over lagged values
                for c1 = a1min:a1max

                    % compute required lag value
                    r1 = (c1 - 21) * 0.05;

                    % create preliminary mask
                    r1m = (a1lag == r1);

                    % iterate over possible finds for lag2
                    for c2 = a2min:a2max

                        % if required
                        if ~isempty(a2lag)

                            % restrict r1
                            r1 = max(-0.2, min(0.8, r1));

                            % compute lag value
                            r2 = (c2 - 21) * 0.05;

                            % update mask and find indices
                            r2m = find(r1m(:) & a2lag(:) == r2);

                            % restrict r2
                            r2 = max(-0.45, min(0.5, r2));

                            % nothing to do for this combination
                            if isempty(r2m)
                                continue;
                            end

                            % generate AR2-matrix
                            if mapopts.arlagv2w
                                v = sparse( ...
                                    [spi; spi(1:end-1); spi(2:end); spi(1:end-2); spi(3:end)], ...
                                    [spi; spi(2:end); spi(1:end-1); spi(3:end); spi(1:end-2)], ...
                                    [spi1; r1 .* spi1; r1 .* spi1(1:end-2); r2 .* spi1; r2 .* spi1(1:end-4)], ...
                                    spin, spin, 5 * spin - 6);

                                % generate weighting matrix and apply to design
                                tD = ne_sqrtm(v) * bc.DesignMatrix;

                            % double r-removal
                            else

                                % subtract from design matrix
                                r12 = 1 / (1 / (r1 * r1) - 1);
                                tD = bc.DesignMatrix(3:end, :) - ...
                                    (r1 * (1 - r2) * (1 + r12)) .* bc.DesignMatrix(2:end-1, :) - ...
                                    (r2 * (1 + r12) - r12) .* bc.DesignMatrix(1:end-2, :);

                                % remove timepoints
                                if ~isempty(r1r)
                                    tD(r2r, :) = [];
                                end
                            end

                        % not required
                        else

                            % take preliminary mask
                            r2m = find(r1m(:));

                            % generate AR1-matrix
                            if mapopts.arlagv2w
                                v = sparse( ...
                                    [spi; spi(1:end-1); spi(2:end)], ...
                                    [spi; spi(2:end); spi(1:end-1)], ...
                                    [spi1; r1 .* spi1; r1 .* spi1(1:end-2)], ...
                                    spin, spin, 3 * spin - 2);

                                % generate weighting matrix and apply to design
                                tD = ne_sqrtm(v) * bc.DesignMatrix;

                            % simple r-removal
                            else

                                % apply lag1 to design matrix
                                tD = bc.DesignMatrix(2:end, :) - r1 .* bc.DesignMatrix(1:end-1, :);

                                % remove timepoints (at study boundaries)
                                if ~isempty(r1r)
                                    tD(r1r, :) = [];
                                end
                            end
                        end

                        % then invert
                        itD = inv(tD' * tD);

                        % compute statistic (for matched indices)
                        mapc.Map(cc).VMPData(r2m) = glmtstat(c(cci, cc)', ...
                            double(bmaps(r2m, :)), itD(cci, cci), semap(r2m));
                    end
                end

            % calculate t stats (simple computation for no lags)
            else
                if istransio(bc.GLMData.BetaMaps)
                    cci = find(c(:, cc) ~= 0);
                    mapc.Map(cc).VMPData = reshape(glmtstat(c(cci, cc)', ...
                        double(bc.GLMData.BetaMaps(:, :, :, cci)), ...
                        iXX(cci, cci), semap), szmap);
                else
                    mapc.Map(cc).VMPData = reshape(glmtstat(c(:, cc)', ...
                        bc.GLMData.BetaMaps, iXX, semap), szmap);
                end
            end

            % set name
            mapc.Map(cc).Name = mapopts.names{cc};
        end

    % MTCs
    case 2

        % create SMP
        map = xff('new:smp');
        mapc = map.C;

        % make initial setting first
        mapc.NrOfVertices = bc.NrOfVertices;

        % calculation of SE map
        semap = sqrt((1 - (bc.GLMData.MultipleRegressionR .^ 2)) .* bc.GLMData.MCorrSS / nval);
        semap(semap == 0) = Inf;

        % generate tmaps
        for cc = 1:size(c, 2)

            % map options
            mapc.Map(cc).DF1 = nval;
            mapc.Map(cc).DF2 = 0;
            mapc.Map(cc).BonferroniValue = bc.NrOfVoxelsForBonfCorrection;
            mapc.Map(cc).UseRGBColor = 0;

            % calculate t stats
            if istransio(bc.GLMData.BetaMaps)
                cci = find(c(:, cc) ~= 0);
                mapc.Map(cc).SMPData = reshape(glmtstat(c(cci, cc)', ...
                    bc.GLMData.BetaMaps(:, cci), iXX(cci, cci), semap), szmap);
            else
                mapc.Map(cc).SMPData = reshape(glmtstat(c(:, cc)', ...
                    bc.GLMData.BetaMaps, iXX, semap), szmap);
            end
            % set name
            mapc.Map(cc).Name = mapopts.names{cc};
        end

     % unknown types
    otherwise
        error('neuroelf:xff:unknownOption', 'Unknown GLM ProjectType: %d', bc.ProjectType);
end

% set back to memory
map.C = mapc;

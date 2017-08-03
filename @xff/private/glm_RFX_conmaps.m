function conmaps = glm_RFX_conmaps(xo, c, mapopts)
% GLM::RFX_conmaps  - return contrast maps (sums of betas)
%
% FORMAT:       maps = glm.RFX_conmaps(c, [, mapopts])
%
% Input fields:
%
%       c           NxC contrast vector
%       mapopts     structure with optional fields
%        .autobal   auto-balance contrast if pos and abs(neg) == 1 (true)
%        .meanr     boolean flag, remove mean from map (added as cov)
%        .meanrmsk  mask to get mean from (object or XxYxZ logical)
%        .subsel    subject selection (otherwise all subjects)
%
% Output fields:
%
%       maps        either 4D or 2D data array (VTC/MTC)
%
% Using: findfirst, lsqueeze.

% Version:  v1.1
% Build:    16040621
% Date:     Apr-06 2016, 9:47 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2012, 2014, 2016, Jochen Weber
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
findfirst = ne_methods.findfirst;

% main argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
   ~isa(c, 'double') || isempty(c) || any(isnan(c(:)) | isinf(c(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get content handle and check predictor separation
bc = xo.C;
if bc.SeparatePredictors ~= 2
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% is not supported for FMR-based data
if ~any(bc.ProjectType == [1, 2])
    error('neuroelf:xff:unsupported', 'RFX map extraction of FMRs are not yet supported.');
end

% get subject predictors
ffxspred = glm_SubjectPredictors(xo);

% RFX
isrfx = (bc.ProjectTypeRFX > 0);
if isrfx
    numsubs = numel(bc.GLMData.Subject);
    numspred = size(bc.GLMData.Subject(1).BetaMaps, ndims(bc.GLMData.Subject(1).BetaMaps));

% FFX
else
    ffxpred = bc.Predictor;
    ffxpred = {ffxpred(:).Name2};
    ffxpred = ffxpred(:);
    ffxsubs = glm_Subjects(xo);
    numsubs = numel(ffxsubs);
    numspred = numel(ffxspred) + 1;
end

% does contrast vector need transpose
if any(numel(c) == [numspred, numspred - 1])
    c = c(:)';
end

% check row length
if ~any(size(c, 2) == [numspred, numspred - 1])
    error('neuroelf:xff:badArgument', 'Invalid first-level contrast spec.');
end

% for VTC
if bc.ProjectType == 1
    if isrfx
        msz = size(bc.GLMData.RFXGlobalMap);
    else
        msz = size(bc.GLMData.MCorrSS);
    end

% for MTC
else

    % use only ONE dimension (numel)
    if isrfx
        msz = numel(bc.GLMData.RFXGlobalMap);
    else
        msz = numel(bc.GLMData.MCorrSS);
    end
end

% options
if nargin < 3 || ~isstruct(mapopts) || numel(mapopts) ~= 1
    mapopts = struct;
end

% auto-balancing
if ~isfield(mapopts, 'autobal') || ~islogical(mapopts.autobal) || numel(mapopts.autobal) ~= 1
    mapopts.autobal = true;
end

% remove mean
if ~isfield(mapopts, 'meanr') || ~islogical(mapopts.meanr) || numel(mapopts.meanr) ~= 1
    mapopts.meanr = false;
end

% mean-removal mask
if isfield(mapopts, 'meanrmsk') && numel(mapopts.meanrmsk) == 1 && ...
    xffisobject(mapopts.meanrmsk, true, 'msk')
    mbc = mapopts.meanrmsk.C;
    if numel(mbc.Mask) == prod(msz)
        mapopts.meanrmsk = ne_methods.lsqueeze(mbc.Mask > 0);
    else
        error('neuroelf:xff:sizeMismatch', 'Invalid mask supplied.');
    end
elseif isfield(mapopts, 'meanrmsk') && islogical(mapopts.meanrmsk) && ...
    numel(mapopts.meanrmsk) == prod(msz)
    mapopts.meanrmsk = ne_methods.lsqueeze(mapopts.meanrmsk);
else
    mapopts.meanrmsk = [];
end

% further restrict to available betamaps (for RFX)
if isempty(mapopts.meanrmsk) && mapopts.meanr && isrfx
    mapopts.meanrmsk = all(bc.GLMData.Subject(1).BetaMaps ~= 0, ...
        ndims(bc.GLMData.Subject(1).BetaMaps));
    for sc = 1:numsubs
        mapopts.meanrmsk = (mapopts.meanrmsk & ...
            all(bc.GLMData.Subject(1).BetaMaps ~= 0, ...
            ndims(bc.GLMData.Subject(1).BetaMaps)));
    end
    mapopts.meanrmsk = ne_methods.lsqueeze(mapopts.meanrmsk);
else
    mapopts.meanrmsk = false;
    mapopts.meanr = false;
end

% final check
meanrmsk = mapopts.meanrmsk;
if ~any(meanrmsk)
    mapopts.meanr = false;
end

% compute 1 / sum(mask(:))
meanrsum = 1 ./ sum(meanrmsk(:));

% subject selection (for this function, only double supported)
if ~isfield(mapopts, 'subsel') || ~isa(mapopts.subsel, 'double') || isempty(mapopts.subsel) || ...
    any(isinf(mapopts.subsel(:)) | isnan(mapopts.subsel(:))) || ...
    numel(unique(round(mapopts.subsel(:)))) ~= numel(mapopts.subsel) || ...
    any(mapopts.subsel(:) < 1 | mapopts.subsel(:) > numsubs)
    mapopts.subsel = 1:numsubs;
else
    mapopts.subsel = round(mapopts.subsel(:)');
end
subsel = mapopts.subsel;
numsubs = numel(subsel);

% transpose contrast, and padd with (mean-confound) zero if necessary
c = c';
if size(c, 1) == (numspred - 1)
    c(end + 1, :) = 0;
end

% number of contrasts
nummaps = size(c, 2);

% access and storage argument
if bc.ProjectType == 1
    subsa = {':', ':', ':', [], []};
    subsr = {':', ':', ':', []};
    mdim = 3;
else
    subsa = {':', [], []};
    subsr = {':', []};
    mdim = 1;
end

% extraction
conmaps = zeros([msz, numsubs, nummaps]);
for cc = 1:nummaps

    % access correct map (storage)
    subsa{end} = cc;

    % fill contrast maps
    csum = zeros([ones(1, mdim), numsubs]);
    csumc = sum(c(:, cc));
    if abs(csumc) < sqrt(eps) && abs(sum(c(c(:, cc) > 0, cc)) - 1) < sqrt(eps) && ...
        abs(1 + sum(c(c(:, cc) < 0, cc))) < sqrt(eps) && mapopts.autobal
        autobal = true;
        possum = zeros(numsubs, 1);
        negsum = zeros(numsubs, 1);
        pconmaps = zeros([msz, numsubs]);
        nconmaps = zeros([msz, numsubs]);
    else
        autobal = false;
    end
    for pc = 1:numspred
        if c(pc, cc) ~= 0

            % RFX
            if isrfx
                subsr{end} = pc;
                for sc = 1:numsubs
                    subsa{end-1} = sc;
                    conext = bc.GLMData.Subject(subsel(sc)).BetaMaps(subsr{:});
                    if all(isinf(conext(:)) | isnan(conext(:)) | conext(:) == 0)
                        continue;
                    end
                    if mapopts.meanr
                        conext(meanrmsk) = conext(meanrmsk) - ...
                            meanrsum .* sum(conext(meanrmsk));
                    end
                    if autobal
                        if c(pc, cc) > 0
                            possum(sc) = possum(sc) + c(pc, cc);
                            pconmaps(subsa{1:end-1}) = pconmaps(subsa{1:end-1}) + ...
                                c(pc, cc) .* conext;
                        else
                            negsum(sc) = negsum(sc) + c(pc, cc);
                            nconmaps(subsa{1:end-1}) = nconmaps(subsa{1:end-1}) + ...
                                c(pc, cc) .* conext;
                        end
                    else
                        conmaps(subsa{:}) = conmaps(subsa{:}) + c(pc, cc) .* conext;
                        csum(sc) = csum(sc) + c(pc, cc);
                    end
                end

            % FFX
            else

                % find index
                for sc = 1:numsubs
                    subsa{end-1} = sc;
                    keepsubi = findfirst(~cellfun('isempty', regexpi(ffxpred, ...
                        sprintf('^subject\\s+%s:\\s*%s', ...
                        ffxsubs{subsel(sc)}, ffxspred{pc}))));
                    if ~isempty(keepsubi)
                        subsr{end} = keepsubi;
                        conext = bc.GLMData.BetaMaps(subsr{:});
                        if all(isinf(conext(:)) | isnan(conext(:)) | conext(:) == 0)
                            continue;
                        end
                        if mapopts.meanr
                            conext(meanrmsk) = conext(meanrmsk) - ...
                                meanrsum .* sum(conext(meanrmsk));
                        end
                        if autobal
                            if c(pc, cc) > 0
                                possum(sc) = possum(sc) + c(pc, cc);
                                pconmaps(subsa{1:end-1}) = pconmaps(subsa{1:end-1}) + ...
                                    c(pc, cc) .* conext;
                            else
                                negsum(sc) = negsum(sc) + c(pc, cc);
                                nconmaps(subsa{1:end-1}) = nconmaps(subsa{1:end-1}) + ...
                                    c(pc, cc) .* conext;
                            end
                        else
                            conmaps(subsa{:}) = conmaps(subsa{:}) + c(pc, cc) .* conext;
                            csum(sc) = csum(sc) + c(pc, cc);
                        end
                    end
                end
            end
        end
    end
    if autobal
        if ~all(abs(possum - 1) < sqrt(eps))
            pconmaps = reshape(reshape(pconmaps, prod(msz), numsubs) * diag(1 ./ possum), [msz, numsubs]);
        end
        if ~all(abs(1 + negsum) < sqrt(eps))
            nconmaps = reshape(reshape(nconmaps, prod(msz), numsubs) * diag(-1 ./ negsum), [msz, numsubs]);
        end
        subsa{end-1} = ':';
        conmaps(subsa{:}) = pconmaps + nconmaps;
        badsum = (possum == 0 | negsum == 0);
        if any(badsum)
            subsa{end-1} = badsum;
            conmaps(subsa{:}) = NaN;
        end
    else
        csumb = find(csum(:) ~= csumc);
        if ~isempty(csumb)
            subsa{end-1} = csumb;
            subsr{end} = csumb;
            if all(c(:, cc) >= 0) || all(c(:, cc) <= 0)
                if ~all(csum == csumc)
                    conmaps(subsa{:}) = conmaps(subsa{:}) .* ...
                        repmat(csumc ./ csum(subsr{:}), [msz, 1]);
                end
            else
                conmaps(subsa{:}) = NaN;
            end
        end
    end
end

% replace invalid maps with 0
if ~isrfx
    conmaps(isnan(conmaps)) = 0;
end

function [out, vbv] = glm_SingleTrialSVM(xo, opts)
% GLM::SingleTrialSVM  - perform single-trial GLM based SVM classification
%
% FORMAT:       output = glm.SingleTrialSVM(opts);
%
% Input fields:
%
%       opts        mandatory settings:
%        .condsel   1xN cell array with condition selection (N >= 2)
%        .holdout   1x1 number of trials for holdout (per condition)
%        .nperm     number of permutations (without label scrambling)
%        .type      either of 'searchlight' or 'voi'
%
%                   general optional settings:
%        .sperm     number of permutations (with label scrambling, 0)
%        .ztrans    z-transform within class across trials (default: false)
%
%                   for type := 'voi', mandatory options
%        .voi       VOI object
%
%                   and optional
%        .voiidx    index into ROI object (default: all)
%
% Output fields:
%
%       output      output, dependent on type
%                   for 'searchlight', VMP with per-subject maps + RFX map
%                   for 'voi', SxVx2/3 subject-by-VOIidx array with mean
%                   accuracy and error measure (SD or LB/UB from scrambling)
%
% Using: applyspmsnc, multimatch, ne_svm*.

% Version:  v1.1
% Build:    16020412
% Date:     Feb-04 2016, 12:56 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2015, 2016, Jochen Weber
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
% applyspmsnc   = ne_methods.applyspmsnc;
multimatch    = ne_methods.multimatch;
ne_svmtrain   = ne_methods.ne_svmtrain;
ne_svmpredict = ne_methods.ne_svmpredict;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
    numel(opts) ~= 1 || ~isstruct(opts) || ~isfield(opts, 'condsel') || ...
   ~iscell(opts.condsel) || numel(opts.condsel) < 2 || ...
    any(cellfun('isempty', opts.condsel(:))) || ...
   ~isfield(opts, 'holdout') || ~isa(opts.holdout, 'double') || numel(opts.holdout) ~= 1 || ...
    isinf(opts.holdout) || isnan(opts.holdout) || opts.holdout < 1 || ...
   ~isfield(opts, 'nperm') || ~isa(opts.nperm, 'double') || numel(opts.nperm) || ...
    isinf(opts.nperm) || isnan(opts.nperm) || opts.nperm < 1 || ...
   ~isfield(opts, 'type') || ~ischar(opts.type) || isempty(opts.type) || ...
   ~any(lower(opts.type(1)) == 'sv')
    error('neuroelf:xff:badArgument', 'Invalid call to GLM::SingleTrialSVM.');
end

% searchlight needs implementation
if lower(opts.type(1)) == 's'
    error('neuroelf:xff:notYetImplemented', 'Searchlight not yet implemented.');
end

% optional settings
if ~isfield(opts, 'sperm') || ~isa(opts.sperm, 'double') || numel(opts.sperm) ~= 1 || ...
    isinf(opts.sperm) || isnan(opts.sperm) || opts.sperm < 40
    opts.sperm = 0;
end
if ~isfield(opts, 'ztrans') || numel(opts.ztrans) ~= 1 || ...
   (~isa(opts.ztrans, 'double') && ~islogical(opts.ztrans))
    opts.ztrans = false;
else
    opts.ztrans = logical(opts.ztrans);
end

% get content
bc = xo.C;
if bc.ProjectType ~= 1 || bc.ProjectTypeRFX ~= 1
    error('neuroelf:xff:badArgument', 'Only valid for VTC based RFX-GLM files.');
end

% number of classify cross-validation folds
nfolds = opts.nfold;

% number of hold-out samples (per class)
nholdout = opts.holdout;

% number of subjects
subs = glm_Subjects(xo);
nsubs = numel(subs);

% get subject predictors
sps = glm_SubjectPredictors(xo);

% get indices of condition entries
condsel = opts.condsel(:)';
for cc = 1:numel(condsel)
    if ischar(condsel{cc})
        condsel{cc} = find(~cellfun('isempty', regexpi(sps, condsel{cc}(:)')));
    elseif iscell(condsel{cc})
        condsel{cc} = find(multimatch(sps, condsel{cc}(:), true) > 0);
    elseif ~isa(condsel{cc}, 'double') || isempty(condsel{cc}) || ...
        any(isinf(condsel{cc}(:)) | isnan(condsel{cc}(:)) | condsel{cc}(:) < 1 | condsel{cc}(:) > numel(sps))
        error('neuroelf:xff:badOption', 'Invalid condition selection.');
    else
        condsel{cc} = unique(round(condsel{cc}(:)));
    end
    condsel{cc} = condsel{cc}(:);
    for cci = 1:(cc-1)
        if ~isempty(intersect(condsel{cc}, condsel{cci}))
            error('neuroelf:xff:badOption', 'Invalid condition selection.');
        end
    end
end
condidx = condsel;
condnum = 0;
for cc = 1:numel(condidx)
    condidx{cc} = (condnum+1):(condnum+numel(condidx{cc}));
    condnum = condnum + numel(condidx{cc});
end

% for searchlight
if lower(opts.type(1)) == 's'
    
    % preset out and vbv
    out = [];
    vbv = [];
    
% for VOI
elseif lower(opts.type(1)) == 'v'
    
    % required settings
    if ~isfield(opts, 'voi') || numel(opts.voi) ~= 1 || ~xffisobject(opts.voi, true, 'voi')
        error('neuroelf:xff:badOption', 'Bad or missing VOI object in call.');
    end
    voic = opts.voi.C;

    % optional settings
    if ~isfield(opts, 'voiidx') || ~isa(opts.voiidx, 'double') || isempty(opts.voiidx) || ...
        any(isinf(opts.voiidx(:)) | isnan(opts.voiidx(:)) | opts.voiidx(:) < 1)
        opts.voiidx = 1:numel(voic.VOI);
    else
        opts.voiidx = unique(min(numel(voic.VOI), round(opts.voiidx(:))))';
    end
    nvoi = numel(opts.voiidx);
    
    % access (individual voxel) raw data
    [vb, vbv] = glm_VOIBetas(xo, opts.voi, struct('vl', opts.voiidx));

    % generate matrix holding the N-fold accuracy values
    out = NaN .* zeros(nsubs, nfolds, nvoi);
    
    % iterate over VOIs
    for vc = 1:nvoi

        % repeat over subjects
        for sc = 1:nsubs
            
            % copy indices
            condidxc = condidx;
            
            % collect voxel values
            voxvals = cat(2, vbv{sc, cat(1, condsel{:}), vc});
            
            % remove bad conditions
            badconds = all(isinf(voxvals) | isnan(voxvals) | voxvals == 0, 1);
            for cc = 1:numel(condidx)
                condidxc{cc}(badconds(condidx{cc})) = [];
            end
            condidxa = cat(2, condidxc{:});
            condidxn = min(cellfun('prodofsize', condidxc));
            
            % remove bad voxels
            badvox = any(isinf(voxvals(:, condidxa)) | isnan(voxvals(:, condidxa)), 2);
            voxvals(badvox, :) = [];

            % number of learning set and prediction indices
            learnto = condidxn - nholdout;
            predictfrom = learnto + 1;

            % training labels
            tl = [ones(f1, 1); 2 .* ones(f2, 1)];
            tsl = [ones(nholdout, 1); 2 .* ones(nholdout, 1)];

            % iterate over subjects
            for sc = 1:nsubs

                % show progress
                fprintf('Working on subject %s...\n', subs{sc});
                pause(0.01);

                % create random vectors
                [rvd, rvi1] = sort(randn(e1, nfolds));
                [rvd, rvi2] = sort(randn(e2, nfolds));

                % iterate over folds
                for fc = 1:nfolds

                    % concatenate data over selected folds
                    d1 = double(cat(2, vbv{sc, sensneg(rvi1(1:f1, fc))})');
                    d2 = double(cat(2, vbv{sc, resneg(rvi2(1:f1, fc))})');

                    % classifier training
                    trained = ne_svmtrain(tl, [d1; d2]);

                    % concatenate test data
                    testdata = cat(1, ...
                        double(cat(2, vbv{sc, sensneg(rvi1(t1:end, fc))})'), ...
                        double(cat(2, vbv{sc, resneg(rvi2(t2:end, fc))})'));

                    % check accuracy
                    [plab, acc] = ne_svmpredict(tsl, testdata, trained);

                    % store accuracy
                    class_accuracy(sc, fc) = 0.01 * acc(1);
                end
            end
        end
    end
end

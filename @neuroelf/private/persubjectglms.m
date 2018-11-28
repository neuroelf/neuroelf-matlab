function varargout = persubjectglms(mdm, opts)
% persubjectglms  - compute one GLM per subject
%
% FORMAT:       [glm1, glm2] = persubjectglms(mdm [, opts])
%
% Input fields:
%
%       mdm         MDM object (or filename) with at least 3 subjects
%       opts        optional settings
%        .cmbffx    combine outputs to one FFX GLM (default: false)
%        .cmbrfx    combine outputs to one RFX GLM (default: true)
%        .loadglms  load existing per-subject GLMs (default: false)
%        .outpatt   outfile pattern, default: '%M_%S_FFX.glm'
%        .prtpnorm  normalize parameters in PRTs (default: false)
%        .sglobsigs per-subject global signals
%        .spmdir    if set (to a valid folder), use SPM for regression
%        .subsel    subject selection (cell array with IDs)
%
% Output fields:
%
%       glm1, glm2  combined GLMs (order is FFX first if both are used)
%
% Note: additionally all fields to MDM::ComputeGLM are supported in opts.

% Version:  v1.1
% Build:    16061611
% Date:     Jun-16 2016, 11:21 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2015, 2016, Jochen Weber
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
if nargin < 1 || isempty(mdm) || (~ischar(mdm) && (numel(mdm) ~= 1 || ~isxff(mdm, 'mdm')))
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
if ischar(mdm)
    mdmf = mdm(:)';
    mdm = [];
else
    mdmf = '';
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if isfield(opts, 'cmbrfx') && ~isfield(opts, 'cmbffx') && islogical(opts.cmbrfx) && ...
    numel(opts.cmbrfx) == 1 && ~opts.cmbrfx
    opts.cmbffx = true;
end
if ~isfield(opts, 'cmbffx') || ~islogical(opts.cmbffx) || numel(opts.cmbffx) ~= 1
    opts.cmbffx = false;
end
if ~isfield(opts, 'cmbrfx') || ~islogical(opts.cmbrfx) || numel(opts.cmbrfx) ~= 1
    opts.cmbrfx = true;
end
if ~isfield(opts, 'loadglms') || ~islogical(opts.loadglms) || numel(opts.loadglms) ~= 1
    opts.loadglms = false;
end
if ~isfield(opts, 'outpatt') || ~ischar(opts.outpatt)
    opts.outpatt = '%M_%S_FFX.glm';
elseif numel(opts.outpatt) < 6 || isempty(strfind(opts.outpatt(:)', '%S')) || ...
   ~strcmpi(lsqueeze(opts.outpatt(end-3:end))', '.glm') || ...
    sum(opts.outpatt(:)' == '%') > 2
    error('neuroelf:general:badArgument', 'Invalid output filename pattern.');
end
opts.outpatt = opts.outpatt(:)';
if ~isfield(opts, 'prtpnorm') || ~islogical(opts.prtpnorm) || numel(opts.prtpnorm) ~= 1
    opts.prtpnorm = false;
end
if ~isfield(opts, 'subsel') || ~iscell(opts.subsel) || isempty(opts.subsel) || ...
   ~ischar(opts.subsel{1}) || isempty(opts.subsel{1})
    opts.subsel = [];
end
if ~isfield(opts, 'sglobsigs') || ~iscell(opts.sglobsigs)
    opts.sglobsigs = {};
end
if ~isfield(opts, 'spmdir') || isempty(opts.spmdir) || ~ischar(opts.spmdir) || ...
    exist(opts.spmdir(:)', 'dir') ~= 7 || ...
    exist([opts.spmdir(:)', filesep, 'spm.m'], 'file') ~= 2
    usespm = false;
else
    spmdir = opts.spmdir(:)';
    usespm = ~isempty(which('spm'));
    if usespm
        try
            if ~any(strcmpi(spm('ver'), {'spm8', 'spm12b'}))
                usespm = false;
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            usespm = false;
        end
    end
end

% load MDM
if ~isempty(mdmf)
    try
        mdm = xff(mdmf);
        if ~isxff(mdm, 'mdm')
            error('neuroelf:general:badArgument', 'Not an MDM file: %s.', mdmf);
        end
    catch ne_eo;
        if isxff(mdm)
            mdm.ClearObject;
        end
        rethrow(ne_eo);
    end
end
mdmc = getcont(mdm);

% patch pattern again
if ~isempty(strfind(opts.outpatt, '%M'))
    [mdmfp, mdmfn] = fileparts(mdm.FilenameOnDisk);
    if isempty(mdmfn)
        mdmfn = 'UNSAVED';
    elseif ~isempty(mdmfp)
        mdmfn = strrep(strrep([mdmfp '/' mdmfn], '\', '/'), '//', '/');
    end
    opts.outpatt = strrep(opts.outpatt, '%M', mdmfn);
    if sum(opts.outpatt == '%') ~= 1
        if ~isempty(mdmf)
            mdm.ClearObject;
        else
            setcont(mdm, mdmc);
        end
        error('neuroelf:general:badArgument', 'Invalid output filename pattern.');
    end
end

% subject selection
ssids = mdm.Subjects;
fsids = mdm.Subjects(true);
if ~isempty(opts.subsel)
    ksub = (multimatch(ssids, opts.subsel(:)) > 0);
    ksubi = zeros(numel(ssids), 1);
    ksubi(ksub) = 1:sum(ksub);
    ksubs = (multimatch(fsids, opts.subsel(:)) > 0);
    if isfield(opts, 'motpars') && iscell(opts.motpars) && ...
        numel(opts.motpars) == size(mdm.XTC_RTC, 1)
        mdm.RunTimeVars.MotionParameters = lsqueeze(opts.motpars(ksubs));
        opts.motpars = true;
    elseif isfield(mdm.RunTimeVars, 'MotionParameters') && ...
        iscell(mdm.RunTimeVars.MotionParameters) && ...
        numel(mdm.RunTimeVars.MotionParameters) == size(mdm.XTC_RTC, 1)
        mdm.RunTimeVars.MotionParameters = lsqueeze(mdm.RunTimeVars.MotionParameters(ksubs));
    end
    if isfield(mdm.RunTimeVars, 'CovariatesData') && ...
        size(mdm.RunTimeVars.CovariatesData, 1) == numel(ssids)
        mdm.RunTimeVars.CovariatesData = mdm.RunTimeVars.CovariatesData(ksub, :);
    end
    if isfield(mdm.RunTimeVars, 'Groups') && iscell(mdm.RunTimeVars.Groups) && ...
       ~isempty(mdm.RunTimeVars.Groups) && size(mdm.RunTimeVars.Groups, 2) == 2
        for gc = size(mdm.RunTimeVars.Groups, 1):-1:1
            if ~isempty(mdm.RunTimeVars.Groups{gc, 2})
                ksubg = false(numel(ssids), 1);
                ksubg(mdm.RunTimeVars.Groups{gc, 2}) = true;
                mdm.RunTimeVars.Groups{gc, 2} = ksubi(ksub & ksubg);
                if isempty(mdm.RunTimeVars.Groups{gc, 2})
                    mdm.RunTimeVars.Groups(gc, :) = [];
                end
            end
        end
    end
    mdm.XTC_RTC = mdm.XTC_RTC(ksubs, :);
    fsids = mdm.Subjects(true);
else
    if isfield(opts, 'motpars') && iscell(opts.motpars) && ...
        numel(opts.motpars) == size(mdm.XTC_RTC, 1)
        mdm.RunTimeVars.MotionParameters = opts.motpars(:);
        opts.motpars = true;
    end
    if isfield(opts, 'xconfound') && iscell(opts.xconfound) && ...
        numel(opts.xconfound) == size(mdm.xconfound, 1)
        mdm.RunTimeVars.XConfounds = opts.xconfound(:);
        opts.xconfound = true;
    end
end

% find subject IDs
sids = mdm.Subjects;
if numel(sids) < 2
    if ~isempty(mdmf)
        mdm.ClearObject;
    else
        setcont(mdm, mdmc);
    end
    error('neuroelf:general:badArgument', 'MDM requires at least 2 subjects.');
end

% per-subject global signals must match
if numel(sids) ~= size(opts.sglobsigs, 1)
    if isequal(opts.sglobsigs, {'SubjectGlobSigs'}) && ...
        isfield(mdmc.RunTimeVars, 'SubjectGlobSigs') && ...
        isstruct(mdmc.RunTimeVars.SubjectGlobSigs) && ...
        numel(mdmc.RunTimeVars.SubjectGlobSigs) == 1
        opts.sglobsigs = cell(numel(sids), 1);
        for sc = 1:numel(sids)
            if isfield(mdmc.RunTimeVars.SubjectGlobSigs, sids{sc})
                opts.sglobsigs(sc, 1:numel(mdmc.RunTimeVars.SubjectGlobSigs.(sids{sc}))) = ...
                    mdmc.RunTimeVars.SubjectGlobSigs.(sids{sc});
            end
        end
    else
        opts.sglobsigs = {};
    end
end

% settings
mdm.RFX_GLM = 0;
mdm.SeparatePredictors = 2;
mdm.NrOfStudies = size(mdm.XTC_RTC, 1);
nstudy = mdm.NrOfStudies;

% iterate over subjects
glms = cell(numel(sids), 1);
for sc = 1:numel(sids)

    % load only
    if opts.loadglms
        glms{sc} = strrep(opts.outpatt, '%S', sids{sc});
        if exist(glms{sc}, 'file') ~= 2
            if ~isempty(mdmf)
                mdm.ClearObject;
            else
                setcont(mdm, mdmc);
            end
            error( ...
                'neuroelf:BadArgument', ...
                'Cannot load GLM for subject %s.', ...
                sids{sc} ...
            );
        end
        continue;
    end

    % make a copy
    smdm = mdm.CopyObject;

    % match ID to list
    kidx = find(strcmpi(ssids, sids{sc}));
    sidx = find(strcmpi(fsids, sids{sc}));

    % keep those entries
    smdm.XTC_RTC = mdm.XTC_RTC(sidx, :);
    if isfield(mdm.RunTimeVars, 'CovariatesData') && ...
        size(mdm.RunTimeVars.CovariatesData, 1) == numel(sids)
        smdm.RunTimeVars.CovariatesData = mdm.RunTimeVars.CovariatesData(kidx, :);
    end
    if isfield(mdm.RunTimeVars, 'MotionParameters') && ...
        numel(mdm.RunTimeVars.MotionParameters) == nstudy && ...
        islogical(opts.motpars) && opts.motpars
        smdm.RunTimeVars.MotionParameters = mdm.RunTimeVars.MotionParameters(sidx);
    end
    if isfield(mdm.RunTimeVars, 'XConfounds') && ...
        numel(mdm.RunTimeVars.XConfounds) == nstudy && ...
        islogical(opts.xconfound) && opts.xconfound
        smdm.RunTimeVars.XConfounds = mdm.RunTimeVars.XConfounds(sidx);
    end
    
    % use SPM?
    if usespm
        try
            subdir = [spmdir filesep sids{sc}];
            if exist(subdir, 'dir') == 0
                mkdir(spmdir, sids{sc});
            else
                mdelete(findfiles(subdir, '*.*', '-d1'));
            end
            
            % create first level design across runs
            
        catch ne_eo;
            fprintf('Error using SPM: %s.\n', ne_eo.message);
            neuroelf_lasterr(ne_eo);
            usespm = false;
        end
    end

    % compute and save GLM
    try
        glm = [];
        if ~isempty(opts.sglobsigs)
            opts.globsigs = opts.sglobsigs(sc, :)';
            opts.globsigs(isempty(opts.globsigs)) = [];
            for gsc = 1:numel(opts.globsigs)
                if ischar(opts.globsigs{gsc})
                    opts.globsigs{gsc} = xff(opts.globsigs{gsc});
                end
            end
        end
        glm = smdm.ComputeGLM(opts);
        if ~isempty(opts.sglobsigs)
            clearxffobjects(opts.globsigs);
        end
        glm.SaveAs(strrep(opts.outpatt, '%S', sids{sc}));
        glm.SaveRunTimeVars;
        glms{sc} = glm.FilenameOnDisk;
        glm.ClearObject;
        smdm.ClearObject;
    catch ne_eo;
        if isxff(glm)
            glm.ClearObject;
        end
        smdm.ClearObject;
        if ~isempty(mdmf)
            mdm.ClearObject;
        else
            setcont(mdm, mdmc);
        end
        rethrow(ne_eo);
    end
end

% combine GLMs
if opts.cmbffx || opts.cmbrfx

    % load GLMs (with transio)
    for sc = 1:numel(glms)
        try
            glms{sc} = xff(glms{sc}, 't');
            if ~isxff(glms{sc}, 'glm') || glms{sc}.ProjectTypeRFX > 0
                error('neuroelf:error:badGLMFile', ...
                    'Bad GLM file for subject %s.', sids{sc});
            end
            if opts.cmbffx && isfield(glms{sc}, 'SubjectSPMsn') && ...
                isstruct(glms{sc}.SubjectSPMsn) && numel(glms{sc}.SubjectSPMsn) == 1 && ...
               ~isempty(fieldnames(glms{sc}.SubjectSPMsn))
                error('neuroelf:error:badCombination', ...
                    'Subject-GLMs with SPM-normalization cannot be FFX-combined.');
            end
        catch ne_eo;
            clearxffobjects(glms);
            if ~isempty(mdmf)
                mdm.ClearObject;
            else
                setcont(mdm, mdmc);
            end
            rethrow(ne_eo);
        end
    end
end

% RFX
nss = numel(glms);
if opts.cmbrfx

    % try/catch
    try

        % make a copy of the first GLM
        rfx = glms{1}.CopyObject;

        % collect some data
        nrtp = 0;
        nrct = 0;
        nrcf = 0;
        msk = zeros(size(glms{1}.GLMData.MCorrSS));
        sts = cell(1, nss);
        glmf = repmat(struct('SourceFile', '', 'SubjectID', '', ...
            'Predictors', [], 'iXX', [], 'DF1', -1, 'SEMap', [], 'ARLag', []), [1, nss]);
        prs = glms{1}.Predictor;
        pns = lsqueeze({prs(1:(glms{1}.NrOfPredictors - glms{1}.NrOfConfounds)).Name2});
        pns = regexprep(pns, '^Subject\s+(\S+)\:\s*', '');
        pnc = zeros(12, 0);
        for sc = 1:nss
            nrtp = nrtp + glms{sc}.NrOfTimePoints;
            nrct = nrct + glms{sc}.NrOfConfounds;
            nrcf = nrcf + numel(glms{sc}.NrOfConfoundsPerStudy);
            sts{sc} = glms{sc}.Study;
            prs = glms{sc}.Predictor;
            pnn = lsqueeze({prs(1:(glms{sc}.NrOfPredictors - glms{sc}.NrOfConfounds)).Name2});
            pnn = regexprep(pnn, '^Subject\s+(\S+)\:\s*', '');
            pns = uunion(pns, pnn);
            if size(pnc, 2) < numel(pns)
                for pc = (size(pnc, 2)+1):numel(pns)
                    pnc(:, pc) = prs(findfirst(strcmp(pnn, pns{pc}))).RGB(:);
                end
            end
            msk = msk + double(glms{sc}.GLMData.MCorrSS > 0);
            glmf(sc).SourceFile = glms{sc}.FilenameOnDisk;
            glmf(sc).iXX = glms{sc}.iXX;

            % calculation of SE map (from @xff/private/glm_FFX_tMap.m)
            nval = glms{sc}.NrOfTimePoints - glms{sc}.NrOfPredictors;
            semap = sqrt( ...
                (1 - (double(glms{sc}.GLMData.MultipleRegressionR) .^ 2)) .* ...
                double(glms{sc}.GLMData.MCorrSS) / nval);
            semap(semap == 0) = Inf;
            if ~isempty(glms{sc}.GLMData.ARLag)
                semap = sqrt(nval / ...
                    (nval - (glms{sc}.SerialCorrelation * glms{sc}.NrOfStudies))) .* semap;
                nval = nval - glms{sc}.SerialCorrelation * glms{sc}.NrOfStudies;
                if sum(glms{sc}.DesignMatrix(:) ~= 0) < (0.5 * numel(glms{sc}.DesignMatrix))
                    glmf(sc).DesignMatrix = sparse(glms{sc}.DesignMatrix);
                else
                    glmf(sc).DesignMatrix = single(glms{sc}.DesignMatrix);
                end
            end
            glmf(sc).Predictors = glms{sc}.Predictor;
            glmf(sc).DF1 = nval;
            glmf(sc).SEMap = semap;
            glmf(sc).ARLag = glms{sc}.GLMData.ARLag;
            if istransio(glmf(sc).ARLag)
                glmf(sc).ARLag = resolve(glmf(sc).ARLag);
            end
        end
        sts = catstruct(sts{:});
        msk = (msk >= (0.75 * nss));
        ntp = (numel(pns) + 1) * nss;
        rfxspmsn = struct;
        rfxtrfpl = struct;

        % restructure
        rfx.ProjectTypeRFX = 1;
        rfx.NrOfSubjects = nss;
        rfx.NrOfSubjectPredictors = numel(pns) + 1;
        rfx.NrOfTimePoints = nrtp;
        rfx.NrOfPredictors = ntp;
        rfx.NrOfConfounds = nrct;
        rfx.NrOfStudies = numel(sts);
        rfx.NrOfStudiesWithConfounds = nrcf;
        rfx.NrOfConfoundsPerStudy = zeros(1, nrcf);
        rfx.SeparatePredictors = 2;
        rfx.NrOfVoxelsForBonfCorrection = sum(msk(:));
        rfx.CortexBasedStatisticsMaskFile = '';
        rfx.Study = sts;
        rfxPredictor = rfx.Predictor;
        rfxPredictor(ntp).Name1 = '';
        rfxPredictor = rfxPredictor(:);
        rfx.DesignMatrix = [];
        rfx.iXX = [];
        rfxGLMData.MultipleRegressionR = [];
        rfxGLMData.MCorrSS = [];
        rfxGLMData.BetaMaps = [];
        rfxGLMData.XY = [];
        rfxGLMData.TimeCourseMean = [];
        rfxGLMData.RFXGlobalMap = single(msk);
        rfxGLMData.Subject = repmat(struct('BetaMaps', ...
            repmat(single(0), [size(msk), numel(pns) + 1])), 1, nss);
        rfxGLMData.RunTimeVars.Subject = repmat( ...
            struct('ARLag', [], 'iXX', [], 'SEMap', []), 1, nss);
        rfx.RunTimeVars.MotionParameters = cell(0, 1);
        rfx.RunTimeVars.PerSubjectGLMs = glmf;

        % generate predictor array
        preds = emptystruct({'Name1', 'Name2', 'RGB'}, [numel(pns), 1]);
        for pc = 1:numel(pns)
            preds(pc).RGB = reshape(pnc(:, pc), 4, 3);
        end

        % iterate over subjects
        ts = 1;
        tp = 1;
        for sc = 1:nss

            % get subject ID
            ssts = glms{sc}.Study;
            sid = ssts(1).NameOfAnalyzedFile;
            [sidf, sid] = fileparts(sid);
            sid = regexprep(sid, '^([^_]+)_.*$', '$1');
            tsid = makelabel(sid);
            rfx.RunTimeVars.PerSubjectGLMs(sc).SubjectID = tsid;

            % RunTimeVars
            rtv = glms{sc}.RunTimeVars;

            % subject-related fields
            if isfield(rtv, 'MotionParameters') && ...
                iscell(rtv.MotionParameters) && ...
                numel(rtv.MotionParameters) == glms{sc}.NrOfStudies
                rfx.RunTimeVars.MotionParameters = [ ...
                    rfx.RunTimeVars.MotionParameters(:); rtv.MotionParameters(:)];
            end
            if isfield(rtv, 'SubjectSPMsn') && ...
                isstruct(rtv.SubjectSPMsn) && ...
                numel(rtv.SubjectSPMsn) == 1 && ...
                isfield(rtv.SubjectSPMsn, tsid)
                rfxspmsn.(tsid) = rtv.SubjectSPMsn.(tsid);
            end
            if isfield(rtv, 'SubjectTrfPlus') && ...
                isstruct(rtv.SubjectTrfPlus) && ...
                numel(rtv.SubjectTrfPlus) == 1 && ...
                isfield(rtv.SubjectTrfPlus, tsid)
                rfxtrfpl.(tsid) = rtv.SubjectTrfPlus.(tsid);
            end

            % fill in header fields
            nst = numel(ssts);
            rfx.NrOfConfoundsPerStudy(ts:ts+nst-1) = glms{sc}.NrOfConfoundsPerStudy;
            prs = glms{sc}.Predictor;
            pnn = lsqueeze({prs([1:(glms{sc}.NrOfPredictors - glms{sc}.NrOfConfounds), end]).Name2});
            pnn = regexprep(pnn, '^Subject\s+(\S+)\:\s*', '');
            pni = multimatch(pns, pnn);
            rfxPredictor(tp:tp+numel(pns)-1) = preds;
            siXX = NaN(numel(pns) + 1);
            sipc = zeros(1, numel(pns) + 1);
            for pc = 1:numel(pns)
                rfxPredictor(tp).Name1 = sprintf('Predictor: %d', tp);
                rfxPredictor(tp).Name2 = sprintf('Subject %s: %s', sid, pns{pc});

                % also fill in beta map?
                if pni(pc) > 0
                    rfxGLMData.Subject(sc).BetaMaps(:, :, :, pc) = ...
                        glms{sc}.GLMData.BetaMaps(:, :, :, pni(pc));
                    sipc(pc) = pni(pc);
                end

                % increase counter
                tp = tp + 1;
            end
            siXX(sipc > 0, sipc > 0) = glmf(sc).iXX(sipc(sipc > 0), sipc(sipc > 0));
            
            % store some additional data
            rfxGLMData.RunTimeVars.Subject(sc) = struct( ...
                'ARLag', glmf(sc).ARLag, ...
                'iXX',   siXX, ...
                'SEMap', glmf(sc).SEMap);

            % add constant predictor info
            rfxPredictor(ntp+sc-nss) = glms{sc}.Predictor(end);
            rfxPredictor(ntp+sc-nss).Name1 = sprintf('Predictor: %d', ntp+sc-nss);
            rfxPredictor(ntp+sc-nss).Name2 = sprintf('Subject %s: Constant', sid);

            % combine constant maps
            cprs = ~cellfun('isempty', regexpi(lsqueeze({prs.Name2}), '\s+constant$'));
            cmaps = glms{sc}.GLMData.BetaMaps(:, :, :, cprs);
            for stc = 1:numel(ssts)
                cmaps(:, :, :, stc) = ssts(stc).NrOfTimePoints .* cmaps(:, :, :, stc);
            end
            rfxGLMData.Subject(sc).BetaMaps(:, :, :, end) = ...
                (1 / sum(cat(1, ssts.NrOfTimePoints))) .* sum(cmaps, 4);

            % increase counter
            ts = ts + nst;
        end

        % store output
        rfx.Predictor = rfxPredictor;
        rfx.GLMData = rfxGLMData;
        if numel(rfx.RunTimeVars.MotionParameters) ~= numel(rfx.Study)
            rfx.RunTimeVars.MotionParameters = cell(0, 1);
        end
        if ~isempty(fieldnames(rfxspmsn))
            rfx.RunTimeVars.SubjectSPMsn = rfxspmsn;
        end
        if ~isempty(fieldnames(rfxtrfpl))
            rfx.RunTimeVars.SubjectTrfPlus = rfxtrfpl;
        end

    % handle errors
    catch ne_eo;
        if ~isempty(mdmf)
            mdm.ClearObject;
        else
            setcont(mdm, mdmc);
        end
        clearxffobjects(glms);
        rfx.ClearObject;
        rethrow(ne_eo);
    end

    % where in output?
    if opts.cmbffx
        varargout{1} = [];
        varargout{2} = rfx;
    else
        varargout{1} = rfx;
    end
end

% FFX combination
if opts.cmbffx

    % try/catch
    try

        % always combine two at a time
        while numel(glms) > 1

            % iterate from the end
            for sc = numel(glms):-2:2

                % combine two GLMs
                nglm = glms{sc-1}.JoinFFX(glms{sc});

                % drop two old GLMs
                clearxffobjects(glms(sc-1:sc));

                % replace first one
                glms{sc-1} = nglm;

                % and remove second from list
                glms(sc) = [];
            end
        end

    % handle errors
    catch ne_eo;
        if ~isempty(mdmf)
            mdm.ClearObject;
        else
            setcont(mdm, mdmc);
        end
        clearxffobjects(glms);
        rethrow(ne_eo);
    end

    % set in output
    varargout{1} = glms{1}.CopyObject;
end

% clean up
if ~isempty(mdmf)
    mdm.ClearObject;
else
    setcont(mdm, mdmc);
end
clearxffobjects(glms);

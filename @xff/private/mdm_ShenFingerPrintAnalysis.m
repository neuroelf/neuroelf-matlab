function [stats, fpm, h] = mdm_ShenFingerPrintAnalysis(xo, atype, opts)
% MDM::ShenFingerFingerPrintAnalysis  - perform finger print analysis
%
% Format:       [stats, fpm, h] = mdm.ShenFingerPrintAnalysis(atype, opts)
%
% Input fields:
%
%       atype       analysis type, one of 'behav' or 'identity', where
%                   - 'behav' performs an correlational analysis with
%                      a behavioral covariate (e.g. task performance or RT)
%                   - 'identity' performs an analysis comparing fingerprints
%                      of different runs of participants to one another
%       opts        optional settings (required for 'behav')
%        .cov       behavioral covariate (if named, must be in RunTimeVars)
%        .fromrtv   use existing FingerPrint from VTC RunTimeVars (true)
%        .graphs    create graphs in new figure (true)
%        .overlap   only use edges from all leave-one-out rounds (true)
%        .prep4conn run VTC::PrepForConnectivity on VTCs (default: false)
%                   the opts argument is passed into this function!
%        .progress  either {true} or a 1x1 xfigure::progress or xprogress
%        .robust    perform robust regression (default: false)
%        .savefpm   target filename to save FingerPrint matrices ('')
%        .subsel    cell array with subject IDs to work on
%        .thresh    threshold for edge inclusion in predictive model (0.01)
%        .updrtv    update FingerPrint in VTC RunTimeVars (true)
%        .validcov  validation covariate (from separate dataset!)
%        .validfpm  validation finger-print data (from separate dataset!)
%                   unless both are given, use within-dataset validation
%
% Output fields:
%
%       stats       struct with statistics output fields
%        .fullmodel binary or weight matrix for edge selection (two-tailed)
%        .negmodel  binary or weight matrix for edge selection (neg-tail)
%        .posmodel  binary or weight matrix for edge selection (pos-tail)
%       fpm         finger print matrices (NxNxSubjects or NxNxRuns)
%       h           graphics handles (struct)
%        .fig       figure handle
%        .ax1, ...  axes handles
%
% Note: opts is also passed into AFT::ShenFingerPrint, so these optional
%       settings can be given as well.
%
% Using: calcbetas, correlpvalue, cov_nd, findfirst, fisherr2z,
%        fitrobustbisquare, multimatch, robcorrcoef.

% Version:  v1.1
% Build:    16061617
% Date:     Jun-16 2016, 5:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% NeuroElf library
global ne_methods;
calcbetas = ne_methods.calcbetas;
correlpvalue = ne_methods.correlpvalue;
cov_nd = ne_methods.cov_nd;
findfirst = ne_methods.findfirst;
fisherr2z = ne_methods.fisherr2z;
fitrobustbisquare = ne_methods.fitrobustbisquare;
multimatch = ne_methods.multimatch;
robcorrcoef = ne_methods.robcorrcoef;

% check input
if nargin < 2 || ~xffisobject(xo, true, 'mdm') || ~ischar(atype) || isempty(atype) || ...
   ~any(strcmpi(atype(:)', {'behav'; 'identity'}))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
atype = lower(atype(1));
bc = xo.C;
hc = xo.H;
h = struct;
try
    if ~isfield(hc, 'FilesChecked') || ~islogical(hc.FilesChecked) || ~hc.FilesChecked
        mdm_CheckFiles(xo, struct('autofind', true, 'silent', true));
        bc = xo.C;
    end
catch xfferror
    rethrow(xfferror);
end
mdmsubs = mdm_Subjects(xo);
rtv = bc.RunTimeVars;
vtcs = bc.XTC_RTC(:, 1);
if any(cellfun('isempty', regexpi(vtcs, '\.vtc$')))
    error('neuroelf:xff:badArgument', 'MDM must be VTC based.');
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    if atype == 'b'
        error('neuroelf:xff:badArgument', 'Behavioral analysis requires opts argument.');
    end
    opts = struct;
end
if atype == 'b' && (~isfield(opts, 'cov') || ...
   (~isa(opts.cov, 'double') && ~ischar(opts.cov)) || isempty(opts.cov))
    error('neuroelf:xff:badArgument', 'Behavioral analysis requires opts.cov option.');
elseif atype == 'b'
    bname = 'behavioral var';
    if isa(opts.cov, 'double')
        opts.cov = opts.cov(:);
        if numel(opts.cov) > numel(mdmsubs) || any(isinf(opts.cov))
            error('neuroelf:xff:badArgument', 'Invalid opts.cov option.');
        end
    else
        if ~isfield(rtv, 'CovariatesNames') || ~iscell(rtv.CovariatesNames) || ...
           ~any(strcmpi(rtv.CovariatesNames(:), opts.cov(:)')) || ...
           ~isfield(rtv, 'CovariatesData') || ~isa(rtv.CovariatesData, 'double') || ...
            isempty(rtv.CovariatesData) || size(rtv.CovariatesData, 2) ~= numel(rtv.CovariatesNames)
            error('neuroelf:xff:badArgument', 'Invalid or missing opts.cov name.');
        end
        bname = opts.cov(:)';
        opts.cov = rtv.CovariatesData(:, findfirst(strcmpi(rtv.CovariatesNames, bname)));
    end
end
if ~isfield(opts, 'fromrtv') || ~islogical(opts.fromrtv) || numel(opts.fromrtv) ~= 1
    opts.fromrtv = true;
end
if ~isfield(opts, 'networks') || ~islogical(opts.networks) || numel(opts.networks) ~= 1
    opts.networks = false;
end
if opts.networks
    n_node = 8;
else
    n_node = 268;
end
if ~isfield(opts, 'graphs') || ~islogical(opts.graphs) || numel(opts.graphs) ~= 1
    opts.graphs = true;
end
if ~isfield(opts, 'overlap') || ~islogical(opts.overlap) || numel(opts.overlap) ~= 1
    opts.overlap = true;
end
if ~isfield(opts, 'prep4conn') || ~islogical(opts.prep4conn) || numel(opts.prep4conn) ~= 1
    opts.prep4conn = false;
end
if ~isfield(opts, 'progress') || numel(opts.progress) ~= 1 || ...
   (~islogical(opts.progress) && ~isa(opts.progress, 'xprogress') && ~isxfigure(opts.progress, true))
    opts.progress = true;
end
if ~isfield(opts, 'robcorr') || ~islogical(opts.robcorr) || numel(opts.robcorr) ~= 1
    opts.robcorr = false;
end
if ~isfield(opts, 'robust') || ~islogical(opts.robust) || numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'savefpm') || ~ischar(opts.savefpm) || isempty(opts.savefpm)
    opts.savefpm = '';
else
    opts.savefpm = opts.savefpm(:)';
end
if ~isfield(opts, 'subsel') || ~iscell(opts.subsel)
    opts.subsel = mdmsubs;
else
    opts.subsel = opts.subsel(:);
    for ssc = numel(opts.subsel):-1:1
        if ~ischar(opts.subsel{ssc}) || isempty(opts.subsel{ssc}) || ...
           ~any(strcmpi(opts.subsel{ssc}(:)', mdmsubs))
            opts.subsel(ssc) = [];
        else
            opts.subsel{ssc} = opts.subsel{ssc}(:)';
        end
    end
end
subjids = opts.subsel;
if isempty(subjids)
    error('neuroelf:xff:badArgument', 'Invalid or missing subject IDs supplied.');
end
if any(strcmp(subjids, ''))
    error('neuroelf:xff:invalidObject', 'Invalid subject IDs for some subjects.');
end
if ~isfield(opts, 'thresh') || ~isa(opts.thresh, 'double') || numel(opts.thresh) ~= 1 || ...
    isinf(opts.thresh) || isnan(opts.thresh) || opts.thresh <= 0 || opts.thresh >= 0.5
    opts.thresh = 0.01;
end
if ~isfield(opts, 'updrtv') || ~islogical(opts.updrtv) || numel(opts.updrtv) ~= 1
    opts.updrtv = true;
end
if ~isfield(opts, 'validcov') || ~isa(opts.validcov, 'double') || isempty(opts.validcov) || ...
    any(isinf(opts.validcov(:)) | isnan(opts.validcov(:)))
    opts.validcov = [];
    opts.validfpm = [];
else
    opts.validcov = opts.validcov(:);
end
if ~isfield(opts, 'validfpm') || ~isa(opts.validfpm, 'double') || ndims(opts.validfpm) ~= 3 || ...
    isempty(opts.validfpm) || size(opts.validfpm, 1) ~= n_node || size(opts.validfpm, 2) ~= n_node || ...
    size(opts.validfpm, 3) ~= numel(opts.validcov)
    opts.validcov = [];
    opts.validfpm = [];
end
vtcsub = mdm_Subjects(xo, true);
vtcmatch = multimatch(vtcsub, subjids);
vtcidx = find(vtcmatch > 0);
numsubs = numel(subjids);
numvtcs = numel(vtcidx);

% determine progress bar capabilities
try
    closepbar = false;
    pbarstep = 1 / 86400;
    pbarnext = now + pbarstep; 
    if islogical(opts.progress)
        if opts.progress
            pbar = xprogress;
            xprogress(pbar, 'settitle', 'Computing finger-print matrices...');
            xprogress(pbar, 0, 'Processing VTCs...', 'visible', 0, 1);
            closepbar = true;
        else
            pbar = [];
        end
    else
        pbar = opts.progress;
        pbarvis = pbar.Visible;
        pbar.Progress(0, 'Processing VTCs...');
        pbar.Visible = 'on';
    end
catch xfferror
    pbar = [];
    neuroelf_lasterr(xfferror);
end

% prepare finger-print matrices (and covariate)
if atype == 'b'
    fpm = zeros(n_node, n_node, numsubs);

    % limit covariate to selected subjects
    if numel(opts.cov) == numel(mdmsubs) && numsubs < numel(mdmsubs)
        opts.cov = opts.cov(multimatch(mdmsubs, subjids) > 0);
    end
else
    fpm = zeros(n_node, n_node, numvtcs);
end

% generate/collect fingerprints
vtc = {[], []};
for vc = 1:numel(vtcidx)

    % allow for errors
    try

        % progress
        vtcfile = vtcs{vtcidx(vc)};
        [vtcpath, vtcfilesh] = fileparts(vtcfile);
        if numel(vtcfilesh) > 32
            vtcfilesh = [vtcfilesh(1:14) '...' vtcfilesh(end-13:end)];
        end
        
        % load VTC
        vtc{1} = xff(vtcfile, 't');
        vtc{2} = [];

        % use from RunTimeVars
        if opts.fromrtv && isfield(vtc{1}.C.RunTimeVars, 'FingerPrintParcels') && ...
            isa(vtc{1}.C.RunTimeVars.FingerPrintNetworks, 'double') && ...
            all(size(vtc{1}.C.RunTimeVars.FingerPrintNetworks) == 8) && ...
            isa(vtc{1}.C.RunTimeVars.FingerPrintParcels, 'double') && ...
            all(size(vtc{1}.C.RunTimeVars.FingerPrintParcels) == 268)

            % progress
            if ~isempty(pbar) && now >= pbarnext
                pbar.Progress((vc - 1) / numel(vtcidx), ['Reading ' vtcfilesh '.vtc']);
                pbarnext = now + pbarstep;
            end

            % get from RunTimeVars
            if opts.networks
                fp = vtc{1}.C.RunTimeVars.FingerPrintNetworks;
            else
                fp = vtc{1}.C.RunTimeVars.FingerPrintParcels;
            end

        % generate
        else

            % progress
            if ~isempty(pbar)
                pbar.Progress((vc - 1) / numel(vtcidx), ['Processing ' vtcfilesh '.vtc']);
                pbarnext = now;
            end

            % prepare
            if istransio(vtc{1}.C.VTCData)
                vtc{1}.C.VTCData = resolve(vtc{1}.C.VTCData);
            end
            if opts.prep4conn
                vtc{2} = vtc{1};
                vtc{1} = [];
                vtc{1} = vtc_PrepForConnectivity(vtc{2}, opts);
            end

            % compute fingerprint matrix
            fp = aft_ShenFingerPrint(vtc{1}, opts);

            % store
            if opts.updrtv
                try
                    if isempty(vtc{2})
                        aft_SaveRunTimeVars(vtc{1});
                    else
                        vtc{2}.C.RunTimeVars.FingerPrintNetworks = ...
                            vtc{1}.C.RunTimeVars.FingerPrintNetworks;
                        vtc{2}.C.RunTimeVars.FingerPrintNetworksTC = ...
                            vtc{1}.C.RunTimeVars.FingerPrintNetworksTC;
                        vtc{2}.C.RunTimeVars.FingerPrintParcels = ...
                            vtc{1}.C.RunTimeVars.FingerPrintParcels;
                        vtc{2}.C.RunTimeVars.FingerPrintParcelsTC = ...
                            vtc{1}.C.RunTimeVars.FingerPrintParcelsTC;
                        aft_SaveRunTimeVars(vtc{2});
                    end
                catch xfferror
                    fprintf('Error updating RunTimeVars for %s: %s.\n', vtcfile, xfferror.message);
                end
            end
        end
        
        % fisherize
        nanfp = isnan(fp);
        fp(nanfp) = 0;
        fp = fisherr2z(fp);
        fp(nanfp) = NaN;

        % sum to storage
        if atype == 'b'
            fpm(:, :, vtcmatch(vtcidx(vc))) = fpm(:, :, vtcmatch(vtcidx(vc))) + fp;
        else
            fpm(:, :, vc) = fp;
        end

        % clear object
        clearxffobjects(vtc);

    % deal with errors
    catch xfferror;
        clearxffobjects(vtc);
        if ~isempty(pbar)
            if closepbar
                closebar(pbar);
            else
                pbar.Visible = pbarvis;
            end
        end
        rethrow(xfferror);
    end
end

% normalize (average)
if atype == 'b'
    usub = unique(vtcmatch(vtcidx));
    for sc = 1:numel(usub)
        fpm(:, :, usub(sc)) = (1 / sum(vtcmatch == usub(sc))) .* fpm(:, :, usub(sc));
    end
end

% save matrices
if ~isempty(opts.savefpm)
    try
        save(opts.savefpm, fpm);
    catch xfferror
        fprintf('Error saving FingerPrint matrices: %s.\n', xfferror.message);
    end
end

% new progress bar
if ~isempty(pbar) && closepbar
    closebar(pbar);
    pbar = xprogress;
    if atype == 'i'
        xprogress(pbar, 'settitle', 'Performing Fingerprint identity analysis...');
    else
        xprogress(pbar, 'settitle', 'Generating Fingerprint->Behavior model...');
    end
    xprogress(pbar, 0, 'Regressing fingerprints...', 'visible', 0, 1);
end

% for identity
if atype == 'i'

    % compute full cross correlation martix of FPMs
    fpmi = true(n_node * n_node, 1);
    fpmi(1:(n_node+1):end) = false;
    fpmm = reshape(fpm, numel(fpmi), numvtcs);
    if opts.robust
        ccm = robcorrcoef(fpmm(fpmi, :));
    else
        ccm = corrcoef(fpmm(fpmi, :));
    end

    % for each subject
    for sc = 1:numsubs

        % compute average within-subject matrix correlation
    end

% for behavioral analysis
% adapted from _files/shenparcel/external_validation_share.m
% Copyright 2015 Monica Rosenberg, Emily Finn, and Dustin Scheinost
else

    % training data 
    train_mats = fpm;
    behav = opts.cov;
    n_sub = size(train_mats,3);
    n_train_sub = n_sub - 1;

    % validation data
    if isempty(opts.validcov)
        valid_behav = opts.cov;
        valid_mats = fpm;
    else
        valid_behav = opts.validcov;
        valid_mats = opts.validfpm;
    end
    n_valid_sub = size(valid_mats, 3);

    aa = ones(n_node, n_node);
    if opts.robcorr
        aa_upp = aa;
        aa_upp(1:(n_node+1):end) = 0;
    else
        aa_upp = triu(aa, 1);
    end
    upp_id = find(aa_upp);
    n_edge = numel(upp_id);

    % with leave-one-out
    if opts.overlap

        pos_mask_all = false(n_node, n_node, n_sub);
        neg_mask_all = false(n_node, n_node, n_sub);
        upp_X = [ones(n_train_sub, 1), zeros(n_train_sub, 1)];

        % iterate over leave-one-out subjects
        for excl_sub = 1:n_sub

            % progress
            if ~isempty(pbar)
                pbar.Progress((excl_sub - 1) / n_sub, ...
                    ['Leave-one-out behavioral model (' subjids{excl_sub} ')...']);
            end

            % exclude data from left-out subject
            train_mats_tmp = train_mats;
            train_mats_tmp(:, :, excl_sub) = [];
            train_behav = behav;
            train_behav(excl_sub) = [];

            % create n_train_sub x n_edge matrix
            train_vect = reshape(train_mats_tmp, n_node * n_node, n_train_sub)';
            upp_vect = train_vect(:, upp_id);

            % relate behavior to edge strength across training subjects
            cp = zeros(n_edge, 1);
            cr = zeros(n_edge, 1);

            % for robust
            if opts.robust
                for edgec = 1:n_edge
                    if any(isinf(upp_vect(:, edgec)) | isnan(upp_vect(:, edgec)))
                        continue;
                    end
                    upp_X(:, 2) = upp_vect(:, edgec);
                    [rb, rr, rw, rstats] = fitrobustbisquare(upp_X, train_behav);
                    cp(edgec) = rstats.p(2);
                    cr(edgec) = rstats.t(2);
                end

                % convert t to r
                cr = sign(cr) .* sqrt((cr .* cr ./ (n_train_sub - 2)) ./ ...
                    (1 + (cr .* cr ./ (n_train_sub - 2))));

            % for OLS (corrcoef)
            else
                upp_vect(:, any(isinf(upp_vect) | isnan(upp_vect))) = 0;
                [cv, cr] = cov_nd(upp_vect', repmat(train_behav(:)', n_edge, 1));
                cr(all(upp_vect == 0, 1)) = 0;
                cp = correlpvalue(cr, n_train_sub);
            end

            % select edges based on threshold
            pos_edge = false(1, n_edge);
            neg_edge = false(1, n_edge);

            pos_edge(cp < opts.thresh & cr > 0) = true;
            neg_edge(cp < opts.thresh & cr < 0) = true;

            pos_mask = false(n_node, n_node);
            neg_mask = false(n_node, n_node);

            pos_mask(upp_id) = pos_edge;
            neg_mask(upp_id) = neg_edge;

            % save masks from every iteration of the leave-one-out loop
            pos_mask_all(:, :, excl_sub) = pos_mask;
            neg_mask_all(:, :, excl_sub) = neg_mask;
        end

        % select edges that appear in all iteraitons of leave-one-out
        pos_overlap = double(all(pos_mask_all, 3));
        neg_overlap = double(all(neg_mask_all, 3));

    % all-in-one (not leave-one-out)
    else

        % create n_sub x n_edge matrix
        train_vect = reshape(train_mats, n_node * n_node, n_sub)';
        upp_vect = train_vect(:, upp_id);
        upp_X = [ones(n_sub, 1), zeros(n_sub, 1)];

        % relate behavior to edge strength across ALL subjects in training set
        cp = zeros(n_edge, 1);
        cr = zeros(n_edge, 1);

        if opts.robust
            for edgec = 1:n_edge
                if any(isinf(upp_vect(:, edgec)) | isnan(upp_vect(:, edgec)))
                    continue;
                end
                upp_X(:, 2) = upp_vect(:, edgec);
                [rb, rr, rw, rstats] = fitrobustbisquare(upp_X, behav);
                cp(edgec) = rstats.p(2);
                cr(edgec) = rstats.t(2);
            end
            cr = sign(cr) .* sqrt((cr .* cr ./ (n_sub - 2)) ./ ...
                (1 + (cr .* cr ./ (n_sub - 2))));
        else
            upp_vect(:, any(isinf(upp_vect) | isnan(upp_vect))) = 0;
            [cv, cr] = cov_nd(upp_vect', repmat(behav(:)', n_edge, 1));
            cr(all(upp_vect == 0, 1)) = 0;
            cp = correlpvalue(cr, n_sub);
        end

        % select edges based on threshold
        pos_edge = false(1, n_edge);
        neg_edge = false(1, n_edge);

        pos_edge(cp < opts.thresh & cr > 0) = true;
        neg_edge(cp < opts.thresh & cr < 0) = true;

        pos_mask = false(n_node, n_node);
        neg_mask = false(n_node, n_node);

        pos_mask(upp_id) = pos_edge;
        neg_mask(upp_id) = neg_edge;

        pos_overlap = double(pos_mask);
        neg_overlap = double(neg_mask);
    end

    % store in stats
    stats = struct('fullmodel', pos_overlap + neg_overlap, ...
        'posmodel', pos_overlap, 'negmodel', neg_overlap);

    % sum edges for all subjects in the training set
    train_pos_sum = zeros(n_sub, 1);
    train_neg_sum = zeros(n_sub, 1);

    % remove Inf/NaNs from matrices to allow summing
    badmats = any(isinf(train_mats) | isnan(train_mats) | train_mats == 0, 3) | ...
        any(isinf(valid_mats) | isnan(valid_mats) | valid_mats == 0, 3);
    train_mats(repmat(badmats, [1, 1, n_sub])) = 0;
    for k = 1:n_sub
        train_pos_sum(k) = sum(sum(pos_overlap .* train_mats(:, :, k)));
        train_neg_sum(k) = sum(sum(neg_overlap .* train_mats(:, :, k)));
    end
    stats.train_behav = behav(:);
    stats.train_fpm = train_mats;
    stats.full_reg = [ones(n_sub, 1), train_pos_sum, train_neg_sum];
    stats.pos_reg = [ones(n_sub, 1), train_pos_sum];
    stats.neg_reg = [ones(n_sub, 1), train_neg_sum];

    % build model with training data
    if opts.robust
        glm_fit = fitrobustbisquare(stats.full_reg, behav(:));
        b_pos = fitrobustbisquare(stats.pos_reg, behav(:));
        b_neg = fitrobustbisquare(stats.neg_reg, behav(:));
    else
        glm_fit = calcbetas(stats.full_reg, behav(:));
        b_pos = calcbetas(stats.pos_reg, behav(:));
        b_neg = calcbetas(stats.neg_reg, behav(:));
    end
    stats.full_betas = glm_fit;
    stats.pos_betas = b_pos;
    stats.neg_betas = b_neg;

    % generate predictions for validation set
    pred_pos = zeros(n_valid_sub, 1);
    pred_neg = zeros(n_valid_sub, 1);
    pred_glm = zeros(n_valid_sub, 1);

    valid_pos_sum = zeros(n_valid_sub, 1);
    valid_neg_sum = zeros(n_valid_sub, 1);

    valid_mats(repmat(badmats, [1, 1, n_valid_sub])) = 0;
    for vs = 1:n_valid_sub
        valid_pos_sum(vs) = sum(sum(pos_overlap .* valid_mats(:, :, vs)));
        valid_neg_sum(vs) = sum(sum(neg_overlap .* valid_mats(:, :, vs)));
        pred_glm(vs) = glm_fit(1) + glm_fit(2) * valid_pos_sum(vs) + ...
            glm_fit(3)*valid_neg_sum(vs);
        pred_pos(vs) = (b_pos(2) * valid_pos_sum(vs)) + b_pos(1);
        pred_neg(vs) = (b_neg(2) * valid_neg_sum(vs)) + b_neg(1);
    end

    % correlate predicted and observed behavior
    stats.valid_behav = valid_behav(:);
    stats.valid_fpm = valid_mats;
    [stats.valid_r_full, stats.valid_p_full] = corr(valid_behav(:), pred_glm);
    [stats.valid_r_pos, stats.valid_p_pos] = corr(valid_behav(:), pred_pos);
    [stats.valid_r_neg, stats.valid_p_neg] = corr(valid_behav(:), pred_neg);

    % graphing
    if opts.graphs
        h.fig = figure;
        h.ax1 = subplot(2,2,1);
        scatter(valid_behav, pred_pos);
        title(['pos r = ', num2str(round(stats.valid_r_pos * 100) / 100), ...
            ', p = ' num2str(stats.valid_p_pos)]);
        xlabel(['Observed ' bname]);
        ylabel(['Predicted ' bname]);
        h.ax2 = subplot(2,2,2);
        scatter(valid_behav, pred_neg);
        title(['neg r = ', num2str(round(stats.valid_r_neg * 100) / 100), ...
            ', p = ' num2str(stats.valid_p_neg)]);
        xlabel(['Observed ' bname]);
        ylabel(['Predicted ' bname]);
        h.ax3 = subplot(2,2,3);
        scatter(valid_behav, pred_glm);
        title(['glm r = ', num2str(round(stats.valid_r_full * 100) / 100), ...
            ', p = ' num2str(stats.valid_p_full)]);
        xlabel(['Observed ' bname]);
        ylabel(['Predicted ' bname]);
    end
end

% progress bar
if ~isempty(pbar)
    if closepbar
        closebar(pbar);
    else
        pbar.Visible = pbarvis;
    end
end

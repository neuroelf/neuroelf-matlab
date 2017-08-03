function [betas, irtc, ptc, se] = vvd_CalcROIBetas(xo, mdm, options)
% VVD::CalcROIBetas  - calculate beta values for ROIs
%
% FORMAT:       [betas, irtc, ptc, se] = vvd.CalcROIBetas(mdm [, options]);
%
% Input fields:
%
%       mdm         mdm xff object (must match with VVD's VTCs)
%       options     1x1 struct with optional parameters
%        .autoreg   auto-regression factor ({0}, 1 .. 16)
%        .tctrans   time course transformation ('none', {'psc'}, 'z')
%
% Output fields:
%
%       betas       1x1 struct with ROIs as fieldnames and Px1xS betas
%       irtc        equal struct with PxPxS inverse matrices
%       ptc         equal struct with predicted time course (as in VVD)
%       se          equal struct with 1x1xS standard errors
%                   where each P := number of predictors and
%                   S := number of studies/subjects
%
% Using: makelabel, orthvec, psctrans, ztrans.

% to-do
%
%        .seppred   separate predictors ('no', 'subject', {'study'})

% Version:  v1.1
% Build:    16020314
% Date:     Feb-03 2016, 2:51 PM EST
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

% neuroelf library and singleton factory
global ne_methods xffsngl;

% argument check
if nargin < 2 || numel(xo) ~= 1 || numel(mdm) ~= 1 || ...
   ~xffisobject(xo, true, 'vvd') || ~xffisobject(mdm, true, 'mdm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
mdmc = mdm.C;

% check options
if nargin < 3 || ~isstruct(options) || numel(options) ~= 1
    options = struct;
end

% auto-regression
autoreg = 0;
if isfield(options, 'autoreg') && isa(options.autoreg, 'double') && ...
    numel(options.autoreg) == 1 && ~isinf(options.autoreg) && ...
   ~isnan(options.autoreg) && options.autoreg > 0
    autoreg = min(16, round(options.autoreg));
end

% separate predictors
seppred = 2;
if isfield(options, 'seppred') && ischar(options.seppred) && ~isempty(options.seppred)
    switch (lower(options.seppred(:)'))
        case 'no'
            seppred = 0;
        case 'subject'
            seppred = 1;
        case 'study'
            seppred = 2;
        otherwise
            error('neuroelf:xff:invalidOption', ...
                'Invalid option for seppred field: %s.', options.seppred(:)');
    end
end

% transformation
tctrans = 1;
if isfield(options, 'tctrans') && ischar(options.tctrans) && ~isempty(options.tctrans)
    switch (lower(options.tctrans(:)'))
        case 'none'
            tctrans = 0;
        case 'psc'
            tctrans = 1;
            psctrans = ne_methods.psctrans;
        case 'z'
            tctrans = 2;
            ztrans = ne_methods.ztrans;
        otherwise
            error('neuroelf:xff:invalidOption', ...
                'Invalid option for tctrans field: %s.', options.tctrans(:)');
    end
end

% prepare ROI names
rois = struct;
vnames = bc.VOINames;
makelabel = ne_methods.makelabel;
for rc = 1:numel(vnames)
    rname = makelabel(vnames{rc});
    if isfield(rois, rname)
        warning('neuroelf:xff:ambiguousData', 'Ambiguous VOI name: ''%s''.', vnames{rc});
        rname = makelabel(sprintf('%s_%06d', ...
            rname(1:min(24, numel(rname))), fix(1e6 * rand(1, 1))));
    end
    rois.(rname) = rc;
end
roil = fieldnames(rois);
nroi = numel(roil);

% check VTCs
vvd_vtc = bc.VTC;
vvd_vtcs = {vvd_vtc(:).FileName};
nvtc = numel(vvd_vtcs);
mdm_vtcs = mdmc.XTC_RTC(:, 1)';
for vc = 1:numel(mdm_vtcs)
    [vtcf{1:3}] = fileparts(vvd_vtcs{vc});
    vvd_vtcs{vc} = [vtcf{2}, vtcf{3}];
    [vtcf{1:3}] = fileparts(mdm_vtcs{vc});
    mdm_vtcs{vc} = [vtcf{2}, vtcf{3}];
end
mdm_rtcs = mdmc.XTC_RTC(:, end)';
vtcl = struct;
vtcr = struct;
rls = xffsngl.CONF.reloadsame;
xffsngl.CONF.reloadsame = true;
ntpt = 0;
for vc = 1:nvtc
    vtcn = vvd_vtcs{vc};
    vtcm = find(strcmpi(vtcn, mdm_vtcs));
    if isempty(vtcm)
        xffsngl.CONF.reloadsame = rls;
        error('neuroelf:xff:argumentMismatch', ...
            'Required VTC entry not in MDM definition file: ''%s''.', vtcn);
    end
    [vtcf{1:2}] = fileparts(vtcn);
    vtcl.(makelabel(vtcf{2})) = vtcm(1);
    rrtc = mdm_rtcs{vtcm(1)};
    [rtcf{1:3}] = fileparts(rrtc);
    if exist([rtcf{2}, rtcf{3}], 'file') == 2
        rrtc = [rtcf{2}, rtcf{3}];
    end
    if exist(rrtc, 'file') ~= 2
        xffsngl.CONF.reloadsame = rls;
        error('neuroelf:xff:fileNotFound', 'Required SDM file not found: ''%s''.', rrtc);
    end
    rtcl = makelabel(rtcf{2});
    if ~any(strcmpi(rtcl, fieldnames(vtcr)))
        try
            trtc = xff(rrtc);
            vtcr.(rtcl) = getcont(trtc);
            xff(0, 'clearobj', trtc.L);
            rrtc = vtcr.(rtcl);
            if bc.VTC(vc).NrOfVolumes ~= rrtc.NrOfDataPoints
                error('neuroelf:xff:sizeMismatch', ...
                    'Timecourse/SDM length mismatch for ''%s''.', vtcn);
            end
            rtcn = rrtc.SDMMatrix;
            rtcm = mean(rtcn);
            rtcv = var(rtcn);
            if ~any(rtcm ~= 0 & rtcv == 0)
                rtcn(:, end+1) = 1;
                vtcr.(rtcl).SDMMatrix = rtcn;
            end
            rtct = rtcn';
            rtci = pinv(rtct * rtcn);
            vtcr.(rtcl).InverseMatrix = rtci;
            vtcr.(rtcl).ProductMatrix = rtci * rtct;
        catch xfferror
            xffsngl.CONF.reloadsame = rls;
            error('neuroelf:xff:errorReadingFile', ...
                'SDM could not be read/parsed: ''%s'' (%s).', rrtc, xfferror.message);
        end
    end
    ntpt = ntpt + bc.VTC(vc).NrOfVolumes;
end
xffsngl.CONF.reloadsame = rls;

% check SDMs
rtcl = fieldnames(vtcr);
for rc = 1:numel(rtcl)
    if rc == 1
        np = size(vtcr.(rtcl{rc}).SDMMatrix, 2);
    elseif np ~= size(vtcr.(rtcl{rc}).SDMMatrix, 2)
        error('neuroelf:xff:badFileContents', ...
            'The given SDMs must match in number of predictors.');
    end
end

% calculate statistics per ROI
glmr = cell(nvtc, nroi);
tcf = 1;
tcs = 0;
orthvec = ne_methods.orthvec;
for rc = 1:nroi

    % iterate over VTCs
    for vc = 1:nvtc

        % get correct SDM
        [vtcf{1:2}] = fileparts(vvd_vtcs{vc});
        vtcm = vtcl.(makelabel(vtcf{2}));
        rrtc = mdm_rtcs{vtcm};
        [rtcf{1:2}] = fileparts(rrtc);
        rrtc = vtcr.(makelabel(rtcf{2}));
        irtc = rrtc.InverseMatrix;
        prtc = rrtc.ProductMatrix;
        rrtc = rrtc.SDMMatrix;

        % get time course
        tc = bc.VTC(vc).Values(:, rc);
        df1 = size(rrtc, 1) - 1;
        df2 = df1 + 1 - size(rrtc, 2) - autoreg;

        % transform tc
        switch (tctrans)

            % PSC
            case 1
                [tc, tcf] = psctrans(tc);

            % z
            case 2
                [tc, tcf, tcs] = ztrans(tc);
        end

        % calculate betas
        betas = prtc * tc(:);

        % calculate predicted timecourse
        ptc = rrtc * betas;

        % calculate residuals
        resid = tc - ptc;

        % auto-regression
        for ac = 1:autoreg
            resid((ac + 1):end) = orthvec(resid((ac + 1):end), resid(1:end - ac));
        end

        % calculate standard error
        se = std(resid) * sqrt(df1/df2);

        % re-transform ptc
        if tctrans > 0
            ptc = ptc / tcf + tcs;
        end
        glmr{vc, rc} = {betas, irtc, ptc, se};
    end
end

% prepare output
betas = struct;
irtc  = struct;
ptc   = struct;
se    = struct;

% which separation mode
switch (seppred)

    % no separation at all
    case 0

    % subject separation
    case 1

    % study separation
    case 2

        % put results into output directly
        for rc = 1:nroi

            % get roi label and prepare structs
            lroi = roil{rc};
            betas.(lroi) = zeros(np,  1, nvtc);
            irtc.(lroi)  = zeros(np, np, nvtc);
            ptc.(lroi)   = zeros(ntpt, 1);
            se.(lroi)    = zeros( 1,  1, nvtc);

            % put results in output
            ttpt = 1;
            for vc = 1:nvtc
                betas.(lroi)(:, 1, vc) = glmr{vc, rc}{1};
                irtc.(lroi)(:, :, vc) = glmr{vc, rc}{2};
                pptc = glmr{vc, rc}{3};
                lptc = numel(pptc);
                ptc.(lroi)(ttpt:(ttpt + lptc - 1)) = pptc;
                ttpt = ttpt + lptc;
                se.(lroi)(vc) = glmr{vc, rc}{4};
            end
        end
end

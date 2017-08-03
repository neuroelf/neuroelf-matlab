function vtcc = vtccondaverage(ovtcs, opts)
% vtccondaverage  - create group average of per-subject average VTCs
%
% FORMAT:       gvtc = vtccondaverage(svtcs [, opts])
%
% Input fields:
%
%       svtcs       per-subject average VTCs (cell array)
%       opts        optional settings
%        .grouping  either of {'ffx'}, 'mixed', 'rfx', 'wrfx'
%        .mask      use a mask other than the one in the VTCs
%        .robust    perform robust 2nd level stats (RFX only)
%
% Output fields:
%
%       gvtc        group-average VTC

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:31 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012 - 2014, 2016, Jochen Weber
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
if nargin < 1 || ...
   ~iscell(ovtcs) || ...
    isempty(ovtcs) || ...
   ~all(cellfun(@ischar, ovtcs(:))) || ...
    any(cellfun('isempty', regexpi(ovtcs(:), '\.vtc$')))
    error( ...
        'neuroelf:BadArgument', ...
        'First argument must be a list of VTC filenames.' ...
    );
end
ovtcs = ovtcs(:);
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'grouping') || ...
   ~ischar(opts.grouping) || ...
    isempty(opts.grouping) || ...
   ~any(lower(opts.grouping(1)) == 'fmrw')
    opts.grouping = 'f';
else
    opts.grouping = lower(opts.grouping(1));
end
if ~isfield(opts, 'groupstc') || ...
   ~ischar(opts.groupstc) || ...
    numel(opts.groupstc) ~= 1
    opts.groupstc = 'd';
end
if ~isfield(opts, 'mask') || ...
    isempty(opts.mask) || ...
   ~islogical(opts.mask) || ...
    ndims(opts.mask) ~= 3
    opts.mask = [];
end
if ~isfield(opts, 'robust') || ...
   ~islogical(opts.robust) || ...
    numel(opts.robust) ~= 1
    opts.robust = false;
end

% then re-load the per-subject VTCs
vtco = ovtcs;
vtcn = numel(vtco);
try
    for sc = 1:vtcn
        vtco{sc} = xff(ovtcs{sc}, 't');
        if sc == 1
            vtcc = vtco{sc}.CopyObject;
            vtcc.VTCData = vtcc.VTCData(:, :, :, :);
            if isempty(opts.mask)
                opts.mask = squeeze(any(vtcc.VTCData ~= 0, 1));
            end
        end
        vtcod = getcont(vtco{sc});
        vtco{sc}.ClearObject;
        vtco{sc} = vtcod;
    end
catch ne_eo;
    clearxffobjects(vtco);
    rethrow(ne_eo);
end
mask12 = squeeze(sum(sum(opts.mask, 1), 2));
cmask12 = [0; cumsum(mask12(:))];

% get values
vsz = size(vtcc.VTCData);
vsz(1) = [];
vsz12 = vsz(1) * vsz(2);
dimt = vtcc.NrOfVolumes;
dima = vtcc.RunTimeVars.NrOfVolumesPerTC;
dimc = vtcc.RunTimeVars.NrOfConditions;
dimtc = vtcc.RunTimeVars.NrOfTCsPerCondition;

% generate new progress bar
try
    atotal = sum(mask12) * vtcn;
    acount = 0;
    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 640, 36]);
    xprogress(pbar, 'settitle', 'Group-averaging subject VTCs...');
    xprogress(pbar, 0, 'Averaging...', 'visible', 0, atotal);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    pbar = [];
end

% create data to hold the per-subjects data (data, error, number of onsets)
dimac = dima * dimc;

% indices
cis = dimtc * dima;
tci = lsqueeze(repmat((1:dima)', 1, dimc) + ones(dima, 1) * (0:cis:(dimt-1)));
sci = tci + dima;
wci = sci + dima;
rma = [1, 1, vtcn];

% work along 3rd dimension
for slc = 1:vsz(3)
%for slc = 34

    % progress
    if ~isempty(pbar)
        xprogress(pbar, acount, sprintf('Averaging slice %d/%d...', slc, vsz(3)));
    end

    % nothing to do
    if cmask12(slc) == cmask12(slc+1)
        continue;
    end

    % scrap data
    smask = opts.mask(:, :, slc);
    smask = smask(:);
    ssmask = sum(smask);
    sdata = zeros([dimac, ssmask, vtcn]);
    edata = sdata;
    wdata = sdata;
    onscount = zeros(1, dimc);
    tcoweights = zeros(dimt, 1);

    % get data for subjects
    for sc = 1:vtcn

        % get VTC object
        svtcc = vtco{sc};

        % add to count of onsets
        onscount = onscount + svtcc.RunTimeVars.NrOfConditionOnsets;

        % get data
        slicew = svtcc.VTCData(:, :, :, slc);
        sliced = reshape(slicew(tci, :, :), dimac, vsz12);
        slicee = reshape(slicew(sci, :, :), dimac, vsz12);
        tcoweights = tcoweights + svtcc.RunTimeVars.TCOnsetWeights;
        if opts.robust
            slicew = repmat(svtcc.RunTimeVars.TCOnsetWeights(wci), 1, vsz12) .* ...
                reshape(slicew(wci, :, :), dimac, vsz12);
        else
            slicew = repmat(svtcc.RunTimeVars.TCOnsetWeights(tci), 1, vsz12);
        end

        % add to array
        sdata(:, :, sc) = sliced(:, smask);
        edata(:, :, sc) = slicee(:, smask);
        wdata(:, :, sc) = slicew(:, smask);
        if ~isempty(pbar)
            acount = acount + mask12(slc);
            xprogress(pbar, acount);
        end
    end

    % keep track of total weights (for d.f.)
    swdata = sum(wdata, 3);

    % fixed effects
    if opts.grouping == 'f'

        % normalize weights (to sum to 1)
        nwdata = wdata ./ repmat(swdata, rma);

        % weighted mean
        mdata = sum(nwdata .* sdata, 3);

        % standard deviation of weighted data
        sddata = sqrt(sum(nwdata .* edata .* edata, 3));

        % compute error metric (if not SD)
        if opts.groupstc ~= 'd'

            % standard error (defined by sqrt of sum of weights)
            sddata = sddata ./ sqrt(swdata);
        end

    % for midex effects
    elseif opts.grouping == 'm'

        % compute standard error (reciprocal as weight)
        wedata = sqrt(wdata) ./ edata;

        % normalize weights to 1
        nwdata = wedata ./ repmat(sum(wedata, 3), rma);

        % weighted mean
        mdata = sum(nwdata .* sdata, 3);

        % deviances (with re-normalized weights)
        rnwdata = vtcn .* nwdata;
        ddata = sdata - repmat(mdata, rma);

        % standard deviation
        sddata = sqrt(sum(rnwdata .* ddata .* ddata, 3) ./ (max(0, sum(rnwdata, 3) - 1)));

        % compute error metric
        if opts.groupstc ~= 'd'

            % standard error
            sddata = sddata ./ sqrt(sum(rnwdata, 3));
        end

    % RFX
    elseif opts.grouping == 'r'

        % compute plain mean
        mdata = mean(sdata, 3);

        % and standard deviation
        sddata = std(sdata, [], 3);

        % compute error metric
        if opts.groupstc ~= 'd'

            % standard error
            sddata = (1 / sqrt(vtcn)) .* sddata;
        end

    % weighted RFX
    else

        % normalize weights (to max of 1)
        nwdata = wdata ./ repmat(max(wdata, [], 3), rma);

        % weighted mean
        mdata = sum(nwdata .* sdata, 3);

        % deviances
        ddata = sdata - repmat(mdata, rma);

        % standard deviation
        sddata = sqrt(sum(nwdata .* ddata .* ddata, 3) ./ (sum(nwdata, 3) - 1));

        % compute error metric
        if opts.groupstc ~= 'd'

            % standard error
            sddata = sddata ./ sqrt(sum(nwdata, 3));
        end
    end

    % convert to t-score
    if opts.groupstc == 't'
        sddata = mdata ./ sddata;
        sddata(isinf(sddata) | isnan(sddata)) = 0;
    end

    % re-normalize to one
    if opts.robust
        swdata = swdata ./ (tcoweights(wci) * ones(1, ssmask));
    else
        swdata = swdata ./ (tcoweights(tci) * ones(1, ssmask));
    end

    % store in full-size data
    mmdata = zeros(dimac, vsz12);
    mmdata(:, smask) = mdata;
    msddata = zeros(dimac, vsz12);
    msddata(:, smask) = sddata;

    % store in VTC
    vtcc.NrOfConditionOnsets = onscount;
    vtcc.VTCData(tci, :, :, slc) = reshape(mmdata, [dimac, vsz(1:2)]);
    vtcc.VTCData(sci, :, :, slc) = reshape(msddata, [dimac, vsz(1:2)]);
    if opts.robust
        mswdata = zeros(dimac, vsz12);
        mswdata(:, smask) = swdata;
        vtcc.VTCData(wci, :, :, slc) = reshape(mswdata, [dimac, vsz(1:2)]);
    end
    vtcc.RunTimeVars.TCOnsetWeights = tcoweights;
end

% clear bar
if ~isempty(pbar)
    closebar(pbar);
end

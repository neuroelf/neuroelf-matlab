function vmp = vtc_AveragePeak(xo, opts)
% VTC::AveragePeakTime  - compute peak and time-to-peak values (voxel-wise)
%
% FORMAT:       vmp = vtc.AveragePeak([opts])
%
% Input fields:
%
%       opts        optional struct with fields
%        .conds     conditions (default: all available)
%        .resamptr  resampling TR (to infer maximum value, default: 100)
%        .tthresh   t-threshold which is used to mask voxels (default: 3)
%
% Output fields:
%
%       vmp         VMP object with peak and time-to-peak maps
%
% Using: flexinterpn_method, limitrangec, multimatch.

% Version:  v1.1
% Build:    16021412
% Date:     Feb-14 2016, 12:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2013, 2014, 2016, Jochen Weber
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
flexinterpn_method = ne_methods.flexinterpn_method;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc')
    error('neuroelf:xff:badArgument', 'Invalid object in call.');
end
bc = xo.C;
if ~isfield(bc.RunTimeVars, 'AvgVTC') || ~islogical(bc.RunTimeVars.AvgVTC) || ...
    numel(bc.RunTimeVars.AvgVTC) ~= 1 || ~bc.RunTimeVars.AvgVTC
    error('neuroelf:xff:badArgument', 'Only valid for Average-VTCs.');
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'conds') || ~iscell(opts.conds) || isempty(opts.conds)
    opts.conds = bc.RunTimeVars.ConditionNames;
else
    opts.conds = opts.conds(:);
    for cc = numel(opts.conds):-1:1
        if ~ischar(opts.conds{cc}) || isempty(opts.conds{cc})
            opts.conds(cc) = [];
        else
            opts.conds{cc} = opts.conds{cc}(:)';
        end
    end
    if numel(unique(opts.conds)) ~= numel(opts.conds)
        error('neuroelf:xff:badArgument', 'Conditions must be uniquely specified.');
    end
end
if ~isfield(opts, 'resamptr') || ~isa(opts.resamptr, 'double') || numel(opts.resamptr) ~= 1 || ...
    isinf(opts.resamptr) || isnan(opts.resamptr) || opts.resamptr <= 0
    opts.resamptr = 100;
else
    opts.resamptr = min(1000, ceil(opts.resamptr));
end
if ~isfield(opts, 'tthresh') || ~isa(opts.tthresh, 'double') || numel(opts.tthresh) ~= 1 || ...
    isinf(opts.tthresh) || isnan(opts.tthresh) || opts.tthresh < 0
    opts.tthresh = 3;
end
trrat = opts.resamptr / bc.TR;

% which conditions
condi = ne_methods.multimatch(opts.conds, bc.RunTimeVars.ConditionNames);
if any(condi < 1)
    error('neuroelf:xff:badArgument', 'Unknown condition name given.');
end

% data
vptc = bc.RunTimeVars.NrOfVolumesPerTC;
tcpc = bc.RunTimeVars.NrOfTCsPerCondition;
ncnd = bc.RunTimeVars.NrOfConditions;
sdse = lower(bc.RunTimeVars.TCNames{2});
tce = bc.RunTimeVars.AvgWindowTo;
tcs = 0:opts.resamptr:tce;

% get data
vtcd = bc.VTCData;
if size(vtcd, 1) ~= (vptc * tcpc * ncnd)
    error('neuroelf:xff:badObject', 'Invalid AverageVTC object.');
end
if istransio(vtcd)
    vtcd = resolve(vtcd);
end
vtcsz = size(vtcd);
vtcsz(1) = [];

% create VMP
vmp = xff('new:vmp');
vmpc = vmp.C;

% adapt settings
vmpc.Resolution = bc.Resolution;
vmpc.XStart = bc.XStart;
vmpc.XEnd = bc.XEnd;
vmpc.YStart = bc.YStart;
vmpc.YEnd = bc.YEnd;
vmpc.ZStart = bc.ZStart;
vmpc.ZEnd = bc.ZEnd;

% set correct number of maps
vmpc.Map(1).VMPData = single(zeros(vtcsz));
vmpc.Map = repmat(vmpc.Map(1), 1, 3 * numel(condi));
vmpc.NrOfMaps = numel(vmpc.Map);

% squeeze data
vtcd = reshape(vtcd, size(vtcd, 1), prod(vtcsz));

% iterate over conditions
for cc = 1:numel(condi)

    % set maps data
    mapc = 3 * cc - 2;
    vmpc.Map(mapc).Type = 15;
    vmpc.Map(mapc).Name = sprintf('%s peak value', opts.conds{cc});
    vmpc.Map(mapc).LowerThreshold = 0.2;
    vmpc.Map(mapc).UpperThreshold = 1.5;
    mapc = mapc + 1;
    vmpc.Map(mapc).Type = 15;
    vmpc.Map(mapc).Name = sprintf('%s time-to-peak value (s)', opts.conds{cc});
    vmpc.Map(mapc).LowerThreshold = 0;
    vmpc.Map(mapc).UpperThreshold = 0.001 * (tcs(end) + tcs(2) - tcs(1));
    mapc = mapc + 1;
    vmpc.Map(mapc).Type = 15;
    vmpc.Map(mapc).Name = sprintf('%s weighted AUC', opts.conds{cc});
    vmpc.Map(mapc).LowerThreshold = 1;
    vmpc.Map(mapc).UpperThreshold = 2;

    % get data indices
    cci = condi(cc);
    di = (1:vptc) + ((cci - 1) * vptc * tcpc);

    % get actual data
    mdata = vtcd(di, :);

    % mask
    mmask = any(mdata ~= 0, 1);
    maski = find(mmask);
    mdata = mdata(:, mmask);

    % errors (and weights)
    edata = vtcd(di + vptc, mmask);
    wdata = bc.RunTimeVars.NrOfConditionOnsets(cci);
    if tcpc > 2 && strcmp(sdse, 'sd')
        wdata = wdata .* vtcd(di + 2 * vptc, mmask);
    end

    % compute error metric
    if strcmp(sdse, 'sd')
        edata = edata ./ sqrt(wdata);
        edata(isinf(edata) | isnan(edata)) = Inf;
    end
    if sdse(1) == 's'
        edata = mdata ./ edata;
        edata(isinf(edata) | isnan(edata)) = 0;
    end

    % decide on direction, mean t must be > threshold
    emask = ne_methods.limitrangec((1 / opts.tthresh) .* mean(edata), -1, 1, 0);
    emask(abs(emask) ~= 1) = 0;

    % re-mask
    mmask(maski(emask == 0)) = false;
    maski(emask == 0) = [];
    maskn = numel(maski);
    if maskn == 0
        continue;
    end
    mdata = mdata(:, emask ~= 0);
    edata = edata(:, emask ~= 0);
    emask(emask == 0) = [];
    pmask = emask > 0;
    nmask = ~pmask;

    % peak value and time to peak
    peakv = zeros(1, maskn, 1);
    peakt = zeros(1, maskn, 1);
    aucv = zeros(1, maskn, 1);
    aucw = zeros(1, maskn, 1);

    % iterate over time
    for tc = 1:numel(tcs)

        % sample value
        scoord = [Inf, Inf; 1 + (tc - 1)* trrat, 1; 2 * vptc, 1; size(mdata)];
        smdata = flexinterpn_method(mdata, scoord, 'cubic');
        sedata = flexinterpn_method(edata, scoord, 'cubic');

        % higher than threshold
        semask = abs(sedata) >= opts.tthresh;

        % for those compare
        sehigher = semask & pmask & (smdata > peakv);
        selower  = semask & nmask & (smdata < peakv);

        % update values and time
        peakv(sehigher) = smdata(sehigher);
        peakv(selower)  = smdata(selower);
        peakt(sehigher | selower) = 0.001 * tcs(tc);

        % update weighted AUC
        sedata = 1 + log(1 + abs(sedata));
        aucv = aucv + (0.001 * tcs(tc)) .* smdata .* sedata;
        aucw = aucw + sedata;
    end

    % weighted AUC
    aucv = aucv ./ aucw;
    aucv(isinf(aucv) | isnan(aucv)) = 0;

    % set maps data
    mapc = 3 * cc - 2;
    vmpc.Map(mapc).VMPData(mmask) = peakv;
    mapc = mapc + 1;
    vmpc.Map(mapc).UpperThreshold = max(peakt);
    vmpc.Map(mapc).VMPData(mmask) = peakt;
    mapc = mapc + 1;
    vmpc.Map(mapc).UpperThreshold = max(aucv);
    vmpc.Map(mapc).VMPData(mmask) = aucv;
end

% set back
vmp.C = vmpc;

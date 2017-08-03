function vmp = vtc_AverageMovieVMP(xo, opts)
% VTC::AverageMovieVMP  - create VMP with movie maps from averaged VTC
%
% FORMAT:       vmp = vtc.AverageMovieVMP([opts])
%
% Input fields:
%
%       opts        optional struct with fields
%        .conds     conditions (default: all available)
%        .deedge    scale alpha at the edge of the brain by 0.5 (false)
%        .resamptr  resampling TR (to infer maximum value, default: 100)
%        .rgbcolor  4x3 colors for lower/upper positive/negative tails
%        .scaleto1  scale each voxel (activity) to maximum (default: false)
%        .scaletoz  if scaleto1 is true, use fisherr2z (default: 0)
%        .signalmax maximum signal value (for color scale, default: 2)
%        .tthresh   t-threshold used to mask voxels (default: 0.5*log(N))
%        .tthreshu  upper t-threshold (full visibility, default: 4*tthresh)
%        .winlength window length (default: from VTC)
%
% Output fields:
%
%       vmp         VMP object with peak and time-to-peak maps
%
% Note: conds can also contain a difference, in which case, the particles
%       will be weighed and averaged
%
% Using: ddeblank, erode3d, fisherr2z, flexinterpn_method, multimatch,
%        splittocellc.

% Version:  v1.1
% Build:    16021412
% Date:     Feb-14 2016, 12:46 PM EST
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
ddeblank           = ne_methods.ddeblank;
fisherr2z          = ne_methods.fisherr2z;
flexinterpn_method = ne_methods.flexinterpn_method;
multimatch         = ne_methods.multimatch;
splittocellc       = ne_methods.splittocellc;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc')
    error('neuroelf:xff:badArgument', 'Invalid object in call.');
end
bc = xo.C;
if ~isfield(bc.RunTimeVars, 'AvgVTC') || ~islogical(bc.RunTimeVars.AvgVTC) || ...
    numel(bc.RunTimeVars.AvgVTC) ~= 1 || ~bc.RunTimeVars.AvgVTC
    error('neuroelf:xff:badArgument', 'Only valid for Average-VTCs.');
end
cnames = ddeblank(bc.RunTimeVars.ConditionNames);
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'conds') || ~iscell(opts.conds) || isempty(opts.conds)
    opts.conds = cnames;
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
if ~isfield(opts, 'deedge') || ~islogical(opts.deedge) || numel(opts.deedge) ~= 1
    opts.deedge = false;
end
if ~isfield(opts, 'resamptr') || ~isa(opts.resamptr, 'double') || numel(opts.resamptr) ~= 1 || ...
    isinf(opts.resamptr) || isnan(opts.resamptr) || opts.resamptr <= 0
    opts.resamptr = 100;
else
    opts.resamptr = max(50/3, min(1000, opts.resamptr));
end
if ~isfield(opts, 'rgbcolor') || ~isa(opts.rgbcolor, 'double') || ~isequal(size(opts.rgbcolor), [4, 3]) || ...
    any(isinf(opts.rgbcolor(:)) | isnan(opts.rgbcolor(:)) | opts.rgbcolor(:) < 0 | opts.rgbcolor(:) > 255)
    opts.rgbcolor = [224, 16, 0; 255, 224, 16; 0, 32, 160; 0, 160, 192];
else
    if all(opts.rgbcolor(:) <= 1)
        opts.rgbcolor = 255 .* opts.rgbcolor;
    end
    opts.rgbcolor = round(opts.rgbcolor);
end
if ~isfield(opts, 'scaleto1') || ~islogical(opts.scaleto1) || numel(opts.scaleto1) ~= 1
    opts.scaleto1 = false;
end
if ~isfield(opts, 'scaletoz') || ~isa(opts.scaletoz, 'double') || numel(opts.scaletoz) ~= 1 || ...
    isinf(opts.scaletoz) || isnan(opts.scaletoz) || opts.scaletoz < 0
    opts.scaletoz = 0;
end
scaletoz = opts.scaletoz;
if scaletoz > 0
    scaletoz = fisherr2z(min(5, scaletoz), true);
end
if ~isfield(opts, 'signalmax') || ~isa(opts.signalmax) || numel(opts.signalmax) ~= 1 || ...
    isinf(opts.signalmax) || isnan(opts.signalmax) || opts.signalmax <= 0
    opts.signalmax = 2;
end
if ~isfield(opts, 'tthresh') || ~isa(opts.tthresh, 'double') || numel(opts.tthresh) ~= 1 || ...
    isinf(opts.tthresh) || isnan(opts.tthresh) || opts.tthresh < 0
    opts.tthresh = -1;
end
if ~isfield(opts, 'tthreshu') || ~isa(opts.tthreshu, 'double') || numel(opts.tthreshu) ~= 1 || ...
    isinf(opts.tthreshu) || isnan(opts.tthreshu) || opts.tthreshu < opts.tthresh
    if opts.tthresh > 0
        opts.tthreshu = 4 * opts.tthresh;
    else
        opts.tthreshu = 4;
    end
end
trrat = opts.resamptr / bc.TR;
if ~isfield(opts, 'winlength') || ~isa(opts.winlength, 'double') || numel(opts.winlength) ~= 1 || ...
    isinf(opts.winlength) || isnan(opts.winlength) || opts.winlength <= 0
    opts.winlength = [];
end

% resolved conditions array
conds = zeros(numel(opts.conds), numel(cnames));

% first, replace a "-" or ">" in condition name to _MINUS_ and _GT_
for cc = 1:numel(cnames)
    cname = strrep(strrep(cnames{cc}, '>', '_GT_'), '-', '_MINUS_');
    if ~strcmp(cname, cnames{cc})
        for cci = 1:numel(opts.conds)
            opts.conds{cci} = strrep(opts.conds{cci}, cnames{cc}, cname);
        end
        cnames{cc} = cname;
    end
end
for cc = 1:size(conds, 1)

    % get particles
    if any(opts.conds{cc} == '>') || any(opts.conds{cc} == '-')
        cparts = ddeblank(splittocellc(opts.conds{cc}, '->', false, true));
        if numel(cparts) ~= 2
            error('neuroelf:xff:badArgument', 'Invalid contrast specification.');
        end
        cplus = cparts{1};
        cminus = cparts{2};
    else
        cplus = ddeblank(opts.conds{cc});
        cminus = '';
    end

    % for positive and negative, find terms
    cplus = ddeblank(splittocellc(cplus, ' +,;', false, true));
    cminus = ddeblank(splittocellc(cminus, ' +,;', false, true));

    % empty plus set
    if isempty(cplus)
        error('neuroelf:xff:badArgument', 'Invalid contrast specification.');
    end
    cplusi = multimatch(unique(lower(cplus(:))), lower(cnames(:)));
    if any(cplusi) == 0
        error('neuroelf:xff:badArgument', 'Unknown condition in contrast...');
    end
    conds(cc, cplusi) = 1;

    % negative set
    if ~isempty(cminus)
        cminusi = multimatch(unique(lower(cminus(:))), lower(cnames(:)));
        if any(cminusi) == 0 || ~isempty(intersect(cplusi(:), cminusi(:)))
            error('neuroelf:xff:badArgument', 'Unknown or invalid condition in contrast...');
        end
        conds(cc, cminusi) = -1;
    end
end

% data
vptc = bc.RunTimeVars.NrOfVolumesPerTC;
tcpc = bc.RunTimeVars.NrOfTCsPerCondition;
ncnd = bc.RunTimeVars.NrOfConditions;
tce = bc.RunTimeVars.AvgWindowTo;
if ~isempty(opts.winlength)
    tce = min(tce, opts.winlength);
end
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
vmpc.Map(1).RunTimeVars = struct;
vmpc.Map = [vmpc.Map(1), vmpc.Map(1)];
vmpc.Map(1).Type = 15;
vmpc.Map(1).RGBLowerThreshPos = opts.rgbcolor(1, :);
vmpc.Map(1).RGBUpperThreshPos = opts.rgbcolor(2, :);
vmpc.Map(1).RGBLowerThreshNeg = opts.rgbcolor(3, :);
vmpc.Map(1).RGBUpperThreshNeg = opts.rgbcolor(4, :);
vmpc.Map(1).UseRGBColor = 1;
vmpc.Map(1).LowerThreshold = 0.025;
if opts.scaleto1
    if opts.scaletoz > 0
        vmpc.Map(1).UpperThreshold = opts.scaletoz;
    else
        vmpc.Map(1).UpperThreshold = 1;
    end
else
    vmpc.Map(1).UpperThreshold = opts.signalmax;
end
vmpc.Map(1).RunTimeVars.AlphaMap = 2;
vmpc.Map(2).Type = 31;
vmpc.Map(2).LowerThreshold = opts.tthresh;
vmpc.Map(2).UpperThreshold = opts.tthreshu;
vmpc.Map = repmat(vmpc.Map, 1, size(conds, 1) * numel(tcs));
vmpc.NrOfMaps = numel(vmpc.Map);

% deedging
if opts.deedge
    vtcmask = squeeze(mean(abs(vtcd), 1) > 0.25 * mean(abs(vtcd(:))));
    vtcmask = 0.5 .* double(vtcmask) + 0.5 .* double(ne_methods.erode3d(vtcmask));
    vtcmask = vtcmask(:)';
end

% squeeze data
vtcd = reshape(vtcd, size(vtcd, 1), prod(vtcsz));

% mask
mmask = any(vtcd ~= 0, 1);
nmask = sum(mmask(:));

% iterate over conditions
for cc = 1:size(conds, 1)

    % get data indices
    ccii = conds(cc, :);

    % only one condition (which then must be positive!)
    cci = find(ccii ~= 0);
    if numel(cci) == 1

        % data index
        di = (1:vptc) + ((cci - 1) * vptc * tcpc);

        % get actual data
        mdata = vtcd(di, mmask);

        % errors (and weights)
        edata = vtcd(di + vptc, mmask);
        ndata = bc.RunTimeVars.NrOfConditionOnsets(cci);
        if tcpc > 2
            wdata = ndata .* vtcd(di + 2 * vptc, mmask);
        else
            wdata = ndata .* repmat(bc.RunTimeVars.TCOnsetWeights(di), 1, nmask);
        end

        % compute error stat for alpha (from SD -> SE -> t-stat)
        edata = edata ./ sqrt(wdata);
        edata(isinf(edata) | isnan(edata)) = Inf;
        edata = mdata ./ edata;
        edata(isinf(edata) | isnan(edata)) = 0;

        % threshold
        if opts.tthresh < 0
            tthresh = 0.5 * log(ndata);
            tthreshu = opts.tthreshu * tthresh;
        end

    % multiple (collapse and contrast)
    else

        % first take care of positive conditions
        cci = find(ccii > 0);
        pmdata = zeros(vptc, nmask);
        pedata = zeros(vptc, nmask);
        pwdata = zeros(vptc, nmask);
        pndata = 0;
        for ccic = 1:numel(cci)
            di = (1:vptc) + ((cci(ccic) - 1) * vptc * tcpc);
            mdata = vtcd(di, mmask);
            edata = vtcd(di + vptc, mmask);
            ndata = bc.RunTimeVars.NrOfConditionOnsets(cci(ccic));
            if tcpc > 2
                wdata = ndata .* vtcd(di + 2 * vptc, mmask);
            else
                wdata = ndata .* repmat(bc.RunTimeVars.TCOnsetWeights(di), 1, nmask);
            end
            pmdata = pmdata + wdata .* mdata;
            pedata = pedata + wdata .* edata .* edata;
            pwdata = pwdata + wdata;
            pndata = pndata + ndata;
        end

        % average
        pmdata = pmdata ./ pwdata;
        pedata = sqrt(pedata ./ pwdata);

        % find negative contrast particles
        cci = find(ccii < 0);
        if ~isempty(cci)
            nmdata = zeros(vptc, nmask);
            nedata = zeros(vptc, nmask);
            nwdata = zeros(vptc, nmask);
            nndata = 0;
            for ccic = 1:numel(cci)
                di = (1:vptc) + ((cci(ccic) - 1) * vptc * tcpc);
                mdata = vtcd(di, mmask);
                edata = vtcd(di + vptc, mmask);
                ndata = bc.RunTimeVars.NrOfConditionOnsets(cci(ccic));
                if tcpc > 2
                    wdata = ndata .* vtcd(di + 2 * vptc, mmask);
                else
                    wdata = ndata .* repmat(bc.RunTimeVars.TCOnsetWeights(di), 1, nmask);
                end
                nmdata = nmdata + wdata .* mdata;
                nedata = nedata + wdata .* edata .* edata;
                nwdata = nwdata + wdata;
                nndata = nndata + ndata;
            end

            % combine negative and positive data
            nmdata = nmdata ./ nwdata;
            nedata = sqrt(nedata ./ nwdata);
            pmdata = pmdata - nmdata;
            pedata = pedata + nedata;
            pwdata = pwdata + nwdata;
            pndata = pndata + nndata;
        end

        % copy over
        mdata = pmdata;

        % error estimate
        pedata = pedata;
    end

    % error is absolute (for alpha maps!)
    edata = abs(edata);

    % scaling
    if opts.scaleto1
        scdata = zeros(1, size(mdata, 2));
    end

    % iterate over time
    for tc = 1:numel(tcs)

        % sample value
        scoord = [Inf, Inf; 1 + (tc - 1) * trrat, 1; 2 * vptc, 1; size(mdata)];
        smdata = flexinterpn_method(mdata, scoord, 'cubic');
        sedata = flexinterpn_method(edata, scoord, 'cubic');

        % keep tally
        if opts.scaleto1
            scdata = max(scdata, abs(smdata));
        end

        % set maps data
        mapc = 1 + 2 * (numel(tcs) * (cc - 1) + (tc - 1));
        vmpc.Map(mapc).Name = sprintf('%s signal average at %6.3fs', opts.conds{cc}, 0.001 * tcs(tc));
        vmpc.Map(mapc).VMPData(mmask) = single(smdata);
        vmpc.Map(mapc).RunTimeVars.AlphaMap = mapc + 1;
        mapc = mapc + 1;
        vmpc.Map(mapc).Name = sprintf('%s signal alpha at %6.3fs', opts.conds{cc}, 0.001 * tcs(tc));
        if opts.tthresh < 0
            vmpc.Map(mapc).LowerThreshold = tthresh;
            vmpc.Map(mapc).UpperThreshold = tthreshu;
        end
        if opts.deedge
            vmpc.Map(mapc).VMPData(mmask) = single(sedata .* vtcmask(mmask));
        else
            vmpc.Map(mapc).VMPData(mmask) = single(sedata);
        end
    end

    % rescale?
    if opts.scaleto1
        scdata(scdata == 0) = Inf;
        if scaletoz > 0
            scdata = scaletoz ./ scdata;
        else
            scdata = 1 ./ scdata;
        end
        for tc = 1:numel(tcs)
            mapc = 1 + 2 * (numel(tcs) * (cc - 1) + (tc - 1));
            if scaletoz > 0
                vmpc.Map(mapc).VMPData(mmask) = single(fisherr2z( ...
                    scdata .* double(vmpc.Map(mapc).VMPData(mmask))));
            else
                vmpc.Map(mapc).VMPData(mmask) = single(scdata .* vmpc.Map(mapc).VMPData(mmask));
            end
        end
    end
end

% set back
vmp.C = vmpc;

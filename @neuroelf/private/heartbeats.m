function [btpos, bs, bf, bv, cpos, cgood, cdata, hrv] = heartbeats(sig, opts)
% heartbeats  - detect heart beats and frequency in physio data
%
% FORMAT:       [bp, bs, bf, bv, cp, wgd, wd, hrv] = heartbeats(sig [, opts])
%
% Input fields:
%
%       sig         Sx1 numeric signal
%       opts        optional settings
%        .badt      t-threshold for detecing irregularity (default: 25)
%        .bppre     pre-defined positions (no detection, only inspection)
%        .calc      preprocessing calculcation, one of
%                   {'none'} - don't to anything (default)
%                   {'absdiff'} - use abs(diff(data))
%                       if .detlength is not given, will be set to 0.02
%                       if .skewdt is not given, will be set to 0.1
%                   {'diffsq'} - square the diff of the data
%                       if .detlength is not given, will be set to 0.01
%                       if .skewdt is not given, will be set to 0.05
%                   {'fourthz'} - fourth power of the z-transformed data
%                       if .detlength is not given, will be set to 0.02
%                       if .skewdt is not given, will be set to 0.04
%                   {'squarez'} - square the z-transformed data
%                       if .detlength is not given, will be set to 0.03
%                       if .skewdt is not given, will be set to 0.1
%                   {'thirdz'} - third power of the abs z-transformed data
%                       if .detlength is not given, will be set to 0.02
%                       if .skewdt is not given, will be set to 0.06
%        .cleanup   interactive cleanup (default: false)
%        .correct   correct for unlikely beats (default: true)
%        .detlength detection length threshold in seconds (default: 0.05)
%        .freq      data frequency in Hz (default: 1000)
%        .hrvrfreq  resampling frequency for computehrv (default: 10)
%        .limrange  1x2 double, limitrange on input data (default: [])
%        .pflength  pre-filter length in seconds (default: 0.025)
%        .pfreps    pre-filter repetitions (default: 2)
%        .plot      plot mean +/- std estimate of signal (default: false)
%        .plotfreq  samples per second to plot (default: 50)
%        .plotwin   plot window size in seconds (default: 6)
%        .poshlp    position help (new beats), either {'max'}, 'min', 'off'
%        .resfreq   resample data prior to detection (default: [])
%        .segsize   segmentation size in seconds (default: 5)
%        .segstep   stepping (window shift) in seconds (default: 1)
%        .skewdt    skewness detection threshold multiplier (default: 0.5)
%        .title     string added to the title window (default: '')
%        .winsor    winsorizing threshold in std's (default: 3)
%
% Output fields:
%
%       bp          beat positions
%       bs          beat positions (in seconds)
%       bf          frequency estimate for each beat
%       bv          values at peak
%       cp          estimate of centers of windows (one value less!)
%       wgd         guess whether window is good or not
%       wd          windowed data (in 100Hz resolution, interpolated)
%       hrv         output of computehrv(bp)
%
% Note: this function is still preliminary, other options passed on to
%       computehrv (if 8th output is requested), with .hrvrfreq being
%       set to .resfreq

% Version:  v1.1
% Build:    16062212
% Date:     Jun-22 2016, 12:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% argument check
if nargin < 1 || ...
   ~isnumeric(sig) || ...
    numel(sig) ~= max(size(sig))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid electrode data provided.' ...
    );
end
sig = sig(:);
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'badt') || ...
   ~isa(opts.badt, 'double') || ...
    numel(opts.badt) ~= 1 || ...
    isinf(opts.badt) || ...
    isnan(opts.badt) || ...
    opts.badt <= 0
    opts.badt = 25;
else
    opts.badt = min(60, max(5, opts.badt));
end
if ~isfield(opts, 'bppre') || ...
    isempty(opts.bppre) || ...
   ~isnumeric(opts.bppre)
    opts.bppre = [];
end
opts.bppre = double(opts.bppre(:));
opts.bppre(isinf(opts.bppre) | isnan(opts.bppre)) = [];
if ~isfield(opts, 'calc') || ...
   ~ischar(opts.calc) || ...
    isempty(opts.calc) || ...
   ~any(lower(opts.calc(1)) == 'adfnst')
    opts.calc = 'n';
else
    opts.calc = lower(opts.calc(1));
end
switch (opts.calc)
    case {'a'}
        if ~isfield(opts, 'detlength')
            opts.detlength = 0.02;
        end
        if ~isfield(opts, 'skewdt')
            opts.skewdt = 0.1;
        end
    case {'d'}
        if ~isfield(opts, 'detlength')
            opts.detlength = 0.01;
        end
        if ~isfield(opts, 'skewdt')
            opts.skewdt = 0.05;
        end
    case {'f'}
        if ~isfield(opts, 'detlength')
            opts.detlength = 0.02;
        end
        if ~isfield(opts, 'skewdt')
            opts.skewdt = 0.04;
        end
    case {'s'}
        if ~isfield(opts, 'detlength')
            opts.detlength = 0.03;
        end
        if ~isfield(opts, 'skewdt')
            opts.skewdt = 0.2;
        end
    case {'t'}
        if ~isfield(opts, 'detlength')
            opts.detlength = 0.02;
        end
        if ~isfield(opts, 'skewdt')
            opts.skewdt = 0.08;
        end
end
if ~isfield(opts, 'cleanup') || ...
   ~islogical(opts.cleanup) || ...
    numel(opts.cleanup) ~= 1
    opts.cleanup = false;
end
if ~isfield(opts, 'correct') || ...
   ~islogical(opts.correct) || ...
    numel(opts.correct) ~= 1
    opts.correct = true;
end
if ~isfield(opts, 'detlength') || ...
   ~isa(opts.detlength, 'double') || ...
    numel(opts.detlength) ~= 1 || ...
    isinf(opts.detlength) || ...
    isnan(opts.detlength) || ...
    opts.detlength <= 0.001 || ...
    opts.detlength > 0.5
    opts.detlength = 0.05;
end
if ~isfield(opts, 'freq') || ...
   ~isa(opts.freq, 'double') || ...
    numel(opts.freq) ~= 1 || ...
    isinf(opts.freq) || ...
    isnan(opts.freq) || ...
    opts.freq <= 0
    opts.freq = 1000;
end
if ~isfield(opts, 'hrvrfreq') || ...
   ~isa(opts.hrvrfreq, 'double') || ...
    numel(opts.hrvrfreq) ~= 1 || ...
    isinf(opts.hrvrfreq) || ...
    isnan(opts.hrvrfreq) || ...
    opts.hrvrfreq < 1
    opts.hrvrfreq = 10;
end
if ~isfield(opts, 'limrange') || ...
   ~isa(opts.limrange, 'double') || ...
    numel(opts.limrange) ~= 2 || ...
    any(isnan(opts.limrange)) || ...
    opts.limrange(2) <= opts.limrange(1)
    opts.limrange = [];
end
if ~isfield(opts, 'pflength') || ...
   ~isa(opts.pflength, 'double') || ...
    numel(opts.pflength) ~= 1 || ...
    isinf(opts.pflength) || ...
    isnan(opts.pflength) || ...
    opts.pflength < 0.001 || ...
    opts.pflength > 0.1
    opts.pflength = 0.025;
end
if ~isfield(opts, 'pfreps') || ...
   ~isa(opts.pfreps, 'double') || ...
    numel(opts.pfreps) ~= 1 || ...
    isinf(opts.pfreps) || ...
    isnan(opts.pfreps) || ...
    opts.pfreps < 0 || ...
    opts.pfreps > 12
    opts.pfreps = 2;
else
    opts.pfreps = round(opts.pfreps);
end
if ~isfield(opts, 'plot') || ...
   ~islogical(opts.plot) || ...
    numel(opts.plot) ~= 1
    opts.plot = false;
end
if ~isfield(opts, 'plotfreq') || ...
   ~isa(opts.plotfreq, 'double') || ...
    numel(opts.plotfreq) ~= 1 || ...
    isinf(opts.plotfreq) || ...
    isnan(opts.plotfreq) || ...
    opts.plotfreq <= 0
    opts.plotfreq = 50;
end
opts.plotfreq = min(1000, opts.plotfreq);
if ~isfield(opts, 'plotwin') || ...
   ~isa(opts.plotwin, 'double') || ...
    numel(opts.plotwin) ~= 1 || ...
    isinf(opts.plotwin) || ...
    isnan(opts.plotwin) || ...
    opts.plotwin <= 0
    opts.plotwin = 6;
else
    opts.plotwin = min(20, max(2, opts.plotwin));
end
if ~isfield(opts, 'poshlp') || ...
   ~ischar(opts.poshlp') || ...
    isempty(opts.poshlp) || ...
   ~any(strcmpi(opts.poshlp(:)', {'max', 'min', 'off'}))
    opts.poshlp = 'max';
else
    opts.poshlp = lower(opts.poshlp(:)');
end
if ~isfield(opts, 'resfreq') || ...
   ~isa(opts.resfreq, 'double') || ...
    numel(opts.resfreq) ~= 1 || ...
    isinf(opts.resfreq) || ...
    isnan(opts.resfreq) || ...
    opts.resfreq <= 0 || ...
    opts.resfreq == opts.freq
    opts.resfreq = [];
end
if ~isfield(opts, 'segsize') || ...
   ~isa(opts.segsize, 'double') || ...
    numel(opts.segsize) ~= 1 || ...
    isinf(opts.segsize) || ...
    isnan(opts.segsize) || ...
    opts.segsize < 2 || ...
    opts.segsize > 20
    opts.segsize = 5;
end
if ~isfield(opts, 'segstep') || ...
   ~isa(opts.segstep, 'double') || ...
    numel(opts.segstep) ~= 1 || ...
    isinf(opts.segstep) || ...
    isnan(opts.segstep) || ...
    opts.segstep < 0.1 || ...
    opts.segstep > (0.5 * opts.segsize)
    opts.segstep = 1;
end
if ~isfield(opts, 'skewdt') || ...
   ~isa(opts.skewdt, 'double') || ...
    numel(opts.skewdt) ~= 1 || ...
    isinf(opts.skewdt) || ...
    isnan(opts.skewdt) || ...
    opts.skewdt <= 0 || ...
    opts.skewdt > 2
    opts.skewdt = 0.5;
end
if ~isfield(opts, 'title') || ...
   ~ischar(opts.title) || ...
    isempty(opts.title)
    opts.title = '';
else
    if numel(opts.title) < 64
        opts.title = [opts.title(:)' ' - '];
    else
        opts.title = opts.title(:)';
        opts.title = [opts.title(1:29) ' ... ' opts.title(end-28:end) ' - '];
    end
end
if ~isfield(opts, 'winsor') || ...
   ~isa(opts.winsor, 'double') || ...
    numel(opts.winsor) ~= 1 || ...
    isinf(opts.winsor) || ...
    isnan(opts.winsor) || ...
    opts.winsor < 0
    opts.winsor = 3;
else
    opts.winsor = min(8, opts.winsor);
end

% preprocessing -> limit range
if ~isempty(opts.limrange)

    % any inf
    if any(isinf(opts.limrange))

        % use one-sided version (max or min)
        if ~isinf(opts.limrange(1))
            sig = max(opts.limrange(1), sig);
        elseif ~isinf(opts.limrange(2))
            sig = min(opts.limrange(2), sig);
        end

    % two numbers
    else
        sig = limitrangec(sig, ...
            opts.limrange(1), opts.limrange(2), 0.5 * sum(opts.limrange));
    end
end

% preprocessing -> resampling
if ~isempty(opts.resfreq)
    sig = resampleaa(sig, opts.freq / opts.resfreq);
    opts.freq = opts.resfreq;
    opts.resfreq = [];
end
if isfield(opts, 'hrvrfreq')
    opts.resfreq = opts.hrvrfreq;
end

% set nans to zero
sig(isinf(sig) | isnan(sig)) = 0;

% -> calculations
switch (opts.calc)
    case {'a'}
        sig = abs(diff(sig));
    case {'d'}
        sig = diff(sig) .^ 2;
    case {'f'}
        sig = ztrans(sig) .^ 4;
    case {'s'}
        sig = ztrans(sig) .^ 2;
    case {'t'}
        sig = abs(ztrans(sig)) .^ 3;
end

% dimensions
numsig = numel(sig);

% btpos not defined by opts.bppre?
if isempty(opts.bppre)

    % filtering
    filtsize = smoothkern(opts.pflength * opts.freq);
    segmsize = round(opts.segsize * opts.freq);
    segmstep = round(opts.freq);
    segmmins = opts.detlength * opts.freq;
    filtered = true;

    % create beats data
    beats = uint8(0);
    beats(numsig, 1) = 0;

    % windowed z-trans
    fsig = zeros(size(sig));
    for so = 1:segmstep:(numsig + 2 * segmstep - 2 * segmsize)
        sof = min(numsig, so + 2 * segmsize - 1);
        beats(so:sof) = beats(so:sof) + 1;
        zsig = ztrans(sig(so:sof));
        zsig(isnan(zsig)) = 0;
        fsig(so:sof) = fsig(so:sof) + zsig;
    end
    sig = fsig ./ double(beats);

    % fix signal to contain only values within +/- 3std range
    if opts.winsor > 0
        sig(sig > opts.winsor) = opts.winsor;
        sig(sig < -opts.winsor) = -opts.winsor;
    end

    % filter data
    fsig = sig;
    if numel(filtsize) > 1 && ...
        opts.pfreps > 0
        for fc = 1:opts.pfreps
            fsig = flexinterpn(fsig, [Inf; 1; 1; numel(fsig)], filtsize, 1, 0);
        end
    end

    % depending on signal
    ssig = opts.skewdt * skew(sig);
    if skew(sig) < 0
        fsig = -fsig;
        ssig = -ssig;
    end

    % for each segment, find beats
    beats(:) = 0;
    for so = 1:segmstep:(numsig - segmsize)

        % get segment
        sseg = sig(so:(so + segmsize - 1));
        fseg = fsig(so:(so + segmsize - 1));

        % find sub-segments > mean
        ssm = diff(double(fseg > ssig));
        ssp = find(ssm > 0);
        ssn = find(ssm < 0);
        if isempty(ssp) || ...
            isempty(ssn)
            continue;
        end
        if ssp(1) > ssn(1)
            ssp = [1; ssp(:)];
        end
        numseg = min(numel(ssp), numel(ssn));

        % for each sub-segment
        for sc = 1:numseg

            % if segment too small, continue
            if (ssn(sc) - ssp(sc)) < segmmins
                continue;
            end

            % determine maximum
            [mv, mp] = max(sseg(ssp(sc):ssn(sc)));

            % if not first value in segment
            if mp > 1

                % at maximum ?
                if mv == opts.winsor

                    % replace with mean position
                    mp = round(mean(find(sseg(ssp(sc):ssn(sc)) == opts.winsor)));
                end

                % add as detected peak
                mp = mp + so + ssp(sc) - 2;
                beats(mp) = beats(mp) + 1;
            end
        end
    end

    % identify beats if found in at least 2 segments
    btpos = find(beats >= uint8(2));

    % (almost) no beats? we can't do much but simply use a poor man's guess
    if numel(beats) < (0.1 * numsig / opts.freq)

        % consider rough estimates for [75, 80, 85, 90, 95, and 100] BPM
        acstep = unique(round(60 .* opts.freq ./ [75, 80, 85, 90, 95, 100]));
        acsdif = mean(diff(acstep));
        acr = zeros(1, numel(acstep));
        for acc = 1:numel(acr)
            [acv, acr(acc)] = cov_nd(fsig', fsig', acstep(acc));
        end

        % and find best match
        acstep = acstep(maxpos(acr));

        % then look if shifting a bit still improves the result
        while ~isempty(acsdif)
            acsdif = round(0.5 * acsdif);
            acr = [0, max(acr), 0];
            [acv, acr(1)] = cov_nd(fsig', fsig', acstep - acsdif);
            [acv, acr(3)] = cov_nd(fsig', fsig', acstep + acsdif);
            if acr(1) > acr(2) && ...
                acr(1) > acr(3)
                acstep = acstep - acsdif;
            elseif acr(3) > acr(2)
                acstep = acstep + acsdif;
            else
                acsdif = [];
            end
        end

        % simply find max overall value in each period
        btstep = 1:acstep:numsig;
        btpos = zeros(numel(btstep) - 1, 1);
        for btspos = 1:numel(btpos)
            btpos(btspos) = btstep(btspos) + ...
                maxpos(sig(btstep(btspos):btstep(btspos+1)-1)) - 1;
        end
    end

% btpos given
else
    btpos = opts.bppre;
    fsig = sig;
    filtered = false;
end

% didn't work?
if numel(btpos) < 3
    warning( ...
        'neuroelf:DetectionFailed', ...
        'Detection of beats failed with given settings.' ...
    );
    bs = [];
    bf = [];
    bv = [];
    cpos = [];
    cgood = [];
    cdata = [];
    hrv = struct;
    return;
end

% remove unlikely beats
if opts.correct
    bdiff = diff(btpos);
    mbdiff = robustmean(bdiff);
    btpos(1 + find(bdiff < (0.333 * mbdiff))) = [];

    % and insert where gap too wide
    bdiff = diff(btpos);
    mbdiff = robustmean(bdiff);
    cbdiff = ceil(mbdiff);
    gaps = find(bdiff > 3 * mbdiff);
    for gapc = 1:numel(gaps)
        gapf = round(btpos(gaps(gapc))+0.5*mbdiff:mbdiff:btpos(gaps(gapc)+1));
        for gapfc = 1:numel(gapf)
            gapf(gapfc) = ...
                gapf(gapfc) + maxpos(sig(gapf(gapfc):gapf(gapfc)+cbdiff)) - 1;
        end
        btpos = [btpos; gapf(:)];
    end
    if ~isempty(gaps)
        btpos = sort(btpos(:));
    end
end

% compute range values
cpos = floor(0.5 .* (btpos(1:end-1) + btpos(2:end)));
cdif = diff(cpos);
cdif(cdif < 1) = 1;

% create list of periods and centers to sample from
[cdfm, w] = robustmean(cdif);
cdfs = std(cdif(w >= 0.5));
cdfl = cdfm - 2.5 * cdfs;
cdfu = cdfm + 2.5 * cdfs;
cgood = ((cdif >= cdfl) & (cdif <= cdfu));

% cleanup requested and necessary
if opts.cleanup

    % build "cheap" mean of good signals
    cdata = zeros(101, numel(cdif));
    cstp = 0.01 .* cdif;
    for cc = 1:size(cdata, 2)
        cdata(:, cc) = sig(round(cpos(cc):cstp(cc):cpos(cc+1) + 0.5 * cstp(cc)));
    end
    mdata = ztrans(mean(cdata(:, cgood), 2));

    % replace estimate with how well the signal fits the mean
    [badcv, badr] = cov_nd(cdata', mdata(:, ones(1, numel(cdif)))');
    badr(isnan(badr) | isinf(badr)) = 0;
    badt = correltstat(badr, 100);

    % force first and last!
    badt([1, end]) = 0;
    bad = find(badt < opts.badt);

    % create figure and axes
    hFig = neuroelf_file('f', 'heartbeats_cleanup');
    hTag = hFig.TagStruct;
    oax = hTag.AX_heartbeats_overview;
    oaxm = oax.MLHandle;

    % compute overview signal
    osig = zeros(oax.Position(3) - oax.Position(1), 1);
    osigf = numsig / numel(osig);
    osigx = (1 / opts.freq) .* ((0.5:numel(osig)) .* osigf)';
    osigf = round(1 + (0:numel(osig)) .* osigf);
    osigf(end) = min(numsig, osigf(end));
    for osc = 1:numel(osig)
        osig(osc) = std(diff(sig(osigf(osc):osigf(osc+1))));
    end
    osig = log(osig);
    osig = max(mean(osig) - std(osig), osig);
    osig = osig - min(osig);
    hold(oaxm, 'on');
    oaxyl = ceil(max(osig));
    badl = zeros(numel(bad), 1);
    for badc = 1:numel(bad)
        bppos = btpos(bad(badc)) / opts.freq;
        badl(badc) = plot(oaxm, [bppos; bppos], [0; oaxyl]);
    end
    set(plot(oaxm, osigx, osig), 'Color', [1, 0, 0]);
    set(oaxm, 'XLim', [0, numsig / opts.freq], 'YTick', []);

    % set colors
    hTag.AX_heartbeats_plot.Color = [1, 1, 1];
    set(oaxm, 'Color', [1, 1, 1]);

    % setup persistent storage in UserData
    hFig.UserData = struct( ...
        'ax',     hTag.AX_heartbeats_plot.MLHandle, ...
        'axpos',  hTag.AX_heartbeats_plot.Position, ...
        'bad',    bad, ...
        'break',  false, ...
        'btnum',  1, ...
        'btpos',  btpos, ...
        'cpos',   cpos, ...
        'freq',   opts.freq, ...
        'fsig',   fsig, ...
        'hFig',   hFig, ...
        'hTag',   hTag, ...
        'numsig', numsig, ...
        'oax',    oaxm, ...
        'oaxl',   [], ...
        'oaxpos', oax.Position, ...
        'oaxyl',  oaxyl, ...
        'osig',   osig, ...
        'osigx',  osigx, ...
        'plwin',  ceil(0.5 * opts.freq * opts.plotwin), ...
        'poshlp', opts.poshlp, ...
        'sig',    sig, ...
        'title',  opts.title, ...
        'xfigure_UDTagStruct', hTag);
    hFig.CloseRequestFcn = @hb_cancel;
    hFig.WindowButtonDownFcn = @hb_btdown;
    hFig.HandleVisibility = 'callback';
    hTag.LB_heartbeats_beats.Callback = @hb_select;
    hTag.BT_heartbeats_remove.Callback = @hb_remove;
    hTag.BT_heartbeats_removea.Callback = {@hb_remove, 1};
    hTag.BT_heartbeats_accept.Callback = {@hb_accept, 1};
    hTag.BT_heartbeats_acceptb.Callback = {@hb_accept, -1};
    if filtered
        hTag.CB_heartbeats_fsig.Callback = @hb_fsig;
    else
        hTag.CB_heartbeats_fsig.Enable = 'off';
    end

    % reposition overview if available
    delete(findobj('type', 'figure', 'Tag', 'MFIG_neuroelf_heartbeats_shape'));
    ovf = findobj('type', 'figure', 'Tag', 'MFIG_neuroelf_heartbeats_overview');
    if ~isempty(ovf)
        pos = hFig.Position;
        pos(2) = pos(2) - round(0.5 * pos(4));
        hFig.Position = pos;
        pos(2) = pos(2) + pos(4) + 20;
        pos(4) = pos(4) - 40;
        set(ovf, 'Position', pos);
        if numel(ovf) > 1
            delete(ovf(2:end));
        end
    end

    % for each bad spot
    while hFig.UserData.btnum <= numel(bad)

        % show signal
        hb_showbeats(hFig);

        % wait for figure to close
        uiwait(hFig.MLHandle);

        % break?
        if hFig.UserData.break
            break;
        end
    end

    % updated btpos
    btpos = hFig.UserData.btpos;

    % delete figure
    hFig.Delete;
    drawnow;

% cleanup recommended and obviously not scripted?
elseif sum(cgood) ~= numel(cgood) && ...
    nargout > 1

    % recommend then
    warning( ...
        'neuroelf:Suggestion', ...
        'Not all detected beats are in good shape. Cleanup suggested.' ...
    );
end

% remove unlikely beats
bdiff = diff(btpos);
btpos(1 + find(bdiff < (0.2 * opts.freq))) = [];

% compute range values
cpos = floor(0.5 .* (btpos(1:end-1) + btpos(2:end)));
cdif = diff(cpos);

% didn't work?
if numel(cdif) < 3
    warning( ...
        'neuroelf:DetectionFailed', ...
        'Detection of beats failed with given settings.' ...
    );
    return;

% otherwise
else

    % create list of periods and centers to sample from
    [cdfm, w] = robustmean(cdif);
    cdfs = std(cdif(w >= 0.5));
    cdfl = cdfm - 2.5 * cdfs;
    cdfu = cdfm + 2.5 * cdfs;
    cgood = ((cdif >= cdfl) & (cdif <= cdfu));
end

% compute heart-rate
bf = (60 .* opts.freq) ./ cdif;
bf = [bf(1); bf(:); bf(end)];

% get values
bv = sig(btpos);

% compute second indices
bs = btpos ./ opts.freq;

% create 2-dim array and resample each period between centers
cdata = zeros(101, numel(cdif));
cstp = 0.01 .* cdif;
rsfac = 0.01 * cdfm;
if rsfac >= 1
    kern = smoothkern(rsfac);
    kern(kern < ((1 / rsfac) .^ 3)) = [];
    if isempty(kern)
        kern = [0;0;1;0;0];
    end
else
    kern = [0;0;1;0;0];
end
kern = flexinterpn_method(kern, [Inf; 1; 1/4096; numel(kern)], 'cubic');
for cc = 1:size(cdata, 2)
    cdata(:, cc) = flexinterpn(sig, ...
        [Inf; cpos(cc); cstp(cc); cpos(cc+1) + 0.5 * cstp(cc)], kern, 4096);
end

% computehrv?
if nargout > 7
    hrv = computehrv(btpos, opts);
end

% plot ?
if opts.plot

    % delete old figures
    delete(findobj('type', 'figure', 'Tag', 'MFIG_neuroelf_heartbeats_overview'));
    delete(findobj('type', 'figure', 'Tag', 'MFIG_neuroelf_heartbeats_shape'));

    % compute mean
    mdata = mean(cdata(:, cgood), 2);

    % but generate std from all data!
    sdata = std(cdata, [], 2);

    % create second figure for mean +/- std plot
    f = figure;
    set(f, 'NumberTitle', 'off', ...
        'Name', ['NeuroElf - ' opts.title 'heartbeats: shape mean +/- std'], ...
        'Tag', 'MFIG_neuroelf_heartbeats_shape');
    ax = axes('Parent', f);

    % then plot
    tcplot(ax, (-0.5:0.01:0.5)', mdata, sdata, sdata, struct( ...
        'color',   [.5, .5, 0.1], ...
        'lwidth',  2, ...
        'scolor',  [0.4, 0.4, 0.9], ...
        'spalpha', 0.4));

    % create figure for simplified plot
    f = figure;
    set(f, 'NumberTitle', 'off', ...
        'Name', ['NeuroElf - heartbeats - ' opts.title 'timecourse overview'], ...
        'Tag', 'MFIG_neuroelf_heartbeats_overview');
    ax = axes('Parent', f);

    % plot ztrans'ed long, resampled signal
    rsfac = (opts.freq/opts.plotfreq);
    sig = resampleaa(ztrans(sig), rsfac, 1, 0);
    spos = (rsfac / opts.freq) .* (0:numel(sig)-1)';
    plot(ax, spos, sig);

    % add scattered heartrate
    hold(ax, 'on');
    scatter(bs, 2.5 + (1 / 60) .* bf);
end



% sub function (taken from convones.m)
function co = hb_convones(f, cn, uw)
fn = numel(f);
lf = fn + cn - 1;
co = zeros(lf, 1);
cimed = [f(:); zeros(cn - 1, 1)];

% weight function
if uw
    wfunc = hb_convones(ones(fn, 1), cn, false);
end

% continue as long as cnum is not 0 !
cpos = 0;
mfac = fix(cn / 2);
cfac = 1;
while cn > 0
    if mod(cn, 2)
        co = co + [zeros(cpos, 1); cimed(1:lf-cpos)];
        cpos = cpos + cfac;
    end
    if cfac > mfac
        break;
    end
    cimed = [cimed(1:lf-cfac) ; zeros(cfac, 1)] + ...
            [zeros(cfac, 1) ; cimed(1:lf-cfac)];
    cfac = cfac * 2;
    cn = fix(cn / 2);
end

% weight now
if uw
    co = co ./ wfunc;
end


% cancel
function hb_cancel(varargin)

% get figure
hFigMLH = findobj('Tag', 'Wnd_NeuroElf_heartbeats');

% update UserData
ud = get(hFigMLH, 'UserData');
ud.break = true;
set(hFigMLH, 'UserData', ud);

% make uiwait -> resume
uiresume;


% add beat
function hb_btdown(varargin)

% get figure
hFigMLH = findobj('Tag', 'Wnd_NeuroElf_heartbeats');

% check figure
if nargin < 1 || ...
    isempty(varargin{1}) || ...
    varargin{1}(1) ~= hFigMLH
    return;
end

% more settings
ud = get(hFigMLH, 'UserData');
hFig = ud.hFig;

% get current point and compare to main axes
cp = hFig.CurrentPoint;
ocp = cp - ud.oaxpos(1:2);
cp = cp - ud.axpos(1:2);
if all(cp >= 0) && ...
    all(cp <= ud.axpos(3:4))

    % compute position for new beat
    xl = get(ud.ax, 'XLim');
    cp = cp ./ ud.axpos(3:4);
    xd = xl(2) - xl(1);
    xp = ud.freq * (xl(1) + xd * cp(1));

    % check closer vicinity of position
    if ~strcmp(ud.poshlp, 'off')
        xpcfrom = round(max(1, xp - 0.1 * ud.freq));
        xpcto = round(min(ud.numsig, xp + 0.1 * ud.freq));
        if strcmp(ud.poshlp, 'max')
            xp = xpcfrom + maxpos(ud.sig(xpcfrom:xpcto)) - 1;
        else
            xp = xpcfrom + minpos(ud.sig(xpcfrom:xpcto)) - 1;
        end
    end

    % add to btpos
    hFig.UserData.btpos = unique([ud.btpos(:); xp]);

    % re-draw
    hb_showbeats(hFig);

% compare to overview axes
elseif all(ocp >= 0) && ...
    all(ocp <= ud.oaxpos(3:4))

    % get clicked area
    xp = ocp(1) * ud.numsig / ud.oaxpos(3);

    % find nearest window and choose that as the next "bad" trial
    xs = ud.cpos(ud.bad);
    hFig.UserData.btnum = minpos(abs(xs - xp));

    % re-draw
    hb_showbeats(hFig);
end


% select beat(s)
function hb_select(varargin)

% get figure and settings
hFigMLH = findobj('Tag', 'Wnd_NeuroElf_heartbeats');
ud = get(hFigMLH, 'UserData');
hTag = ud.hTag;

% get listbox selection
ri = hTag.LB_heartbeats_beats.Value;
bl = ud.blines;

% set all colors to green
set(bl, 'Color', [0, 0.9, 0.25], 'LineWidth', 2, 'LineStyle', '-');

% set selected colors to red
set(bl(ri), 'Color', [1, 0, 0], 'LineWidth', 2, 'LineStyle', '--');
drawnow;


% remove beat(s)
function hb_remove(varargin)

% get figure and settings
hFigMLH = findobj('Tag', 'Wnd_NeuroElf_heartbeats');
ud = get(hFigMLH, 'UserData');
hFig = ud.hFig;
hTag = ud.hTag;

% check listbox value
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1 || ...
   ~isequal(varargin{3}, 1)
    ri = hTag.LB_heartbeats_beats.Value;
else
    rs = hTag.LB_heartbeats_beats.String;
    if ~iscell(rs)
        rs = cellstr(rs);
    end
    ri = 1:numel(rs);
end
if isempty(ri)
    return;
end

% remove selected beats
ud.btpos = setdiff(ud.btpos, ud.bplot(ri));

% update
hFig.UserData = ud;

% re-draw
hb_showbeats(hFig);


% accept this window
function hb_accept(varargin)

% check nargin
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1
    return;
end

% get figure and settings
hFigMLH = findobj('Tag', 'Wnd_NeuroElf_heartbeats');
ud = get(hFigMLH, 'UserData');

% what to do
if abs(varargin{3}) == 1
    ud.btnum = ud.btnum + varargin{3};
    if ud.btnum < 1
        ud.btnum = numel(ud.bad);
    end
else
    if varargin{3} < 0
        ud.btnum = 1;
    else
        ud.btnum = numel(ud.bad);
    end
end

% update
set(hFigMLH, 'UserData', ud);

% leave ?
if ud.btnum > numel(ud.bad)
    uiresume;
    return;
end

% otherwise show
hb_showbeats(ud.hFig);


% switch filtering
function hb_fsig(varargin)

% pass on to hb_showbeats
hFigMLH = findobj('Tag', 'Wnd_NeuroElf_heartbeats');
ud = get(hFigMLH, 'UserData');
hb_showbeats(ud.hFig);


% show curve
function hb_showbeats(hFig)

% get tags
hTag = hFig.TagStruct;
ud = hFig.UserData;
ax = ud.ax;
bad = ud.bad;
btnum = ud.btnum;
btpos = ud.btpos;
cpos = ud.cpos;
oax = ud.oax;

% delete axes' children
delete(get(ax, 'Children'));
hold(ax, 'on');
delete(ud.oaxl);

% get position to dislay
dfrom = max(1, cpos(bad(btnum)) - ud.plwin);
dto = min(ud.numsig, cpos(bad(btnum)+1) + ud.plwin);
dpos = round(0.5 * (dfrom + dto));

% plot
if hTag.CB_heartbeats_fsig.Value == 0
    plot(ax, (1 / ud.freq) .* (dfrom:dto)', ztrans(ud.sig(dfrom:dto)));
else
    plot(ax, (1 / ud.freq) .* (dfrom:dto)', ztrans(ud.fsig(dfrom:dto)));
end

% plot cursor
ud.oaxl = plot(oax, [dpos; dpos] ./ ud.freq, [0; ud.oaxyl]);
set(ud.oaxl, 'Color', [0, 0.8, 0], 'LineWidth', 3);

% plot beats in window
bplot = btpos(btpos >= dfrom & btpos <= dto);
bppos = zeros(numel(bplot), 1);
bpnames = cell(numel(bplot), 1);
blines = zeros(numel(bplot), 1);
for lc = 1:numel(bplot)
    bppos(lc) = bplot(lc) / ud.freq;
    bpnames{lc} = sprintf('beat at %.3fs', bppos(lc));
    blines(lc) = plot(ax, [bppos(lc); bppos(lc)], [-3; 3]);
end
set(blines, 'Color', [0, 0.9, 0.25], 'LineWidth', 2);
ud.blines = blines;
ud.bplot = bplot;
ud.bppos = bppos;
hFig.UserData = ud;

% set string in listbox
hTag.LB_heartbeats_beats.String = bpnames;
hTag.LB_heartbeats_beats.Max = hTag.LB_heartbeats_beats.Min + 2;
hTag.LB_heartbeats_beats.Value = [];

% make figure visible
hFig.Name = sprintf('NeuroElf - heartbeats cleanup - %swindow %d of %d', ...
    ud.title, ud.btnum, numel(ud.bad));
hFig.Visible = 'on';
drawnow;

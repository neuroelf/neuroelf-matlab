function [sel, varc, spot, sets] = plotcurves(hfile, opts)
% plotcurves  - plot different curves and select those of interest
%
% FORMAT:       [sel, varc, spot, sets] = plotcurves(obj [, opts]]);
%
% Input fields:
%
%       obj         xff object with plotable data or SxC double data
%       opts        optional settings
%        .cuediff   difference between cue and onset (only useful if
%                   fixed interval!; default: -2)
%        .curves    Cx2 cell array with names and a 1x3 double array
%                   containing [channelnumber, onset, offset]
%                   - for curves in a set, the mean/std. error is computed
%                   - onset and offset are given in seconds
%        .dchannel  data channel (for onset detection, default: 1)
%        .dfilt     detection filter length (in seconds, default: 0.5)
%        .dminlat   minimum detection latency (after one point, def: 0.5)
%        .freq      data sampling frequency (if not in object, def: 100)
%        .ochannel  onset channel number (will add curves, default: [])
%        .odchannel original data channel (allows to filter effect)
%        .onsets    Ox1 onset vector in seconds
%        .owin      1x2 onset window in seconds (default: [-2, 18])
%        .sets      Sx2 cell array with names and lists of curve indices
%        .spot      Sx1 cell array with Cx1 X or Cx2 x/y time-points
%        .spotnames Sx1 cell array with strings
%        .spottype  Sx1 cell array with either of
%                   'cue', {'free'}, 'max', 'min', 'minmax', 'onset'
%        .var       1xV struct with fields (up to 6 variables)
%          .calc    either of 'dx', {'dy'}, 'mean', 'std', 'var', 'x', 'y'
%          .name    variable name
%          .spot    1x1 or 1x2 index into S
%          .trans   apply additional transformation to each value, one of
%                   {'none'}, 'log', 'log+1', 'sqrt'
%
% Output fields:
%
%       sel         selection (Cx1 boolean array)
%       varc        computed variables
%       spot        updated spots
%       sets        updates sets (new indices)
%
% Note: where appropriate, the options are taken from obj.RunTimeVars
%       if not specified in the call!

% Version:  v0.9d
% Build:    14061709
% Date:     Jun-17 2014, 9:50 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% variable for UI stuff
global ne_ui;

% allow double data to be used
delobj = false;
if nargin > 0 && ...
    isa(hfile, 'double')
    delobj = true;
    ntt = xff('new:ntt');
    ntt.Data = hfile;
    hfile = ntt;
end

% argument check
if nargin < 1 || ...
    numel(hfile) ~= 1 || ...
   ~isxff(hfile, {'acq', 'ntt'})
    error( ...
        'neuroelf:InvalidCall', ...
        'Invalid call to plotcurves.' ...
    );
end
ft = lower(hfile.Filetype);
switch (ft)
    case {'acq'}
        nc = hfile.NrOfChannels;
        freq = 1000 / hfile.MillisecsPerSample;
    case {'ntt'}
        nc = size(hfile.Data, 2);
        freq = 100;
end
if nargin < 2 || ...
    numel(opts) ~= 1 || ...
   ~isstruct(opts)
    opts = struct;
end
if isfield(opts, 'freq') && ...
    isa(opts.freq, 'double') && ...
    numel(opts.freq) == 1 && ...
   ~isinf(opts.freq) && ...
   ~isnan(opts.freq) && ...
    opts.freq >= 0.01
    freq = opts.freq;
end
if ~isfield(opts, 'cuediff') || ...
   ~isa(opts.cuediff, 'double') || ...
    numel(opts.cuediff) ~= 1 || ...
    isinf(opts.cuediff) || ...
    isnan(opts.cuediff)
    opts.cuediff = -2;
end
if ~isfield(opts, 'curves') || ...
   ~iscell(opts.curves) || ...
    size(opts.curves, 2) ~= 2 || ...
    ndims(opts.curves) > 2
    opts.curves = cell(0, 2);
end
for cc = 1:size(opts.curves, 1)
    if ~ischar(opts.curves{cc, 1}) || ...
        isempty(opts.curves{cc, 1}) || ...
       ~isa(opts.curves{cc, 2}, 'double') || ...
        size(opts.curves{cc, 2}, 2) ~= 3 || ...
        any(isinf(opts.curves{cc, 2}(:)) | isnan(opts.curves{cc, 2}(:))) || ...
        any(any(opts.curves{cc, 2} < 1, 2) | opts.curves{cc, 2}(:, 1) > nc)
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid curve selection for curve %d.', ...
            cc ...
        );
    end
end
if ~isfield(opts, 'dchannel') || ...
   ~isa(opts.dchannel, 'double') || ...
    numel(opts.dchannel) ~= 1 || ...
    isinf(opts.dchannel) || ...
    isnan(opts.dchannel) || ...
    opts.dchannel < 1 || ...
    opts.dchannel ~= fix(opts.dchannel) || ...
    opts.dchannel > nc
    opts.dchannel = 1;
end
if ~isfield(opts, 'dfilt') || ...
   ~isa(opts.dfilt, 'double') || ...
    numel(opts.dfilt) ~= 1 || ...
    isinf(opts.dfilt) || ...
    isnan(opts.dfilt) || ...
    opts.dfilt < 0 || ...
    opts.dfilt > 30
    opts.dfilt = 0.5;
end
if opts.dfilt == 0;
    dfiltk = [];
else
    dfiltk = smoothkern(freq * opts.dfilt, 1.42e-5 / (freq * opts.dfilt));
    dfiltk = flexinterpn_method(dfiltk, [Inf; 1; 1 / 4096; numel(dfiltk)], 'cubic');
end
if ~isfield(opts, 'dminlat') || ...
   ~isa(opts.dminlat, 'double') || ...
    numel(opts.dminlat) ~= 1 || ...
    isinf(opts.dminlat) || ...
    isnan(opts.dminlat) || ...
    opts.dminlat < 0
    opts.dminlat = 0.5;
end
if ~isfield(opts, 'owin') || ...
   ~isa(opts.owin, 'double') || ...
    isempty(opts.owin) || ...
    numel(opts.owin) > 2 || ...
    any(isinf(opts.owin) | isnan(opts.owin)) || ...
    (opts.owin(end) - opts.owin(1)) < 0 || ...
    (opts.owin(end) - opts.owin(1)) > 60
    opts.owin = [-2, 18];
elseif numel(opts.owin) == 1
    opts.owin = opts.owin + [0, 20];
end
if ~isfield(opts, 'ochannel') || ...
   ~isa(opts.ochannel, 'double') || ...
    numel(opts.ochannel) ~= 1 || ...
    isinf(opts.ochannel) || ...
    isnan(opts.ochannel) || ...
    opts.ochannel < 1 || ...
    opts.ochannel ~= fix(opts.ochannel) || ...
    opts.ochannel > nc
    opts.ochannel = [];
    onsets = [];
else
    och = hfile.ChannelData(opts.ochannel);
    mmoch = minmaxmean(och);
    mmoch = 0.5 * (mmoch(1) + mmoch(2));
    och = (och >= mmoch);
    onsets = (1 / freq) .* find(~och(1:end-1) & och(2:end));
end
if ~isfield(opts, 'odchannel') || ...
   ~isa(opts.odchannel, 'double') || ...
    numel(opts.odchannel) ~= 1 || ...
    isinf(opts.odchannel) || ...
    isnan(opts.odchannel) || ...
    opts.odchannel < 1 || ...
    opts.odchannel ~= fix(opts.odchannel) || ...
    opts.odchannel > nc
    opts.odchannel = [];
end
if isfield(opts, 'onsets') && ...
    isa(opts.onsets, 'double') && ...
    numel(opts.onsets) == max(size(opts.onsets)) && ...
   ~any(isinf(opts.onsets) | isnan(opts.onsets) | opts.onsets < 0)
    onsets = opts.onsets(:);
end
if ~isempty(onsets)
    ncr = size(opts.curves, 1);
    opts.curves = [opts.curves; cell(numel(onsets), 2)];
    for cc = 1:numel(onsets)
        opts.curves{ncr + cc, 1} = ...
            sprintf('Onset %d, t=%.3f', cc, onsets(cc));
        opts.curves{ncr + cc, 2} = [opts.dchannel, onsets(cc) + opts.owin];
    end
end
if ~isfield(opts, 'sets') || ...
   ~iscell(opts.sets) || ...
    isempty(opts.sets) || ...
    size(opts.sets, 2) ~= 2 || ...
    ndims(opts.sets) > 2
    if nc == (size(opts.curves, 1) + numel(onsets)) || ...
       ~isfield(opts, 'sets') || ...
       ~iscell(opts.sets) || ...
        size(opts.sets, 2) ~= 2 || ...
        ndims(opts.sets) > 2
        opts.sets = {'Overall mean', 1:nc};
        if ~isempty(onsets)
            opts.sets(2, :) = {'Onsets mean', nc+1:nc+numel(onsets)};
        end
    end
end
if ~isfield(opts, 'spot') || ...
   ~iscell(opts.spot) || ...
    all(cellfun('isempty', opts.spot))
    if ~isfield(hfile.RunTimeVars, 'plotcurves') || ...
       ~isstruct(hfile.RunTimeVars.plotcurves) || ...
       ~isfield(hfile.RunTimeVars.plotcurves, 'spot')
        if ~isfield(opts, 'spot') || ...
           ~iscell(opts.spot)
            opts.spot = {};
        end
    else
        opts.spot = hfile.RunTimeVars.plotcurves.spot;
    end
else
    opts.spot = opts.spot(:);
end
if isempty(opts.spot)
    opts.spot = {cat(1, opts.curves{:, 2})};
    opts.spot{1} = opts.spot{1}(:, 2);
end
if ~isfield(opts, 'spotnames') || ...
   ~iscell(opts.spotnames) || ...
    numel(opts.spotnames) ~= numel(opts.spot)
    opts.spotnames = cell(1, numel(opts.spot));
else
    opts.spotnames = opts.spotnames(:);
end
if ~isfield(opts, 'spottype') || ...
   ~iscell(opts.spottype) || ...
    numel(opts.spottype) ~= numel(opts.spot)
    opts.spottype = cell(1, numel(opts.spot));
else
    opts.spottype = opts.spottype(:);
end
for sc = numel(opts.spot):-1:1
    if ~isa(opts.spot{sc}, 'double') || ...
        ndims(opts.spot{sc}) ~= 2 || ...
       (~isempty(opts.spot{sc}) && ...
        size(opts.spot{sc}, 1) ~= size(opts.curves, 1)) || ...
       ~any(size(opts.spot{sc}, 2) == [1, 2]) || ...
        any(isinf(opts.spot{sc}(:)) | isnan(opts.spot{sc}(:))) || ...
        any(opts.spot{sc}(:, 1) < 0)
        opts.spot(sc) = [];
        opts.spotnames(sc) = [];
        opts.spottype(sc) = [];
        continue;
    end
    if ~ischar(opts.spotnames{sc}) || ...
        isempty(opts.spotnames{sc})
        opts.spotnames{sc} = sprintf('Spot %d', sc);
    else
        opts.spotnames{sc} = opts.spotnames{sc}(:)';
    end
    if ~ischar(opts.spottype{sc}) || ...
        isempty(opts.spottype{sc}) || ...
       ~any(strcmpi(opts.spottype{sc}(:)', ...
            {'c', 'cue', 'f', 'free', 'max', 'min', 'minmax', 'n', 'none', 'o', 'onset'}))
        opts.spottype{sc} = 'free';
    else
        opts.spottype{sc} = lower(opts.spottype{sc}(:)');
    end
end
findspot = false;
for sc = 1:numel(opts.spot)
    if isempty(opts.spot{sc})
        findspot = true;
        opts.spot{sc} = zeros(size(opts.curves, 1), 1);
        for ssc = 1:size(opts.curves, 1)
            spcrv = opts.curves{ssc, 2};
            if opts.spottype{sc}(1) == 'm'
                spotc = diff(hfile.SampleData( ...
                    spcrv(1), spcrv(2), spcrv(3), ...
                    ceil(freq * (spcrv(3) - spcrv(2))), freq));
                if ~isempty(dfiltk)
                    spotc = flexinterpn(spotc, ...
                        [Inf; 1; 1; numel(spotc)], dfiltk, 4096, 0);
                end
                spotc = sign(spotc);
                spotc(find(spotc(1:end-1) ~= spotc(2:end))) = 0;
                if sc == 1
                    spots = ceil(1 - freq * opts.owin(1));
                else
                    spots = min(numel(spotc), ...
                        ceil(2 + freq * (opts.dminlat + ...
                        opts.spot{sc-1}(ssc, 1) - spcrv(2))));
                end
            end
            switch (opts.spottype{sc})
                case {'cue'}
                    spott = ceil(1 - freq * (opts.owin(1) + opts.cuediff));
                case {'f', 'free', 'n', 'none', 'o', 'onset'}
                    spott = ceil(1 - freq * opts.owin(1));
                case {'max'}
                    if spotc(spots <= 0)
                        spots = findfirst(spotc > 0, spots);
                        if isempty(spots)
                            spots = numel(spotc);
                        end
                    end
                    spott = findfirst(spotc < 0, spots);
                case {'min'}
                    if spotc(spots >= 0)
                        spots = findfirst(spotc < 0, spots);
                        if isempty(spots)
                            spots = numel(spotc);
                        end
                    end
                    spott = findfirst(spotc > 0, spots);
                case {'minmax'}
                    spott = findfirst(sign(spotc) == -sign(spotc(1)), spots);
            end
            if isempty(spott)
                spott = numel(spotc);
            end
            opts.spot{sc}(ssc, 1) = spcrv(2) + (spott - 1) / freq;
        end
    end
end
if ~isfield(opts, 'var') || ...
   ~isstruct(opts.var) || ...
    isempty(opts.var) || ...
   ~isfield(opts.var, 'calc') || ...
   ~isfield(opts.var, 'spot') || ...
   ~isfield(opts.var, 'trans')
    opts.var = struct('calc', 'y', 'spot', 1, 'trans', 'none');
end
opts.var = opts.var(:);
for vc = numel(opts.var):-1:1
    if ~isa(opts.var(vc).spot, 'double') || ...
       ~any(numel(opts.var(vc).spot) == [1, 2]) || ...
        any(isinf(opts.var(vc).spot) | isnan(opts.var(vc).spot) | ...
            opts.var(vc).spot < 1 | opts.var(vc).spot > numel(opts.spot) | ...
            opts.var(vc).spot ~= fix(opts.var(vc).spot)) || ...
        numel(unique(opts.var(vc).spot)) ~= numel(opts.var(vc).spot)
        opts.var(vc) = [];
        continue;
    end
    if ~ischar(opts.var(vc).calc) || ...
        isempty(opts.var(vc).calc) || ...
       ~any(strcmpi(opts.var(vc).calc(:)', ...
            {'dx', 'dy', 'm', 'mean', 's', 'std', 'v', 'var', 'x', 'y'}))
        if numel(opts.var(vc).spot) > 1
            opts.var(vc).calc = 'dy';
        else
            opts.var(vc).calc = 'y';
        end
    else
        opts.var(vc).calc = lower(opts.var(vc).calc(:)');
    end
    if any(strcmp(opts.var(vc).calc, ...
            {'dx', 'dy', 'm', 'mean', 's', 'std', 'v', 'var'})) && ...
        numel(opts.var(vc).spot) < 2
        error( ...
            'neuroelf:BadArgument', ...
            'The %d calculation requires 2 spots.', ...
            opts.var(vc).calc ...
        );
    end
    if ~ischar(opts.var(vc).trans) || ...
        isempty(opts.var(vc).trans) || ...
       ~any(strcmpi(opts.var(vc).trans(:)', ...
            {'n', 'none', 'l', 'log', 'log+1', 'r', 's', 'sqrt'}))
        opts.var(vc).trans = 'none';
    else
        opts.var(vc).trans = lower(opts.var(vc).trans(:)');
    end
    if ~isfield(opts.var(vc), 'name') || ...
       ~ischar(opts.var(vc).name) || ...
        isempty(opts.var(vc).name)
        opts.var(vc).name = sprintf('Var %d', vc);
    end
end
if numel(opts.var) > 6
    opts.var = opts.var(1:6);
end
if isempty(opts.curves)
    opts.curves = cell(nc, 2);
    for cc = 1:nc
        opts.curves{cc, 1} = sprintf('channel %d', cc);
        opts.curves{cc, 2} = [cc, 0, numel(hfile.ChannelData(cc)) / freq];
    end
end

% fill in initial spot y values
for sc = 1:numel(opts.spot)
    if size(opts.spot{sc}, 2) == 1
        opts.spot{sc}(:, 2) = 0;
        for cc = 1:size(opts.curves, 1)
            opts.spot{sc}(cc, 2) = hfile.SampleData(opts.curves{cc, 2}(1),  ...
                opts.spot{sc}(cc, 1), opts.spot{sc}(cc, 1), 1, freq);
        end
    end
end

% prepare
ne_ui.plotcurves = struct( ...
    'ax',    -1, ...
    'axlx',  [], ...
    'axly',  [], ...
    'axpos', [-1, -1], ...
    'btdwn', false, ...
    'cpbtd', [-1, -1], ...
    'cpos',  [-1, -1], ...
    'freq',  freq, ...
    'hFig',  [], ...
    'hFigM', -1, ...
    'hTag',  struct, ...
    'hevnt', false, ...
    'mctrl', false, ...
    'mshft', false, ...
    'obj',   hfile, ...
    'opts',  opts, ...
    'ovch',  0, ...
    'sel',   [], ...
    'spot',  [], ...
    'spota', [], ...
    'spotx', [], ...
    'varc',  zeros(size(opts.curves, 1), numel(opts.var)), ...
    'zoomp', []);

% pre-compute values
pc_compute;

% load figure
try
    hFig = neuroelf_file('f', 'plotcurves');
    hTag = hFig.TagStruct;

catch ne_eo;
    error( ...
        'neuroelf:xfigureError', ...
        'Error creating UI for colorpicker: %s.', ...
        ne_eo.message ...
    );
end

% make settings
if isfield(hfile.RunTimeVars, 'plotcurves') && ...
    isstruct(hfile.RunTimeVars.plotcurves) && ...
    isfield(hfile.RunTimeVars.plotcurves, 'sel') && ...
    numel(hfile.RunTimeVars.plotcurves.sel) == size(opts.curves, 1)
    sel = hfile.RunTimeVars.plotcurves.sel;
else
    if findspot && ...
       ~isempty(opts.var) && ...
       ~isempty(opts.var(1).name) && ...
        lower(opts.var(1).name(1)) == 'l'
        ovarc = ne_ui.plotcurves.varc;
        [mlat, wlat] = robustmean(ovarc(:, 1));
        mvarc = mean(ovarc(wlat >= 0.5, :));
        svarc = std(ovarc(wlat >= 0.5, :));
        lvarc = repmat(mvarc - 1.5 * svarc, size(ovarc, 1), 1);
        uvarc = repmat(mvarc + 1.5 * svarc, size(ovarc, 1), 1);
        sel = (wlat >= 0.5) & ...
            ~(sum(ovarc < lvarc | ovarc > uvarc, 2) >= (1 + 0.25 * size(ovarc, 2)));
    else
        sel = true(size(opts.curves, 1), 1);
    end
end
cnames = [opts.sets(:, 1); opts.curves(:, 1)];
for cnc = size(opts.sets, 1)+1:numel(cnames)
    if sel(cnc - size(opts.sets, 1))
        cnames{cnc} = sprintf('(+) %s', cnames{cnc});
    else
        cnames{cnc} = sprintf('(-) %s', cnames{cnc});
    end
end
hTag.CB_plotcurves_focurve.Value = 1;
hTag.LB_plotcurves_curves.String = cnames;
hTag.LB_plotcurves_curves.Value = 1 + size(opts.sets, 1);
set(hTag.AX_plotcurves_zoombar.MLHandle, ...
    'Color',  [.6, .6, .6], ...
    'XColor', [.6, .6, .6], ...
    'YColor', [.6, .6, .6], ...
    'XLim',   [0, 1], ...
    'YLim',   [0, 1], ...
    'XTick',  [], ...
    'YTick',  [], ...
    'ZTick',  []);
for vc = 1:numel(opts.var)
    hFig.SetGroupVisible(sprintf('GVar%d', vc));
    hTag.(sprintf('TX_plotcurves_var%d', vc)).String = opts.var(vc).name;
end
chax = hTag.AX_plotcurves_plot.MLHandle;
set(chax, 'Units', 'pixels');
if isempty(opts.odchannel)
    hTag.Plot = plot(chax, [0; 1], [0; 1]);
else
    hTag.Plot = plot(chax, [0; 1], [0, 0; 1, 1]);
end
set(hTag.Plot(1), 'LineWidth', 2);
hTag.CursorX = line([0.5; 0.5], [0.001; 0.999], 'Color', [0, 0, 0], 'Parent', chax);
hTag.CursorY = line([0; 0.999], [0.5; 0.5], 'Color', [0, 0, 0], 'Parent', chax);
axpos = get(chax, 'Position');
hold(hTag.AX_plotcurves_plot.MLHandle, 'on');

% put in global structure
ne_ui.plotcurves.ax = chax;
ne_ui.plotcurves.axpos = axpos;
ne_ui.plotcurves.cpos = round(0.5 * hTag.AX_plotcurves_plot.Position(3:4));
ne_ui.plotcurves.hFig = hFig;
ne_ui.plotcurves.hFigM = hFig.MLHandle;
ne_ui.plotcurves.hTag = hTag;
ne_ui.plotcurves.sel = sel;

% set callbacks
hTag.DD_plotcurves_filter.Callback = @pc_setfilter;
hTag.LB_plotcurves_curves.Callback = @pc_showcurve;
hTag.CB_plotcurves_usecurve.Callback = @pc_toggleusecurve;
hTag.ED_plotcurves_var1.Callback = {@pc_updatevar, 1};
hTag.ED_plotcurves_var2.Callback = {@pc_updatevar, 2};
hTag.ED_plotcurves_var3.Callback = {@pc_updatevar, 3};
hTag.ED_plotcurves_var4.Callback = {@pc_updatevar, 4};
hTag.ED_plotcurves_var5.Callback = {@pc_updatevar, 5};
hTag.ED_plotcurves_var6.Callback = {@pc_updatevar, 6};
hFig.KeyPressFcn = @pc_keypress;
hFig.KeyReleaseFcn = @pc_keyrelease;
hFig.WindowButtonDownFcn = @pc_btdown;
hFig.WindowButtonMotionFcn = @pc_btmove;
hFig.WindowButtonUpFcn = @pc_btup;

% show current selection
pc_showcursor;
pc_showcurve;

% wait for dialog to finish
hFig.HandleVisibility = 'callback';
hFig.Visible = 'on';
uiwait(hFig.MLHandle);

% get selection
sel = ne_ui.plotcurves.sel;

% apply transformation
varc = ne_ui.plotcurves.varc;
for vc = 1:numel(opts.var)
    switch (opts.var(vc).trans)
        case {'l', 'ln', 'log'}
            ne_ui.plotcurves.varc(:, vc) = log(varc(:, vc));
        case {'log+1'}
            ne_ui.plotcurves.varc(:, vc) = log(1 + varc(:, vc));
        case {'log10'}
            ne_ui.plotcurves.varc(:, vc) = log10(varc(:, vc));
        case {'log2'}
            ne_ui.plotcurves.varc(:, vc) = log2(varc(:, vc));
        case {'r', 'root', 'sqrt'}
            ne_ui.plotcurves.varc(:, vc) = sqrt(varc(:, vc));
    end
end

% delete object if necessary
if delobj
    hfile.ClearObject;

% otherwise...
else

    % add output to RunTimeVars of object
    hfile.RunTimeVars.plotcurves = struct( ...
        'sel',       sel, ...
        'sets',      {ne_ui.plotcurves.opts.sets}, ...
        'spot',      {ne_ui.plotcurves.opts.spot}, ...
        'spotnames', {opts.spotnames}, ...
        'spottype',  {opts.spottype}, ...
        'var',       {opts.var}, ...
        'varc',      ne_ui.plotcurves.varc);

    % try to save
    try
        if ~isempty(hfile.FilenameOnDisk)
            hfile.SaveRunTimeVars;
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        warning( ...
            'neuroelf:SaveMATFailed', ...
            'Couldn''t save RunTimeVars to disk; please store manually.' ...
        );
    end
end

% output arguments
if nargout > 2
    varc = ne_ui.plotcurves.varc;
    if nargout > 3
        spot = ne_ui.plotcurves.spot;
        if nargout > 4
            sets = ne_ui.plotcurves.sets;
        end
    end
end

% remove field
ne_ui = rmfield(ne_ui, 'plotcurves');



% functions



% compute variables
function pc_compute(varargin)
global ne_ui;
c = ne_ui.plotcurves;

% curve selection
if nargin > 0 && ...
    isa(varargin{1}, 'double')
    cis = varargin{1}(:)';
else
    cis = 1:size(c.opts.curves, 1);
end

% get conf
curve = c.opts.curves(:, 2);
obj = c.obj;
spot = c.opts.spot;
var = c.opts.var;

% create new varc
varc = zeros(numel(cis), numel(var));

% iterate over curves
for cc = 1:numel(cis)

    % curve index
    ci = cis(cc);

    % iterate over vars
    for vc = 1:numel(var)

        % get snippet of data that we need
        if numel(var(vc).spot) > 1 && ...
            any(var(vc).calc(1) == 'msv')
            snip = obj.ChannelData(curve{ci}(1), floor(x1):ceil(x2));
        end

        % depending on calculation type
        switch (var(vc).calc)
            case {'dx'}
                val = spot{var(vc).spot(2)}(ci, 1) - spot{var(vc).spot(1)}(ci, 1);
            case {'dy'}
                val = spot{var(vc).spot(2)}(ci, 2) - spot{var(vc).spot(1)}(ci, 2);
            case {'m', 'mean'}
                val = sum(snip) / numel(snip);
            case {'s', 'std'}
                if numel(snip) > 1
                    val = std(snip);
                else
                    val = 0;
                end
            case {'v', 'var'}
                if numel(snip) > 1
                    val = varc(snip);
                else
                    val = 0;
                end
            case {'x'}
                val = spot{var(vc).spot(1)}(ci, 1);
            case {'y'}
                val = spot{var(vc).spot(1)}(ci, 2);
        end

        % store in matrix
        varc(cc, vc) = val;
    end
end

% store in global matrix
ne_ui.plotcurves.varc(cis, :) = varc;


% update var from UI
function pc_updatevar(varargin)
global ne_ui;
c = ne_ui.plotcurves;

% argument check
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1 || ...
   ~any(varargin{3} == 1:6)
    return;
end
vars = varargin{3};

% get current selection
csel = c.hTag.LB_plotcurves_curves.Value;
if isempty(csel) || ...
    csel <= size(c.opts.sets, 1)
    return;
end
csel = csel - size(c.opts.sets, 1);

% get vars
varc = c.varc;

% get current value
curval = varc(csel, vars);

% get string
uic = sprintf('ED_plotcurves_var%d', vars);
newval = c.hTag.(uic).String;

% check
try
    newval = str2double(newval);
    if numel(newval) ~= 1 || ...
        isinf(newval) || ...
        isnan(newval)
        newval = curval;
    end
    ne_ui.plotcurves.varc(csel, vars) = newval;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    newval = curval;
end

% set new string
c.hTag.(uic).String = sprintf('%g', newval);


% hit mouse
function pc_btdown(varargin)
global ne_ui;
c = ne_ui.plotcurves;

% check if within axes
f = c.hFig;

% get current point and compare to axes
cp = f.CurrentPoint - c.axpos(1:2);
if ~all(cp >= 0) || ...
   ~all(cp <= c.axpos(3:4))
    return;
end

% set button down true
ne_ui.plotcurves.btdwn = true;
ne_ui.plotcurves.cpbtd = cp;


% release mouse
function pc_btup(varargin)
global ne_ui;

% set button down false
ne_ui.plotcurves.btdwn = false;

% now, get config
c = ne_ui.plotcurves;

% clear potential zoom patch
if ~isempty(c.zoomp)
    delete(c.zoomp);
    ne_ui.plotcurves.zoomp = [];
    pc_handlezoom;
end


% move cursor
function pc_btmove(varargin)
global ne_ui;
c = ne_ui.plotcurves;

% check figure
if nargin < 1 || ...
    varargin{1} ~= c.hFigM
    return;
end

% more settings
f = c.hFig;

% get current point and compare to axes
cp = f.CurrentPoint - c.axpos(1:2);
if (~all(cp >= 0) || ...
    ~all(cp <= c.axpos(3:4)))
    return;
end

% set into global structure and show
ne_ui.plotcurves.cpos = cp;
pc_showcursor;

% the button was pressed
if c.btdwn
    pc_handlemove;
end


% handle move
function pc_handlemove
global ne_ui;
c = ne_ui.plotcurves;

% only handle one event
if c.hevnt
    return;
end
ne_ui.plotcurves.hevnt = true;

% spot
if ~isempty(c.spota)

    % get curve number
    csel = c.hTag.LB_plotcurves_curves.Value;
    if isempty(csel) || ...
        csel <= size(c.opts.sets, 1)
        ne_ui.plotcurves.hevnt = false;
        return;
    end
    csel = csel - size(c.opts.sets, 1);

    % current spot positions
    spotp = zeros(numel(c.opts.spot), 2);
    for spc = 1:size(spotp, 1)
        spotp(spc, :) = c.opts.spot{spc}(csel, 1:2);
    end

    % what position to put the spot to
    spos = pc_axescp(c.cpos);
    if ~c.mctrl
        if c.spota > 1
            spos(1) = max(spotp(c.spota - 1, 1), spos(1));
        end
        if c.spota < size(spotp, 1)
            spos(1) = min(spotp(c.spota + 1, 1), spos(1));
        end
    end

    % force to curve
    if c.hTag.CB_plotcurves_focurve.Value > 0

        % which data
        if c.mshft && ...
           ~isempty(c.opts.odchannel)
            chan = c.opts.odchannel;
        else
            chan = c.opts.curves{csel, 2}(1);
        end

        % sample
        spos(2) = c.obj.SampleData(chan, spos(1), spos(1), 1, c.freq);
    end

    % update spot
    ne_ui.plotcurves.opts.spot{c.spota}(csel, :) = spos;

    % replot spots
    pc_showspots(csel);

    % recompute variables
    pc_compute(csel);

    % update variables
    varc = ne_ui.plotcurves.varc(csel, :);
    for vc = 1:numel(varc)
        c.hTag.(sprintf('ED_plotcurves_var%d', vc)).String = ...
            sprintf('%g', varc(vc));
    end

% zoom
else

    % create zoom patch
    if isempty(c.zoomp)

        % get position
        zpos = pc_axescp(c.cpos);
        ne_ui.plotcurves.zoomp = ...
            plot(c.ax, zpos(ones(1, 5), 1), zpos(ones(1, 5), 2));
        set(ne_ui.plotcurves.zoomp, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
        ne_ui.plotcurves.hevnt = false;
        return;
    end

    % give correct zoomp position
    zpos = pc_axescp(c.cpos);
    xd = get(c.zoomp, 'XData');
    yd = get(c.zoomp, 'YData');
    xd(2:3) = zpos(1);
    yd(3:4) = zpos(2);
    set(c.zoomp, 'XData', xd, 'YData', yd);
end
ne_ui.plotcurves.hevnt = false;


% handle zoom
function pc_handlezoom
global ne_ui;
c = ne_ui.plotcurves;

% only handle one event
if c.hevnt
    return;
end
ne_ui.plotcurves.hevnt = true;

% get two positions
zpos = [pc_axescp(c.cpbtd); pc_axescp(c.cpos)];
minzp = min(zpos);
maxzp = max(zpos);

% set new axes limits
pc_showcurve(([minzp; maxzp])');

% unblock event
ne_ui.plotcurves.hevnt = false;


% handle keyboard input
function pc_keypress(src, ke, varargin)
global ne_ui;

% check src
if numel(src) ~= 1
    return;
end

% configuration
c = ne_ui.plotcurves;
axlx = c.axlx;
axly = c.axly;
axdx = axlx(2) - axlx(1);
axdy = axly(2) - axly(1);
hTag = c.hTag;

% get Key and Modifier from keyboard event (see Matlab docu!)
kk = ke.Key;
mn = ke.Modifier;

% determine which modifiers are pressed
km = false(1, 4);
if ~isempty(mn)
    try
        km = [ ...
            any(strcmpi('alt', mn)), ...
            any(strcmpi('control', mn)), ...
            any(strcmpi('shift', mn)), ...
            any(strcmpi('command', mn))];
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% what to do
if ~any(km)
    switch lower(kk)
        case {'downarrow'}
            hTag.LB_plotcurves_curves.Value = ...
                min(numel(hTag.LB_plotcurves_curves.String), ...
                hTag.LB_plotcurves_curves.Value + 1);
            pc_showcurve;
        case {'uparrow'}
            hTag.LB_plotcurves_curves.Value = ...
                max(1, hTag.LB_plotcurves_curves.Value - 1);
            pc_showcurve;
        case {'d'}
            axlx = [axlx(1) - 0.25 * axdx, axlx(2) + 0.25 * axdx];
            axly = [axly(1) - 0.25 * axdy, axly(2) + 0.25 * axdy];
            pc_showcurve([axlx; axly]);
        case {'f'}
            hTag.CB_plotcurves_focurve.Value = ...
                double(hTag.CB_plotcurves_focurve.Value == 0);
        case {'i'}
            if strcmpi(hTag.CB_plotcurves_usecurve.Enable, 'on')
                hTag.CB_plotcurves_usecurve.Value = 1;
                pc_toggleusecurve;
            end
        case {'m'}
            axlx = [axlx(1) + 0.2 * axdx, axlx(2) - 0.2 * axdx];
            axly = [axly(1) + 0.2 * axdy, axly(2) - 0.2 * axdy];
            pc_showcurve([axlx; axly]);
        case {'o'}
            if strcmpi(hTag.CB_plotcurves_usecurve.Enable, 'on')
                hTag.CB_plotcurves_usecurve.Value = 0;
                pc_toggleusecurve;
            end
        case {'p'}
            if ~isempty(c.spota)
                if c.mshft && ...
                    numel(c.hTag.Plot) > 1
                    pobj = c.hTag.Plot(end);
                else
                    pobj = c.hTag.Plot(1);
                end
                csel = c.hTag.LB_plotcurves_curves.Value;
                if isempty(csel) || ...
                    csel <= size(c.opts.sets, 1)
                    return;
                end
                csel = csel - size(c.opts.sets, 1);
                pobjxd = get(pobj, 'XData');
                pobjyd = get(pobj, 'YData');
                mmyd = minmaxmean(pobjyd);
                switch (c.opts.spottype{c.spota})
                    case {'f', 'free', 'minmax'}
                        nyd = floor(numel(yd) / 3);
                        mbyd = (0.5 / nyd) * ...
                            (sum(pobjyd(1:nyd)) + sum(pobjyd(end-nyd:end)));
                        mcyd = sum(pobj(nyd:2*nyd)) / (nyd + 1);
                        if mbyd > mcyd
                            npx = mmyd(4);
                        else
                            npx = mmyd(5);
                        end
                    case {'min'}
                        npx = mmyd(4);
                    case {'max'}
                        npx = mmyd(5);
                    otherwise
                        return;
                end
                ne_ui.plotcurves.opts.spot{c.spota}(csel, :) = ...
                    [pobjxd(npx), pobjyd(npx)];
                pc_showspots(csel);
            end
        case {'r'}
            pc_showcurve;
        case {'u'}
            if strcmpi(hTag.CB_plotcurves_usecurve.Enable, 'on')
                hTag.CB_plotcurves_usecurve.Value = ...
                    double(hTag.CB_plotcurves_usecurve.Value == 0);
                pc_toggleusecurve;
            end
    end

% with control/shift
else
    if km(2)
        ne_ui.plotcurves.mctrl = true;
    end
    if km(3)
        ne_ui.plotcurves.mshft = true;
        switch lower(kk)
            case {'downarrow'}
                axly = axly - 0.1 * axdy;
                pc_showcurve([axlx; axly]);
            case {'leftarrow'}
                axlx = axlx - 0.1 * axdx;
                pc_showcurve([axlx; axly]);
            case {'rightarrow'}
                axlx = axlx + 0.1 * axdx;
                pc_showcurve([axlx; axly]);
            case {'uparrow'}
                axly = axly + 0.1 * axdy;
                pc_showcurve([axlx; axly]);
        end
    end
end

function pc_keyrelease(varargin)
global ne_ui;
ne_ui.plotcurves.mctrl = false;
ne_ui.plotcurves.mshft = false;


% set filter
function pc_setfilter(varargin)
global ne_ui;
c = ne_ui.plotcurves;

% pass on with correct arguments
pc_showcurve([c.axlx; c.axly]);


% show cursor
function pc_showcursor(varargin)
global ne_ui;
c = ne_ui.plotcurves;
f = c.hFig;

% compute new position for lines
cp = c.cpos;
xl = get(c.ax, 'XLim');
yl = get(c.ax, 'YLim');
cp = cp ./ c.axpos(3:4);
xd = xl(2) - xl(1);
yd = yl(2) - yl(1);
xp = xl(1) + xd * cp(1);
yp = yl(1) + yd * cp(2);
xc = xl(1) + xd .* min(1, max(0, cp(1) + [-0.05, 0.05]));
yc = yl(1) + yd .* min(1, max(0, cp(2) + [-0.07, 0.07]));
set(c.hTag.CursorX, 'XData', [xp, xp], 'YData', yc);
set(c.hTag.CursorY, 'XData', xc, 'YData', [yp, yp]);
set(c.ax, 'XLim', xl, 'YLim', yl);
set(f.MLHandle, 'Name', ...
    sprintf('NeuroElf - plotcurves  -  x = %.3f, y = %.4f', mean(xc), mean(yc)));

% no spots or currently moving
if isempty(c.spot) || ...
    c.btdwn
    return;
end

% set all spots to standard
set(c.spot, 'LineStyle', '-');
ne_ui.plotcurves.spota = [];

% check if a spot is closeby
spotx = [(c.spotx(:, 1) - xl(1)) ./ xd, (c.spotx(:, 2) - yl(1)) ./ yd];
spotxd = sqrt(sum((spotx - cp(ones(1, size(spotx, 1)), :)) .^ 2, 2));
spotxd(~cellfun('isempty', regexpi(c.opts.spottype, '^(o|onset|c|cue)$'))) = Inf;
[minspd, minspdp] = min(spotxd);
if minspd <= 0.08
    set(c.spot(minspdp), 'LineStyle', '--');
    ne_ui.plotcurves.spota = minspdp;
else
    ne_ui.plotcurves.spota = [];
end


% show curve
function pc_showcurve(varargin)
global ne_ui;
c = ne_ui.plotcurves;

% get current selection
csel = c.hTag.LB_plotcurves_curves.Value;
if isempty(csel)
    return;
end

% get handles
curve = c.opts.curves(:, 2);
obj = c.obj;
var = c.opts.var;
varc = c.varc;

% set
if csel <= size(c.opts.sets, 1)

    % update checkbox
    c.hTag.CB_plotcurves_usecurve.Enable = 'off';
    c.hTag.CB_plotcurves_usecurve.Value = 0;

% curve
else

    % what curve
    csel = csel - size(c.opts.sets, 1);

    % limits
    if nargin < 1 || ...
       ~isa(varargin{1}, 'double') || ...
       ~isequal(size(varargin{1}), [2, 2])
        axlx = curve{csel}(2:3);
        axly = [];
    else
        axlx = varargin{1}(1, :);
        axly = varargin{1}(2, :);
    end

    % additional filter
    flt = str2double( ...
        c.hTag.DD_plotcurves_filter.String{c.hTag.DD_plotcurves_filter.Value}(1:4));

    % resample curve
    xft = curve{csel}(1);
    cdata = obj.SampleData(xft, axlx(1), axlx(2), 1000, c.freq, flt);
    xdata = (axlx(1):(axlx(2) - axlx(1))/999.99:axlx(2))';

    % make plot
    set(c.hTag.Plot(1), 'XData', xdata, 'YData', cdata);
    if ~isempty(c.opts.odchannel)
        set(c.hTag.Plot(2), 'XData', xdata, 'YData', obj.SampleData( ...
            c.opts.odchannel, axlx(1), axlx(2), 1000, c.freq, 0, 'linear'));
    end

    % update limits
    if isempty(axly)
        mmms = minmaxmean(cdata, 5);
        mmmv = 0.2 * sqrt(mmms(6) + eps);
        axly = [mmms(1) - mmmv, mmms(2) + mmmv];
    end
    c.hTag.AX_plotcurves_plot.XLim = axlx;
    c.hTag.AX_plotcurves_plot.YLim = axly;

    % update checkbox
    c.hTag.CB_plotcurves_usecurve.Enable = 'on';
    c.hTag.CB_plotcurves_usecurve.Value = double(c.sel(csel));

    % update overview plot if necessary
    if xft ~= c.ovch
        ne_ui.plotcurves.ovch = xft;
        lastsmp = numel(obj.ChannelData(xft)) / c.freq;
        xdata = (0:lastsmp/999.99:lastsmp)';
        plot(c.hTag.AX_plotcurves_overview.MLHandle, xdata, ...
            obj.SampleData(xft, 0, lastsmp, 1000, c.freq));
    end
    ne_ui.plotcurves.axlx = axlx;
    ne_ui.plotcurves.axly = axly;

    % show spots
    pc_showspots(csel);

    % update variables
    for vc = 1:numel(var)
        c.hTag.(sprintf('ED_plotcurves_var%d', vc)).String = ...
            sprintf('%g', varc(csel, vc));
    end
end

% update cursor
pc_showcursor;


% show spots
function pc_showspots(csel)
global ne_ui;
c = ne_ui.plotcurves;
spot = c.opts.spot;
spott = c.opts.spottype;
axlx = c.axlx;
axly = c.axly;

% show spots
if ~isempty(c.spot)
    delete(c.spot);
end
plh = zeros(numel(spot), 1);
spotx = zeros(numel(spot), 2);
for spc = 1:numel(spot)
    spcc = spot{spc}(csel, 1:2);
    spcfx = 0.1 * (axlx(2) - axlx(1));
    spcfy = 0.1 * (axly(2) - axly(1));
    switch (spott{spc})
        case {'c', 'cue'}
            mt = 'c';
            pt = [-0.2, 0.5; 0, 0; 0.2, 0.5];
            lw = 2;
        case {'f', 'free'}
            mt = 'k';
            pt = [0, -0.5; 0, 0.5; -0.3, 0; 0, -0.5; 0.3, 0; 0, 0.5];
            lw = 1;
        case {'max'}
            mt = 'r';
            pt = [-0.3, -0.5; 0.3, -0.5];
            lw = 2;
        case {'min'}
            mt = 'b';
            pt = [-0.3, 0.5; 0.3, 0.5];
            lw = 2;
        case {'minmax'}
            mt = 'g';
            pt = [-0.25, -0.5; 0.25, -0.5; -0.25, 0.5; 0.25, 0.5];
            lw = 1;
        case {'n', 'none'}
            mt = 'k';
            pt = [0, -0.4; -0.2, -0.4; 0.2, -0.4; 0, -0.4; ...
                  0,  0.4; -0.2,  0.4; 0.2,  0.4; 0,  0.4];
            lw = 1;
        case {'o', 'onset'}
            mt = 'm';
            pt = [0, -0.3; 0.3, 0; 0, 0.3];
            lw = 2;
    end
    xpt = [[0, 0]; pt; [0, 0]];
    xpt = [spcc(1) + spcfx .* xpt(:, 1), spcc(2) + spcfy .* xpt(:, 2)];
    plh(spc) = plot(c.ax, xpt(:, 1), xpt(:, 2), [mt '-']);
    set(plh(spc), 'LineWidth', lw);
    spotx(spc, :) = xpt(1, :);
end
ne_ui.plotcurves.spot = plh;
ne_ui.plotcurves.spotx = spotx;


% use curve checkbox
function pc_toggleusecurve(varargin)
global ne_ui;
c = ne_ui.plotcurves;

% which curve
cselr = c.hTag.LB_plotcurves_curves.Value;
if isempty(cselr) || ...
    cselr <= size(c.opts.sets, 1)
    return;
end
csel = cselr - size(c.opts.sets, 1);

% update value
ne_ui.plotcurves.sel(csel) = (c.hTag.CB_plotcurves_usecurve.Value) > 0;

% update string
if ne_ui.plotcurves.sel(csel)
    c.hTag.LB_plotcurves_curves.String{cselr}(2) = '+';
else
    c.hTag.LB_plotcurves_curves.String{cselr}(2) = '-';
end


% compute within-axes position of cursor
function cp = pc_axescp(cp)
global ne_ui;
c = ne_ui.plotcurves;
xl = get(c.ax, 'XLim');
yl = get(c.ax, 'YLim');
cp = cp ./ c.axpos(3:4);
xd = xl(2) - xl(1);
yd = yl(2) - yl(1);
cp = [xl(1) + xd * cp(1), yl(1) + yd * cp(2)];

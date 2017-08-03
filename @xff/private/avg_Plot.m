function [ax, acrvs] = avg_Plot(xo, aplot, splot, ax, vcol, opt)
% AVG::Plot  - plot averaged data according to AVG spec
%
% FORMAT:       [ax, acrvs] = avg.Plot(aplot [, splot, ax, vcol, opt])
%
% Input fields:
%
%       aplot       average data obtained by AVG::Average method
%       splot       optional standard error bars data
%       ax          if given and a valid axes object, use this
%                   if not given, a new figure/axes is created
%       vcol        colors used for multiple VOIs/voxels
%       opt         optional arguments as struct fields
%        .adddata   DxN double array (will be zero padded)
%        .addedcol  Nx3 RGB colors for added data (dim(N) must match)
%        .evtdur    event duration (default: from 1st curve in AVG)
%        .figsize   if given and 1x2 vector, resize figure of axes
%        .linewidth line width of plot lines
%        .quality   quality (only used for JPG snapshots)
%        .snapshot  if given and a filename, write a snapshot of figure
%        .snpclose  close after snapshot
%        .title     chart title
%        .xlabel    label for x axis (default: 'Time (scans)')
%        .ylabel    label for y axis (default: 'fMRI response')
%
% Output fields:
%
%       ax          axes handle where plot has been done
%       acrvs       curves' handles

% Version:  v1.1
% Build:    16041910
% Date:     Apr-19 2016, 10:04 AM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'avg') || ...
   ~isa(aplot, 'double') || any(isinf(aplot(:)))
    error('neuroelf:xff:BadArguments', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if size(aplot, 1) ~= bc.NrOfTimePoints || size(aplot, 3) ~= numel(bc.Curve)
    error('neuroelf:xff:badArgument', 'Invalid aplot argument size.');
end
if nargin > 2 && isa(splot, 'double') && numel(size(splot)) == numel(size(aplot)) && ...
    all(size(splot) == size(aplot)) && ~any(isinf(splot(:)))
    usesplot = true;
else
    usesplot = false;
end
asize = size(aplot);
if numel(asize) < 3
    asize(3) = 1;
end
useax = false;
if nargin > 3 && (isa(ax, 'double') || isa(ax, 'matlab.graphics.axis.Axes')) && numel(ax) == 1
    try
        axp = get(ax);
        if isstruct(axp) && numel(axp) == 1 && isfield(axp, 'Type') && ...
            ischar(axp.Type) && strcmpi(axp.Type, 'axes')
            useax = true;
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end
usevcol = false;
if nargin > 4 && isa(vcol, 'double') && size(vcol, 1) == asize(2) && ...
    size(vcol, 2) == 3 && numel(size(vcol)) == 2 && ...
   ~any(isinf(vcol(:)) | isnan(vcol(:)) | vcol(:) < 0 | vcol(:) > 255)
    usevcol = true;
    if any(vcol(:) > 1)
        vcol = vcol / 255;
    end
end
if nargin < 6 || ~isstruct(opt) || numel(opt) ~= 1
    opt = struct;
end

% create a figure if needed
if ~useax
    figure;
    ax = axes;
end

% make some initial settings to axes object
set(ax, 'Box',    'on');
set(ax, 'Color',  [0, 0, 0]);
set(ax, 'XColor', [1, 1, 1]);
set(ax, 'YColor', [1, 1, 1]);
try
	fh = get(ax, 'Parent');
    set(fh, 'Color', [0, 0, 0], 'NumberTitle', 'off', ...
        'Units', 'pixels', 'Position', [40, 60, 560, 360], ...
        'Name', sprintf('Event-related averaging plot: %s', ...
        aft_FilenameOnDisk(xo)));
catch xfferror
    neuroelf_lasterr(xfferror);
end
hold(ax, 'on');
set(ax, 'Units', 'normalized');
set(ax, 'Position', [0.1, 0.1, 0.85, 0.85]);

% get AVG settings for data
lpre = -bc.PreInterval;
lpos = bc.PostInterval;
if isfield(opt, 'evtdur') && isa(opt.evtdur, 'double') && numel(opt.evtdur) == 1 && ...
   ~isinf(opt.evtdur) && ~isnan(opt.evtdur) && opt.evtdur > 0 && opt.evtdur < lpos
    edur = opt.evtdur;
else
    edur = bc.Curve(1).EventDuration;
end

% prepare aplot data
aplot = reshape(aplot, [asize(1), prod(asize(2:end))]);
pllen = asize(1);

% set more properties on plot
if usesplot
    splot = reshape(splot, size(aplot));
    smin = aplot - splot;
    smax = aplot + splot;
    smin = floor(min(smin(:)));
    smax = ceil(max(smax(:)));
else
    smin = floor(min(aplot(:)));
    smax = ceil(max(aplot(:)));
end

% plot standard error if requested
if usesplot
    scol = cell(1, size(aplot, 2));
    for ac = 1:numel(scol)
        scol{ac} = zeros(asize(1), 1);
        for pc = 1:asize(1)
            xpos = lpre + pc - 1;
            stde = [-splot(pc, ac), splot(pc, ac)];
            lh = line([xpos, xpos], aplot(pc, ac) + stde, 'Parent', ax);
            scol{ac}(pc) = lh;
            set(lh, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 12);
        end
    end
end

% plot curves on top of standard errors
acrvs = plot((lpre:lpos)', aplot);

% added data?
if isfield(opt, 'adddata') && isa(opt.adddata, 'double') && ~isempty(opt.adddata) && ...
    numel(size(opt.adddata)) == 2

    % get data to be added
    adddata = opt.adddata;
    adsz = size(adddata);
    if adsz(1) < pllen
        adddata = [adddata; zeros(adsz(1) - pllen, adsz(2))];
    elseif adsz(1) > pllen
        adddata = adddata(1:pllen, :);
    end

    % plot added data
    addcrvs = plot((lpre:lpos)', adddata);

    % coloring?
    if isfield(opt, 'addedcol') && isa(opt.addedcol, 'double') && ...
        size(opt.addedcol, 1) == adsz(2) && size(opt.addedcol, 2) == 3 && ...
       ~any(isinf(opt.addedcol(:)) | isnan(opt.addedcol(:)) | ...
            opt.addedcol(:) < 0 | opt.addedcol(:) > 255.5)
        addedcol = opt.addedcol;
        if any(addedcol(:) > 1)
            addedcol = round(addedcol) / 255;
        end
        for cc = 1:adsz(2)
            set(addcrvs(cc), 'Color', addedcol(cc, :));
        end
    end
end

% add lines for event
l1 = line([0, 0], [smin, smax], 'Parent', ax);
l2 = line([edur, edur], [smin, smax], 'Parent', ax);
set([l1; l2], 'Color', [1, 1, 1]);

% extend axes *a bit*
set(ax, 'XLim', [lpre - 0.5, lpos + 0.5]);
set(ax, 'YLim', [smin, smax]);
set(ax, 'YTick', smin:smax);
set(ax, 'YGrid', 'on');

% set texts
if isfield(opt, 'xlabel') && ischar(opt.xlabel)
    set(get(ax, 'XLabel'), 'String', opt.xlabel(:)');
else
    set(get(ax, 'XLabel'), 'String', 'Time (scans)');
end
if isfield(opt, 'ylabel') && ischar(opt.ylabel)
    set(get(ax, 'YLabel'), 'String', opt.ylabel(:)');
else
    set(get(ax, 'YLabel'), 'String', 'fMRI response');
end
if isfield(opt, 'title') && ischar(opt.title)
    set(get(ax, 'Title'), 'String', opt.title(:)');
else
    set(get(ax, 'Title'), 'String', 'Event-related averaging plot');
end
set([get(ax, 'XLabel'); get(ax, 'YLabel'); get(ax, 'Title')], ...
    'Color', [1, 1, 1], 'FontWeight', 'bold');

% set width property on curves
if isfield(opt, 'linewidth') && isa(opt.linewidth, 'double') && numel(opt.linewidth) == 1
    set(acrvs, 'LineWidth', opt.linewidth);
else
    set(acrvs, 'LineWidth', 3);
end

% set color for single curve
if asize(3) == 1

    % set vcolor if given
    if usevcol
        for vc = 1:asize(2)
            set(acrvs(vc), 'Color', vcol(vc, :));
            if usesplot
                set(scol{vc}, 'Color', vcol(vc, :));
            end
        end
    end

% set color for multiple curves
else

    % iterate over curves
    for cc = 1:asize(3)

        % single voi
        if asize(2) == 1
            set(acrvs(cc), 'Color', bc.Curve(cc).TimeCourseColor1 / 255);

        % multiple vois
        else
            as2 = asize(2);
            for vc = 1:as2
                set(acrvs(vc + (cc - 1) * as2), 'Color', vcol(vc, :) / 2 + ...
                    bc.Curve(cc).TimeCourseColor1 / 510);
            end
        end
    end
end

% resize figure
if isfield(opt, 'figsize') && isa(opt.figsize, 'double') && numel(opt.figsize) == 2 && ...
   ~any(isinf(opt.figsize) | isnan(opt.figsize) | opt.figsize < 64 | ...
        opt.figsize > 3072)
    try
        ssize = get(0, 'ScreenSize');
        set(fh, 'Units', 'pixels');
        fpos = get(fh, 'Position');
        mpos = min(ssize(3:4), fpos(1:2) + opt.figsize(:)') - [16, 80];
        set(fh, 'Position', [mpos - opt.figsize(:)', opt.figsize(:)']);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

% snapshot
if isfield(opt, 'snapshot') && ischar(opt.snapshot) && ~isempty(opt.snapshot)
    snapshot = opt.snapshot(:)';
    [spath{1:3}] = fileparts(snapshot);
    pause(0.1);
    refresh(fh);
    pause(0.1);
    fr = getframe(fh);
    switch (lower(spath{3}))

        % use defaults
        case {'.bmp', '.gif', '.png'}
            try
                imwrite(fr.cdata, snapshot);
            catch xfferror
                warning('neuroelf:xff:callError', ...
                    'Error calling imwrite for snapshot: %s.', xfferror.message);
            end
        case '.jpg'
            if isfield(opt, 'quality') && isa(opt.quality, 'double') && ...
                numel(opt.quality) == 1 && ~isnan(opt.quality) && ...
                opt.quality > 0 && opt.quality <= 100
                jpgqual = opt.quality;
            else
                jpgqual = 90;
            end
            try
                imwrite(fr.cdata, snapshot, 'jpg', 'Quality', jpgqual);
            catch xfferror
                warning('neuroelf:xff:CallError', ...
                    'Error calling imwrite for snapshot: %s.', xfferror.message);
            end
        otherwise
            warning('neuroelf:xff:badArgument', ...
                'Invalid image format selected for snapshot.');
    end
    if isfield(opt, 'snpclose') && ~isempty(opt.snpclose) && opt.snpclose(1)
        delete(fh);
        pause(0.1);
    end
end

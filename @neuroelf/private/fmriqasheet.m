function [qs, qsc, qstxt] = fmriqasheet(q, odt)
% fmriqasheet  - create figure with output of fmriquality
%
% FORMAT:       [qs, qsc, qstxt = ] fmriqasheet(q [, odt])
%
% Input fields:
%
%       q           return structure from fmriquality call
%       odt         outlier detection threshold (nr of criteria, default: 3)
%
% Output fields:
%
%       qs          figure handle (1x1 double)
%       qsc         figure content (as image)
%       qstxt       quality sheet output as text structure
%        .gfgTC     global time course (foreground voxels)
%        .mpars     motion parameters
%        .sSNR      spatial SNR
%        .tSNR      temporal SNR

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:19 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2015, 2016, Jochen Weber
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
if nargin > 0 && ...
    ischar(q) && ...
    numel(q) == size(q, 2)
    if any(strcmpi(q, {'*.*', '*.fq'}))
        [qfile, qpath] = uigetfile( ...
            {'*.fq', 'fMRIQuality files (*.fq)'}, ...
            'Please select an fMRIquality file...', 'MultiSelect', 'off');
        if isequal(qfile, 0) || ...
            isequal(qpath, 0) || ...
            isempty(qfile)
            qs = [];
            qsc = [];
            return;
        end
        if isempty(qpath)
            qpath = pwd;
        end
        q = [qpath '/' qfile];
    end
    if exist(q, 'file') == 2
        try
            q = load(q, '-mat');
        catch ne_eo;
             neuroelf_lasterr(ne_eo);
        end
    else
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid filename or pattern provided.' ...
        );
    end
end
if nargin > 0 && ...
    isstruct(q) && ...
    numel(q) == 1 && ...
    numel(fieldnames(q)) == 1
    if isfield(q, 'q')
        q = q.q;
    elseif isfield(q, 'fq')
        q = q.fq;
    end
end
if nargin < 1 || ...
    numel(q) ~= 1 || ...
   ~isstruct(q) || ...
   ~isfield(q, 'Dims') || ...
   ~isfield(q, 'Filename') || ...
   ~isfield(q, 'Masks') || ...
   ~isfield(q, 'Raw') || ...
   ~isfield(q, 'Quality')
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument to fmriqasheet.' ...
    );
end
if nargin < 2 || ...
   ~isa(odt, 'double') || ...
    numel(odt) ~= 1 || ...
    isinf(odt) || ...
    isnan(odt) || ...
    odt < 1 || ...
    odt > 4
    odt = 3;
else
    odt = floor(odt);
end

% motion detection run?
if isfield(q, 'MotCorr') && ...
    isstruct(q.MotCorr) && ...
    isfield(q.MotCorr, 'Params')
    m = true;
else
    m = false;
end

% create figure
qs = figure;
qsc = [];
qstxt = struct( ...
    'gfgTC',   '', ...
    'gfgslTC', '', ...
    'mpars',   '', ...
    'imHist',  [], ...
    'imMPar',  [], ...
    'imQPar',  [], ...
    'imsSNR',  [], ...
    'imslTC',  [], ...
    'imtSNR',  [], ...
    'maxRot',  '?.?deg ?.?deg ?.?deg', ...
    'maxTrns', '?.??mm ?.??mm ?.??mm', ...
    'qatxt',   {{}}, ...
    'resFFT',  '0.0', ...
    'sSNR',    '0.0', ...
    'tSNR',    '0.0');

% get ROOT object properties
try
    rp = get(0);
    rs = rp.ScreenSize;
    
% if that fails
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    
    % assume LARGE screen size (for high-res graphics)
    rs = [1, 1, 2560, 1440];
end
rc = 0.5 .* rs(3:4);
if (rs(3) / rs(4)) >= sqrt(0.5)
    rs = 2 * round(0.45 * [sqrt(0.5) * rs(4), rs(4)]);
else
    rs = 2 * round(0.45 * [rs(3), sqrt(0.5) * rs(3)]);
end

% make figure settings
set(qs, ...
    'NumberTitle', 'off', ...
    'Name', ['fMRI data quality sheet: ' q.Filename], ...
    'Position', [rc - 0.5 * rs, rs]);
figure(qs);

% add subplots for output
cols = 2;
targets = [1:2:6, 2:2:6];

% get some numbers
nvol = size(q.TC.Slices, 1);
nslc = size(q.TC.Slices, 2);
volsize = q.Dims(1:3);
fgvox = sum(q.Masks.Foreground(:));
fgslvox = sum(reshape(q.Masks.Foreground, [volsize(1) * volsize(2), volsize(3)]));

% plot PSC time courses or slices with outliers marked
if m
    tcp = subplot(3, cols, targets(1));
else
    tcp = subplot(3, cols, targets([1, 4]));
end
qstxt.gfgTC = sprintf('%.3f\n', ...
    (sum(diag(fgslvox(:)) * q.TC.TF_ForeSlicesWeighted') ./ fgvox)');
qstxt.gfgslTC = sprintf([repmat('%.3f\t', 1, nslc) '\n'], q.TC.TF_ForeSlicesWeighted');
tcdata = psctrans(q.TC.TF_ForeSlicesWeighted);
tcdata(:, all(tcdata == 0)) = 100;
tcdata = tcdata - 100;
plot(2 * repmat(0:(nslc - 1), nvol, 1) + tcdata);
hold(tcp, 'on');
olm = uint8(zeros(2 * nslc, nvol, 3));
olms = olm;
olm(:, :, 1) = 255;
olms(:, :, 3) = 160;
olh = image(olm, 'Parent', tcp);
olhs = image(olms, 'Parent', tcp);
olm = repmat(0.25 * log(1 + q.Quality.Outliers.Volumes'), 2 * nslc, 1);
olms = zeros(2 * nslc, nvol);
set(olh, 'AlphaData', olm, 'YData', [0, 2 * nslc - 1]);
set(olhs, 'AlphaData', olms, 'YData', [-0.5, 2 * nslc - 1.5]);
set(tcp, 'XLim', [0, nvol + 1], 'YLim', [-1, 2 * nslc]);
xlabel(tcp, 'Volumes (TRs)');
ylabel(tcp, 'foreground slices mean signal (% change)');

% show overall global signal-to-noise ratio info
subplot(3, cols, targets(2));
id = repmat(uint8(floor(packmosaic(scaledata(double(q.Quality.GlobalSNRImage))))), [1, 1, 3]);
image(id);
qstxt.imsSNR = id;

% SNR over time
if isfield(q.Masks, 'ForegroundClipped')
    snrimage = repmat(uint8(floor(packmosaic(scaledata( ...
        double(q.Quality.LocalSNRImage) .* double(~q.Masks.ForegroundClipped))))), [1, 1, 3]);
    fclipped = find(lsqueeze(packmosaic(uint8(q.Masks.ForegroundClipped))) > 0);
    if ~isempty(fclipped)
        snrimage(fclipped) = 255;
        snrimage(fclipped + size(snrimage, 1) * size(snrimage, 2)) = 64;
        snrimage(fclipped + 2 * size(snrimage, 1) * size(snrimage, 2)) = 64;
    end
else
    snrimage = repmat(uint8(floor(packmosaic(scaledata(double(q.Quality.LocalSNRImage))))), [1, 1, 3]);
end
subplot(3, cols, targets(5));
image(snrimage);
qstxt.imtSNR = snrimage;

% and histogram of SNR
if m
    ihp = subplot(3, cols, targets(3));
else
    ihp = subplot(3, cols, targets([3, 6]));
end
if isfield(q.Masks, 'ForegroundClipped')
    hist(q.Quality.GlobalSNRImage( ...
        q.Masks.Foreground & (~q.Masks.ForegroundClipped)), 250);
else
    hist(q.Quality.GlobalSNRImage(q.Masks.Foreground), 250);
end
xlabel(ihp, 'Intensity (a.u.)');
ylabel(ihp, 'absolute frequency (number of voxels)');

% detect outliers
dout = find(q.Quality.Outliers.Volumes >= odt);
if isfield(q.TC, 'TF_ForeSlicesNoiseSpike')
    nout = numel(q.TC.TF_ForeSlicesNoiseSpike);
    sout = find(q.TC.TF_ForeSlicesNoiseSpike(:) >= 0.2);
else
    nout = 1;
    sout = [];
end

% print out some info
line = repmat('-', 1, 72);
qstxt.qatxt{end+1} = 'fMRI Quality statistics:';
qstxt.qatxt{end+1} = line;
if isfield(q, 'Res')
    qstxt.qatxt{end+1} = sprintf( ...
        ' - data dimensions:   %d x %d x %d (%.1fmm x %.1fmm x %.1fmm)', ...
        volsize, volsize .* q.Res);
else
    qstxt.qatxt{end+1} = sprintf(' - data dimensions:   %d x %d x %d', volsize);
end
qstxt.qatxt{end+1} = sprintf(' - number of volumes: %d', nvol);
qstxt.qatxt{end+1} = sprintf(' - outlier volumes:   %d (%.1f%%) [ %s]', ...
    sum(q.Quality.Outliers.Volumes >= odt), 100 * q.Quality.Outliers.VolumeRatio, ...
    sprintf('%d ', dout));
qstxt.qatxt{end+1} = sprintf(' - outlier slices:    %d (%.1f%%)', ...
    numel(sout), 100 * numel(sout) / nout);
qstxt.qatxt{end+1} = line;
if isfield(q, 'Raw') && ...
    isfield(q.Raw, 'MeanImage95Pct')
    mi95p = q.Raw.MeanImage95Pct;
else
    mi95p = sort(q.Raw.MeanImage(q.Masks.Foreground));
    mi95p = mi95p([floor(1 + 0.025 * numel(mi95p)), ceil(0.975 * numel(mi95p))]);
end
qstxt.qatxt{end+1} = sprintf( ...
    ' - 95%% central intensity interval:  [%-6.2f .. %-6.2f]', mi95p(1), mi95p(2));
if isfield(q.Quality, 'GlobalSNRMean') && ...
    isfield(q.Quality, 'GlobalSNRMedian')
    gsmean = q.Quality.GlobalSNRMean;
    gsmedian = q.Quality.GlobalSNRMedian;
else
    gsmean = mean(q.Quality.GlobalSNRImage(q.Masks.Foreground));
    gsmedian = median(q.Quality.GlobalSNRImage(q.Masks.Foreground));
end
qstxt.qatxt{end+1} = sprintf( ...
    ' - average / median spatial  SNR:   %-6.2f / %-6.2f', gsmean, gsmedian);
qstxt.sSNR = sprintf('%.2f', gsmean);
if isfield(q.Quality, 'LocalSNRMean') && ...
    isfield(q.Quality, 'LocalSNRMean')
    lsmean = q.Quality.LocalSNRMean;
    lsmedian = q.Quality.LocalSNRMedian;
else
    lsmean = mean(q.Quality.LocalSNRImage(q.Masks.Foreground));
    lsmedian = median(q.Quality.LocalSNRImage(q.Masks.Foreground));
end
qstxt.qatxt{end+1} = sprintf( ...
    ' - average / median temporal SNR:   %-6.2f / %-6.2f', lsmean, lsmedian);
qstxt.tSNR = sprintf('%.2f', lsmean);
if isfield(q.TempFiltered, 'FFTMeanSD') && ...
    isfield(q.TempFiltered, 'FFTMeanSDPowerMeanStd') && ...
   ~isempty(q.TempFiltered.FFTMeanSD)
    qstxt.qatxt{end+1} = sprintf( ...
        ' - average / std residual FFT power:   %-6g / %-6g', ...
        q.TempFiltered.FFTMeanSDPowerMeanStd(1), ...
        q.TempFiltered.FFTMeanSDPowerMeanStd(2));
    qstxt.resFFT = sprintf('%6g', q.TempFiltered.FFTMeanSDPowerMeanStd(1));
end

% motion correction stuff
if m

    % text info
    qstxt.mpars = sprintf([repmat('%.3f\t', 1, 6) '\n'], q.MotCorr.Params');
    maxmot = max(q.MotCorr.Params, [], 1) - min(q.MotCorr.Params, [], 1);
    qstxt.qatxt{end+1} = line;
    qstxt.maxRot = sprintf('%-5.3fdeg %-5.3fdeg %-5.3fdeg', maxmot(4:6));
    qstxt.maxTrns = sprintf('%-5.3fmm  %-5.3fmm  %-5.3fmm', maxmot(1:3));
    qstxt.qatxt{end+1} = sprintf(' - maximal translation: %s', qstxt.maxTrns);
    qstxt.qatxt{end+1} = sprintf(' - maximal rotation:    %s', qstxt.maxRot);

    % parameters
    mpp = subplot(3, cols, targets(4));
    plot(q.MotCorr.Params);
    set(mpp, 'XLim', [0, nvol + 1]);
    xlabel(mpp, 'Volumes (TRs)');
    ylabel(mpp, 'translation (mm) / rotation (degrees)');
    legend(mpp, 'x-trans', 'y-trans', 'z-trans', 'x-yaw', 'y-pitch', 'z-roll', ...
        'Location', 'SouthWest');

    % important time courses
    qtp = subplot(3, cols, targets(6));
    plot(3 * repmat(0:(size(q.TC.Quality, 2) - 1), nvol, 1) + q.TC.Quality);
    set(qtp, 'XLim', [0, nvol + 1], 'YLim', [-10, 14]);
    xlabel(qtp, 'Volumes (TRs)');
    ylabel(qtp, 'quality measures (a.u.)');
    legend(qtp, 'global-TC', 'brain-TC', 'outliers', 'global-filtered', 'smooth-est', ...
        'Location', 'SouthWest');

% otherwise some more stuff
else
    mpp = [];
end

% combine text together
qstxt.qatxt = sprintf('%s\n', qstxt.qatxt{:});
disp(qstxt.qatxt);

% no filename given, we can't go on
if isempty(q.Filename) || ...
    exist(regexprep(q.Filename, '\,\d+$', ''), 'file') ~= 2

    % second output
    if nargout > 1

        % write into tempfile in hi-res
        tempfile = [tempname '.png'];
        set(qs, 'PaperPositionMode', 'auto');
        pause(0.5);
        print(qs, '-r600', '-dpng', tempfile);
        set(qs, 'PaperPositionMode', 'manual');

        % read tempfile and resize to 300DPI
        qsc = imread(tempfile);
        delete(tempfile);
        qsc = image_resize(qsc, ceil(0.5 * size(qsc, 1)), ceil(0.5 * size(qsc, 2)));
        
        % more outputs
        if nargout > 2
            qsz = size(qsc);
            qsz = qsz(1, [2, 1]);
            qsu = get(qs, 'Units');
            set(qs, 'Units', 'normalized');
            axpos = [1, -1, 1, 1] .* get(tcp, 'Position') + [-0.05, 1.02, 0.08, 0.04];
            axpos = 1 + round([qsz .* axpos(1:2), ...
                qsz .* (axpos(1:2) + [axpos(3), -axpos(4)])]);
            qstxt.imslTC = qsc(axpos(4):axpos(2), axpos(1):axpos(3), :);
            axpos = [1, -1, 1, 1] .* get(ihp, 'Position') + [-0.05, 1.02, 0.08, 0.04];
            axpos = 1 + round([qsz .* axpos(1:2), ...
                qsz .* (axpos(1:2) + [axpos(3), -axpos(4)])]);
            qstxt.imHist = qsc(axpos(4):axpos(2), axpos(1):axpos(3), :);
            set(qs, 'Units', qsu);
            if ~isempty(mpp)
                axpos = [1, -1, 1, 1] .* get(mpp, 'Position') + [-0.05, 1.02, 0.08, 0.04];
                axpos = 1 + round([qsz .* axpos(1:2), ...
                    qsz .* (axpos(1:2) + [axpos(3), -axpos(4)])]);
                qstxt.imMPar = qsc(axpos(4):axpos(2), axpos(1):axpos(3), :);
                axpos = [1, -1, 1, 1] .* get(qtp, 'Position') + [-0.05, 1.02, 0.08, 0.04];
                axpos = 1 + round([qsz .* axpos(1:2), ...
                    qsz .* (axpos(1:2) + [axpos(3), -axpos(4)])]);
                qstxt.imQPar = qsc(axpos(4):axpos(2), axpos(1):axpos(3), :);
            end
        end
    end
    return;
end

% try to get a handle to the object data
try
    obj = xff(q.Filename, 't');
    rtv = obj.RunTimeVars;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    return;
end
saved = true;
if ~isfield(rtv, 'Discard') || ...
   ~isnumeric(rtv.Discard) || ...
    isempty(rtv.Discard) || ...
    any(rtv.Discard(:) ~= fix(double(rtv.Discard(:))))
    rtv.Discard = dout;
    if ~isempty(dout)
        saved = false;
    end
else
    rtv.Discard = double(rtv.Discard(:));
end
if ~isfield(rtv, 'SpikeSlices') || ...
   ~islogical(rtv.SpikeSlices) || ...
   ~isequal(size(rtv.SpikeSlices), q.Dims(1, [4, 3]))
    if isfield(q.TC, 'TF_ForeSlicesNoiseSpike')
        rtv.SpikeSlices = (q.TC.TF_ForeSlicesNoiseSpike > 0.2);
    else
        rtv.SpikeSlices = false(q.Dims(1, [4, 3]));
    end
end

% re-set outliers
olm = zeros(1, nvol);
olm(rtv.Discard) = 0.5;
set(olh, 'AlphaData', repmat(olm, 2 * nslc, 1));
olms = 0.5 .* double(rtv.SpikeSlices(:, ceil(0.5:0.5:nslc))');
set(olhs, 'AlphaData', olms);

% set up a few things with the figure
udt = struct( ...
    'btdown', false, ...
    'chndl',  -1, ...
    'dsc',    rtv.Discard, ...
    'mpos',   [-1, -1], ...
    'nslc',   nslc, ...
    'nvol',   nvol, ...
    'object', obj, ...
    'ofile',  obj.FilenameOnDisk, ...
    'ola',    get(olh, 'Parent'), ...
    'olh',    olh, ...
    'olhs',   olhs, ...
    'olm',    olm, ...
    'olms',   olms, ...
    'q',      q, ...
    'set',    true, ...
    'shift',  false, ...
    'rtv',    rtv, ...
    'saved',  saved, ...
    'vmark',  -1);
set(qs, ...
    'CloseRequestFcn',       {@fmriqasheet_close, qs}, ...
    'DeleteFcn',             {@fmriqasheet_delete, qs}, ...
    'HandleVisibility',      'callback', ...
    'Name',                  sprintf( ...
        'fMRI data quality sheet: %s (%d outlier%s)', ...
        udt.ofile, numel(udt.dsc), plurals(numel(udt.dsc))), ...
    'Units',                 'normalized', ...
    'UserData',              udt, ...
    'Visible',               'on', ...
    'WindowButtonDownFcn',   {@fmriqasheet_btdown, qs}, ...
    'WindowButtonMotionFcn', {@fmriqasheet_btmove, qs}, ...
    'WindowButtonUpFcn',     {@fmriqasheet_btup, qs}, ...
    'WindowKeyPressFcn',     {@fmriqasheet_keydown, qs}, ...
    'WindowKeyReleaseFcn',   {@fmriqasheet_keyup, qs});

% second output
if nargout > 1

    % write into tempfile in hi-res
    tempfile = [tempname '.png'];
    set(qs, 'PaperPositionMode', 'auto');
    pause(0.5);
    print(qs, '-r600', '-dpng', tempfile);
    set(qs, 'PaperPositionMode', 'manual');

    % read tempfile and resize to 300DPI
    qsc = imread(tempfile);
    delete(tempfile);
    qsc = image_resize(qsc, ceil(0.5 * size(qsc, 1)), ceil(0.5 * size(qsc, 2)));

    % third output
    if nargout > 2
        qsz = size(qsc);
        qsz = qsz(1, [2, 1]);
        qsu = get(qs, 'Units');
        set(qs, 'Units', 'normalized');
        axpos = [1, -1, 1, 1] .* get(tcp, 'Position') + [-0.05, 1.02, 0.08, 0.04];
        axpos = 1 + round([qsz .* axpos(1:2), ...
            qsz .* (axpos(1:2) + [axpos(3), -axpos(4)])]);
        qstxt.imslTC = qsc(axpos(4):axpos(2), axpos(1):axpos(3), :);
        axpos = [1, -1, 1, 1] .* get(ihp, 'Position') + [-0.05, 1.02, 0.08, 0.04];
        axpos = 1 + round([qsz .* axpos(1:2), ...
            qsz .* (axpos(1:2) + [axpos(3), -axpos(4)])]);
        qstxt.imHist = qsc(axpos(4):axpos(2), axpos(1):axpos(3), :);
        set(qs, 'Units', qsu);
        if ~isempty(mpp)
            axpos = [1, -1, 1, 1] .* get(mpp, 'Position') + [-0.05, 1.02, 0.08, 0.04];
            axpos = 1 + round([qsz .* axpos(1:2), ...
                qsz .* (axpos(1:2) + [axpos(3), -axpos(4)])]);
            qstxt.imMPar = qsc(axpos(4):axpos(2), axpos(1):axpos(3), :);
            axpos = [1, -1, 1, 1] .* get(qtp, 'Position') + [-0.05, 1.02, 0.08, 0.04];
            axpos = 1 + round([qsz .* axpos(1:2), ...
                qsz .* (axpos(1:2) + [axpos(3), -axpos(4)])]);
            qstxt.imQPar = qsc(axpos(4):axpos(2), axpos(1):axpos(3), :);
        end
    end
end



% UI functions



% close request
function fmriqasheet_close(varargin)

% get figure, UserData
f = varargin{3};
udt = get(f, 'UserData');

% not saved
if ~udt.saved

    % request saving
    savertv = questdlg(sprintf('Update RunTimeVars.Discard field in file %s?', ...
        udt.ofile), 'NeuroElf - user request', 'No', 'Yes', 'Yes');
    if ischar(savertv) && ...
        strcmpi(savertv, 'yes')
        try
            udt.object.RunTimeVars = udt.rtv;
            udt.object.SaveRunTimeVars;
        catch ne_eo;
            uiwait(warndlg(['Error updating RunTimeVars: ' ne_eo.message], ...
                'NeuroElf - error', 'modal'));
            neuroelf_lasterr(ne_eo);
        end
    end

    % update udt
    udt.saved = true;
    set(f, 'UserData', udt);
end

% delete
if ishandle(udt.chndl)
    delete(udt.chndl);
end
delete(f);


% figure deletion
function fmriqasheet_delete(varargin)

% get figure, UserData
f = varargin{3};
udt = get(f, 'UserData');

% not saved
if ~udt.saved

    % request saving
    savertv = questdlg(sprintf('Update RunTimeVars.Discard field in file %s?', ...
        udt.ofile), 'NeuroElf - user request', 'No', 'Yes', 'Yes');
    if ischar(savertv) && ...
        strcmpi(savertv, 'yes')
        try
            udt.object.RunTimeVars = udt.rtv;
            udt.object.SaveRunTimeVars;
        catch ne_eo;
            uiwait(warndlg(['Error updating RunTimeVars: ' ne_eo.message], ...
                'NeuroElf - error', 'modal'));
            neuroelf_lasterr(ne_eo);
        end
    end

    % update udt
    udt.saved = true;
    set(f, 'UserData', udt);
end

% clear object
clearxffobjects({udt.object});


% button down
function fmriqasheet_btdown(varargin)

% get figure, UserData
f = varargin{3};
udt = get(f, 'UserData');

% get position
pos = get(f, 'CurrentPoint');

% compare to important axes
olpos = get(udt.ola, 'Position');
olpos(5:6) = olpos(1:2) + olpos(3:4);

% click not within axes
if pos(1) < olpos(1) || ...
    pos(1) > olpos(5) || ...
    pos(2) < olpos(2) || ...
    pos(2) > olpos(6)

    % return early
    return;
end

% compute X/Y position
oll = get(udt.ola, 'XLim');
oly = get(udt.ola, 'YLim');
vmark = max(1, min(udt.nvol, ...
    round((oll(2) - oll(1)) * (pos(1) - olpos(1)) / olpos(3) + oll(1))));
smark = max(1, min(udt.nslc, ...
    round(0.5 * ((oly(2) - oly(1)) * (pos(2) - olpos(2)) / olpos(4) + oly(1) + 2))));

% if shift is not pressed
if isfield(udt.q.TempFiltered, 'Data') && ...
   ~isempty(udt.q.TempFiltered.Data) && ...
   ~udt.shift
    if udt.chndl > 0 && ...
        ishandle(udt.chndl)
        fmriqasheet_showslice(f, vmark, smark, udt.chndl);
    else
        udt.chndl = fmriqasheet_showslice(f, vmark, smark);
    end
    set(f, 'UserData', udt);
    return;
end

% store current settings
udt.btdown = true;
udt.dsc = udt.rtv.Discard;
udt.mpos = pos;

% determine which volume to mark also
udt.vmark = vmark;

% set in discarded volumes for display
if any(udt.dsc == vmark)
    udt.dsc(udt.dsc == vmark) = [];
    udt.set = false;
else
    udt.dsc = sort([udt.dsc(:); vmark]);
    udt.set = true;
end

% update image
olm = zeros(1, numel(udt.olm));
olm(udt.dsc) = 0.5;
set(udt.olh, 'AlphaData', repmat(olm, 2 * udt.nslc, 1));

% update in UserData
set(f, 'UserData', udt, 'Name', sprintf( ...
    'fMRI data quality sheet: %s (%d outlier%s)', ...
    udt.ofile, numel(udt.dsc), plurals(numel(udt.dsc))));


% mouse move
function fmriqasheet_btmove(varargin)

% get figure, UserData
f = varargin{3};
udt = get(f, 'UserData');

% if button not pressed, return
if ~udt.btdown;
    return;
end

% get position
pos = get(f, 'CurrentPoint');

% determine which volume to mark also
oll = get(udt.ola, 'XLim');
olpos = get(udt.ola, 'Position');
vmark = max(1, min(udt.nvol, ...
    round((oll(2) - oll(1)) * (pos(1) - olpos(1)) / olpos(3) + oll(1))));
if vmark < udt.vmark
    vmark = vmark:udt.vmark;
else
    vmark = udt.vmark:vmark;
end

% set or unset
if udt.set
    udt.dsc = union(udt.rtv.Discard(:), vmark(:));
else
    udt.dsc = setdiff(udt.rtv.Discard(:), vmark(:));
end

% update image
olm = zeros(1, numel(udt.olm));
olm(udt.dsc) = 0.5;
set(udt.olh, 'AlphaData', repmat(olm, 2 * udt.nslc, 1));

% update in UserData
set(f, 'UserData', udt, 'Name', sprintf( ...
    'fMRI data quality sheet: %s (%d outlier%s)', ...
    udt.ofile, numel(udt.dsc), plurals(numel(udt.dsc))));


% button up
function fmriqasheet_btup(varargin)

% get figure, UserData
f = varargin{3};
udt = get(f, 'UserData');

% update UserData
udt.btdown = false;
udt.rtv.Discard = udt.dsc;
udt.mpos = [-1, -1];
udt.saved = false;
set(f, 'UserData', udt, 'Name', sprintf( ...
    'fMRI data quality sheet: %s (%d outlier%s)', ...
    udt.ofile, numel(udt.dsc), plurals(numel(udt.dsc))));


% key down
function fmriqasheet_keydown(src, ke, varargin)

% get UserData
udt = get(src, 'UserData');
if ~isstruct(udt) || ...
   ~isfield(udt, 'shift')
    return;
end

% shift pressed
if ~isempty(ke.Modifier) && ...
    any(strcmpi(ke.Modifier, 'shift'))
    udt.shift = true;
else
    udt.shift = false;
end

% set UserData
set(src, 'UserData', udt);


% key up
function fmriqasheet_keyup(src, ke, varargin)

% get UserData
udt = get(src, 'UserData');
if ~isstruct(udt) || ...
   ~isfield(udt, 'shift')
    return;
end

% shift pressed
if ~isempty(ke.Modifier) && ...
    any(strcmpi(ke.Modifier, 'shift'))
    udt.shift = true;
else
    udt.shift = false;
end

% set UserData
set(src, 'UserData', udt);


% show slice image
function chndl = fmriqasheet_showslice(pf, vol, slice, chin)

% no handle given
if nargin < 4 || ...
   ~ishandle(chin)
    chndl = figure;
    figure(chndl);
    ax = axes;
else
    chndl = chin;
    ax = get(chndl, 'Children');
end

% get UserData
udt = get(pf, 'UserData');

% create image
simage = uint8(0);
sdim = udt.q.Dims;
simage(3 * sdim(1), 3 * sdim(2)) = 0;

% compute mean over data
fdata = udt.q.TempFiltered.Data(:, :, slice, :);
mdata = mean(fdata, 4);
fmask = udt.q.Masks.Foreground(:, :, slice);

% set scaled data
simage(1:sdim(1), 1:sdim(2)) = fmask .* floor(scaledata(fdata(:, :, 1, vol - 1)));
simage(1:sdim(1), sdim(2)+1:2*sdim(2)) = fmask .* floor(scaledata(fdata(:, :, 1, vol)));
simage(1:sdim(1), 2*sdim(2)+1:3*sdim(2)) = fmask .* floor(scaledata(fdata(:, :, 1, vol + 1)));

% remove mean
fdata = fdata - repmat(mdata, [1, 1, 1, size(fdata, 4)]);

% set more data
mm = 127.5 / max(abs(lsqueeze(fdata(:, :, 1, vol-1:vol+1))));
simage(sdim(1)+1:2*sdim(1), 1:sdim(2)) = ...
    floor(127.75 + mm .* fdata(:, :, 1, vol - 1));
simage(sdim(1)+1:2*sdim(1), sdim(2)+1:2*sdim(2)) = ...
    floor(127.75 + mm .* fdata(:, :, 1, vol));
simage(sdim(1)+1:2*sdim(1), 2*sdim(2)+1:3*sdim(2)) = ...
    floor(127.75 + mm .* fdata(:, :, 1, vol + 1));

% compute FFT
fdata = squeeze(double(fdata(:, :, :, vol-1:vol+1)));
for vc = 1:3
    fdata(:, :, vc) = abs(fftshift(fftn(fdata(:, :, vc))));
end

% set more data
mm = 255.999 / max(lsqueeze(fdata));
simage(2*sdim(1)+1:3*sdim(1), 1:sdim(2)) = floor(mm .* fdata(:, :, 1));
simage(2*sdim(1)+1:3*sdim(1), sdim(2)+1:2*sdim(2)) = floor(mm .* fdata(:, :, 2));
simage(2*sdim(1)+1:3*sdim(1), 2*sdim(2)+1:3*sdim(2)) = floor(mm .* fdata(:, :, 3));

% set image
image(repmat(simage, [1, 1, 3]), 'Parent', ax);
set(ax, 'YDir', 'reverse', ...
    'XLim', [0.5, size(simage, 2) + 0.5], 'YLim', [0.5, size(simage, 1) + 0.5]);

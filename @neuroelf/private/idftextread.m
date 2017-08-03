function [eyepts, timepts, idfcfg] = idftextread(idffile, options)
% idftextread  - read a text-converted IDF file (iView X)
%
% FORMAT:       [eyepts, timepts, idfcfg] = idftextread(idffile, options)
%
% Input fields:
%
%       idffile     filename of input file
%       options     1x1 struct with optional fields
%        .bfilter   1x1, short for .lfilter = binomial(N,[0:N])./(2^N)
%        .cal       Nx8 matrix for calibration matrix, order:
%                   screen left X/Y, right X/Y, file left X/Y, right X/Y
%        .cols      1x4 double array, columns to take as eye positions
%                   (default: [4, 5, 6, 7])
%        .deblink   false|true, try to remove eye blinks (default: true)
%        .detrend   1x1 double, 0: don't detrend, 1: de-mean, 2: de-linear,
%                   3...N: de-spline with wavelets of length L
%        .freq      1x1, data frequency (default: auto detect)
%        .lfilter   1xN linear filter weights (e.g. [0.25 0.5 0.25],
%                   default: [1])
%        .range     1x2 double, seconds to select from input
%        .smooth    1x3 double, [passes, min diff to use, smoothing factor]
%        .timeunit  1x1 double, factor to calculate seconds from time
%                   (default: 1e-6)
%        .translat  false|true, try value translation into degrees
%                   (only valid if unit is NOT rawpix; default: false)
%        .unit      unit in which eye positions are read, one out of
%                   'degrees', 'mm', 'radiens', 'rawpix', 'screenpix'
%                   (default: 'rawpix')
%
% Output fields:
%
%       eyepts      Nx4 double list with coordinates, order:
%                   left X/Y, right X/Y
%       timepts     Nx1 double list with time indices from file
%       idfcfg      1x1 struct with parameters from IDF text file
%
% See also: asciiread

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 2:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
   ~ischar(idffile) || ...
    isempty(idffile) || ...
    exist(idffile(:)', 'file') ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Missing or invalid idffile argument.' ...
    );
end

% startup
disp('idftextread - IDF text file parser');
idfcfg  = struct;

% parse options
defopt = { ...
    'bfilter',  'double',       'nonempty',     0; ...
    'cal',      'double',       'noinfnan',     []; ...
    'cols',     'double',       'nonempty',     4:7; ...
    'deblink',  'logical',      'nonempty',     true; ...
    'detrend',  'double',       'nonempty',     0; ...
    'freq',     'double',       'nonempty',     0; ...
    'lfilter',  'double',       'nonempty',     1; ...
    'range',    'double',       'nonempty',     [0, Inf]; ...
    'smooth',   'double',       'nonempty',     [0, 0, 0]; ...
    'timeunit', 'double',       'nonempty',     1e-6; ...
    'translat', 'logical',      'nonempty',     false; ...
    'unit',     'char',         'nonempty',     'rawpix' ...
};
if nargin < 2 || ...
  ~isstruct(options) || ...
   numel(options) ~= 1
    options = checkstruct(struct, defopt);
else
    options = checkstruct(options, defopt);
end
defopt = checkstruct(struct, defopt);

bfilter = options.bfilter(1);
if bfilter >= 1 && ...
   ~isnan(bfilter) && ...
   ~isinf(bfilter)
    bfilter = fix(bfilter);
else
    bfilter = defopt.bfilter;
end
if bfilter > 0
    options.lfilter = binomial(bfilter, 0:bfilter) ./ (2^bfilter);
end

cols = fix(options.cols(:)');
if length(cols) < 4 || ...
    length(cols) > 8 || ...
   any(cols < 4 | cols > 24 | isnan(cols) | isinf(cols))
    cols = defopt.cols;
end
cols = cols - 3;

deblink = options.deblink(1);

detrendv = fix(options.detrend(1));
if isnan(detrendv) || ...
    isinf(detrendv) || ...
    detrendv < 0
    detrendv = defopt.detrend;
end

lfilter = options.lfilter(:)';
if any(isinf(lfilter) | isnan(lfilter) | lfilter < 0)
    lfilter = defopt.lfilter;
else
    lfilter = lfilter ./ sum(lfilter);
    lfilter = lfilter(lfilter >= 0.001);
end
mlf = round((length(lfilter)+1)/2);

freq = options.freq(1);
if isnan(freq) || ...
    isinf(freq) || ...
    freq < 1
    freq = defopt.freq;
end

range = options.range(:)';
if length(range) ~= 2 || ...
    any(range < 0 | isnan(range) | isinf(range))
    range = defopt.range;
end

smoothv = options.smooth(:)';
if length(smoothv) < 3 || ...
    any(smoothv < 0 | isnan(smoothv) | isinf(smoothv))
    smoothv = defopt.smooth;
else
    smoothv = smoothv(1:3);
    smoothv(1) = fix(smoothv(1));
    if smoothv(1) > 3
        smoothv(1) = 3;
    end
    if smoothv(2) > 0.1
        smoothv(2) = 0.1;
    end
    if smoothv(3) > 1
        smoothv(3) = 1;
    end
end

timeunit = options.timeunit(1);
if isnan(timeunit) || ...
    isinf(timeunit) || ...
    timeunit < 1e-9 || ...
    timeunit > 1
    timeunit = defopt.timeunit;
end

translat = options.translat(1);

unit = lower(options.unit(:)');
if ~any(strcmp({'degrees', 'mm', 'radiens', 'rawpix', 'screenpix'}, unit))
    unit = defopt.unit;
end

clear defopt;

% load file
disp([' -> reading IDF input file [' idffile ']...']);
[idf, starts, ends, lnum] = asciiread(idffile);

% reading settings
disp(' -> parsing file header...');
lc = 1;
tline = textline(idf, starts, ends, lc);
while ~isempty(tline) && ...
    tline(1) == '#'
    lopt = tline(4:end);
    if any(lopt == ':')
        loptcol = find(lopt == ':');
        loptname = upper(makelabel(lopt(1:loptcol(1)-1)));
        try
            idfcfg.(loptname) = lopt(loptcol(1)+2:end);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
    lc = lc + 1;
    tline = textline(idf, starts, ends, lc);
end
clear lopt loptcol loptname;

for tfield = { ...
    'SAMPLE_RATE', ...
    'CALIBRATION_AREA', ...
    'STIMULUS_DIMENSION_MM_', ...
    'HEAD_DISTANCE_MM_', ...
    'NUMBER_OF_SAMPLES' ...
}
    if isfield(idfcfg, tfield{1})
        try
            idfcfg.(tfield{1}) = eval(idfcfg.(tfield{1}));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
end
clear tfield;

if isfield(idfcfg, 'CALIBRATION_AREA')
    if isfield(idfcfg, 'STIMULUS_DIMENSION_MM_')
        idfcfg.MM_PER_PIX_ = idfcfg.STIMULUS_DIMENSION_MM_ ./ idfcfg.CALIBRATION_AREA;
        idfcfg.PIX_PER_MM_ = idfcfg.CALIBRATION_AREA ./ idfcfg.STIMULUS_DIMENSION_MM_;
    end
end

if isfield(idfcfg, 'CALIBRATION_TYPE')
    try
        cpnum = regexp(idfcfg.CALIBRATION_TYPE, '^(\d+)\-.*$', 'tokens');
        if ~isempty(cpnum)
            idfcfg.CALIBRATION_TYPE = str2double(cpnum{1}{1});
        else
            idfcfg = rmfield(idfcfg, 'CALIBRATION_TYPE');
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        idfcfg = rmfield(idfcfg, 'CALIBRATION_TYPE');
    end
    clear cpnum;
end

numpat  = '([\-\+]?\d*\.?\d+[de]?\d*)';
if isfield(idfcfg, 'CALIBRATION_TYPE') && ...
   isfield(idfcfg, 'CALIBRATION_TARGETS')
    cpnum = idfcfg.CALIBRATION_TYPE;
    cpts = splittocell(lower(idfcfg.CALIBRATION_TARGETS), 'point');
    if ~isempty(cpts) && ...
        isempty(cpts{1})
        cpts(1) = [];
    end
    if ~isempty(cpts) && ...
       ~isempty(cpts{1}) && ...
        cpnum <= length(cpts)
        disp('    -> detecting calibration points...');
        disp('       p#:   Cam-X   Cam-Y   Scr-X   Scr-Y');
        ncpts   = zeros(cpnum, 6);
        numpair = ['\(' numpat ',\s*' numpat '\)'];
        tnumpat = ['^\d+\:\s+raw' numpair '\s*screen' numpair '\s*coefficient' numpair '\s*$'];
        for cc = 1:cpnum
            cptm = regexp(cpts{cc}, tnumpat, 'tokens');
            if ~isempty(cptm) && ...
                length(cptm{1}) == 6
                try
                    ncpts(cc, :) = eval([cptm{1}{1} ',' cptm{1}{2} ',' ...
                        cptm{1}{3} ',' cptm{1}{4} ',' cptm{1}{5} ',' cptm{1}{6}]);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
                disp(sprintf('       %02d: %7.2f %7.2f   %5d   %5d', [cc,ncpts(cc, 1:4)]));
            end
        end
        gcpts = ncpts(ncpts(:, 1) ~= 0, :);
        [cr, cp] = corrcoef(gcpts);
        disp(sprintf('       -> X-correlation: r=%7.5f, p=%7.5f', cr(3,1), cp(3,1)));
        disp(sprintf('       -> Y-correlation: r=%7.5f, p=%7.5f', cr(4,2), cp(4,2)));
        idfcfg.CALIBRATION_POINTS = ncpts;
        if any(size(gcpts)~=size(ncpts))
            disp(sprintf('          (only based on %d points!)', size(gcpts, 1)));
            idfcfg.CALIBRATION_GOODPOINTS = gcpts;
        else
            gcpts = ncpts;
        end
        cc = [cr(3,1), cr(4,2); cp(3,1), cp(4,2)];
        idfcfg.CALIBRATION_CORRCOEF = cc;
        if all(cc(1,:) > 0.9) && ...
            all(cc(2,:) < 0.05)
            [slx, offx, r2x] = rsquare(gcpts(:,3), gcpts(:,1));
            [sly, offy, r2y] = rsquare(gcpts(:,4), gcpts(:,2));
            idfcfg.CALIBRATION_SLOPES = [slx, offx, r2x; sly, offy, r2y];
        end
    end
    clear cc cp cpnum cptm cpts cr gcpts ncpts numpair offx offy r2x r2y slx sly tnumpat;
end

if isfield(idfcfg, 'DATE')
    try
        idfcfg.DATE = datenum(idfcfg.DATE);
    catch ne_eo;
    	neuroelf_lasterr(ne_eo);
    end
end

if isfield(idfcfg, 'SAMPLE_RATE') && ...
    freq == 0
    freq = idfcfg.SAMPLE_RATE;
end

% reading sample points
disp(' -> reading sample points...');
if isfield(idfcfg, 'NUMBER_OF_SAMPLES')
    snum = idfcfg.NUMBER_OF_SAMPLES;
else
    snum = length(starts) - lc;
end
eyepts  = zeros(snum, 4);
timepts = zeros(snum, 1);
mc = 1;
pc = 1;
sc = size(cols, 2);
npc    = {['\s+' numpat]};
eyepat = ['^(\d+)\s+smp\s+\d+' sprintf('%s', npc{ones(1, max(cols))})];
msgpat = '^(\d+)\s+msg\s+\#?\s*([^\:]+)\:\s+(.*)$';
msgs = struct('time', 0, 'msgtype', '', 'message', '');
timenotset = true;
firsttm = 0;
while lc <= lnum && ...
    pc <= snum
    tline = lower(textline(idf, starts, ends, lc));
    epm = regexp(tline, eyepat, 'tokens');
    if ~isempty(epm) && ...
       ~isempty(epm{1})
        if timenotset
            try
                firsttm = str2double(epm{1}{1});
                idfcfg.FIRST_TIME = firsttm;
                timenotset = false;
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
        try
            timepts(pc) = str2double(epm{1}{1}) - firsttm;
            for cc = 1:sc
                eyepts(pc, cc) = str2double(epm{1}{cols(cc)+1});
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
        pc = pc + 1;
    else
        mpm = regexp(tline, msgpat, 'tokens');
        if ~isempty(mpm) && ...
           ~isempty(mpm{1})
            if timenotset
                try
                    firsttm = str2double(mpm{1}{1});
                    idfcfg.FIRST_TIME = firsttm;
                    timenotset = false;
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
            try
                msgs(mc) = struct( ...
                    'time', str2double(mpm{1}{1}) - firsttm, ...
                    'msgtype', mpm{1}{2}, ...
                    'message', mpm{1}{3});
                mc = mc + 1;
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
    lc = lc + 1;
end
if (pc-1) > snum
    disp('    too many samples.');
else
    if pc <= snum
        eyepts  = eyepts(1:(pc-1), :);
        timepts = timepts(1:(pc-1));
    end
end
pc = size(eyepts, 1);
if mc == 1
    msgs(1) = [];
end
idfcfg.MESSAGES = msgs;
clear ans cc cols ends epm eyepat firsttm idf lc lnum mc mpm msgpat msgs npc numpat sc snum starts timenotset tline;

% replacing blink values
if deblink
    disp(' -> removing eye blinks...');
    maxeye = max(eyepts, [], 1);
    mineye = min(eyepts, [], 1);
    mneye  = mean(eyepts, 1);
    for cc = 1:4
        if abs(maxeye(cc) - mneye(cc)) < abs(mineye(cc) - mneye(cc))
            bleye = mineye(cc);
        else
            bleye = mineye(cc);
        end
        bfound = find(eyepts(:, cc) == bleye);
        if ~isempty(bfound)
            for bc = 1:length(bfound)
                br = bfound(bc);
                if br > 1
                    eyepts(br, :) = eyepts(br-1, :);
                end
            end
        end
    end
    clear bc bfound bleye br cc maxeye mineye mneye;
    idfcfg.DEBLINK = true;
else
    idfcfg.DEBLINK = false;
end

% detecting frequency (if needed)
disp(' -> frequency detection...');
mtdiff = mean(diff(timepts));
dfreq = 1 / (mtdiff * timeunit);
disp(sprintf('    detected: %7.2fHz', dfreq));
idfcfg.DETECTED_FREQUENCY = dfreq;
if freq == 0
    freq = dfreq;
end
clear dfreq mtdiff;

% selecting range
if range(1) > 0 || ...
   (~isinf(range(2)) && range(2) < (pc/freq))
    disp(' -> selecting range...');
    idfcfg.RANGE_SELECTION = range;
    rangesel = find(timepts >= (range(1)/timeunit) & timepts <= (range(2)/timeunit));
    eyepts = eyepts(rangesel, :);
    timepts = timepts(rangesel);
    clear rangesel;
end
clear freq pc range;

% applying filter
if mlf > 1
    disp(' -> applying filter...');
    disp(sprintf('    filter weights: %s', sprintf('%5.3f ', lfilter)));
    idfcfg.FILTER_WEIGHTS = lfilter;
    lflt = length(lfilter);
    feyepts = lfilter(1) .* eyepts(1:(end+1-lflt),:);
    for fcc = 2:length(lfilter)
        feyepts = feyepts + lfilter(fcc) .* eyepts(fcc:(end+fcc-lflt),:);
    end
    eyepts = feyepts;
    clear fcc feyepts lflt;
    timepts = timepts(mlf:size(eyepts,1)+mlf-1);
end
timepts = timepts .* timeunit;
clear bfilter lfilter mlf timeunit;

% detrending
if detrendv > 0
    pc = size(eyepts, 1);
    disp(' -> detrending data...');
    if detrendv > pc
        detrendv = 2;
    end
    idfcfg.MEAN_VALUES = mean(eyepts, 1);
    switch fix(detrendv)
      case {1}
        idfcfg.DETRENDING_METHOD = 'constant';
        eyepts = detrend(eyepts, 'constant');
      case {2}
        idfcfg.DETRENDING_METHOD = 'linear';
        eyepts = detrend(eyepts, 'linear');
      otherwise
        idfcfg.DETRENDING_METHOD = sprintf('spline_%d', detrendv);
        eyepts1 = eyepts(1, :);
        medc   = 1;
        medi   = zeros(length(1:detrendv:pc));
        detr05 = fix(detrendv / 2);
        for mc = 1:detrendv:pc
            medi(medc) = (mc+(min(mc+detrendv-1,pc)))/2;
            meds(medc, 1:4) = median(eyepts(max(1,mc-detr05):min(mc+detr05+detrendv-1,pc),:));
            medc = medc + 1;
        end
        seyepts = zeros(size(eyepts));
        for cc = 1:4
            seyepts(:, cc) = interp1(fix(medi), meds(:,cc), 1:pc, 'spline');
        end
        eyepts = eyepts - seyepts;
        eyepts = eyepts + (ones(pc, 1) * (eyepts1 - eyepts(1, :)));
        eyepts = eyepts - ones(pc, 1) * mean(eyepts, 1);
        clear cc detr05 eyepts1 mc medc medi meds seyepts;
    end
    if deblink
        disp('    advanced deblinking...');
        for cc = 1:4
            vreye = var(eyepts(:, cc)) * 0.6;
            [nh, ch] = hist(eyepts(:, cc), 100);
            snh = cumsum(nh) / sum(nh);
            minp = find(snh > 0.05);
            maxp = find(snh < 0.95);
            minp = ch(minp(1));
            maxp = ch(maxp(end));
            eyepts(eyepts(:, cc) < min(-vreye, minp), cc) = 0;
            eyepts(eyepts(:, cc) > max(vreye, maxp), cc) = 0;
        end
        clear cc ch maxp minp nh pc snh vreye;
    end
end
clear deblink;

% smoothing
if all(smoothv > 0)
    smoothf = smoothv(3);
    smoothr = 1 - smoothf;
    for sc = 1:smoothv(1)
        disp(sprintf(' -> smoothing (pass %d)...', sc));
        svr = var(eyepts) .* smoothv(2);
        deyepts = [0,0,0,0;diff(eyepts)];
        for cc = 1:4
            seyepts = find(deyepts(:, cc) < svr(cc));
            for rc = 2:length(seyepts)
                eyepts(seyepts(rc), cc) = smoothf * eyepts(seyepts(rc) - 1, cc) + smoothr * eyepts(seyepts(rc), cc);
            end
            for rc = length(seyepts)-1:-1:1
                eyepts(seyepts(rc), cc) = smoothf * eyepts(seyepts(rc) + 1, cc) + smoothr * eyepts(seyepts(rc), cc);
            end
        end
    end
end
clear cc deyepts rc sc seyepts smoothf smoothr smoothv svr;

% translation
if translat && ...
   ~strcmp(unit, 'rawpix');
    disp(' -> trying translating coordinates into degrees...');
    idfcfg.TRANSLATION = '';
    if strcmp(unit, 'screenpix')
        idfcfg.TRANSLATION = 'screenpix->';
        if detrendv == 0 && ...
            isfield(idfcfg, 'CALIBRATION_AREA')
            eyepts(:, [1,3]) = eyepts(:, [1,3]) - (idfcfg.CALIBRATION_AREA(1) / 2);
            eyepts(:, [2,4]) = eyepts(:, [2,4]) - (idfcfg.CALIBRATION_AREA(2) / 2);
        else
            if detrendv == 0
                warning( ...
                    'neuroelf:BadFileContent', ...
                    'The screen center coordinate is needed to translate pixels.' ...
                );
                translat = false;
            end
        end
        if ~isfield(idfcfg, 'MM_PER_PIX_')
            warning( ...
                'neuroelf:BadFileContent', ...
                'Can''t determine mm per pixel resolution for translation.' ...
            );
            translat = false;
        else
            eyepts(:, [1,3]) = eyepts(:, [1,3]) .* idfcfg.MM_PER_PIX_(1);
            eyepts(:, [2,4]) = eyepts(:, [2,4]) .* idfcfg.MM_PER_PIX_(2);
            unit = 'mm';
        end
    end
    if translat && ...
        strcmp(unit, 'mm')
        idfcfg.TRANSLATION = [idfcfg.TRANSLATION 'mm->'];
        if ~isfield(idfcfg, 'HEAD_DISTANCE_MM_')
            warning( ...
                'neuroelf:BadFileContent', ...
                'The head distance is needed for translation.' ...
            );
            translat = false;
        else
            eyepts(:) = atan(eyepts(:) ./ idfcfg.HEAD_DISTANCE_MM_);
            unit = 'radiens';
        end
    end
    if translat > 1 && ...
        strcmp(unit, 'radiens')
        idfcfg.TRANSLATION = [idfcfg.TRANSLATION 'radiens->'];
        eyepts = eyepts .* (180/pi);
        unit = 'degrees';
    end
    if strcmp(unit, 'degrees')
        idfcfg.TRANSLATION = [idfcfg.TRANSLATION 'degrees'];
    end
end
clear detrendv translat unit;

% set fields in struct
if isfield(options, 'instruct') && ...
    nargout > 2
    idfcfg.EYE_POINTS  = eyepts;
    idfcfg.TIME_POINTS = timepts;
end

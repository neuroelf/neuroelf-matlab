function slidequant(impath, thresh, range, final)
%SLIDEQUANT  Slide (image) quantification.
%   SLIDEQUANT(FOLDER) quantifies the slide images (JPG/JPEG) in the given
%   path FOLDER.
%
%   SLIDEQUANT(FOLDER, THRESH, RANGE) uses threshold THRESH and range
%   RANGE for the quantification. If not given, these values will be
%   requested from the user.

% select folder if necessary
if nargin < 1 || ~ischar(impath) || isempty(impath) || exist(impath, 'dir') ~= 7
    impath = uigetdir(pwd, 'Please select folder containing image files...');
    if isequal(impath, 0) || isempty(impath) || exist(impath, 'dir') ~= 7
        return;
    end
end
imfiles = cat(1, dir([impath filesep '*.jpg']), dir([impath filesep '*.jpeg']));
imfiles(~cellfun('isempty', regexpi({imfiles.name}, '_filtered'))) = [];
imfiles(~cellfun('isempty', regexpi({imfiles.name}, '^\._'))) = [];

% thresh and range
if nargin < 3 || ~isa(thresh, 'double') || ~isa(range, 'double') || numel(thresh) ~= 1 || numel(range) ~= 1
    uivals = inputdlg({'Threshold (hard)', 'Value range (fuzziness)'}, ...
        'Please enter threshold criterion', 1, {'32', '0'});
    if ~iscell(uivals) || numel(uivals) ~= 2 || isempty(uivals{1}) || isempty(uivals{2})
        return;
    end
    thresh = str2double(uivals{1});
    range = str2double(uivals{2});
end
if isnan(thresh) || isinf(thresh)
    thresh = 32;
end
if isnan(range) || isinf(range)
    range = 0;
end

% check for sub-folders
if isempty(imfiles) && (nargin < 4 || ~islogical(final) || numel(final) ~= 1 || ~final)
    vans = questdlg(['Would you like to scan the subfolders in ' impath '?'], 'Input', 'Yes', 'No', 'Yes');
    if ischar(vans) && strcmpi(vans, 'yes')
        imdirs = dir(impath);
        imdirs = imdirs(cat(1, imdirs.isdir));
        if ~isempty(imdirs) && (strcmp(imdirs(1).name, '.') || strcmp(imdirs(1).name, '..'))
            imdirs(1) = [];
        end
        if ~isempty(imdirs) && (strcmp(imdirs(1).name, '.') || strcmp(imdirs(1).name, '..'))
            imdirs(1) = [];
        end
        for dc = 1:numel(imdirs)
            fibrosis_select([impath filesep imdirs(dc).name], thresh, range, true);
        end
    end
    return;
elseif isempty(imfiles)
    return;
end

% some stats
fprintf('%d images found in %s; using thresh=%d, range=%d\n\n', numel(imfiles), impath, thresh, range);

% iterate
imtext = cell(numel(imfiles)+1, 1);
imtext{1} = sprintf('filename\ttotal_pixel\tliver_pixel\tfibrosis_pixel\tfraction (%%)');
for ic = 1:numel(imfiles)

    % load image
    try
        imfile = [impath filesep imfiles(ic).name];
        im = imread(imfile);
    catch eo
        uiwait(warndlg(['Error reading image: ' eo.message], 'Error', 'modal'));
        return;
    end

    % convert to gray-scale if necessary
    im = double(im);
    if size(im, 3) > 1
        im = mean(im, 3);
    end

    % call function
    stats = fs_stats(im, thresh, range);

    % print some stats
    imtext{ic+1} = sprintf('%s\t%d\t%d\t%d\t%.3g', [impath filesep imfiles(ic).name], ...
        numel(im), stats.liver, stats.fib, 100 * stats.fib / stats.liver);
    fprintf('Image %30s: Liver pixels: %8d / Fibrosis pixels: %6d // ratio: %.4f\n', ...
        imfiles(ic).name, stats.liver, stats.fib, stats.fib / stats.liver);

    % and write out filtered image
    [~, imfile, imext] = fileparts(imfile);
    imwrite(stats.rim, [impath filesep imfile '_filtered' imext], 'Quality', 100);
end
tfile = strrep(strrep(datestr(now), ':', ''), ' ', '_');
fid = fopen([impath filesep 'fibrosis_stats_' tfile '.txt'], 'w');
fprintf(fid, '%s\n', imtext{:});
fclose(fid);


% do the work
function stats = fs_select(im)

% init stats
global mystats;

% general structure
mystats = struct( ...
    'Figure', [], ...
    'Images', im, ...
    'LiverPixels', [], ...
    'FibroticPixels', []);

% create figure
f = figure;
f.Name = 'Fibrosis pixel selector';
f.NumberTitle = 'off';
f.Units = 'pixels';
f.Position(3:4) = imszh + [144, 16];

% axes for image
a = axes('Parent', f, 'Units', 'pixels', 'Position', [8, 8, imszh]);

% image

a.XTick = [];
a.YTick = [];

% done button
done = uicontrol(f, 'Style', 'pushbutton', 'String', 'Done', ...
    'Position', [712, 8, 80, 24]);
done.Callback = @fs_done;

% wait
mystats.Figure = f;
mystats.ImageHist = imh;
mystats.ImageMedian = immed;
mystats.ScaledImage = rim;
waitfor(f);

% return stats
stats = mystats;


% stats
function st = fs_stats(im, thresh, range)

% image size and histogram
nim = numel(im);
imsz = size(im);
imszh = ceil(0.5 .* imsz(1, [2, 1]));
imh = histcounts(im(:), -.5:255.5);
imch = cumsum(imh);
imrch = round((255 / nim) .* imch);
immed = find(imch >= (0.5 * nim), 1, 'first');

% 2D smoothing kernel (3 pixel FWHM)
smk = [.0023, .0188, 0.0909, .231, .314, .231, .0909, .0188, .0023];
smk = smk' * smk;

% smoothed image
sim = conv2(im, smk, 'same');

% second version of image for display (mix rescaled smoothed and original)
rim = round((2/3) .* imrch(round(sim+1)) + (1/3) .* im);

% without background detection, liver is everything
backg = false(imsz);
liver = true(imsz);

% try to detect background
c1 = reshape(rim(3:22, 3:22), 400, 1);
c2 = reshape(rim(imsz(1)-21:imsz(1)-2, 3:22), 400, 1);
c3 = reshape(rim(3:22, imsz(2)-21:imsz(2)-2), 400, 1);
c4 = reshape(rim(imsz(1)-21:imsz(1)-2, imsz(2)-21:imsz(2)-2), 400, 1);
mc1 = median(c1);
mc2 = median(c2);
mc3 = median(c3);
mc4 = median(c4);
c1m = (c1 >= (mc1 - 1) & c1 <= (mc1 + 1));
c2m = (c2 >= (mc2 - 1) & c2 <= (mc2 + 1));
c3m = (c3 >= (mc3 - 1) & c3 <= (mc3 + 1));
c4m = (c4 >= (mc4 - 1) & c4 <= (mc4 + 1));
c1b = sum(c1m) / 400;
c2b = sum(c2m) / 400;
c3b = sum(c3m) / 400;
c4b = sum(c4m) / 400;
[mcv, mcp] = max([c1b, c2b, c3b, c4b], [], 2);

% any corner at least 55% detected background
if mcv >= 0.55
    
    % which corner
    switch mcp
        case 1
            [xi, yi] = find(reshape(c1m, 20, 20));
            xi = xi + 2;
            yi = yi + 2;
        case 2
            [xi, yi] = find(reshape(c2m, 20, 20));
            xi = xi + imsz(1) - 22;
            yi = yi + 2;
        case 3
            [xi, yi] = find(reshape(c3m, 20, 20));
            xi = xi + 2;
            yi = yi + imsz(2) - 22;
        case 4
            [xi, yi] = find(reshape(c4m, 20, 20));
            xi = xi + imsz(1) - 22;
            yi = yi + imsz(2) - 22;
    end

    % original background values
    obv = sort(im(xi + imsz(1) .* (yi - 1)));
    mino = obv(ceil(0.125 * numel(obv))) - 1;
    maxo = obv(floor(0.875 * numel(obv))+1);

    % cluster image
    [~, cv] = clustercoordsc(im >= mino & im <= maxo, 4, ceil(0.001 * nim));

    % first pass of foreground (liver) is where cv == 0
    liver = (cv == 0);
    
    % keep largest chunk only
    [lcs, lcv] = clustercoordsc(liver, 4, ceil(0.25 * nim));
    [~, lcs] = max(lcs);
    if ~isempty(lcs)
        liver = (lcv == lcs);

        % smooth this image to only keep things with certain characteristics
        sliver = conv2(double(conv2(double(~liver), smk, 'same') == 1), smk, 'same');

        % only keep clusters with at least 400 pixels
        [ncs, ncv] = clustercoordsc(sliver == 1, 4, ceil(0.0005 * nim));

        % now construct background
        backg = false(imsz);
        for cc = 1:numel(ncs)
            backp = (cv == median(cv(ncv == cc)));
            if sum(backp(:)) < (2 * sum(ncv(:) == cc))
                backg = backg | backp;
            end
        end

        % and only keep one chunk liver
        [lcs, lcv] = clustercoordsc(~backg, 4, ceil(0.25 * nim));
        [~, lcs] = max(lcs);
        liver = (lcv == lcs);
    else
        liver = true(imsz);
    end
end

% get values within liver
livvals = im(liver);

% weigh fiber voxels around thresh
livvals = max(0, min(1, (1 / range) .* (livvals - thresh)));

% pack together
rim = reshape(repmat(uint8(round(rim)), [1, 1, 3]), nim, 3);
rim(liver, 1) = uint8(min(255, im(liver) + round(80 .* livvals)));
rim(liver, 2) = uint8(min(255, im(liver) + round(128 .* (1 - livvals))));
rim(backg, 3) = uint8(min(255, im(backg) + 48));
rim = reshape(rim, [imsz, 3]);
st = struct('rim', rim, 'backg', sum(backg(:)), 'liver', sum(liver(:)), 'fib', round(sum(livvals)));


% done
function varargout = fs_done(varargin)

% global stats
global mystats;

% do some work
varargout = cell(1, nargout);

% delete figure
delete(mystats.Figure);

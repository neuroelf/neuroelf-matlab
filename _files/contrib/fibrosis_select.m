function fibrosis_select(impath, thresh, range, final)
%FIBROSIS_SELECT  Write out RGB images indicating fibrotic tissue.
%   FIBROSIS_SELECT without any arguments will ask for a folder to operate
%   on. All images in that folder (and if no images found, optionally
%   sub-folders) with the extension of either .jpg or .jpeg will be
%   processed. Without arguments, the threshold and fuzziness range will
%   be requested from the user per UI dialog.
%
%   FIBROSIS_SELECT(IMFOLDER) operates on the specified folder.
%
%   FIBROSIS_SELECT(IMFOLDER, THRESH, RANGE) also specifies the threshold
%   and fuzziness range (better for scripting).
%
%   FIBROSIS_SELECT(IMFOLDER, THRESH, RANGE, FINALFLAG) will not scan
%   subfolders if FINALFLAG is set to true.
%
%   The processed (RGB) images will be written in parallel to the source
%   images (same sub-folder) with a name of FILENAME_filtered.EXT and a
%   JPG quality setting of 100.
%
%   See also IMWRITE.

% last edit: 04/02/2016 (C) JW.

% select folder if necessary
if nargin < 1 || ~ischar(impath) || isempty(impath) || exist(impath, 'file') < 2
    impath = uigetdir(pwd, 'Please select folder containing image files...');
    if isequal(impath, 0) || isempty(impath) || exist(impath, 'dir') ~= 7
        return;
    end
end
if exist(impath, 'file') == 2
    [impath, imfiles, imext] = fileparts(impath);
    if isempty(impath)
        impath = pwd;
    end
    imfiles = struct('name', [imfiles imext], 'isdir', false);
else
    imfiles = cat(1, dir([impath filesep '*.jpg']), dir([impath filesep '*.jpeg']));
end
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

% global handles
global ne_methods;
if isempty(ne_methods)
    n = neuroelf;
end

% image size and histogram
nim = numel(im);
imsz = size(im);
%imszh = ceil(0.5 .* imsz(1, [2, 1]));
imh = histcounts(im(:), -.5:255.5);
imch = cumsum(imh);
imrch = round((255 / nim) .* imch);
%immed = find(imch >= (0.5 * nim), 1, 'first');

% varkerns
% varkern1 = zeros(5, 5);
% varkern1(3, 3) = 1;
% varkern2 = varkern1;
% verkern1 = varkern1 - 1/25;
% varkern2([1,end], :) = -1/16;
% varkern2(:, [1,end]) = -1/16;

% expand image (4 on each side, later truncated again)
im = cat(1, im(1:4, :, :), im, im(end-3:end, :, :));
im = cat(2, im(:, 1:4, :), im, im(:, end-3:end, :));

% 2D smoothing kernel (3 pixel FWHM)
smk = [.0023, .0188, 0.0909, .231, .314, .231, .0909, .0188, .0023];
smk = smk' * smk;

% smoothed image
sim = conv2(im, smk, 'same');

% noise image
nsim = sqrt(log(1 + conv2(conv2(abs(im - sim), smk, 'same'), smk, 'same')));

% background is where noise is low
backg = (nsim < (0.5 * median(nsim(:))));
[lcs, lcv] = ne_methods.clustercoordsc(backg, 4, 0.001 .* nim);
for lcc = 1:numel(lcs)
    backgt = squeeze(ne_methods.dilate3d(permute(lcv == lcc, [3, 1, 2])));
    if sum(backgt(:)) > (1.25 * lcs(lcc))
        lcv(lcv == lcc) = 0;
    end
end
liver = (lcv == 0);
liver(1:4, :) = [];
liver(end-3:end, :) = [];
liver(:, 1:4) = [];
liver(:, end-3:end) = [];
[lcs, lcv] = ne_methods.clustercoordsc(liver, 4, 0.25 .* nim);
if numel(lcs) == 0
    liver = true(imsz);
elseif numel(lcs) == 1
    liver = (lcv > 0);
else
    [lcm, lcmp] = max(lcs);
    liver = (lcv == lcmp);
end
backg = ~liver;
im(1:4, :) = [];
im(end-3:end, :) = [];
im(:, 1:4) = [];
im(:, end-3:end) = [];
sim(1:4, :) = [];
sim(end-3:end, :) = [];
sim(:, 1:4) = [];
sim(:, end-3:end) = [];

% second version of image for display (mix rescaled smoothed and original)
rim = round((2/3) .* imrch(round(sim+1)) + (1/3) .* im);

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

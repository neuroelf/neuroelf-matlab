function count = czi_countcells(imfile, chn)
%CZI_COUNTCELLS Count cells in a CZI image.
%   COUNT = CZI_COUNTCELLS(IMAGE)
%

% requires neuroelf
n = neuroelf;

% input check (image filename)
if nargin < 1 || ~ischar(imfile) || isempty(imfile) || exist(imfile, 'file') ~= 2
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end

% assume channel 4 (if not specified)
if nargin < 2
    chn = 4;
end

[impath, imname] = fileparts(imfile);
if isempty(impath)
    impath = pwd;
end

% load image
try
    fprintf('Reading %s...', imfile);
    czi = xff(imfile);
    pause(0.01);
    
    % check image
    if numel(czi) ~= 1 || ~isxff(czi, 'czi')
        if numel(czi) == 1 && isxff(czi, true)
            czi.ClearObject;
        end
        error('neuroelf:general:badFileContent', 'Bad CZI file content.');
    end
    
    % make sure we have the channels we need
    ch = cell(1, max(4, chn));
    chc = 1;
    for cc = 1:numel(czi.RawContent)
        if ~isempty(czi.RawContent(cc).FileType) && any(strcmpi(czi.RawContent(cc).FileType, {'dvimage', 'dvraw'})) && ...
            isstruct(czi.RawContent(cc).FileCooked) && isfield(czi.RawContent(cc).FileCooked, 'Image')
            ch{chc} = czi.RawContent(cc).FileCooked.Image;
            chc = chc + 1;

            % only up to 
            if chc > max(4, chn)
                break;
            end
        end
    end
    czi.ClearObject;
    
    % make sure enough channels were available
    if chc <= chn
        error('neuroelf:general:badFileContent', 'Image must have at least 4 channels.');
    end
catch ne_eo;
    rethrow(ne_eo);
end

% get channel data as counting channel
countc = double(ch{chn});

% get smoothed version
scountc = n.smoothimg(countc, [11.5, 11.5]);

% resample by factor 4
countc = n.resampleaa(n.resampleaa(countc, 0.25, 1), 0.25, 2);
scountc = n.resampleaa(n.resampleaa(scountc, 0.25, 1), 0.25, 2);

% remove last three samples
countc(:, end-2:end) = [];
countc(end-2:end, :) = [];
scountc(:, end-2:end) = [];
scountc(end-2:end, :) = [];

% difference (high-frequency image)
hcountc = countc - scountc;

% smooth an estimate of background as lack of high-frequency stuff
hback = n.smoothimg(double(abs(hcountc) <= 3 * median(abs(hcountc(:)))), [8, 8]) >= 0.7;
[cs, cv] = n.clustercoordsc(hback, 4);
hback = n.smoothimg(double(cv ~= n.maxpos(cs)), [8, 8]) <= 0.3;

% resampling (flexinterpn) grid
szc = size(countc);

% compute gradient
[gradc, grady] = gradient(countc);
gradc = sqrt(gradc .* gradc + grady .* grady);

% find background
cmf = 0.5;
gmf = 2;
background = (countc < (cmf * mean(countc(gradc >= 1))) & gradc <= gmf) | hback | countc == 0;
[cs, cv] = n.clustercoordsc(background, 4);
foreground = (cv ~= n.maxpos(cs));
mfore = mean(countc(foreground));

% copy the channel we need
count = zeros(size(countc));

% set background to -1
count(~foreground) = -1;

% begin thresholding
tc = 1;
tcn = 100;
keeplooking = true;
minsize = 192;
while keeplooking

    % cluster foreground (given the resampling, >= 48 pixels only!)
    n.scaleimage(double(foreground));
    [cs, cv, cl] = n.clustercoordsc(foreground, 4, minsize);
    if isempty(cs)
        break;
    end
    [~, cli] = sort(cl(:, 4));
    cl = cl(cli, :);
    
    % for each cluster
    cli = 1;
    cst = 4 * median(cs) + 2 * std(cs);
    for cc = 1:numel(cs)
        
        % coordinate list
        clist = cl(cli:cli+cs(cc)-1, :);
        cli = cli + cs(cc);
        
        % cluster too large (> median+std) or small(< 0.25 * median)
        if cs(cc) > cst
            continue;
        end
        
        % for small enough clusters, mark in output
        if cs(cc) <= 5000
            count(clist(:, 1) + szc(1) .* (clist(:, 2) - 1)) = tc;
            tc = tc + 1;
        else
            mclist = min(clist(:, 1:2)) - 1;
            [scs, scv, scl, scc] = n.splitclustercoords([ones(cs(cc), 1), clist(:, 1:2) - ones(cs(cc), 1) * mclist], ...
                scountc(clist(:, 1) + szc(1) .* (clist(:, 2) - 1)), 48, [4, 4, 4]);
            for sscs = 1:numel(scc)
                if size(scc{sscs}, 1) >= 384
                    count(scc{sscs}(:, 2) + mclist(1) + szc(1) .* (scc{sscs}(:, 3) + mclist(2) - 1)) = tc;
                    tc = tc + 1;
                end
            end
        end
        if tc > tcn
            fprintf('%d', mod(0.01 * tcn, 10));
            tcn = tcn + 100;
        end
    end
    
    % mask
    cmask = squeeze(n.dilate3d(n.dilate3d(n.dilate3d(n.dilate3d(permute(count > 0, [3, 1, 2]))))));
    countc(cmask) = 0;
    gradc(cmask) = 0;
    
    % find background
    cmf = cmf + 0.15;
    gmf = gmf - 0.1;
    minsize = 2 .* minsize;
    if cmf > 1.2
        break;
    end
    background = (countc < (cmf * mean(countc(gradc >= 1 & countc > 0))) & gradc <= gmf) | hback | countc < (0.5 * mfore);
    [cs, cv] = n.clustercoordsc(background, 4);
    foreground = (cv ~= n.maxpos(cs)) & countc > 0;
end

% sample down again
count = count(1:4:end, 1:4:end);
rgb = uint8(zeros(numel(count), 3));

% for each of the clusters, print out the mean of each of the channels
tab = zeros(tc-1, 8);
for cc = 1:(tc - 1)
    ccc = find(count(:) == cc);
    [ccx, ccy] = ind2sub(size(count), ccc);
    tab(cc, :) = [cc, numel(ccc), mean(ccy), mean(ccx), ...
        mean(ch{1}(ccc)), mean(ch{2}(ccc)), mean(ch{3}(ccc)), mean(ch{4}(ccc))];
    rgb(ccc, :) = uint8(ones(numel(ccc), 1) * round(floor(255.999 .* rand(1, 3))));
end
n.asciiwrite([impath filesep imname sprintf('_%d_cellcounts.txt', cc)], [ ...
    sprintf('Cellnumber\tpixels\txpos\typos\tch1\tch2\tch3\tch4\n'), ...
    sprintf('%03d\t%5d\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\n', tab')]);
imwrite(reshape(rgb, [size(count), 3]), [impath filesep imname '_cellcounts_rgb.jpg'], 'Quality', 85);
fprintf('... %d cells found.\n', cc);
pause(0.01);
count = cc;

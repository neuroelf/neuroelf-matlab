function count = czi_CountCells(xo, opts)
% CZI::CountCells  - count cells in a CZI image
%
% FORMAT:       count = czi.CountCells([opts])
%
% Input fields:
%
%       opts        1x1 struct with optional settings
%        .channel   CZI channel to use (default: 4)
%        .progress  1x1 logical value, show progress (default: false)
%        .rgbmix    1x1 value to mix original with count (default: 2/3)
%        .smkern    general smoothing kernel (pixels, default: 11.5)
%
% Output fields:
%
%       count       1x1 double with count of cells

% Version:  v1.1
% Build:    16080812
% Date:     Aug-08 2016, 12:10 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'czi')
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'channel') || ~isa(opts.channel, 'double') || numel(opts.channel) ~= 1 || ...
    isinf(opts.channel) || isnan(opts.channel) || opts.channel < 1
    opts.channel = 4;
end
chn = opts.channel;
if ~isfield(opts, 'progress') || ~islogical(opts.progress) || numel(opts.progress) ~= 1
    opts.progress = false;
end
sprog = opts.progress;
if ~isfield(opts, 'rgbmix') || ~isa(opts.rgbmix, 'double') || numel(opts.rgbmix) ~= 1 || ...
    isinf(opts.rgbmix) || isnan(opts.rgbmix) || opts.rgbmix < 0 || opts.rgbmix >= 1
    opts.rgbmix = 2/3;
end
if ~isfield(opts, 'smkern') || ~isa(opts.smkern, 'double') || numel(opts.smkern) ~= 1 || ...
    isinf(opts.smkern) || isnan(opts.smkern) || opts.smkern <= 0
    opts.smkern = 11.5;
else
    opts.smkern = min(40, max(0.5, opts.smkern));
end
smk = opts.smkern;

% get filename
[impath, imname] = fileparts(xo.F);
if isempty(impath)
    impath = pwd;
end
if isempty(imname)
    imname = 'CZI_Image';
end

% test image
try
    
    % make sure we have the channels we need
    ch = cell(1, max(4, chn));
    chc = 1;
    for cc = 1:numel(bc.RawContent)
        if ~isempty(bc.RawContent(cc).FileType) && any(strcmpi(bc.RawContent(cc).FileType, {'dvimage', 'dvraw'})) && ...
            isstruct(bc.RawContent(cc).FileCooked) && isfield(bc.RawContent(cc).FileCooked, 'Image')
            ch{chc} = bc.RawContent(cc).FileCooked.Image;
            chc = chc + 1;

            % only up to 
            if chc > numel(ch)
                break;
            end
        end
    end
    
    % make sure enough channels were available
    chc = chc - 1;
    if chc < chn
        error('neuroelf:xff:badFileContent', 'Image must have at least %d channels.', numel(ch));
    end
catch xfferror
    rethrow(xfferror);
end

% get channel data as counting channel
countc = double(ch{chn});
mcountc = max(countc(:));
if mcountc > 255
    countc = (255 / mcountc) .* countc;
end

% get smoothed version
scountc = ne_methods.smoothimg(countc, [smk, smk]);

% resample by factor 4
countc = ne_methods.resampleaa(ne_methods.resampleaa(countc, 0.25, 1), 0.25, 2);
scountc = ne_methods.resampleaa(ne_methods.resampleaa(scountc, 0.25, 1), 0.25, 2);

% remove last three samples
countc(:, end-2:end) = [];
countc(end-2:end, :) = [];
scountc(:, end-2:end) = [];
scountc(end-2:end, :) = [];

% difference (high-frequency image)
hcountc = countc - scountc;

% smooth an estimate of background as lack of high-frequency stuff
hback = ne_methods.smoothimg(double(abs(hcountc) <= 3 * median(abs(hcountc(:)))), [8, 8]) >= 0.7;
[cs, cv] = ne_methods.clustercoordsc(hback, 4);
hback = ne_methods.smoothimg(double(cv ~= ne_methods.maxpos(cs)), [8, 8]) <= 0.3;

% resampling (flexinterpn) grid
szc = size(countc);

% compute gradient
[gradc, grady] = gradient(countc);
gradc = sqrt(gradc .* gradc + grady .* grady);

% find background
cmf = 0.5;
gmf = 2;
background = (countc < (cmf * mean(countc(gradc >= 1))) & gradc <= gmf) | hback | countc == 0;
[cs, cv] = ne_methods.clustercoordsc(background, 4);
foreground = (cv ~= ne_methods.maxpos(cs));
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
    if sprog
        ne_methods.scaleimage(double(foreground));
    end
    [cs, cv, cl] = ne_methods.clustercoordsc(foreground, 4, minsize);
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
            [scs, scv, scl, scc] = ne_methods.splitclustercoords([ones(cs(cc), 1), clist(:, 1:2) - ones(cs(cc), 1) * mclist], ...
                scountc(clist(:, 1) + szc(1) .* (clist(:, 2) - 1)), 48, [4, 4, 4]);
            for sscs = 1:numel(scc)
                if size(scc{sscs}, 1) >= 384
                    count(scc{sscs}(:, 2) + mclist(1) + szc(1) .* (scc{sscs}(:, 3) + mclist(2) - 1)) = tc;
                    tc = tc + 1;
                end
            end
        end
        if tc > tcn && sprog
            fprintf('%d', mod(0.01 * tcn, 10));
            tcn = tcn + 100;
        end
    end
    
    % mask
    cmask = squeeze(ne_methods.dilate3d(ne_methods.dilate3d(ne_methods.dilate3d(ne_methods.dilate3d(permute(count > 0, [3, 1, 2]))))));
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
    [cs, cv] = ne_methods.clustercoordsc(background, 4);
    foreground = (cv ~= ne_methods.maxpos(cs)) & countc > 0;
end

% sample down again
count = count(1:4:end, 1:4:end);
rgb = zeros(numel(count), 3);

% for each of the clusters, print out the mean of each of the channels
tab = zeros(tc - 1, 4 + chc);
for cc = 1:(tc - 1)
    ccc = find(count(:) == cc);
    [ccx, ccy] = ind2sub(size(count), ccc);
    ms = zeros(1, chc);
    for chcc = 1:chc
        ms(chcc) = mean(ch{chcc}(ccc));
    end
    tab(cc, :) = [cc, numel(ccc), mean(ccy), mean(ccx), ms];
    rgb(ccc, :) = ones(numel(ccc), 1) * floor(255 .* rand(1, 3));
end
rgb = reshape(rgb, [size(count), 3]);
rgb = uint8(round(opts.rgbmix .* (255 / double(max(ch{chn}(:)))) .* ...
    double(repmat(ch{chn}, [1, 1, 3])) + (1 - opts.rgbmix) .* rgb));
chline = sprintf('\\tch%d', 1:chc);
dtline = ['%03d\t%5d\t%5.1f\t%5.1f' repmat('\t%5.1f', 1, chc), '\n'];
ne_methods.asciiwrite([impath filesep imname sprintf('_%d_cellcounts.txt', cc)], [ ...
    sprintf(['Cellnumber\tpixels\txpos\typos' chline '\n']), sprintf(dtline, tab')]);
imwrite(rgb, [impath filesep imname '_cellcounts_rgb.jpg'], 'Quality', 95);
if sprog
    fprintf('... %d cells found.\n', cc);
    pause(0.01);
end
count = cc;

function imf = histmatch(imf, opts)
%HISTMATCH  match image by their histograms
%   HISTMATCH(IMAGES) matches the images with the default options.
%
%   HISTMATCH(IMAGES, OPTS) matches images with overridden options.
%
%   IMAGES is a cell array of either filenames or HxWx3 image data.
%
%   OPTS supports the following fields:
%
%    .histbins  number of histogram bins (default: 128)
%    .hsv       match image in HSV space instead [default: false]
%    .matching  approach, one of 'sample', 'shift', or ['spread']
%    .postfix   string added to filenames [default: '_hm']
%    .rgbsep    match each color channel separately [default: true]
%    .smooth    smooth histograms prior to sampling matching [default: 8]
%    .weights   1xNUMEL(IMAGES) vector of weights of shared histogram
%               [default: ones(1, numel(IMAGES))]
%
%   IMAGES = HISTMATCH(IMAGES) returns the image data as cell array.
%
%   The matching approaches are:
%
%   - sample matches the histogram(s) completely (such that values
%     match across the entire spectrum)
%   - shift only shifts the histogram(s) to fit the center of mass
%   - spread shifts and stretches the histogram(s) to fit center of mass
%     and variability

% Version:  v1.1
% Build:    22090910
% Date:     Sep-09 2022, 10:34 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2017, Jochen Weber
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

% check inputs
if nargin < 1 || ~iscell(imf) || isempty(imf)
    error('neuroelf:missingOrBadInput', 'Missing or bad IMAGES argument.');
end
imf = imf(:);
numim = numel(imf);
if any(cellfun('isempty', imf))
    error('neuroelf:emptyInput', 'IMAGES cells must all be non-empty.');
end
imc = cellfun(@class, imf, 'UniformOutput', false);
if ~all(strcmp(imc, imc{1}))
    error('neuroelf:typeMismatch', 'IMAGES cells must all be of the same type.');
end
imc = imc{1};
if ~any(strcmp(imc, {'double', 'uint8', 'char'}))
    error('neuroelf:invalidClass', 'IMAGES cells must either be double, uint8, or char.');
end
if ischar(imf{1})
    imn = imf;
    try
        for c = 1:numim
            imf{c} = imread(imn{c});
        end
        imc = cellfun(@class, imf, 'UniformOutput', false);
        if ~all(strcmp(imc, imc{1}))
            error('neuroelf:typeMismatch', 'IMAGES cells must all be of the same type.');
        end
        imc = imc{1};
        if ~any(strcmp(imc, {'double', 'uint8'}))
            error('neuroelf:invalidClass', 'Loaded images must either be double or uint8.');
        end
    catch ne_eo;
        error('neuroelf:imageReadError', 'Error reading image %d (%s).', ...
            c, ne_eo.message);
    end
else
    imn = cell(size(imf));
end

% options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'histbins') || ~isa(opts.histbins, 'double') || ...
    numel(opts.histbins) ~= 1 || isinf(opts.histbins) || isnan(opts.histbins) || ...
    opts.histbins < 8 || opts.histbins > 256 || opts.histbins ~= round(opts.histbins)
    opts.histbins = 128;
end
if ~isfield(opts, 'hsv') || ~islogical(opts.hsv) || numel(opts.hsv) ~= 1
    opts.hsv = false;
end
if ~isfield(opts, 'matching') || ~ischar(opts.matching) || ...
   ~any(strcmpi(opts.matching(:)', {'sample', 'shift', 'spread'}))
    opts.matching = 'spread';
else
    opts.matching = opts.matching(:)';
end
opts.matching = opts.matching(2);
if ~isfield(opts, 'postfix') || ~ischar(opts.postfix) || isempty(opts.postfix)
    opts.postfix = '_hm';
else
    opts.postfix = opts.postfix(:)';
end
if ~isfield(opts, 'rgbsep') || ~islogical(opts.rgbsep) || numel(opts.rgbsep) ~= 1
    opts.rgbsep = true;
end
if ~isfield(opts, 'smooth') || ~isa(opts.smooth, 'double') || numel(opts.smooth) ~= 1 || ...
    isinf(opts.smooth) || isnan(opts.smooth) || opts.smooth < 0
    opts.smooth = min(8, opts.histbins / 4);
end
if ~isfield(opts, 'weights') || ~isa(opts.weights, 'double') || numel(opts.weights) ~= numim
    opts.weights = ones(numim, 1);
else
    opts.weights = opts.weights(:);
    opts.weights(isinf(opts.weights) | isnan(opts.weights) | opts.weights < 0) = 0;
end
sw = sum(opts.weights);
if sw == 0
    error('neuroelf:badOption', 'Not all weights can be 0.');
end

% for RGB data
hb = opts.histbins;
if ~opts.hsv

    % arguments for histcount
    if imc(1) == 'u'
        hmax = 256;
    else
        hmax = 1;
    end
    hstep = hmax / hb;
    hmax = hmax - 0.5 * hstep;
    harg = {0, hmax, hstep};

    % separate colors
    if opts.rgbsep
        hists = zeros(hb, 3);
    else
        hists = zeros(hb, 1);
    end

% HSV
else
    hstep = 1 / hb;
    hmax = 1 - 0.5 * hstep;
    harg = {0, hmax, hstep};
    
    % separate colors
    if opts.rgbsep
        hists = zeros(hb, 3);
    else
        hists = zeros(hb, 2);
    end
end
hscale = (0:hstep:hmax)';

% iterate over images (first pass)
ims = cell(numim, 1);
h = repmat({hists}, numim, 1);
for c = 1:numim
    
    % get image size
    ims{c} = size(imf{c});
    imn = ims{c}(1) * ims{c}(2);
    
    % RGB space
    if ~opts.hsv
        
        % with separate colors
        if opts.rgbsep
            h{c}(:, 1) = (1 / imn) .* histcount(imf{c}(:, :, 1), harg{:})';
            if opts.weights(c) > 0
                hists(:, 1) = hists(:, 1) + opts.weights(c) .* h{c}(:, 1);
            end
            h{c}(:, 2) = (1 / imn) .* histcount(imf{c}(:, :, 2), harg{:})';
            if opts.weights(c) > 0
                hists(:, 2) = hists(:, 2) + opts.weights(c) .* h{c}(:, 2);
            end
            h{c}(:, 3) = (1 / imn) .* histcount(imf{c}(:, :, 3), harg{:})';
            if opts.weights(c) > 0
                hists(:, 3) = hists(:, 3) + opts.weights(c) .* h{c}(:, 3);
            end
            
        % gray scale (contrast matching)
        else
            h{c}(:) = (1 / imn) .* histcount(mean(imf{c}, 3), harg{:})';
            if opts.weights(c) > 0
                hists = hists + opts.weights(c) .* h{c};
            end
        end
        
    % HSV space
    else
        
        % convert image
        imf{c} = hsvconv(imf{c}, 2);

        % with separate colors
        if opts.rgbsep
            h{c}(:, 1) = (1 / imn) .* histcount(imf{c}(:, :, 1), harg{:})';
            if opts.weights(c) > 0
                hists(:, 1) = hists(:, 1) + opts.weights(c) .* h{c}(:, 1);
            end
            h{c}(:, 2) = (1 / imn) .* histcount(imf{c}(:, :, 2), harg{:})';
            if opts.weights(c) > 0
                hists(:, 2) = hists(:, 2) + opts.weights(c) .* h{c}(:, 2);
            end
            h{c}(:, 3) = (1 / imn) .* histcount(imf{c}(:, :, 3), harg{:})';
            if opts.weights(c) > 0
                hists(:, 3) = hists(:, 3) + opts.weights(c) .* h{c}(:, 3);
            end
            
        % without color, skip HUE slice
        else
            h{c}(:, 1) = (1 / imn) .* histcount(imf{c}(:, :, 2), harg{:})';
            if opts.weights(c) > 0
                hists(:, 1) = hists(:, 1) + opts.weights(c) .* h{c}(:, 1);
            end
            h{c}(:, 2) = (1 / imn) .* histcount(imf{c}(:, :, 3), harg{:})';
            if opts.weights(c) > 0
                hists(:, 2) = hists(:, 2) + opts.weights(c) .* h{c}(:, 2);
            end
        end
    end
    
    % force to double precision
    imf{c} = double(imf{c});
end

% correct for weights
hists = (1 / sw) .* hists;

% find optimal hue as anchor
if opts.hsv && opts.rgbsep
    hn = numel(hscale);
    ht = hists(:, 1);
    ht = ht - mean(ht);
    hm = hmax + 0.5 * hstep;
    ha = sin((2 * pi / hm) .* (0:hstep:hmax) + 0.5 * pi);
    ha = repmat(ha(:), 1, hn);
    for hc = 2:hn
        ha(:, hn + 2 - hc) = ha(1 + mod(hc-1:hc+numel(hscale)-2, numel(hscale)), 1);
    end
    htm = sum(repmat(ht, 1, numel(hscale)) .* ha);
    [~, htm] = max(htm);
    hnh = ceil(0.5 * hn);
    hts = hnh + htm;
    hta = 1 - mod(hscale(htm) + 0.5, 1);
    hists(:, 1) = hists(1+mod(hts-1:hts+hn-2, hn), 1);
    for c = 1:numim
        h{c}(:, 1) = h{c}(1+mod(hts-1:hts+hn-2, hn), 1);
    end
end

% sample
if opts.matching == 'a'
    
% shift/spread
else
    
    % compute center of mass
    com = sum(hscale(:, ones(1, size(hists, 2))) .* hists);
    
    % and spread
    if opts.matching == 'p'
        vod = hscale(:, ones(1, size(hists, 2))) - ones(size(hists, 1), 1) * com;
        som = sqrt(sum(hists .* hists .* vod .* vod));
    end
end

% second pass
for c = 1:numim
    
    % sampling
    if opts.matching == 'a'
        
        imn = ims{c}(1) * ims{c}(2);
        
    % shift/spread
    else
        
        % compute sample mean
        scom = sum(hscale(:, ones(1, size(hists, 2))) .* h{c});
        
        % and spread
        if opts.matching == 'p'
            svod = hscale(:, ones(1, size(h{c}, 2))) - ones(size(h{c}, 1), 1) * scom;
            ssom = sqrt(sum(h{c} .* h{c} .* svod .* svod));
        end
        
        % RGB space
        if ~opts.hsv
            
            % with separate colors
            if ~opts.rgbsep
                com = com([1, 1, 1]);
                scom = scom([1, 1, 1]);
                if opts.matching == 'p'
                    som = som([1, 1, 1]);
                    ssom = ssom([1, 1, 1]);
                end
            end
            
            % iterate over planes
            for s = 1:3
                if opts.matching == 'p'
                    p = imf{c}(:, :, s) - scom(s);
                    p = p .* (som(s) / ssom(s)) + som(s);
                    imf{c}(:, :, s) = min(hmax + 0.499 * hstep, max(0, p + com(s)));
                else
                    imf{c}(:, :, s) = min(hmax + 0.499 * hstep, max(0, ...
                        imf{c}(:, :, s) + (com(s) - scom(s))));
                end
            end
            
        % HSV space
        else
            
            % apply to color as well
            if opts.rgbsep
                for s = 1:3
                    
                    % anchor HUE
                    if s == 1
                        imf{c}(:, :, 1) = mod(imf{c}(:, :, 1) + hta, 1);
                    end
                    if opts.matching == 'p'
                        p = imf{c}(:, :, s) - scom(s);
                        p = p .* (som(s) / ssom(s)) + som(s);
                        imf{c}(:, :, s) = min(hmax + 0.499 * hstep, max(0, p + com(s)));
                    else
                        imf{c}(:, :, s) = min(hmax + 0.499 * hstep, max(0, ...
                            imf{c}(:, :, s) + (com(s) - scom(s))));
                    end
                    if s == 1
                        imf{c}(:, :, 1) = mod(imf{c}(:, :, 1) + (1 - hta), 1);
                    end
                end
            else
                for s = 1:2
                    if opts.matching == 'p'
                        p = imf{c}(:, :, s+1) - scom(s);
                        p = p .* (som(s) / ssom(s)) + som(s);
                        imf{c}(:, :, s+1) = min(hmax + 0.499 * hstep, max(0, p + com(s)));
                    else
                        imf{c}(:, :, s+1) = min(hmax + 0.499 * hstep, max(0, ...
                            imf{c}(:, :, s+1) + (com(s) - scom(s))));
                    end
                end
            end
            
            % convert back
            imf{c} = hsvconv(imf{c}, 1);
        end
    end
    
    % correct class
    if imc(1) == 'u'
        imf{c} = uint8(imf{c});
    end
end

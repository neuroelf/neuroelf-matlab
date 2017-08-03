function [meanimg, zeropart, maxhistp, data] = reccheck(recfile, dsize, wimg)
% reccheck  - perform some preliminary checks on a (PAR/) REC file
%
% FORMAT:       [meanimg, zeropart, maxhistp] = reccheck(recfile, dsize [, wimg])
%
% Input fields:
%
%       recfile     REC filename
%       dsize       data size: [rows, cols, dyns, slices]
%       wimg        write output images (mean and histograms, default: true)
%
% Output fields:
%
%       meanimg     uint8 (grayscale) mean image of the raw data

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
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
if nargin < 2 || ...
   ~ischar(recfile) || ...
    isempty(recfile) || ...
    exist(recfile(:)', 'file') ~= 2 || ...
   ~isa(dsize, 'double') || ...
   ~isequal(size(dsize), [1, 4]) || ...
    any(isinf(dsize) | isnan(dsize) | dsize < 1 | dsize ~= fix(dsize))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
recfile = recfile(:)';
if nargin < 3 || ...
   ~islogical(wimg) || ...
    numel(wimg) ~= 1
    wimg = true;
end

% open the REC file read-only with little-endian byte ordering
try
    rec = 0;
    rec = fopen(recfile, 'r', 'ieee-le');
    if rec < 1
        error( ...
            'neuroelf:FileOpenError', ...
            'Error opening REC file.' ...
        );
    end

    % get filesize
    fseek(rec, 0, 1);
    fsize = ftell(rec);
    fseek(rec, 0, -1);

    % check filesize
    if fsize < (2 * prod(dsize))
        error( ...
            'neuroelf:BadFileSize', ...
            'REC file size too small for 2 * prod(dsize).' ...
        );
    end
catch ne_eo;
    oeo = ne_eo;
    if rec > 0
        try
            fclose(rec);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
    rethrow(oeo);
end

% read the data (Matlab requires to read as a vector)
data = fread(rec, [1, prod(dsize)], '*int16');

% close file
fclose(rec);

% get mean value (where data > 0)
meanvalue = sum(data) ./ sum(data > 0);

% then reshape to 4D array
data = reshape(data, dsize);

% check dyns vs. slices
if dsize(3) < dsize(4)
    data = permute(data, [1, 2, 4, 3]);
    dsize = size(data);
end

% get the max. value (for histogram creation, etc.)
maxvalue = double(max(data(:)));

% compute full histogram
hdata = zeros(1, maxvalue + 1);
for dc = 1:dsize(3)
    hdata = hdata + hist(single(lsqueeze(data(:, :, dc, :))), 0:maxvalue);
end

% store portion of zeros
zeropart = hdata(1) / numel(data);

% then remove
hdata(1) = [];

% find maximum (mode) of values >= mean
[maxhist, maxhistp] = max(hdata(ceil(meanvalue):maxvalue));
maxhistp = maxhistp + ceil(meanvalue - 1);

% plot full histogram
f = figure;
drawnow;
bar(1:maxvalue, hdata, 'histc');

% get screenshot of full histogram
drawnow;
sshot = getframe(f);
fullhist = sshot.cdata;

% get only up to 95th percentile of distribution (to get bulk of data hist)
hdatas = cumsum(hdata) ./ sum(hdata);
hdatau = hdata(hdatas <= 0.95);
perc95 = numel(hdatau);
hdatau99 = hdata(hdatas <= 0.99);

% adapt axes
set(gca, 'XLim', [0, perc95], 'YLim', [0, ceil(1.25 * maxhist)]);

% get second screenshot
drawnow;
sshot = getframe(f);
bulkhist = sshot.cdata;
delete(f);

% prepare mean image
mdata = squeeze(mean(data, 3));
rows = size(mdata, 1);
cols = size(mdata, 2);
slices = size(mdata, 3);
nrows = ceil(sqrt(1.5 * slices));
ncols = ceil(slices / nrows);
meanimg = uint8(zeros(rows*nrows, cols*ncols, 3));

% fill mean image
try
    for cr = 1:nrows
        for cc=1:ncols
            meanimg((cr-1)*rows+1:cr*rows, (cc-1)*cols+1:cc*cols, :) = ...
                repmat(uint8(round(min(255, (255 / numel(hdatau99)) .* ...
                mdata(:, :, (cr-1)*ncols+cc)))), [1, 1, 3]);
        end
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% write images (into current folder)
if wimg
    outfile = recfile(1:end-4);
    imwrite(fullhist, [outfile '_fullhist.png']);
    imwrite(bulkhist, [outfile '_bulkhist.png']);
    imwrite(meanimg, [outfile '_meanimg.png']);
end

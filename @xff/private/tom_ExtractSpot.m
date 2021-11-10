function [cut, im] = tom_ExtractSpot(xo, crd, csz, msz)
% TOM::MarkSpot  - mark and extract a spot around a coordinate
%
% FORMAT:       [cut, im] = tom.ExtractSpot(coord [, cutsize]);
%
% Input fields:
%
%       coord       3-element coordinate around which to cut (tnum, cx, cy)
%       cutsize     optional cut size (default: 256 pixels)
%       marksize    optional mark size (default: 0 = no marking)
%
% Output fields:
%
%       cut         cut-out piece
%       im          full texture image
%
% Using: catstruct.

% Version:  v1.1
% Build:    21102012
% Date:     Oct-20 2021, 12:50 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2021, Jochen Weber
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

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'tom')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
t = xo.C;
tf = xo.F;
td = fileparts(tf);
f = t.Field;
fn = {f.ContentName};
txs = find(strcmpi(fn, 'txtrjpg_') | strcmpi(fn, 'txtrjpga'));
mn = {t.MetaData.Name};
mi = find(strcmpi(mn(:), 'texcameramodelnames'));
modelnames = ne_methods.splittocellc(t.MetaData(mi).Content, char(0));
if numel(mi) ~= 1
    error('neuroelf:xff:badArgument', 'Cannot resolve camera/image names.');
end
if iscell(crd) && numel(crd) > 1
    modelindex = find(strcmpi(modelnames(:), crd{1}));
    if numel(modelindex) ~= 1
        error('neuroelf:xff:badArgument', 'Cannot resolve camera/image names.');
    end
    if numel(crd) == 2 && isa(crd{2}, 'double')
        crd = [modelindex, crd{2}(:)'];
    elseif numel(crd) > 2 && isa(crd{2}, 'double') && numel(crd{2}) == 1 && ...
        isa(crd{3}, 'double') && numel(crd{3}) == 1
        crd = [modelindex, crd{2}, crd{3}];
    end
end
if ~isa(crd, 'double') || numel(crd) ~= 3 || any(isinf(crd) | isnan(crd))
    error('neuroelf:xff:badArgument', 'Bad or missing COORD argument.');
end
if nargin < 3 || ~isa(csz, 'double') || numel(csz) ~= 1 || isinf(csz) || isnan(csz)
    csz = 256;
else
    csz = max(40, abs(round(csz)));
end
if nargin < 4 || ~isa(msz, 'double') || numel(msz) ~= 1 || isinf(msz) || isnan(msz)
    msz = 0;
else
    msz = min(round(0.8 * csz), abs(round(msz)));
end

% image/texture
msel = abs(round(crd(1)));
mcrd = [crd(2), 1 - crd(3)];
if msel > numel(txs)
    error('Invalid texture number.');
end
txfolder = sprintf('%s/textures', td);
if exist(txfolder, 'dir') ~= 7
    mkdir(td, 'textures');
end

% read corresponding image
tname = sprintf('%s/%s.jpg', txfolder, modelnames{msel});
if exist(tname, 'file') ~= 2
    fid = fopen(tname, 'wb');
    if fid < 1
        error('Cannot write temp file.');
    end
    fwrite(fid, txc.Content, 'uint8');
    fclose(fid);
end
im = imread(tname);
isz = size(im);

% mark and cut spot
if any(mcrd > 1)
    mcrd = mcrd ./ isz([2, 1]);
end
row = 1 + floor(isz(1) * mcrd(2));
col = 1 + floor(isz(2) * mcrd(1));
chalf = ceil(0.5 * csz);
if mcrd(1) < 0.5
    fcol = max(1, col - chalf);
    tcol = fcol + csz - 1;
else
    tcol = min(isz(2), col + chalf);
    fcol = tcol + 1 - csz;
end
if mcrd(2) < 0.5
    frow = max(1, row - chalf);
    trow = frow + csz - 1;
else
    trow = min(isz(1), row + chalf);
    frow = trow + 1 - csz;
end
cut = im(frow:trow, fcol:tcol, :);
if msz > 0
    mh = ceil(0.5 * msz);
    crow = row - frow + 1;
    ccol = col - fcol + 1;
    mfrow = max(1, crow - mh);
    mtrow = min(size(cut, 1), crow + mh);
    mfcol = max(1, ccol - mh);
    mtcol = min(size(cut, 2), ccol + mh);
    cut([mfrow, mtrow], mfcol:mtcol, :) = 255;
    cut(mfrow:mtrow, [mfcol, mtcol], :) = 255;
end
im([frow:frow+3, trow-3:trow], fcol:tcol, :) = 255;
im(frow:trow, [fcol:fcol+3, tcol-3:tcol], :) = 255;

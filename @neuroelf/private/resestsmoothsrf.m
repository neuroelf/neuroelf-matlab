function [fwhm, fi] = resestsmoothsrf(res, srf, opts)
% resestsmoothsrf  - estimate SRF-based smoothness from regression residual
%
% FORMAT:       [fwhm, fi] = resestsmoothsrf(res, srf [, opts])
%
% Input fields:
%
%       res         regression residual
%       srf         SRF object (with coordinates and neighbors)
%       opts        optional settings struct
%        .tdim      temporal dimension (default: last)
%
% Output fields:
%
%       fwhm        overall FWHM estimate
%       fi          smoothness per vertex
%
% Note: this function uses ideas taken from a FMRIB webpage at
%       http://www.fmrib.ox.ac.uk/analysis/techrep/tr00df1/tr00df1/node6.html

% Version:  v1.1
% Build:    16031200
% Date:     Mar-12 2016, 12:35 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
   ~isnumeric(res) || ...
    ndims(res) ~= 2 || ...
    numel(srf) ~= 1 || ...
   ~isxff(srf, 'srf') || ...
   ~any(size(res) == srf.NrOfVertices)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
srfh = srf.Handles;
if ~isfield(srfh, 'SmoothingInfo') || ...
   ~isstruct(srfh.SmoothingInfo)
    try
        srf.SetSmoothingInfo;
    catch ne_eo;
        rethrow(ne_eo);
    end
    srfh = srf.Handles;
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'tdim') || ...
   ~isa(opts.tdim, 'double') || ...
    numel(opts.tdim) ~= 1 || ...
    isinf(opts.tdim) || ...
    isnan(opts.tdim) || ...
   ~any(opts.tdim == 1:2)
    opts.tdim = 3 - findfirst(size(res) == srf.NrOfVertices, -1);
end

% get neighbors list
ndist = srfh.SmoothingInfo.NeighborsDist;
nlist = srfh.SmoothingInfo.NeighborsList;
nn = size(nlist, 2);

% z-transform residuals
if ~isa(res, 'double')
    res = double(res);
end
if opts.tdim == 2
    res = res';
end
res = ztrans(res, 1);

% compute correlation
rres = squeeze((1 / (size(res, 1) - 1)) .* ...
    sum(repmat(res, [1, 1, nn]) .* reshape(res(:, nlist(:)), [size(res), nn]), 1));
rres(isinf(rres) | rres > (1 - sqrt(eps))) = NaN;

% compute FWHM estimate
sef = sqrt(8 * log(2));
ft = (sef .* ndist) .* sqrt(-1 ./ (4 .* log(rres)));

% compute outputs
fi = meannoinfnan(ft, 2);
fwhm = median(fi);

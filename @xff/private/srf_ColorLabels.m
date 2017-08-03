function xo = srf_ColorLabels(xo, cols, opts)
% SRF::ColorLabels  - set vertices of labels (networks) to specified colors
%
% FORMAT:       [srf = ] srf.ColorLabels(colors [, opts])
%
% Input fields:
%
%       colors      either Lx3 (RGB) or Lx1 (heatmap) color values
%       opts        optional settings
%        .labels    list of labels to color (default: all found)
%        .lut       color lookup table (for Lx1 values, default: xstandard)
%        .maxval    maximum value (default: largest absolute value)
%        .minval    minimum value (default: smallest absolute value > 0)
%        .upfunc    update function handle (displayed surfaces), default:
%                   {@neuroelf_gui, 'srf_tools', 'update'}

% Version:  v1.1
% Build:    16060813
% Date:     Jun-08 2016, 1:50 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% NeuroElf methods
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'}) || ...
   ~isa(cols, 'double') || isempty(cols) || ~any(size(cols, 2) == [1, 3]) || ...
    any(isinf(cols(:)) | isnan(cols(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
rtv = xo.C.RunTimeVars;
if ~isfield(rtv, 'Labels') || ~isa(rtv.Labels, 'double') || isempty(rtv.Labels) || ...
    size(cols, 1) > numel(rtv.Labels) || ~isfield(rtv, 'LabelMeshVertices') || ...
    numel(rtv.LabelMeshVertices) ~= numel(rtv.Labels)
    error('neuroelf:xff:badObject', 'Surface doesn''t have correct labels.');
end
cslc = cumsum(rtv.LabelMeshVertices(:));
if cslc(end) ~= size(xo.C.VertexCoordinate, 1)
    error('neuroelf:xff:badObject', 'Invalid number of vertices.');
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'labels') || ~isa(opts.labels, 'double') || numel(opts.labels) ~= size(cols, 1)
    opts.labels = rtv.Labels(:);
else
    opts.labels = opts.labels(:);
    if any(isinf(opts.labels) | isnan(opts.labels)) || ~isempty(setdiff(opts.labels, rtv.Labels(:)))
        error('neuroelf:xff:badArgument', 'Unknown label(s) found.');
    end
end
if ~isfield(opts, 'upfunc') || ~iscell(opts.upfunc) || ...
   (~isempty(opts.upfunc) && ~isa(opts.upfunc{1}, 'function_handle'))
    opts.upfunc = {@neuroelf_gui, 'srf_tools', 'update'};
end

% requires LUT coloring
numv = size(cols, 1);
if size(cols, 2) == 1
    acols = abs(cols);
    if ~isfield(opts, 'lut') || ~isa(opts.lut, 'double') || isempty(opts.lut) || ...
        mod(size(opts.lut, 1), 2) ~= 0 || size(opts.lut, 2) ~= 3 || ...
        any(isinf(opts.lut(:)) | isnan(opts.lut(:)) | opts.lut(:) < 0 | opts.lut(:) > 255)
        opts.lut = ([255 * ones(1, 40), zeros(1, 40); ...
            [60:5:255, 60:5:255]; zeros(1, 40), 248, 246, 243, 239, 235:-5:60])';
    else
        opts.lut = round(opts.lut);
    end
    if ~isfield(opts, 'maxval') || ~isa(opts.maxval, 'double') || numel(opts.maxval) ~= 1 || ...
        isinf(opts.maxval) || isnan(opts.maxval) || opts.maxval <= 0
        opts.maxval = max(acols);
    end
    if ~isfield(opts, 'minval') || ~isa(opts.minval, 'double') || numel(opts.minval) ~= 1 || ...
        isinf(opts.minval) || isnan(opts.minval) || opts.minval <= 0
        opts.minval = min(acols(acols ~= 0));
    end
    if opts.maxval <= opts.minval
        opts.maxval = opts.minval + sqrt(eps);
    end
    luts = 0.5 * size(opts.lut, 1);
    colv = 1 + (luts .* (cols < 0)) + (luts - 1) .* ...
        min(1, max(0, (acols - opts.minval) ./ (opts.maxval - opts.minval)));
    colv = [repmat(colv, 3, 1), reshape(ones(numv, 1) * (1:3), 3 * numv, 1)];
    cols = reshape(ne_methods.flexinterpn(opts.lut, colv), numv, 3);
    ncols = find(acols < opts.minval);
    if ~isempty(ncols)
        cols(ncols, :) = repmat(round(255 .* xo.C.ConvexRGBA(1, 1:3)), numel(ncols), 1);
    end
end

% compute start vertices
cslb = 1 + [0; cslc(1:end-1)];

% iterate over labels
ulab = opts.labels;
for lc = 1:numv

    % set colors
    li = find(rtv.Labels == ulab(lc));
    li1 = cslb(li);
    li2 = cslc(li);
    nli = 1 + li2 - li1;
    xo.C.VertexColor(li1:li2, :) = repmat([NaN, cols(lc, :)], nli, 1);
end

% update
if ~isempty(opts.upfunc)
    try
        feval(opts.upfunc{:});
    catch xfferror
        fprintf('Error updating surface: %s\n', xfferror.message);
    end
end

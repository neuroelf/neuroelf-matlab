function c = aft_LabelCoords(xo, opts)
% AFT::LabelCoords  - return a cell array (Nx1) of label coordinates
%
% FORMAT:       c = obj.LabelCoords([opts])
%
% Input fields:
%
%       opts        optional settings
%        .labels    labels (default: all found > 0)
%        .space     coordinate space, one of 'bvi', 'bvs', {'tal'}
%
% Output fields:
%
%       c           Lx1 cell array with Cx3 coordinates
%
% TYPES: HDR, HEAD, VMR

% Version:  v1.1
% Build:    16060711
% Date:     Jun-07 2016, 11:32 AM EST
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

% global methods
global ne_methods;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'hdr', 'head', 'vmr'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
ft = lower(xo.S.Extensions{1});
if strcmp(ft, 'hdr')
    vol = bc.VoxelData(:, :, :, 1);
    ftrf = hdr_CoordinateFrame(xo);
elseif strcmp(ft, 'head')
    vol = bc.Brick(1).Data(:, :, :);
    ftrf = head_CoordinateFrame(xo);
else
    vol = bc.VMRData(:, :, :);
    ftrf = struct('Trf', ne_methods.bvcoordconv(zeros(0, 3), 'bvc2tal', aft_BoundingBox(xo)));
end
vsz = size(vol);
ftrf = ftrf.Trf';
vol = vol(:);
ul = double(unique(vol));
ul(isinf(ul) | isnan(ul) | ul == 0) = [];
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'labels') || ~isa(opts.labels, 'double') || isempty(opts.labels) || ...
    any(isinf(opts.labels(:)) | isnan(opts.labels(:)))
    opts.labels = ul(:);
end
if ~isfield(opts, 'space') || ~ischar(opts.space) || isempty(opts.space) || ...
   ~any(lower(opts.space(end)) == 'ils')
    sp = 'l';
else
    sp = lower(opts.space(end));
end

% generate output
l = unique(opts.labels(:));
c = cell(numel(l), 1);

% iterate over labels
for lc = 1:numel(l)

    % get internal coordinates
    [c1, c2, c3] = ind2sub(vsz, find(vol == l(lc)));

    % create TAL coordinates
    cc = [c1(:), c2(:), c3(:), ones(numel(c1), 1)];
    cc = cc * ftrf;

    % store TAL?
    if sp == 'l'
        c{lc} = cc(:, 1:3);

    % store BVS
    elseif sp == 's'
        c{lc} = min(255, max(0, 129 - round(cc(:, 1:3))));

    % store BVI
    else
        c{lc} = min(255, max(0, 129 - round(cc(:, [2, 3, 1]))));
    end
end

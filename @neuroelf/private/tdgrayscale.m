function tg = tdgrayscale(opts)
% tdgrayscale  - create a gray-scale version of the TD daemon atlas
%
% FORMAT:       tg = tdgrayscale([opts])
%
% Input fields:
%
%       opts        optional settings
%        .csf       CSF value (default: 40)
%        .gm        gray matter value (default: 125)
%        .smooth    post-smoothing kernel
%        .wm        white matter value (default: 175)
%
% Output fields:
%
%       tg          NII object with gray-scale version of TD daemon atlas

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:31 AM EST
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

% argument
if nargin < 1 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'csf') || ...
   ~isa(opts.csf, 'double') || ...
    numel(opts.csf) ~= 1 || ...
    isinf(opts.csf) || ...
    isnan(opts.csf)
    opts.csf = 40;
else
    opts.csf = min(255, max(0, opts.csf));
end
if ~isfield(opts, 'gm') || ...
   ~isa(opts.gm, 'double') || ...
    numel(opts.gm) ~= 1 || ...
    isinf(opts.gm) || ...
    isnan(opts.gm)
    opts.gm = 125;
else
    opts.gm = min(255, max(0, opts.gm));
end
if ~isfield(opts, 'smooth') || ...
   ~isa(opts.smooth, 'double') || ...
    numel(opts.smooth) ~= 1 || ...
    isinf(opts.smooth) || ...
    isnan(opts.smooth) || ...
    opts.smooth < 0
    opts.smooth = 0;
else
    opts.smooth = min(12, opts.smooth);
end
if ~isfield(opts, 'wm') || ...
   ~isa(opts.wm, 'double') || ...
    numel(opts.wm) ~= 1 || ...
    isinf(opts.wm) || ...
    isnan(opts.wm)
    opts.wm = 175;
else
    opts.wm = min(255, max(0, opts.wm));
end

% load file
try
    t = xff(neuroelf_file('t', 'talairach.nii'));
    tg = t.CopyObject;
    t.ClearObject;
catch ne_eo;
    rethrow(ne_eo);
end

% load data
try
    tg.LoadVoxelData;
catch ne_eo;
    tg.ClearObject;
    rethrow(ne_eo);
end

% get labels;
tl = char(tg.IntermedData);
tl = splittocellc(tl(13:end), char([10, 13]), true, true);
tl(1) = [];
tl(end) = [];

% build lookup index
idx = opts.csf .* ones(numel(tl), 1);

% find gray matter
gm = ~cellfun('isempty', regexpi(tl, 'gr[ae]y matter'));
idx(gm) = opts.gm;
wm = ~cellfun('isempty', regexpi(tl, 'white matter'));
wm = wm | ~cellfun('isempty', regexpi(tl, 'brainstem'));
idx(wm) = opts.wm;
mb = ~cellfun('isempty', regexpi(tl, 'midbrain'));
idx(mb) = 0.5 * (opts.gm + opts.wm);

% set into voxeldata
tg.VoxelData(tg.VoxelData ~= 0) = idx(tg.VoxelData(tg.VoxelData ~= 0));

% set new scaling window
tg.SetScalingWindow([0, 255], true);

% smooth data?
if opts.smooth >= 0.3
    tg.SmoothData3D(opts.smooth);
end

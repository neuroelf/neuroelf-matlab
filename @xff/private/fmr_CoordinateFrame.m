function fmrc = fmr_CoordinateFrame(xo, opts)
% FMR::CoordinateFrame  - generate coordinate matrices of FMR
%
% FORMAT:       fmrc = fmr.CoordinateFrame([opts]);
%
% Input fields:
%
%       opts        optional struct with settings
%        .dirs      compute directions (default true)
%        .first     compute coordinate of first voxel (default false)
%        .ndgrid    run ndgrid with FMR resolution (default false)
%        .origin    compute origin (center of STC data, default true)
%        .trans     compute transformation matrix (default true)
%        .transinv  compute inverse matrix (default false)

% Version:  v1.1
% Build:    16021412
% Date:     Feb-14 2016, 12:58 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'dmr', 'fmr'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dirs') || ~islogical(opts.dirs) || isempty(opts.dirs)
    opts.dirs = true;
else
    opts.dirs = opts.dirs(1);
end
if ~isfield(opts, 'first') || ~islogical(opts.first) || isempty(opts.first)
    opts.first = false;
else
    opts.first = opts.first(1);
end
if ~isfield(opts, 'ndgrid') || ~islogical(opts.ndgrid) || isempty(opts.ndgrid)
    opts.ndgrid = false;
else
    opts.ndgrid = opts.ndgrid(1);
end
if ~isfield(opts, 'origin') || ~islogical(opts.origin) || isempty(opts.origin)
    opts.origin = true;
else
    opts.origin = opts.origin(1);
end
if ~isfield(opts, 'trans') || ~islogical(opts.trans) || isempty(opts.trans)
    opts.trans = true;
else
    opts.trans = opts.trans(1);
end
if ~isfield(opts, 'transinv') || ~islogical(opts.transinv) || isempty(opts.transinv)
    opts.transinv = false;
else
    opts.transinv = opts.transinv(1);
end

% get resolution
rsx = bc.InplaneResolutionX;
rsy = bc.InplaneResolutionY;
rsz = bc.SliceThickness + bc.GapThickness;
dmx = bc.ResolutionX;
dmy = bc.ResolutionY;
dmz = bc.NrOfSlices;

% get slice centers
cs1 = [bc.Slice1CenterX; bc.Slice1CenterY; bc.Slice1CenterZ];
csn = [bc.SliceNCenterX; bc.SliceNCenterY; bc.SliceNCenterZ];
css = 0.5 .* (cs1 + csn);
drs = csn - cs1;
if all(drs == 0)
    if dmz > 1
        cs1 = [0; 0; 0.5 * rsz * (1 - dmz)];
        csn = [0; 0; 0.5 * rsz * (dmz - 1)];
    else
        cs1 = [0; 0; -0.5 * rsz];
        csn = [0; 0; 0.5 * rsz];
    end
    css = 0.5 .* (cs1 + csn);
    drs = csn - cs1;
end
drs = drs ./ sqrt(sum(drs .* drs));

% get directions
drx = [bc.RowDirX; bc.RowDirY; bc.RowDirZ];
dry = [bc.ColDirX; bc.ColDirY; bc.ColDirZ];
drx = drx ./ sqrt(sum(drx .* drx));
dry = dry ./ sqrt(sum(dry .* dry));
if any(isinf(drx) | isnan(drx))
    if any(isinf(dry) | isnan(dry))
        drx = [1; 0; 0];
        dry = [0; 1; 0];
    else
        drx = cross(dry, drs);
        drx = drx ./ sqrt(sum(drx .* drx));
    end
elseif any(isinf(dry) | isnan(dry))
    dry = cross(drs, drx);
    dry = dry ./ sqrt(sum(dry .* dry));
end
drz = cross(drx, dry);
if sum(drs .* drz) < 0
    drz = -drz;
end
if abs(sum(drs .* drz) - 1) > 0.01
    error('neuroelf:xff:invalidObject', 'Fields Slice*Center* and ***Dir* don''t match.');
end
if abs(sum(drs .* ((csn - cs1) ./ (rsz * (dmz - 1)))) - 1) > (dmz * 0.01)
    error('neuroelf:xff:invalidObject', ...
        'Fields NrOfSlices, Slice/GapThickness, and directions mismatch.');
end

% compute origin
fst = cs1 - (((dmx + 1) * rsx / 2) * drx + ((dmy + 1) * rsy / 2) * dry + rsz .* drz);

% compute transformation matrix
trf = [rsx, 0, 0, fst(1); 0, rsy, 0, fst(2); 0, 0, rsz, fst(3); 0, 0, 0, 1];

% get direction components
dcp = -trf(1:3, 1:3);
dcp = dcp ./ repmat(sqrt(sum(dcp .* dcp, 2)), [1, 3]);
dcp(4, :) = -cross(dcp(1, :), dcp(2, :));

% build output
fmrc = struct;
fmrc.DimX = dmx;
fmrc.DimY = dmy;
fmrc.DimZ = dmz;
fmrc.DimT = bc.NrOfVolumes;
fmrc.Dimensions = [dmx, dmy, dmz, bc.NrOfVolumes];
fmrc.ResX = rsx;
fmrc.ResY = rsy;
fmrc.ResZ = rsz;
fmrc.Resolution = [rsx, rsy, rsz];
fmrc.Slice1Center = cs1';
fmrc.SliceNCenter = csn';
fmrc.RowDir = dcp(1, :);
fmrc.ColDir = dcp(2, :);
fmrc.SlcDir = dcp(3, :);
fmrc.IsRadiological = (sum(dcp(3, :) .* dcp(4, :)) > 0);
fmrc.Trf = trf;

% set directions
if opts.dirs
    fmrc.DirX = drx;
    fmrc.DirY = dry;
    fmrc.DirZ = drz;
end

% set origin
if opts.first
    fmrc.First = fst;
end

% ndgrid
if opts.ndgrid
    [fmrc.NDgrid{1:3}] = ndgrid(0:(rsx - 1), 0:(rsy - 1), 0:(rsz - 1));
end

% inverse transformation
if opts.transinv
    fmrc.InverseTrans = inv(trf);
end

% set origin
if opts.origin
    fmrc.Origin = css;
end

% set trans
if opts.trans
    fmrc.Trans = trf;
end

function vmrc = vmr_CoordinateFrame(xo, opts)
% VMR::CoordinateFrame  - generate coordinate matrices of VMR
%
% FORMAT:       vmrc = vmr.CoordinateFrame([opts]);
%
% Input fields:
%
%       opts        optional struct with settings
%        .dirs      compute directions (default true)
%        .first     compute coordinate of first voxel (default false)
%        .ndgrid    run ndgrid with VMR resolution (default false)
%        .offsize   get start/end pos (offset + size, default true)
%        .origin    compute origin (center of STC data, default true)
%        .trans     compute transformation matrix (default true)
%        .transinv  compute inverse matrix (default false)

% Version:  v1.1
% Build:    16021110
% Date:     Feb-11 2016, 10:12 AM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr')
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
if ~isfield(opts, 'offsize') || ~islogical(opts.offsize) || isempty(opts.offsize)
    opts.offsize = true;
else
    opts.offsize = opts.offsize(1);
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
rsx = bc.VoxResX;
rsy = bc.VoxResY;
rsz = bc.VoxResZ;
dmx = bc.DimX;
dmy = bc.DimY;
dmz = bc.DimZ;

% get slice centers
cs1 = [bc.Slice1CenterX; bc.Slice1CenterY; bc.Slice1CenterZ];
csn = [bc.SliceNCenterX; bc.SliceNCenterY; bc.SliceNCenterZ];
css = mean([cs1, csn], 2);
drs = csn - cs1;
drs = drs ./ sqrt(sum(drs .* drs));

% get directions
drx = [bc.RowDirX; bc.RowDirY; bc.RowDirZ];
dry = [bc.ColDirX; bc.ColDirY; bc.ColDirZ];
drx = drx ./ sqrt(sum(drx .* drx));
dry = dry ./ sqrt(sum(dry .* dry));
drz = -cross(drx, dry);
cnv = 1;
if sum(drs .* drz) < 0
    drz = -drz;
    cnv = 0;
end
if abs(sum(drs .* drz) - 1) > 0.01
    error('neuroelf:xff:invalidObject', 'Fields Slice*Center* and ***Dir* don''t match.');
end
trf = [[drx, dry, drz, css]; [0, 0, 0, 1]];

% build output
vmrc = struct;

% set resolution
vmrc.DimX = dmx;
vmrc.DimY = dmy;
vmrc.DimZ = dmz;

% set directions
if opts.dirs
    vmrc.DirC = cnv;
    vmrc.DirX = drx;
    vmrc.DirY = dry;
    vmrc.DirZ = drz;
end

% compute origin
if opts.first
    fst = cs1 - ((rsx - 1) * dmx / 2) * drx + ((rsy - 1) * dmy / 2) * dry;
    vmrc.First = fst;
end

% ndgrid
if opts.ndgrid
    [opts.NDgrid{1:3}] = ndgrid(0:(rsx - 1), 0:(rsy - 1), 0:(rsz - 1));
end

% inverse transformation
if opts.transinv
    vmrc.InverseTrans = inv(trf);
end

% get offset and size
if opts.offsize
    if isfield(bc, 'OffsetX')
        vmro = [bc.OffsetX, bc.OffsetY, bc.OffsetZ];
    else
        vmro = [0, 0, 0];
    end
    vmrc.OffSize = [vmro; vmro + [dmx, dmy, dmz] - 1];
end

% set origin
if opts.origin
    vmrc.Origin = css;
end

% set resolution
vmrc.ResX = rsx;
vmrc.ResY = rsy;
vmrc.ResZ = rsz;

% set trans
if opts.trans
    vmrc.Trans = trf;
end

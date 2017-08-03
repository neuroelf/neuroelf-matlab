function varargout = ne_spmmip(varargin)
%ne_spmmip  SPM8-based MIP (maximum intensity projection) of results

% SPM maximum intensity projection
% FORMAT mip = spm_mip(Z,XYZ,M,units)
% Z       - vector point list of SPM values for MIP
% XYZ     - matrix of coordinates of points (mip coordinates)
% M       - voxels - > mip matrix or size of voxels (mm)
% units   - defining space     [default {'mm' 'mm' 'mm'}]
%         - Scalar specifies intensity of grid
%
% mip     - maximum intensity projection
%           if no output, the mip is displayed in current figure.
%__________________________________________________________________________
%
% If the data are 2 dimensional [DIM(3) = 1] the projection is simply an
% image, otherwise:
%
% spm_mip creates and displays a maximum intensity projection of a point
% list of voxel values (Z) and their location (XYZ) in three orthogonal
% views of the brain.  It is assumed voxel locations conform to the space
% defined in the atlas of Talairach and Tournoux (1988); unless the third
% dimension is time.
%
% This routine loads a mip outline from MIP.mat. This is an image with
% contours and grids defining the space of Talairach & Tournoux (1988).
% mip95 corresponds to the Talairach atlas, mip96 to the MNI templates.
% The outline and grid are superimposed at intensity 0.4.
%
% A customised mip outline can be used instead of the default.
%
% A default colormap of 64 levels is assumed. The pointlist image is
% scaled to fit in the interval [0.25,1]*64 for display. Flat images
% are scaled to 1*64.
%
% If M is not specified, it is assumed the XYZ locations are
% in Talairach mm.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston et al.
% $Id: spm_mip.m 6087 2014-07-03 16:14:31Z guillaume $

%-Get units and grid scaling
%--------------------------------------------------------------------------

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% stats var
tail = 1;
if nargin < 4 || numel(varargin{3}) ~= 1 || ~isxff(varargin{3}, 'vmp') || ...
   ~isa(varargin{4}, 'double') || numel(varargin{4}) ~= 1 || ...
    isinf(varargin{4}) || isnan(varargin{4}) || varargin{4} < 1 || varargin{4} > numel(varargin{3}.Map)
    stvar = cc.StatsVar;
    stvix = ch.StatsVarMaps.Value;
    if numel(stvar) ~= 1 || ~isxff(stvar, 'vmp') || numel(stvix) ~= 1
        return;
    end
    if nargin > 2 && ischar(varargin{3}) && strcmp(varargin{3}, '-')
        tail = -1;
    end
else
    stvar = varargin{3};
    stvix = round(varargin{4});
end

% cluster VMP if not yet done
map = stvar.Map(stvix);
if map.EnableClusterCheck > 0 && isempty(map.VMPDataCT)
    stvar.ClusterTable(stvix, []);
    map = stvar.Map(stvix);
end

% find supra-threshold voxels
mapdata = double(map.VMPData);
if map.EnableClusterCheck
    mapdata = mapdata .* double(map.VMPDataCT);
end
if tail < 0
    mapdata = -mapdata;
    tailname = ', negative tail';
else
    tailname = '';
end
statthresh = map.LowerThreshold;
mapvox = find(mapdata > statthresh);
if isempty(mapvox)
    return;
end
switch map.Type
    case 1
        stattype = 't';
    case 2
        stattype = 'r';
    case 3
        stattype = 'CCr';
    case 4
        stattype = 'F';
    otherwise
        stattype = 'raw-value';
end
Z = mapdata(mapvox);
XYZ = bvcoordconv(mapvox, 'bvx2tal', stvar.BoundingBox)';
M = stvar.Resolution .* eye(4);
Grid = 0.4;

%-Scale & offset point list values to fit in [0.25,1]
%==========================================================================
Z = Z - min(Z);
mx = max(Z);
Scal = 8;
Z = (1 + Scal .* Z ./ mx) ./ (Scal + 1);

%-Display format
%==========================================================================
% load various grids, DXYZ, CXYZ, scale (see spm_mip_ui and spm_project)
mipmat = neuroelf_file('p', 'spm8_MIP');
mip = 4 .* mipmat.mip.grid_all + mipmat.mip.mask_all;

% Create maximum intensity projection
%--------------------------------------------------------------------------
mip  = mip ./ max(mip(:));
c = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1] - 0.5;
c = c * M(1:3,1:3);
dim = [(max(c) - min(c)) size(mip)];
try
    d = spm_project(Z,round(XYZ),dim,mipmat.mip.DXYZ,mipmat.mip.CXYZ);
catch ne_eo;
    rethrow(ne_eo);
end
mip = max(d, Grid .* mip);
mip = rot90((1 - mip) .* 64);

%-And display it
%--------------------------------------------------------------------------
f = figure;
varargout{1} = f;
set(f, 'Name', sprintf('SPM8-based MIP: %s (%s > %.3f, one-tailed%s)', ...
    map.Name, stattype, statthresh, tailname), 'NumberTitle', 'off');
set(varargout{1}, 'Units', 'normalized');
ax = axes;
set(ax, 'Position', [0, 0, 1, 1]);
set(varargout{1}, 'Units', 'pixels');
image(mip); axis tight; axis off;
colormap gray

% set keyboard shortcut for screenshot function
set(varargout{1}, 'KeyPressFcn', @ne_spmmip_keypress);



% key-press function
function varargout = ne_spmmip_keypress(src, ke, varargin)
global ne_gcfg;
varargout = cell(1, nargout);
if ~ishandle(src) || ~strcmpi(get(src, 'Type'), 'figure')
    return;
end
kk = ke.Key;
mn = ke.Modifier;
km = false(1, 4);
if ~isempty(mn)
    try
        km = [any(strcmpi('alt', mn)), any(strcmpi('control', mn)), ...
            any(strcmpi('shift', mn)), any(strcmpi('command', mn))];
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% shift-S (screenshot)
if km(3) && strcmpi(kk, 's')
    ne_screenshot(0, 0, src, '', 'high-q');
end

% PUBLIC FUNCTION ne_limitclusters: limit clusters to shape
function varargout = ne_limitclusters(varargin)

% Version:  v1.1
% Build:    20031215
% Date:     Mar-12 2020, 3:06 PM EST
% Author:   Jochen Weber, NeuroElf.net
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2020, Jochen Weber
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get list of clusters
clidx = ch.Clusters.Value(:);

% for an empty list
if isempty(clidx)

    % get voi object handle
    voi = ne_gcfg.voi;

    % ask for a sphere center and radius to add a spherical VOI
    cpos = ne_gcfg.fcfg.cpos;
    srad = ne_gcfg.c.ini.Statistics.SphereSize;
    sphclus = inputdlg({'Sphere center', 'Sphere radius (mm)', 'Name'}, ...
        'NeuroElf - user input', 1, ...
        {sprintf('  %d', cpos), sprintf('  %g', srad), ...
        ['Sphere' sprintf('_%d', cpos) sprintf('_%gmm', srad)]});
    if ~iscell(sphclus) || numel(sphclus) ~= 3 || ...
        isempty(sphclus{1}) || ~ischar(sphclus{1}) || ...
        isempty(sphclus{2}) || ~ischar(sphclus{2}) || ...
        isempty(sphclus{3}) || ~ischar(sphclus{3})
        return;
    end

    % try to convert values
    try
       sphcent = u8str2double(sphclus{1});
       sphrad = str2double(ddeblank(sphclus{2}));
       if ~isa(sphcent, 'double') || numel(sphcent) ~= 3 || ...
           any(isinf(sphcent) | isnan(sphcent) | sphcent < -128 | sphcent > 128) || ...
          ~isa(sphrad, 'double') || numel(sphrad) ~= 1 || ...
           isinf(sphrad) || isnan(sphrad) || sphrad < 0
           return;
       end
       sphrad = min(30, sphrad);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        return;
    end
    ne_gcfg.c.ini.Statistics.SphereSize = sphrad;

    % try to add to current VOI object
    voi.AddSphericalVOI(sphcent, sphrad);
    lastvoi = voi.VOI(end);

    % add string to listbox
    clnames = ch.Clusters.String;
    if ~iscell(clnames)
        clnames = cellstr(clnames);
    end
    clnames{end+1} = sprintf('%s (%d voxels around [%d, %d, %d])', ...
        sphclus{3}, size(lastvoi.Voxels, 1), round(mean(lastvoi.Voxels, 1)));
    ch.Clusters.String = clnames;
    ch.Clusters.Value = [];

    % return already
    return;
end

% try to see if VOI is valid
voi = ne_gcfg.voi;
try
    voi.VOI(clidx);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    warning('neuroelf:GUI:internalError', ...
        'VOI doesn''t contain requested clusters.');
end

% ask for code to use
if nargin < 3 || ~iscell(varargin{3}) || numel(varargin{3}) ~= 2 || ...
   ~isa(varargin{3}{1}, 'double') || numel(varargin{3}{1}) ~= 1 || ...
    isinf(varargin{3}{1}) || isnan(varargin{3}{1}) || varargin{3}{1} < 1 || ...
   ~ischar(varargin{3}{2}) || numel(varargin{3}{2}) ~= 1 || ...
   ~any(lower(varargin{3}{2}) == 'bimsu')
    clim = inputdlg({'Radius:', 'Shape: (b)ox or (s)phere; Function: (i)ntersect, (u)nion, or (m)ap value'}, ...
        'NeuroElf GUI - restrict clusters', 1, ...
        {sprintf('  %d', cc.clim), '  s'});
    if numel(clim) ~= 2 || ~iscell(clim) || isempty(clim{1}) || ~ischar(clim{1}) || ...
        isempty(clim{2}) || ~ischar(clim{2})
        return;
    end
    clim{1}(clim{1} == ' ') = [];
    clim{2}(clim{2} == ' ') = [];
    if isempty(clim{2}) || ~any(lower(clim{2}(1)) == 'bimsu')
        return;
    end
    try
        if clim{1}(1) == '+' || clim{1}(1) == '-'
            climsign = clim{1}(1);
        else
            climsign = ' ';
        end
        clim{1} = str2double(clim{1});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        return;
    end
    if isinf(clim{1}) || isnan(clim{1}) || clim{1} <= 0
        clim{1} = 0;
    end
else
    clim = varargin{3};
end
ne_gcfg.fcfg.clim = clim{1};
if lower(clim{2}(1)) == 'b'
    clim{2} = 'box';
    ccmd = 'restrict';
elseif lower(clim{2}(1)) == 'i'
    clim{2} = '';
    ccmd = 'intersect';
elseif lower(clim{2}(1)) == 'u'
    clim{2} = '';
    ccmd = 'union';
else
    clim{2} = 'sphere';
    ccmd = 'restrict';
end

% old VOI names
onam = ch.Clusters.String;
if ~iscell(onam)
    onam = cellstr(onam);
end

% get current number of clusters
ncl = numel(voi.VOI);

% limit VOIs
opts = struct('rinplace', false, 'rshape', clim{2}, 'rsize', clim{1});
voi.Combine(clidx, ccmd, opts);

% echo
if ne_gcfg.c.echo
    ne_echo('voi', 'Combine', opts);
end

% get new VOIs
nvoi = (ncl+1):numel(voi.VOI);

% intersection
if any(ccmd(1) == 'iu') && numel(nvoi) == 1
    cvox = voi.VOI(end).Voxels;
    if isempty(cvox)
        voi.VOI(end) = [];
        uiwait(warndlg('No voxels remain.', 'NeuroElf - info', 'modal'));
        return;
    end
    ccnt = mean(voi.VOI(end).Voxels);
    vdsq = sum((cvox - ccnt(ones(size(cvox, 1), 1), :)) .^ 2, 2);
    [vdsq, vdsi] = sort(vdsq);
    voi.VOI(end).Voxels = cvox(vdsi, :);
    ccnt = round(ccnt);
    if ccmd(1) == 'i'
        voi.VOI(end).Name = sprintf('ISECT_Cluster_%d_%d_%d_', ccnt);
        onam{end+1} = sprintf('Intersection cluster around [%d,%d,%d] - %d voxel(s)', ...
            ccnt, size(voi.VOI(end).Voxels, 1));
    else
        voi.VOI(end).Name = sprintf('UNION_Cluster_%d_%d_%d_', ccnt);
        onam{end+1} = sprintf('Union cluster around [%d,%d,%d] - %d voxel(s)', ...
            ccnt, size(voi.VOI(end).Voxels, 1));
    end
    ch.Clusters.String = onam;
    ch.Clusters.Value = numel(onam);
    return;
end

% something went wrong...
if numel(nvoi) ~= numel(clidx)
    return;
end

% create new order
[oo, ooi] = sort([1:ncl, clidx(:)' + 0.5]);
voi.VOI = voi.VOI(ooi);

% update names
nnam = cell(numel(voi.VOI), 1);
for clc = 1:numel(nnam)
    if oo(clc) == floor(oo(clc))
        nnam(clc) = onam(oo(clc));
    else
        nnam{clc} = sprintf('%s (rest. %s %.1fmm)', ...
            regexprep(onam{floor(oo(clc))}, '(\s*)\d+\s+voxels?', ...
            sprintf('$1%d voxel', size(voi.VOI(clc).Voxels, 1))), ...
            clim{2}, clim{1});
    end
end
ch.Clusters.String = nnam;
ch.Clusters.Value = clidx(:) + (1:numel(clidx))';

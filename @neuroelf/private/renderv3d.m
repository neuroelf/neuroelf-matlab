function ri  = renderv3d(volume, opts)
% renderv3d  - render a (set of) 3D volume(s) to an image
%
% FORMAT:       ri = renderv3d(vol, opts)
%
% Input fields:
%
%       vol         3D volume of data or 1xV cell array
%       opts        optional settings:
%        .atable    Cx1 alpha-value table, interpolated to match range
%        .backcolor 1x3 background color
%        .ctable    Cx3 color table, also interpolated across value range
%        .imagesize size of the rendered image, defaults to [256, 256]
%        .imax      maximal intensity
%        .imin      minimal intensity
%        .lightvec  light direction (default: [0.67, 0.33, -0.67])
%        .mview     this 4x4 matrix is the viewing matrix, default: eye(4)
%        .rendtype  render type, one of 'bw', 'color', {'mip'}, 'shaded'
%        .warpip    interpolation method used in the warp step
%        .viewvec   view vector (default: [0, 0, 1])
%        .shademat  type of shading (material): 'dull', 'metal', {'shiny'}
%        .update    per-slice update function and arguments ({})
%                   -- additional, optiona parameter to speed up rendering
%        .normals   normalized gradient of the voxel volume; created by
%                   [fy, fx, fz] = gradient(opts.volume);
%                   flength = sqrt(fx .^ 2 + fy .^ 2 + fz .^ 2) + 1e-8;
%                   opts.normals = zeros([size(fx), 3]);
%                   opts.normals(:, :, :, 1) = fx ./ flength;
%                   opts.normals(:, :, :, 2) = fy ./ flength;
%                   opts.normals(:, :, :, 3) = fz ./ flength;
%
% Output fields:
%
%       ri          rendered image
%
% Example:
%
%   % load data
%   vmr = xff([neuroelf_path('colin') '/colin_brain_ICBMnorm.vmr']);
%   v = double(vmr.VMRData);
%   vmr.ClearObject;
%
%   % type of rendering
%   opts.rendtype = 'shaded';
%
%   % color and alpha table
%   opts.atable = [0, 0.1, 0.16, 0.25, 0.8, 0.95, 1];
%   opts.ctable = opts.atable * ones(1, 3);
%
%   % viewer matrix
%   opts.mview = renderv3dvm([0, 0, 0], [0.25, 0.25, 0.25], [0, 0, 0]);
%
%   % render and show image
%   figure;
%   im = renderv3d(v, opts);
%   image(im);
%
% Function is written by D.Kroon University of Twente (April 2009)

% Version:  v0.9d
% Build:    14062716
% Date:     Jun-27 2014, 4:09 PM EST
% Author:   Dirk-Jan Kroon, University of Twente, Enschede, NL
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://www.mathworks.com/matlabcentral/fileexchange/21993-viewer3d

% Copyright (c) 2008, 2014, Dirk-Jan Kroon
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
if nargin < 1 || ...
    isempty(volume) || ...
   (~isnumeric(volume) && ...
    (~iscell(volume) || ...
     ~isnumeric(volume{1}) || ...
      isempty(volume{1})))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if istransio(volume)
    volume = double(volume);
end
if isnumeric(volume)
    volume = {volume};
end
if ndims(volume{1}) < 3 || ...
    size(volume{1}, 4) ~= 1
    error( ...
        'neuroelf:BadArgument', ...
        'Data volume must be 3D.' ...
    );
end
defaultopts = struct( ...
    'afactor',   1, ...
    'atable',    [0.00, 0.00, 0.04, 0.10, 0.25, 0.45, 0.70, 1.00, 1.00], ...
    'avolume',   [], ...
    'backcolor', [0, 0, 0], ...
    'colslice',  false, ...
    'ctable',    [0, 0, 0; 1, 1, 1], ...
    'imagesize', [256, 256], ...
    'imax',      [], ...
    'imin',      [], ...
    'lightvec',  [0.67, 0.33, -0.67], ...
    'material',  [0.7, 0.6, 0.9, 20], ...
    'mview',     eye(4), ...
    'normals',   [], ...
    'rendtype',  'mip', ...
    'shademat',  'shiny', ...
    'slicesel',  1, ...
    'update',    {{}}, ...
    'viewvec',   [0, 0, 1], ...
    'volume',    zeros(3,3,3), ...
    'warpip',    'linear', ...
    'xfromy',    [], ...
    'xtoy',      [], ...
    'xfromz',    [], ...
    'xtoz',      [], ...
    'yfromx',    [], ...
    'ytox',      [], ...
    'yfromz',    [], ...
    'ytoz',      [], ...
    'zfromx',    [], ...
    'ztox',      [], ...
    'zfromy',    [], ...
    'ztoy',      []);
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = defaultopts;
else
    tags = fieldnames(defaultopts);
    for fc = 1:numel(tags)
         if ~isfield(opts, tags{fc})
             opts.(tags{fc}) = defaultopts.(tags{fc});
         end
    end
    if numel(tags) ~= numel(fieldnames(opts))
        warning( ...
            'neuroelf:UnknownOption', ...
            'Unknown fields in opts struct found.' ...
        );
    end
end

% make the data structure from the opts structure
data = opts;
data.volume = volume;
nvols = numel(volume);
if ~iscell(data.avolume) || ...
    isempty(data.avolume)
    data.avolume = cell(1, nvols);
elseif numel(data.avolume) < nvols
    if numel(data.avolume) == 1
        data.avolume = repmat(data.avolume, nvols, 1);
    else
        data.avolume = data.avolume(:);
        data.avolume(end+1:nvols) = data.avolume(1);
    end
end
if isempty(data.afactor) || ...
   ~isnumeric(data.afactor)
    data.afactor = 1;
end

% end points
dsize = size(data.volume{1});
if numel(dsize) > 3
    dsize(4:end) = [];
end
dsiz3 = size(data.volume{1}, 3);
if numel(data.xfromy) ~= dsize(1)
    data.xfromy = ones(dsize(1), 1);
end
if numel(data.xtoy) ~= dsize(1)
    data.xtoy = dsize(2) .* ones(dsize(1), 1);
end
if numel(data.xfromz) ~= dsize(1)
    data.xfromz = ones(dsize(1), 1);
end
if numel(data.xtoz) ~= dsize(1)
    data.xtoz = dsiz3 .* ones(dsize(1), 1);
end
if numel(data.yfromx) ~= dsize(2)
    data.yfromx = ones(dsize(2), 1);
end
if numel(data.ytox) ~= dsize(2)
    data.ytox = dsize(1) .* ones(dsize(2), 1);
end
if numel(data.yfromz) ~= dsize(2)
    data.yfromz = ones(dsize(2), 1);
end
if numel(data.ytoz) ~= dsize(2)
    data.ytoz = dsiz3 .* ones(dsize(2), 1);
end
if numel(data.zfromx) ~= dsiz3
    data.zfromx = ones(dsiz3, 1);
end
if numel(data.ztox) ~= dsiz3
    data.ztox = dsize(1) .* ones(dsiz3, 1);
end
if numel(data.zfromy) ~= dsiz3
    data.zfromy = ones(dsiz3, 1);
end
if numel(data.ztoy) ~= dsiz3
    data.ztoy = dsize(2) .* ones(dsiz3, 1);
end

% alpha, colors, and material per volume
if ~iscell(data.atable)
    data.atable = repmat({data.atable}, [1, nvols]);
end
if ~iscell(data.ctable)
    data.ctable = repmat({data.ctable}, [1, nvols]);
end
if ~iscell(data.material)
    data.material = repmat({data.material}, [1, nvols]);
end

% size check
if numel(data.atable) ~= nvols || ...
    numel(data.ctable) ~= nvols || ...
    numel(data.material) ~= nvols
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid size for alpha, color, or material setting.' ...
    );
end

% settings
for vc = 1:nvols

    % size must match
    tsize = size(data.volume{vc});
    if ~isequal(dsize, tsize(1:min(3, numel(tsize)))) || ...
        size(data.volume{vc}, 4) ~= 1
        error( ...
            'neuroelf:BadData', ...
            'Multiple volumes must match in size.' ...
        );
    elseif ~isempty(data.avolume{vc}) && ...
       ~isequal(dsize, size(data.avolume{vc}))
        error( ...
            'neuroelf:BadData', ...
            'Alpha volumes must match in size.' ...
        );
    end

    % min/max/range
    if numel(data.imax) < numel(data.volume) || ...
        numel(data.imin) < numel(data.volume)
        mmm = minmaxmean(data.volume{vc}, 4);
        if numel(data.imax) < numel(data.volume)
            data.imax(vc) = mmm(2);
        end
        if numel(data.imin) < numel(data.volume)
            data.imin(vc) = mmm(1);
        end
    end
    data.imaxmin(vc) = 1e-8 + data.imax(vc) - data.imin(vc);

    % split ctable into R, G, and B
    if ~isempty(data.ctable{vc})
        if size(data.ctable{vc}, 1) < 16
            if size(data.ctable{vc}, 1) < 2
                data.ctable{vc} = [0, 0, 0; 1, 1, 1];
            end
            data.ctable{vc} = flexinterpn_method(data.ctable{vc}, ...
            [Inf, Inf; 1, 1; max(1, size(data.ctable{vc}, 1) - 1) / 255, 1; ...
            size(data.ctable{vc}, 1), 3], 'linear');
        end
        data.ctable{vc}(data.ctable{vc} < 0) = 0;
        if any(data.ctable{vc}(:) > 2)
            data.ctable{vc} = (1 / 255) .* data.ctable{vc};
        end
        data.ctable{vc}(data.ctable{vc} > 1) = 1;
        data.ctable_r{vc} = data.ctable{vc}(:, 1);
        data.ctable_g{vc} = data.ctable{vc}(:, 2);
        data.ctable_b{vc} = data.ctable{vc}(:, 3);
    end

    % shading type, set phong values (for all volumes)
    switch lower(data.shademat)
        case {'shiny'}
            data.material{vc} = [0.7, 0.6, 0.9, 15];
        case {'dull'}
            data.material{vc} = [0.7, 0.8, 0.0, 10];
        case {'metal'}
            data.material{vc} = [0.7, 0.3, 1.0, 20];
    end
end

% calculate the shear and warp matrices
sizes = size(data.volume{1});
[data.Mshear, data.Mwarp2D, data.c] = makeShearWarpMatrix(data.mview, sizes(1:3));
data.Mwarp2Dinv = inv(double(data.Mwarp2D));
data.Mshearinv = inv(data.Mshear);

% store volume sizes as separate values
data.Iin_sizex = sizes(1);
data.Iin_sizey = sizes(2);
data.Iin_sizez = sizes(3);

% create shear (intimidate) buffer
data.Ibuffer_sizex = 2 * max(sizes);
data.Ibuffer_sizey = data.Ibuffer_sizex;
if any(strcmpi(data.rendtype, {'bw', 'mip'}))
    data.Ibuffer = zeros([data.Ibuffer_sizex, data.Ibuffer_sizey]);
else
    data.Ibuffer = ...
        reshape(ones(data.Ibuffer_sizex * data.Ibuffer_sizey, 1) * data.backcolor, ...
        [data.Ibuffer_sizex, data.Ibuffer_sizey, 3]);
end

% normalize light and viewer vectors
data.lightvec = [data.lightvec(:); 0];
data.lightvec = data.lightvec ./ sqrt(sum(data.lightvec(1:3) .^ 2));
data.viewvec = [data.viewvec(:); 0];
data.viewvec = data.viewvec ./ sqrt(sum(data.viewvec(1:3) .^ 2));

% adjust alpha table by voxel length because of rotation and volume size
lengthcor = 0.01 * sqrt(1 + data.Mshearinv(1, 3) ^ 2 + data.Mshearinv(2, 3) ^ 2) * mean(sizes(1:3));
for vc = 1:nvols
    if numel(data.atable{vc}) < 32
        if numel(data.atable{vc}) < 2
            data.atable{vc} = [0; 1];
        end
        data.atable{vc} = flexinterpn_method(data.atable{vc}(:), ...
            [Inf; 1; max(numel(data.atable{vc}) - 1, 1) / 255; numel(data.atable{vc})], 'linear');
    end
    data.atable{vc} = 1 - (1 - data.atable{vc}(:)) .^ (1 / lengthcor);
    data.atable{vc}(data.atable{vc} < 0) = 0;
    data.atable{vc}(data.atable{vc} > 1) = 1;
end

% pre-set fields
data.afactor = min(1, max(0, data.afactor(1)));
data.alpha_loc = [];
data.intensity_loc = [];
data.intensity_rgb = [];

% shear Rendering
data = shear(data);
data = warp(data);
ri = data.Iout;
if any(strcmpi(data.rendtype, {'color', 'shaded'}))
    ri = uint8(round(255.999 .* ri - 0.49995));
end


% internal sub-functions


% shear function
function data = shear(data)

% depending on the main viewing direction
switch (data.c)

    % along first axis +
    case 1

        % iterate over slices
        for z = 0:(data.Iin_sizex - 1)

            % skip slice
            if isnan(data.xfromy(z+1))
                continue;
            end

            % offset calculation
            xd = (-data.Ibuffer_sizex / 2) + data.Mshearinv(1, 3) * (z - data.Iin_sizex / 2) + data.Iin_sizey / 2;
            yd = (-data.Ibuffer_sizey / 2) + data.Mshearinv(2, 3) * (z - data.Iin_sizex / 2) + data.Iin_sizez / 2;
            xdfloor = floor(xd);
            ydfloor = floor(yd);

            % calculate the coordinates on which a image slice starts and
            % ends in the temporary shear image (buffer)
            pxstart = max(data.xfromy(z+1) - xdfloor, -xdfloor);
            pystart = max(data.xfromz(z+1) - ydfloor, -ydfloor);
            pxend = min(data.Ibuffer_sizex, data.xtoy(z+1) - xdfloor);
            pyend = min(data.Ibuffer_sizey, data.xtoz(z+1) - ydfloor);
            data.px = (pxstart+1:pxend-1);
            data.py = (pystart+1:pyend-1);
            if isempty(data.px) || ...
                isempty(data.py)
                continue;
            end

            % determine x and y coordinates of pixel(s) which will be come current pixel
            xBas = data.px + xdfloor;
            yBas = data.py + ydfloor;

            % linear interpolation constants (percentages)
            xCom = xd - floor(xd);
            yCom = yd - floor(yd);
            perc = [(1 - xCom) * (1 - yCom), (1 - xCom) * yCom, xCom * (1 - yCom), xCom * yCom];

            % get the intensities
            for vc = 1:numel(data.volume)
                slice = double(squeeze(data.volume{vc}(z+1, :, :, :, :)));
                if isempty(data.avolume{vc}) || ...
                    data.rendtype(1) == 'm'
                    aslice = [];
                else
                    aslice = double(squeeze(data.avolume{vc}(z+1, :, :, 1)));
                end
                data = extractfromslice(data, vc, z, slice, aslice, xBas, yBas, perc);
            end
        end

    % along second axis +
    case 2

        % iterate over slices
        for z = 0:(data.Iin_sizey - 1)

            % similar path
            if isnan(data.yfromx(z+1))
                continue;
            end
            xd = (-data.Ibuffer_sizex / 2) + data.Mshearinv(1, 3) * (z - data.Iin_sizey / 2) + data.Iin_sizez / 2;
            yd = (-data.Ibuffer_sizey / 2) + data.Mshearinv(2, 3) * (z - data.Iin_sizey / 2) + data.Iin_sizex / 2;
            xdfloor = floor(xd);
            ydfloor = floor(yd);
            pxstart = max(data.yfromz(z+1) - xdfloor, -xdfloor);
            pystart = max(data.yfromx(z+1) - ydfloor, -ydfloor);
            pxend = min(data.Ibuffer_sizex, data.ytoz(z+1) - xdfloor);
            pyend = min(data.Ibuffer_sizey, data.ytox(z+1) - ydfloor);
            data.px = (pxstart+1:pxend-1);
            data.py = (pystart+1:pyend-1);
            if isempty(data.px) || ...
                isempty(data.py)
                continue;
            end
            xBas = data.px + xdfloor;
            yBas = data.py + ydfloor;
            xCom = xd - floor(xd);
            yCom = yd - floor(yd);
            perc = [(1 - xCom) * (1 - yCom), (1 - xCom) * yCom, xCom * (1 - yCom), xCom * yCom];
            for vc = 1:numel(data.volume)
                slice = double(squeeze(data.volume{vc}(:, z+1, :, :, :)));
                if isempty(data.avolume{vc}) || ...
                    data.rendtype(1) == 'm'
                    aslice = [];
                else
                    aslice = double(squeeze(data.avolume{vc}(:, z+1, :, 1)));
                end
                data = extractfromslice(data, vc, z, slice, aslice, yBas, xBas, perc);
            end
        end

    % third axis +
    case 3

        % iterate over slices
        for z = 0:(data.Iin_sizez - 1)

            % similar path
            if isnan(data.zfromx(z+1))
                continue;
            end
            xd = (-data.Ibuffer_sizex / 2) + data.Mshearinv(1, 3) * (z - data.Iin_sizez / 2) + data.Iin_sizex / 2;
            yd = (-data.Ibuffer_sizey / 2) + data.Mshearinv(2, 3) * (z - data.Iin_sizez / 2) + data.Iin_sizey / 2;
            xdfloor = floor(xd);
            ydfloor = floor(yd);
            pxstart = max(data.zfromx(z+1) - xdfloor, -xdfloor);
            pystart = max(data.zfromy(z+1) - ydfloor, -ydfloor);
            pxend = min(data.Ibuffer_sizex, data.ztox(z+1) - xdfloor);
            pyend = min(data.Ibuffer_sizey, data.ztoy(z+1) - ydfloor);
            data.px = (pxstart+1:pxend-1);
            data.py = (pystart+1:pyend-1);
            if isempty(data.px) || ...
                isempty(data.py)
                continue;
            end
            xBas = data.px + xdfloor;
            yBas = data.py + ydfloor;
            xCom = xd - floor(xd);
            yCom = yd - floor(yd);
            perc = [(1 - xCom) * (1 - yCom), (1 - xCom) * yCom, xCom * (1 - yCom), xCom * yCom];
            for vc = 1:numel(data.volume)
                slice = double(squeeze(data.volume{vc}(:, :, z+1, :, :)));
                if isempty(data.avolume{vc}) || ...
                    data.rendtype(1) == 'm'
                    aslice = [];
                else
                    aslice = double(data.avolume{vc}(:, :, z+1, 1));
                end
                data = extractfromslice(data, vc, z, slice, aslice, xBas, yBas, perc);
            end
        end

    % first axis -
    case 4

        % iterate over slices
        for z = (data.Iin_sizex-1):-1:0

            % similar path
            if isnan(data.xfromy(z+1))
                continue;
            end
            xd = (-data.Ibuffer_sizex / 2) + data.Mshearinv(1, 3) * (z - data.Iin_sizex / 2) + data.Iin_sizey / 2;
            yd = (-data.Ibuffer_sizey / 2) + data.Mshearinv(2, 3) * (z - data.Iin_sizex / 2) + data.Iin_sizez / 2;
            xdfloor = floor(xd);
            ydfloor = floor(yd);
            pxstart = max(data.xfromy(z+1) - xdfloor, -xdfloor);
            pystart = max(data.xfromz(z+1) - ydfloor, -ydfloor);
            pxend = min(data.Ibuffer_sizex, data.xtoy(z+1) - xdfloor);
            pyend = min(data.Ibuffer_sizey, data.xtoz(z+1) - ydfloor);
            data.px = (pxstart+1:pxend-1);
            data.py = (pystart+1:pyend-1);
            if isempty(data.px) || ...
                isempty(data.py)
                continue;
            end
            xBas = data.px + xdfloor;
            yBas = data.py + ydfloor;
            xCom = xd - floor(xd);
            yCom = yd - floor(yd);
            perc = [(1 - xCom) * (1 - yCom), (1 - xCom) * yCom, xCom * (1 - yCom), xCom * yCom];
            for vc = 1:numel(data.volume)
                slice = double(squeeze(data.volume{vc}(z+1, :, :, :, :)));
                if isempty(data.avolume{vc}) || ...
                    data.rendtype(1) == 'm'
                    aslice = [];
                else
                    aslice = double(squeeze(data.avolume{vc}(z+1, :, :, 1)));
                end
                data = extractfromslice(data, vc, z, slice, aslice, xBas, yBas, perc);
            end
        end

    % second axis -
    case 5

        % iterate over slices
        for z = (data.Iin_sizey-1):-1:0

            % similar path
            if isnan(data.yfromx(z+1))
                continue;
            end
            xd = (-data.Ibuffer_sizex / 2) + data.Mshearinv(1, 3) * (z - data.Iin_sizey / 2) + data.Iin_sizez / 2;
            yd = (-data.Ibuffer_sizey / 2) + data.Mshearinv(2, 3) * (z - data.Iin_sizey / 2) + data.Iin_sizex / 2;
            xdfloor = floor(xd);
            ydfloor = floor(yd);
            pxstart = max(data.yfromz(z+1) - xdfloor, -xdfloor);
            pystart = max(data.yfromx(z+1) - ydfloor, -ydfloor);
            pxend = min(data.Ibuffer_sizex, data.ytoz(z+1) - xdfloor);
            pyend = min(data.Ibuffer_sizey, data.ytox(z+1) - ydfloor);
            data.px = (pxstart+1:pxend-1);
            data.py = (pystart+1:pyend-1);
            if isempty(data.px) || ...
                isempty(data.py)
                continue;
            end
            xBas = data.px + xdfloor;
            yBas = data.py + ydfloor;
            xCom = xd - floor(xd);
            yCom = yd - floor(yd);
            perc = [(1 - xCom) * (1 - yCom), (1 - xCom) * yCom, xCom * (1 - yCom), xCom * yCom];
            for vc = 1:numel(data.volume)
                slice = double(squeeze(data.volume{vc}(:, z+1, :, :, :)));
                if isempty(data.avolume{vc}) || ...
                    data.rendtype(1) == 'm'
                    aslice = [];
                else
                    aslice = double(squeeze(data.avolume{vc}(:, z+1, :, 1)));
                end
                data = extractfromslice(data, vc, z, slice, aslice, yBas, xBas, perc);
            end
        end

    % third axis -
    case 6

        % iterate over slices
        for z = (data.Iin_sizez-1):-1:0

            % similar path
            if isnan(data.zfromx(z+1))
                continue;
            end
            xd = (-data.Ibuffer_sizex / 2) + data.Mshearinv(1, 3) * (z - data.Iin_sizez / 2) + data.Iin_sizex / 2;
            yd = (-data.Ibuffer_sizey / 2) + data.Mshearinv(2, 3) * (z - data.Iin_sizez / 2) + data.Iin_sizey / 2;
            xdfloor = floor(xd);
            ydfloor = floor(yd);
            pxstart = max(data.zfromx(z+1) - xdfloor, -xdfloor);
            pystart = max(data.zfromy(z+1) - ydfloor, -ydfloor);
            pxend = min(data.Ibuffer_sizex, data.ztox(z+1) - xdfloor);
            pyend = min(data.Ibuffer_sizey, data.ztoy(z+1) - ydfloor);
            data.px = (pxstart+1:pxend-1);
            data.py = (pystart+1:pyend-1);
            if isempty(data.px) || ...
                isempty(data.py)
                continue;
            end
            xBas = data.px + xdfloor;
            yBas = data.py + ydfloor;
            xCom = xd - floor(xd);
            yCom = yd - floor(yd);
            perc = [(1 - xCom) * (1 - yCom), (1 - xCom) * yCom, xCom * (1 - yCom), xCom * yCom];
            for vc = 1:numel(data.volume)
                slice = double(squeeze(data.volume{vc}(:, :, z+1, :, :)));
                if isempty(data.avolume{vc}) || ...
                    data.rendtype(1) == 'm'
                    aslice = [];
                else
                    aslice = double(data.avolume{vc}(:, :, z+1, 1));
                end
                data = extractfromslice(data, vc, z, slice, aslice, xBas, yBas, perc);
            end
        end
end

% ensure data limits
if strcmp(data.rendtype, 'mip')
    data.Ibuffer = limitrangec(data.Ibuffer, 0, 1, 0);
end


% extract pixel information from slice
function data = extractfromslice(data, vc, z, slice, aslice, x1, y1, perc)

% use MEX file and store accordingly
if size(slice, 3) == 1
    [data.intensity_loc, data.alpha_loc] = renderv3dxia(slice, aslice, x1, y1, perc, data.c);
    data.intensity_rgb = [];
else
    [data.intensity_rgb, data.alpha_loc] = renderv3dxia(slice, aslice, x1, y1, perc, data.c);
end

% update
switch (data.rendtype)
    case {'mip'}
        data = updatebuffer_MIP(data, vc);
    case {'color'}
        renderv3dub(data, vc);
    case {'bw'}
        data = updatebuffer_BW(data, vc);
    case {'shaded'}
        switch (data.c)
            case 1
                data = returnnormal(z+1, x1, y1, data, vc);
            case 2
                data = returnnormal(x1, z+1, y1, data, vc);
            case 3
                data = returnnormal(x1, y1, z+1, data, vc);
            case 4
                data = returnnormal(z+1, x1, y1, data, vc);
            case 5
                data = returnnormal(x1, z+1, y1, data, vc);
            case 6
                data = returnnormal(x1, y1, z+1, data, vc);
        end
        data = updatebuffer_SHADED(data, vc);
end

% update
if ~isempty(data.update) && ...
    iscell(data.update) && ...
    vc == numel(data.volume) && ...
    strcmpi(class(data.update{1}), 'function_handle')
    for uc = 1:size(data.update, 1)
        ua = data.update(uc, :);
        uag = findfirst(~cellfun('isempty', ua(:)), -1);
        ccs = find(cellfun(@ischar, ua));
        if ~isempty(ccs)
            pix = ccs(strcmpi(ua(ccs), '%p'));
        else
            pix = [];
        end
        if uc == 1 && ...
            numel(pix) ~= 1
            return;
        end
        if numel(pix) == 1
            data = warp(data);
            ri = data.Iout;
            if any(strcmpi(data.rendtype, {'color', 'shaded'}))
                ri = uint8(round(255.999 .* ri - 0.49995));
            end
            ua{pix} = ri;
            uag = max(uag, pix);
        end
        feval(ua{1, 1:uag});
    end
end


% compute normals (for shading)
function data = returnnormal(x, y, z, data, vc)
% Calculate the normals for a certain pixel / slice or volume.
% The normals are calculated by normalizing the voxel volume gradient.

% the central pixel positions
x1 = x;
y1 = y;
z1 = z;

% if the gradients are not delivered by the user
if isempty(data.normals)

    % the forward pixel positions
    x2 = x1 + 1;
    y2 = y1 + 1;
    z2 = z1 + 1;

    % everything inside the boundaries
    checkx = (x2 > size(data.volume{vc}, 1));
    checky = (y2 > size(data.volume{vc}, 2));
    checkz = (z2 > size(data.volume{vc}, 3));
    if any(checkx)
        x1(checkx) = x1 - 1;
        x2(checkx) = size(data.volume{vc}, 1);
    end
    if any(checky)
        y1(checky) = y1 - 1;
        y2(checky) = size(data.volume{vc}, 2);
    end
    if any(checkz)
        z1(checkz) = z1 - 1;
        z2(checkz) = size(data.volume{vc}, 3);
    end

    % calculate the forward gradient
    N = cat(4, ...
        double(squeeze(data.volume{vc}(x2, y1, z1) - data.volume{vc}(x1, y1, z1))), ...
        double(squeeze(data.volume{vc}(x1, y2, z1) - data.volume{vc}(x1, y1, z1))), ...
        double(squeeze(data.volume{vc}(x1, y1, z2) - data.volume{vc}(x1, y1, z1))));

    % normalize the gradient data
    nlength = sqrt(N(:, :, :, 1) .^ 2 + N(:, :, :, 2) .^ 2 + N(:, :, :, 3) .^ 2) + 1e-8;
    N = N ./ repmat(nlength, [1, 1, 1, 3]);

% pre-computed normals
else
    N = cat(4, ...
        data.normals{vc}(x1, y1, z1, 1), ...
        data.normals{vc}(x1, y1, z1, 2), ...
        data.normals{vc}(x1, y1, z1, 3));
end

% squeeze
sN = size(N);
N = reshape(N, [sN(1), sN(2), prod(sN(3:end))]);

% rotate the data in case of certain views
if data.c == 2 || ...
    data.c == 5
    N = permute(N, [2, 1, 3]);
end

% return the normals in data
data.N = N;


% update MIP
function data = updatebuffer_MIP(data, vc)

% color volume
if isempty(data.intesity_loc)
    data.intensity_loc = max(data.intensity_rgb, [], 3);
end

% update the current pixel in the shear image buffer
data.intensity_loc = (1 / data.imaxmin(vc)) .* (data.intensity_loc - data.imin(vc));
check = double(data.intensity_loc > data.Ibuffer(data.px, data.py));
data.Ibuffer(data.px, data.py) = ...
    check .* data.intensity_loc + (1 - check) .* data.Ibuffer(data.px, data.py);


% update buffer for gray scale
function data = updatebuffer_BW(data, vc)

% color volume
if isempty(data.intesity_loc)
    data.intensity_loc = ...
        0.212 .* data.intensity_rgb(:, :, 1) + ...
        0.701 .* data.intensity_rgb(:, :, 2) + ...
        0.087 .* data.intensity_rgb(:, :, 3);
end

% calculate index in alpha transparency look up table
if data.imin(vc) ~= 0
    data.intensity_loc = data.intensity_loc - data.imin(vc);
end
if data.imaxmin(vc) ~= 1
    data.intensity_loc = data.intensity_loc ./ data.imaxmin(vc);
end
data.intensity_loc = limitrangec(data.intensity_loc, 0, 1, 0);

% compute index into table
if isempty(data.alpha_loc)
    betaA = numel(data.atable{vc}) - 1;
    indexAlpha = round(data.intensity_loc * betaA) + 1;

    % calculate current alphaimage
    alphaimage = reshape(data.atable{vc}(indexAlpha), size(data.intensity_loc));
else
    alphaimage = data.alpha_loc;
end

% multiply
if data.afactor ~= 1
    alphaimage = data.afactor .* alphaimage;
end

% invert (for addition)
alphaimage_inv = 1 - alphaimage;

% update the current pixel in the shear image buffer
data.Ibuffer(data.px, data.py) = ...
    alphaimage_inv .* data.Ibuffer(data.px, data.py) + alphaimage .* data.intensity_loc;


% update buffer for shaded image
function data = updatebuffer_SHADED(data, vc)

% alpha and color scaling
if data.imin(vc) ~= 0
    data.intensity_loc = data.intensity_loc - data.imin(vc);
end
data.intensity_loc = limitrangec(data.intensity_loc, 0, data.imaxmin(vc), 0);
betaA = (numel(data.atable{vc}) - 1) / data.imaxmin(vc);
betaC = (length(data.ctable_r{vc}) - 1) / data.imaxmin(vc);

% no alpha-information for slice
if isempty(data.alpha_loc)

    % calculate index in alpha transparency look up table
    indexAlpha = round(data.intensity_loc * betaA) + 1;

    % calculate index in color look up table
    if betaA ~= betaC
        indexColor = round(data.intensity_loc * betaC) + 1;
    else
        indexColor = indexAlpha;
    end

    % calculate current alphaimage
    alphaimage = data.atable{vc}(indexAlpha);

% alpha information given
else
    indexColor = round(data.intensity_loc * betaC) + 1;
    alphaimage = data.alpha_loc;
end

% Rotate the light and view vector
data.lightvec2 = data.mview \ data.lightvec;
data.lightvec2 = data.lightvec2 ./ sqrt(sum(data.lightvec2(1:3) .^ 2));
data.viewvec2 = data.mview \ data.viewvec;
data.viewvec2 = data.viewvec2 ./ sqrt(sum(data.viewvec2(1:3) .^ 2));

Ia = 1;
Id = ...
    data.lightvec2(1) .* data.N(:, :, 1) + ...
    data.lightvec2(2) .* data.N(:, :, 2) + ...
    data.lightvec2(3) .* data.N(:, :, 3);

% R = 2.0*dot(N,L)*N - L;
R(:, :, 1) = 2 .* Id .* data.N(:, :, 1) - data.lightvec2(1);
R(:, :, 2) = 2 .* Id .* data.N(:, :, 2) - data.lightvec2(2);
R(:, :, 3) = 2 .* Id .* data.N(:, :, 3) - data.lightvec2(3);

%Is = max(pow(dot(R,V),3),0);
Is = -( ...
    data.viewvec2(1) .* R(:, :, 1) + ...
    data.viewvec2(2) .* R(:, :, 2) + ...
    data.viewvec2(3) .* R(:, :, 3));

% No spectacular highlights on "shadow" part
Is(Id < 0) = 0;

% Specular exponent
Is = Is .^ data.material{vc}(4);

% Phong shading values
Ipar = zeros([size(Id) 2]);
Ipar(:, :, 1) = data.material{vc}(1) .* Ia + data.material{vc}(2) .* Id;
Ipar(:, :, 2) = data.material{vc}(3) .* Is;

% multiply
if data.afactor ~= 1
    alphaimage = data.afactor .* alphaimage;
end

% invert (for addition)
alphaimage_inv = (1 - alphaimage);

% Update the current pixel in the shear image buffer
data.Ibuffer(data.px, data.py, 1) = ...
    alphaimage_inv .* data.Ibuffer(data.px, data.py, 1) + ...
    alphaimage .* (data.ctable_r{vc}(indexColor) .* Ipar(:, :, 1) + Ipar(:, :, 2));
data.Ibuffer(data.px, data.py, 2) = ...
    alphaimage_inv .* data.Ibuffer(data.px, data.py, 2) + ...
    alphaimage .* (data.ctable_g{vc}(indexColor) .* Ipar(:, :, 1) + Ipar(:, :, 2));
data.Ibuffer(data.px, data.py, 3) = ...
    alphaimage_inv .* data.Ibuffer(data.px, data.py, 3) + ...
    alphaimage .* (data.ctable_b{vc}(indexColor) .* Ipar(:, :, 1) + Ipar(:, :, 2));

% warp the shear rendered buffer image
function data = warp(data)

% make affine matrix
M = eye(4);
M(2:3, 2:4) = data.Mwarp2Dinv(1:2, 1:3);
M(2:3, 4) = M(2:3, 4) + data.Mshearinv(1:2, 4);
M(2:3, 4) = M(2:3, 4) + M(2:3, 2:3) * (-0.5 .* (data.imagesize(:) + 1)) + 0.5 .* (1 + ([data.Ibuffer_sizex; data.Ibuffer_sizey]));

% perform the affine transformation
if size(data.Ibuffer, 3) == 1
    data.Iout = reshape(flexinterpn_method(reshape( ...
        data.Ibuffer, [1, data.Ibuffer_sizex, data.Ibuffer_sizey]), ...
        [Inf, Inf, Inf; 1, 1, 1; 1, 1, 1; 1, data.imagesize], ...
        0, M, data.warpip), data.imagesize);
else
    data.Iout = zeros([data.imagesize, 3]);
    for cc = 1:3
        data.Iout(:, :, cc) = reshape(flexinterpn_method(reshape( ...
            data.Ibuffer(:, :, cc), [1, data.Ibuffer_sizex, data.Ibuffer_sizey]), ...
            [Inf, Inf, Inf; 1, 1, 1; 1, 1, 1; 1, data.imagesize], ...
            0, M, data.warpip), data.imagesize);
    end
end

% splits the view matrix in to shear and warp matrices for rendering
function [Mshear, Mwarp2D, c] = makeShearWarpMatrix(mview, sizes)

% find the principal viewing axis
Vo = [mview(1, 2) * mview(2, 3) - mview(2, 2) * mview(1, 3); ...
      mview(2, 1) * mview(1, 3) - mview(1, 1) * mview(2, 3); ...
      mview(1, 1) * mview(2, 2) - mview(2, 1) * mview(1, 2)];
[maxv, c] = max(abs(Vo));

% choose the corresponding Permutation matrix P
switch (c)
    case 1, %yzx
        P=[0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
    case 2, % zxy
        P=[0 0 1 0; 1 0 0 0; 0 1 0 0; 0 0 0 1];
    case 3, % xyz
        P=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
end

% compute the permuted view matrix from mview and P
mview_p = mview /P;

% 180 degrees rotate detection
if mview_p(3, 3) < 0
    c = c + 3;
end

% compute the shear coeficients from the permuted view matrix
Si = (mview_p(2, 2) * mview_p(1, 3) - mview_p(1, 2) * mview_p(2, 3)) / ...
     (mview_p(1, 1) * mview_p(2, 2) - mview_p(2, 1) * mview_p(1, 2));
Sj = (mview_p(1, 1) * mview_p(2, 3) - mview_p(2, 1) * mview_p(1, 3)) / ...
     (mview_p(1, 1) * mview_p(2, 2) - mview_p(2, 1) * mview_p(1, 2));

% compute the translation between the orgins of standard object coordinates
% and intermdiate image coordinates
if c == 1 || ...
    c == 4
    kmax = sizes(1) - 1;
end
if c == 2 || ...
    c == 5
    kmax = sizes(2) - 1;
end
if c == 3 || ...
    c == 6
    kmax = sizes(3) - 1;
end
if Si >= 0 && ...
    Sj >= 0
    Ti = 0;
    Tj = 0;
end
if Si >= 0 && ...
    Sj < 0
    Ti = 0;
    Tj = -Sj * kmax;
end
if Si < 0 && ...
    Sj >= 0
    Ti = -Si * kmax;
    Tj = 0;
end
if Si < 0 && ...
    Sj < 0
    Ti = -Si * kmax;
    Tj = -Sj * kmax;
end

% compute the shear matrix and a 2Dwarp matrix
Mshear = [ 1  0 Si Ti; 0  1 Sj Tj; 0  0  1  0; 0  0  0  1];
Mwarp2D = [mview_p(1, 1), mview_p(1, 2), (mview_p(1, 4) - Ti * mview_p(1, 1) - Tj * mview_p(1, 2)); ...
    mview_p(2, 1), mview_p(2, 2), (mview_p(2, 4) - Ti * mview_p(2, 1) - Tj * mview_p(2, 2)); 0, 0, 1];

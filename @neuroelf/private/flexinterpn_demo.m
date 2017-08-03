function flexinterpn_demo(varargin)
% flexinterpn_demo  - demonstration of flexinterpn
%
% FORMAT:       flexinterpn_demo([options])
%
% Input fields:
%
%       options     list of char options among (with {default})
%         data      either {"colin"} or "rand"
%         disp      either {"image"} or "surf"
%         dtype     a valid numeric datatype/class name of matlab
%                   (matlab only supports single and double !)
%         func      either {"flex"} or "matlab"
%         method    either of "cubic", "lanczos2", "lanczos3", {"linear"},
%                   "nearest"
%         size      either "128", "256", {"512"} or "1024"
%         what      either {"tilt"} or "zoom"
%
%       tiltorzoom  show a "tilt" demo or a "zoom" demo (default: "tilt")
%       ml          if given and true, use Matlab interpn instead
%
% No output fields.

% Version:  v1.1
% Build:    16012321
% Date:     Jan-23 2016, 9:43 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% needed for settings
global xffsngl;

% argument check
data = 'colin';
disp = 'image';
dtype = 'uint8';
func = 'flex';
gaussk = {};
method = 'linear';
size = 512;
what = 'tilt';
for nc = 1:nargin
    if ischar(varargin{nc}) && ...
       ~isempty(varargin{nc})
        if any(strcmpi(varargin{nc}(:)', {'colin', 'rand'}))
            data = lower(varargin{nc}(:)');
            continue;
        end
        if any(strcmpi(varargin{nc}(:)', {'image', 'surf'}))
            disp = lower(varargin{nc}(:)');
            continue;
        end
        if any(strcmpi(varargin{nc}(:)', ...
               {'double', 'single', 'int8', 'int16', 'int32', 'uint8', 'uint16', 'uint32'}))
            dtype = lower(varargin{nc}(:)');
            continue;
        end
        if any(strcmpi(varargin{nc}(:)', {'f', 'flex', 'm', 'matlab', 'ml'}))
            func = lower(varargin{nc}(1));
            continue;
        end
        if any(strcmpi(varargin{nc}(:)', ...
                {'cubic', 'gauss2', 'gauss3', ...
                 'lanczos2', 'lanczos3', 'lanczos4', 'lanczos5', ...
                 'linear', 'nearest'}))
            method = lower(varargin{nc}(:)');
            if strcmp(method(1:5), 'gauss')
                gaussk = {str2double(method(end))};
                method = 'gauss';
            end
            continue;
        end
        if any(strcmpi(varargin{nc}(:)', {'1024', '128', '256', '512'}))
            size = str2double(varargin{nc}(:)');
            continue;
        end
        if any(strcmpi(varargin{nc}(:)', {'tilt', 'zoom'}))
            what = lower(varargin{nc}(:)');
            continue;
        end
    end
end
useimg = false;
if strcmp(disp, 'image')
    useimg = true;
end

% Matlab doesn't have Lanczos interpolation
useml = false;
if func == 'm'
    useml = true;
    if ~strcmp(dtype, 'single')
        dtype = 'single';
    end
    if any(strcmp(method, {'gauss', 'lanczos2', 'lanczos3'}))
        method = 'cubic';
    end
    try
        [xyzc{1:3}] = ndgrid(single(1:256), single(1:256), single(1:256));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        xyzc = {};
        warning( ...
            'neuroelf:OutOfMemory', ...
            'Error creating coordinate grids for Matlab interpn.' ...
        );
    end
end

% load VMR
if strcmp(data, 'colin')
    try
        v = neuroelf_file('c', 'colin.vmr');
        vd = v.VMRData;
        v.ClearObject;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:FileNotFound', ...
            'Could not find/load the COLIN template VMR.' ...
        );
    end
else
    try
        vd = min(max(128 + 32 * randn(256, 256, 256), 0), 255);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:MemoryError', ...
            'Error allocating enough memory for data.' ...
        );
    end
end

% convert to data of type
switch (dtype)
    case {'double'}
        vd = double(vd);
    case {'int16'}
        vd = int16(vd);
    case {'int32'}
        vd = int32(vd);
    case {'int8'}
        vd = int8(0.5 * double(vd));
    case {'single'}
        vd = single(vd);
    case {'uint16'}
        vd = uint16(vd);
    case {'uint32'}
        vd = uint32(vd);
    case {'uint8'}
        vd = uint8(vd);
end

% show blank image
if useimg
    h = scaleimage(zeros(size, size));
    hp = get(h, 'Parent');
    set(get(hp, 'Parent'), ...
        'Units', 'pixels', 'Position', [120, 120, size, size]);
    set(hp, 'Units', 'pixels', 'Position', [1, 1, size, size]);
    prop = 'CData';
    issurf = false;
else
    if size > 512
        size = 512;
    end
    [x, y] = meshgrid(1:size, 1:size);
    h = surf(x, y, zeros(size, size));
    set(h, 'EdgeColor', 'none', 'LineStyle', 'none');
    prop = 'ZData';
    if ispc && ...
        ~xffsngl.CONF.settings.OpenGL.HardwareAccelOnWindows
        try
            opengl('software', true);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            warning( ...
                'neuroelf:OpenGLError', ...
                'Error setting OpenGL to SoftwareAcceleration for PCs.' ...
            );
        end
    end
    colormap((0:1/255:1)' * ones(1, 3));
    hp = get(h, 'Parent');
    set(hp, 'View', [166, 84]);
    set(get(hp, 'Parent'), ...
        'Units', 'pixels', 'Position', [120, 120, 640, 640]);
    set(hp, 'Units', 'pixels', 'Position', [-50, -40, 710, 720]);
    issurf = true;
end
set(hp, 'Color', [0, 0, 0]);

% for tilt
step = 256 / size;
hstep = 0.5 * step;
from = 0.5 + hstep;
to = 256.5 - hstep;
hrange = 0.5 * (to - from);
if strcmp(what, 'tilt')

    % for matlab
    if useml
        % create coordinate array
        [x, y] = ndgrid(to:-step:from, from:step:to);
        x = [x(:), y(:), 128 * ones(numel(x), 1), ones(numel(x), 1)];

    % for flexinterpn
    else

        % set default range
        r = [Inf, Inf, Inf; to, from, 128; -step, step, 1; from, to, 128];
    end

    % and trf matrices
    tf = eye(4);
    tt = tf;
    tf(1:3, 4) = -128;
    tt(1:3, 4) = 128;

    % loop from 0 to 2 * pi
    for a = 0:pi/90:(5*pi+sqrt(eps))

        % compute transformation matrix
        t = (tt * bvtrf([0, 0, 0], 90 * [sin(a), cos(0.3 * a), -sin(0.2 * a)]) * tf);

        % for matlab
        if useml

            % transform coordinates
            xt = x * t';
        end

        % get interpolated data
        if useml
            set(h, prop, reshape( ...
                interpn(xyzc{:}, vd, xt(:, 1), xt(:, 2), xt(:, 3), method), ...
                size, size)');
        else
            set(h, prop, flexinterpn_method(vd, r, 0, t, method, gaussk{:})');
        end

        % then show image
        if issurf
            set(hp, 'ZLim', [0, 255]);
        end
        drawnow;
    end

% for zoom
else

    % set default range
    r = [Inf, Inf, Inf; from, from, 128; step, step, 1; to, to, 128];

    if useml
        z3d = 128 * ones(size, size);
    end

    % loop from full size to 1/4 size in the middle slice
    for z = [0.75:0.01:0.99, ...
             1.00:0.02:2.00, ...
             2.03:0.03:2.99, ...
             3.03:0.04:3.99, ...
             4.02:0.06:6.00, ...
             6.08:0.08:8.00, ...
             8.10:0.10:10.0, ...
             10.2:0.20:15.0]

        % compute range
        rr = [128.5 - hrange / z; step / z; 128.5 + hrange / z];

        % for matlab
        if useml
            [x, y] = ndgrid(rr(1):rr(2):rr(3), rr(1), rr(2), rr(3));
            set(h, prop, interpn(xyzc{:}, vd, x, y, z3d, method)');

        else

            % use flex
            r(2:4, 1:2) = rr(:, ones(1, 2));
            set(h, prop, flexinterpn_method(vd, r, method, gaussk{:})');
        end

        % then show image
        if issurf
            set(hp, 'ZLim', [0, 255]);
        end
        drawnow;
    end
end

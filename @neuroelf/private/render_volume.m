function render_volume(VD, isovalue, limits, mycolor, scale)
% render_vol  - 3-D rendering of an spm_vol struct
%
% this is a (rather poor) approach for a 3D rendering of a
% stack of slices from a 3-D volume
%
% FORMAT:       render_vol(VD, ISOvalue, limits, color [,scale])
%
% Input fields:
%
%       VD          volume data
%       ISOvalue    threshold for "air"
%       limits      3-D limitations (e.g. [1,X, 1,Y, 1,Z])
%       color       surface color for mesh
%       scale       scale factor for input data
%
% See also spm_vol, spm_read_vols.

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% enough arguments ?
if nargin < 2 || ...
   ~isnumeric(VD) || ...
    length(size(VD)) ~= 3 || ...
    isempty(VD) || ...
   ~isnumeric(isovalue) || ...
    numel(isovalue) ~= 1
    error( ...
        'neuroelf:BadArgument',...
        'Too few or bad arguments. Try ''help %s''.',...
        mfilename ...
    );
end
if nargin < 3
    limits = NaN * ones(1, 6);
end
if nargin < 4
    mycolor = [0.2, 0.5, 0.9];
end

% get subvolume
[x y z D] = subvolume(VD, limits);

% while too big
while numel(D) > 1e7
    D = D(1:2:end, 1:2:end, 1:2:end);
    [x y z D] = subvolume(D, NaN * ones(1, 6));
end

% scale values
if nargin > 4
    D = D .* scale;
end

p = patch(isosurface(x,y,z,D,isovalue), ...
    'FaceColor', mycolor, 'EdgeColor', 'none');
patch(isocaps(x,y,z,D,isovalue), 'FaceColor', 'interp', 'EdgeColor', 'none');
isonormals(x,y,z,D,p);

% format view
view(3);
axis tight;
daspect([1 1 .4]);
axis equal;
colormap(gray(256));

camva(6);
box on;
camlight(40, 40);
camlight(-20,-10);
lighting gouraud;

function varargout = ne_srfupdatecoords(varargin)
% ne_srfupdatecoords  - update surface coordinates in patch
%
% FORMAT:       ne_srfupdatecoords([SRC, EVT, srf, patch, props])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       srf         SRF object (must be loaded, if not given all loaded!)
%       patch       1x1 patch UI object (if not given, use from Handles)
%       props       if not given, use srf.Handles.SurfProps
%
% No output fields. (will be set to [])
%
% Example:
%
%     ne_srfupdatecoords(0, 0, srf, srf.Handles.
%
%     this sets the surface viewpoint in satellite window with ID
%     'BS123456' to the standard position (left-hemisphere) with a
%     zoom factor of 1.2 and a morph of 1/3 along the way between
%     the surfaces and the morphing targets

% Version:  v1.1
% Build:    16031911
% Date:     Mar-19 2016, 11:30 AM EST
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

% global config
global ne_gcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% no input SRF
if nargin < 3 || numel(varargin{3}) ~= 1 || ~isxff(varargin{3}, {'fsbf', 'srf', 'tom'})

    % do for all
    srfs = ne_gcfg.h.Scenery.UserData(:, 4);
    for sc = 1:numel(srfs)
        try
            ne_srfupdatecoords(0, 0, srfs{sc});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
    return;
end

% alternative patch given
srf = varargin{3};
srfh = handles(srf);
if nargin > 3 && numel(varargin{4}) == 1 && ...
   (isa(varargin{4}, 'double') || isa(varargin{4}, 'matlab.graphics.primitive.Patch'))
    hp = varargin{4};
    ht = get(hp, 'Parent');
else
    hp = srfh.Surface;
    ht = srfh.SurfaceTransform;
end

% props
if nargin > 4 && iscell(varargin{5}) && numel(varargin{5}) > 2
    pp = varargin{5};
else
    pp = srfh.SurfProps;
end

% get and set coordinates
[p, pn] = btc_meshcn(srf, ne_gcfg.fcfg.srfcfg, ...
    ~strcmpi(get(hp, 'FaceColor'), 'none') || ~strcmpi(ne_gcfg.fcfg.renderer, 'opengl'));
if ~isempty(srf.TriangleVertex)
    set(hp, 'Vertices', p, 'VertexNormals', pn);
else
    set(hp, 'Vertices', p);
end

% general transform
set(ht, 'Matrix', btc_meshtrf(pp));

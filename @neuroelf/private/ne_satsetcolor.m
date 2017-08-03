function varargout = ne_satsetcolor(varargin)
% ne_satsetcolor  - set background color for satellite
%
% FORMAT:       ne_satsetcolor(SRC, EVT, window, bgcolor)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       window      window specifier (used to target satellite windows)
%       bgcolor     1x3 double color (0 .. 255 RGB)
%
% No output fields.
%
% Example:
%
%       ne_satsetcolor(0, 0, 'BS123456', [255, 255, 255]);

% Version:  v1.1
% Build:    16052817
% Date:     May-28 2016, 5:42 PM EST
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

% global variable
global ne_gcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% no window given
if nargin < 3 || ~ischar(varargin{3}) || isempty(varargin{3}) || ~isfield(ne_gcfg.cc, varargin{3}(:)')
    return;
end
tsat = varargin{3}(:)';
ch = ne_gcfg.cc.(tsat);

% for now, only for surfaces
if ~strcmpi(ch.Config.sattype, 'surf')
    return;
end

% size given
if nargin > 3 && isa(varargin{4}, 'double') && numel(varargin{4}) == 3 && ...
   ~any(isinf(varargin{4}) | isnan(varargin{4}) | varargin{4} < 0 | varargin{4} > 255)
    newcolor = varargin{4}(:)';
    if all(newcolor == fix(newcolor)) && any(newcolor > 1)
        newcolor = (1 / 255) .* newcolor;
    end
else
    newcolor = (1 / 255) .* colorpicker(round(255 .* ch.Config.SurfBackColor, {'Axes background'}));
    if numel(newcolor) ~= 3
        return;
    end
end

% with error handling
try

    % set color
    newcolor = min(1, max(0, newcolor));
    set(ch.SatelliteMLH, 'Color', newcolor);
    set(ch.Surface, 'Color', newcolor);
    if ne_gcfg.c.mlversion >= 900
        set(ch.Surface, 'XColor', newcolor', 'YColor', newcolor, 'ZColor', newcolor);
    end

    % and in config
    ne_gcfg.cc.(tsat).Config.SurfBackColor = get(ch.Surface, 'Color');
    if sum(newcolor) >= 1.5
        set(ch.SurfaceStatsText, 'Color', [0, 0, 0]);
    else
        set(ch.SurfaceStatsText, 'Color', [1, 1, 1]);
    end
    ne_setcsrfstatbars(0, 0, tsat);

% error?
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

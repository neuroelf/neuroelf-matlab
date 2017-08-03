function varargout = ne_sceneryselect(varargin)
% ne_sceneryselect  - select visible surfaces and re-set coordinates those
%
% FORMAT:       ne_sceneryselect(SRC, EVT [, srfidx [, uitag]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       srfidx      default: take from control, or double or {srf,...} list
%       uitag       re-select objects (from copied surfaces) in undocked UI
%
% No output fields.
%
% Examples:
%
%     ne_sceneryselect(0, 0, 1:3); % select first three surfaces
%     ne_sceneryselect(0, 0, {SRF}, 'BS123456'); % select SRF in BS123456

% Version:  v1.1
% Build:    16060912
% Date:     Jun-09 2016, 12:30 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% external window
if nargin > 3 && ischar(varargin{4}) && ~isempty(varargin{4}) && ...
    isfield(ne_gcfg.cc, varargin{4}(:)') && isstruct(ne_gcfg.cc.(varargin{4}(:)')) && ...
    isfield(ne_gcfg.cc.(varargin{4}(:)'), 'Config') && isstruct(ne_gcfg.cc.(varargin{4}(:)').Config) && ...
    isfield(ne_gcfg.cc.(varargin{4}(:)').Config, 'sattype') && ...
    ischar(ne_gcfg.cc.(varargin{4}(:)').Config.sattype) && ...
    strcmp(ne_gcfg.cc.(varargin{4}(:)').Config.sattype, 'surf')
    twin = {varargin{4}(:)'};
    ch = ne_gcfg.cc.(varargin{4}(:)');
else
    twin = {};
end

% scenery selection
scu = ch.Scenery.UserData;
if isempty(scu)
    return;
end
if nargin < 3 || isempty(varargin{3}) || (~iscell(varargin{3}) && ~isa(varargin{3}, 'double'))
    sci = ch.Scenery.Value;
elseif iscell(varargin{3})
    sci = zeros(numel(varargin{3}), 1);
    for sc = 1:numel(sci)
        if numel(varargin{3}{sc}) ~= 1 || ~isxff(varargin{3}{sc}, {'fsbf', 'srf'})
            fprintf('Requested element not a valid surface.\n');
            return;
        end
        for sct = 1:size(scu, 1)
            if varargin{3}{sc} == scu{sct, 4}
                sci(sc) = sct;
                break;
            end
        end
    end
    if any(sci == 0)
        fprintf('Requested surface not scenery list.\n');
        return;
    end
    sci = unique(sci(:));
    ch.Scenery.Value = sci;
elseif isa(varargin{3}, 'double')
    sci = varargin{3}(:);
    if any(isinf(sci) | isnan(sci) | sci < 1 | sci > size(scu, 1) | sci ~= fix(sci))
        fprintf('Invalid scenery index selection.');
        return;
    end
    sci = unique(sci(:));
    ch.Scenery.Value = sci;
end
sca = scu(:, 4);

% show surfaces if necessary
if isempty(twin) && cc.page ~= 3
    ne_showpage(0, 0, 3);
end

% apply changes
for sc = 1:numel(sca)

    % get surface and handles
    srf = sca{sc};
    sch = srf.Handles;

    % from root (no additional patch and map information)
    if size(scu, 2) == 4

        % visible
        if any(sci == sc)

            % set visible
            set(sch.Surface, 'Visible', 'on');
            set(sch.SurfaceTransform, 'Visible', 'on');

            % then also update coordinates
            ne_srfupdatecoords(0, 0, srf, sch.Surface, sch.SurfProps);

        % hide (don't process further)
        else
            set(sch.Surface, 'Visible', 'off');
            set(sch.SurfaceTransform, 'Visible', 'off');
        end

    % for specific window
    else
        % visible
        if any(sci == sc)

            % set visible and process
            set(scu{sc, 5}, 'Visible', 'on');
            btc_meshcolor(scu{sc, 4}, true, twin{1}, scu{sc, 5}, true);
            ne_srfupdatecoords(0, 0, srf, scu{sc, 5}, sch.SurfProps);
        else
            set(scu{sc, 5}, 'Visible', 'off');
        end
    end
end

% update position generally
if isempty(twin)
    ne_setsurfpos(0, 0, 1);
else
    ne_setsurfpos(0, 0, twin{:}, 'upshape');
end

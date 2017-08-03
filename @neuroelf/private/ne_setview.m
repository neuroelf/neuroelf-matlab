% FUNCTION ne_setview: set page and sub-page (view)
function varargout = ne_setview(varargin)

% Version:  v0.9b
% Build:    11051121
% Date:     Apr-10 2011, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% not enough input
if nargin < 4
    return;
end

% get configuration
cc = ne_gcfg.fcfg;

% only show page
if varargin{3} > 0
    ne_showpage(0, 0, varargin{3});

    % and update position accordingly
    if any((1:2) == varargin{3})
        ne_gcfg.fcfg.zoom = varargin{4};
        ne_setslicepos;
    elseif varargin{3} == 3
        ne_setsurfpos;
    end

% request new configuration
else

    % rotation, scaling and translation
    if nargin < 5 || ...
       ~iscell(varargin{5}) || ...
        numel(varargin{5}) ~= 3 || ...
       ~isa(varargin{5}{1}, 'double') || ...
        numel(varargin{5}{1}) ~= 3 || ...
       ~isa(varargin{5}{2}, 'double') || ...
        numel(varargin{5}{2}) ~= 3 || ...
       ~isa(varargin{5}{3}, 'double') || ...
        numel(varargin{5}{3})

        rot = cc.strrot;
        scl = cc.strscl;
        tra = cc.strtra;
        newv = inputdlg({'Translation (X, Y, Z):', 'Rotation (degrees):', ...
            'Scaling factor:'}, 'NeuroElf GUI - apply TRF to display', 1, ...
            {sprintf('  %g', tra), sprintf('  %g', rot), sprintf('  %.3f', scl)});

        % cancelled
        if isempty(newv)
            ne_setslicepos(0, 0, cc.cpos);
            return;
        end
    else
        newv = varargin{5};
    end

    % process new settings
    try
        if ischar(newv{1})
            tra = eval(['[' newv{1} ']']);
        else
            tra = newv{1}(:)';
        end
        if ischar(newv{2})
            rot = eval(['[' newv{2} ']']);
        else
            rot = newv{2}(:)';
        end
        if ischar(newv{3})
            scl = eval(['[' newv{3} ']']);
        else
            scl = newv{3}(:)';
        end
        if numel(tra) ~= 3 || ...
            any(isinf(tra) | isnan(tra) | tra < -256 | tra > 256) || ...
            numel(rot) ~= 3 || ...
            any(isinf(rot) | isnan(rot) | rot < -256 | rot > 256) || ...
            numel(scl) ~= 3 || ...
            any(isinf(scl) | isnan(scl) | abs(scl) < 0.1 | abs(scl) > 8)
            return;
        end

        % if valid, make permanent
        ne_gcfg.fcfg.strrot = rot(:)';
        ne_gcfg.fcfg.strtra = tra(:)';
        ne_gcfg.fcfg.strscl = scl(:)';
        ne_gcfg.fcfg.strans = spmtrf(tra(:)', (pi / 180) .* rot(:)', 1 ./ scl(:)');

        % and update screen
        ne_setslicepos(0, 0, cc.cpos);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        return;
    end
end

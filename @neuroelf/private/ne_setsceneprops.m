function varargout = ne_setsceneprops(varargin)
% ne_setsceneprops  - set properties of scenery content
%
% FORMAT:       ne_setsceneprops(SRC, EVT [, surface [, surfprops]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       surface     which surface to set properties for (current if empty)
%       surfprops   properties to set (default: UI-based)
%
% No output fields.
%
% Example:
%
%       neuroelf_gui('setsceneprops', lh, {trans, rot, scale, alpha, fow});
%       % trans = translation
%       % rot = rotation
%       % scale = scaling factor
%       % alpha = 1 - transparency
%       % fow = faces ('f') or wireframe ('w') mode

% Version:  v1.1
% Build:    16052716
% Date:     May-27 2016, 4:11 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% scenery selection
ch = ne_gcfg.h;
scu = ch.Scenery.UserData;
if nargin < 3 || numel(varargin{3}) ~= 1 || ~isxff(varargin{3}, {'fsbf', 'srf'})
    sci = ch.Scenery.Value;
else
    scif = false;
    for sci = 1:size(scu, 1)
        if varargin{3} == scu{sci, 4}
            scif = true;
            break;
        end
    end
    if ~scif
        fprintf('Requested surface not scenery list.\n');
        return;
    end
end
if numel(sci) ~= 1 || ~isxff(scu{sci, 4}, {'fsbf', 'srf'})
    return;
end
sco = scu{sci, 4};
sch = handles(sco);
scp = sch.SurfProps;

% translation, rotation, scaling, and alpha
if nargin < 4 || ~iscell(varargin{4}) || numel(varargin{4}) ~= 5
    tra = scp{1};
    rot = (180 / pi) .* scp{2};
    scl = scp{3};
    alp = scp{4};
    fow = scp{5};
    if scp{7}(1) == 'n'
        dvp = 'no';
    else
        dvp = 'yes';
    end
    iactive = true;
else
    tra = varargin{4}{1};
    rot = varargin{4}{2};
    scl = varargin{4}{3};
    alp = varargin{4}{4};
    fow = varargin{4}{5};
    if ~isa(tra, 'double') || ~isa(rot, 'double') || ~isa(scl, 'double') || ~isa(alp, 'double') || ...
        numel(tra) ~= 3 || any(isinf(tra) | isnan(tra) | abs(tra) > 256) || ...
        numel(rot) ~= 3 || any(isinf(rot) | isnan(rot) | abs(rot) > 360) || ...
        numel(scl) ~= 1 || isinf(scl) || isnan(scl) || scl < 0.2 || scl > 5 || ...
        numel(alp) ~= 1 || isinf(alp) || isnan(alp) || alp < 0 || alp > 1 || ...
       ~ischar(fow) || isempty(fow) || ~any(strcmpi(fow, {'f', 'w'}))
        fprintf('Requested properties invalid.\n');
        return;
    end
    scl = scl([1, 1, 1]);
    fow = lower(fow);
    if scp{7}(1) == 'n'
        dvp = 'no';
    else
        dvp = 'yes';
    end
    iactive = false;
end

% request updated values
if iactive
    newv = inputdlg({'Translation (X, Y, Z):', 'Rotation (degrees):', ...
        'Scaling factors (X, Y, Z):', 'Alpha (transparency)', ...
        'Display surface as (f)aces or (w)ireframe', ...
        'Display vertex coordinate points (yes/no)'}, ...
        'NeuroElf GUI - set surface properties within scene', 1, ...
        {sprintf('  %g', tra), sprintf('  %g', rot), sprintf('  %.3f', scl), ...
         sprintf('  %.3f', alp), ['  ' fow], ['  ' dvp]});

    % cancelled
    if isempty(newv) || ~iscell(newv) || numel(newv) ~= 6;
        return;
    end

    % bad values (overall)
    try
        tra = eval(['[' newv{1} ']']);
        rot = eval(['[' newv{2} ']']);
        scl = eval(['[' newv{3} ']']);
        alp = str2double(newv{4});
        fow = newv{5};
        fow(fow == ' ') = '';
        dvp = newv{6};
        dvp(dvp == ' ') = '';
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        return;
    end
end

% process new settings
try
    if numel(tra) ~= 3 || any(isinf(tra) | isnan(tra) | tra < -256 | tra > 256) || ...
        numel(rot) ~= 3 || any(isinf(rot) | isnan(rot) | rot < -360 | rot > 360) || ...
        numel(scl) ~= 3 || any(isinf(scl) | isnan(scl) | abs(scl) < 0.2 | abs(scl) > 5) || ...
        numel(alp) ~= 1 || isinf(alp) || isnan(alp) || alp < 0 || alp > 1 || ...
        numel(fow) ~= 1 || ~any('fw' == lower(fow)) || ...
        isempty(dvp) || ~any('ny' == lower(dvp(1)))
        return;
    end

    % adapt rotation to radiens
    rot = (pi / 180) .* rot(:)';

    % if valid, make permanent
    if lower(dvp(1)) == 'y'
        dvp = 'flat';
    else
        dvp = 'none';
    end
    sco.SetHandle('SurfProps', {tra(:)', rot, scl(:)', alp, lower(fow), [], dvp});

    % and update screen
    ne_srfupdatecoords(0, 0, sco, sco.Handles.Surface);
    ne_setsurfpos(0, 0, true);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

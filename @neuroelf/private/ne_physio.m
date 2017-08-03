% FUNCTION ne_physio: open physio dialog and handle callbacks
function varargout = ne_physio(varargin)

% Version:  v0.9d
% Build:    14062015
% Date:     Jun-20 2014, 3:35 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% output assignment
if nargout > 0
    varargout = cell(1, nargout);
end

% argument check
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmpi(varargin{3}(:)', {'ecg', 'gsr'}))
    return;
end
kind = lower(varargin{3}(:)');

% bring up dialog
try
    fPhysio = neuroelf_file('f', 'ne_physio');
catch ne_eo;
    rethrow(ne_eo);
end

% make settings
hTag = fPhysio.TagStruct;
if strcmp(kind, 'ecg')
    hTag.DD_physio_datakind.Value = 1;
else
    hTag.DD_physio_datakind.Value = 2;
end
hTag.ED_physio_sf.Callback = @ne_physio_sfedit;
hTag.BT_physio_sf_browse.Callback = @ne_physio_browse;
hTag.BT_physio_sf_load.Callback = @ne_physio_load;
hTag.DD_physio_datakind.Callback = @ne_physio_datakind;
hTag.ED_physio_sfreq.Callback = @ne_physio_newsfreq;
hTag.BT_physio_resample.Callback = @ne_physio_resample;
hTag.CB_physio_ecgrwchan.Callback = @ne_physio_togglerwc;
hTag.BT_physio_cancel.Callback = 'set(gcf, ''Visible'', ''off'');';
hTag.BT_physio_process.Callback = @ne_physio_process;
fPhysio.CloseRequestFcn = 'set(gcf, ''Visible'', ''off'');';

% store settings
fPhysio.UserData = struct( ...
    'hTag', hTag, ...
    'kind', kind, ...
    'obj',  [], ...
    'out',  struct);

% make dialog visible
fPhysio.HandleVisibility = 'callback';
fPhysio.Visible = 'on';

% wait til done
waitfor(fPhysio.MLHandle, 'Visible', 'off');

% get new settings
ud = fPhysio.UserData;

% delete figure
fPhysio.Delete;

% save RunTimeVars of object (for GSR data)
if strcmp(ud.kind, 'gsr') && ...
   ~isempty(ud.obj) && ...
   ~isempty(fieldnames(ud.obj.RunTimeVars))
    try
        ud.obj.SaveRunTimeVars;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        warning( ...
            'neuroelf:ErrorSavingRTV', ...
            'Error saving RunTimeVars with onsets for data file.' ...
        );
    end
end

% then clear object!
if ~isempty(ud.obj)
    ud.obj.ClearObject;
end

% for each output
of = fieldnames(ud.out);
for fc = 1:numel(of)

    % assing in BASE WS
    try
        assignin('base', of{fc}, ud.out.(of{fc}));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        warning( ...
            'neuroelf:ErrorCopyingOutputs', ...
            'Error assigning %s in BASE workspace.', ...
            of{fc} ...
        );
    end
end

% FUNCTION ne_console: handle console callbacks
function varargout = ne_console(varargin)
% ne_console  - handle console UI callbacks
%
% FORMAT:       [output = ] ne_console(src, evt, command [, arguments])
%
% Input fields:
%
%       src         handle of UI object issuing the call (0 for GUI)
%       evt         event data (currently unused)
%       command     internally used string to identify requested action
%
% Output fields:
%
%       output      optionally provided output

% Version:  v1.1
% Build:    16052500
% Date:     May-25 2016, 12:14 AM EST
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
cc = ne_gcfg.fcfg.Console;
ch = ne_gcfg.h.Console;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only with valid input
if nargin < 3 || ~ischar(varargin{3}) || ~any(strcmpi(varargin{3}, ...
    {'close', 'command', 'open'}))
    return;
end
action = lower(varargin{3}(:)');

% without open console, only allow open call
if (~isstruct(ch) || numel(ch) ~= 1 || ~isfield(ch, 'hFig')) && ~strcmp(action, 'open')
    return;
end

% action
switch (action)

    % close UI
    case 'close'

        % delete window
        ch.hFig.Delete;

        % store history
        if ~isempty(cc.hist)
            asciiwrite([neuroelf_path('config') filesep 'history.txt'], ...
                gluetostring(cc.hist, char(10)));
        end

        % remove diary
        diary('off');
        if exist(cc.diaryfile, 'file') > 0
            delete(cc.diaryfile);
        end

        % and config, etc.
        ne_gcfg.fcfg.Console = [];
        ne_gcfg.h.Console = [];

    % open UI
    case 'open'

        % only allow one instance
        if isstruct(ch) && numel(ch) == 1 && isfield(ch, 'hFig') && isxfigure(ch.hFig, true)

            % bring forward
            figure(ch.hFig.MLHandle);
            drawnow;
            return;
        end

        % load figure
        try
            hFig = neuroelf_file('f', 'console');
            if ~isxfigure(hFig, true)
                error('neuroelf:xfigure:errorCreatingTFG', ...
                    'Couldn''t create console figure.');
            end
        catch ne_eo;
            uiwait(warndlg(ne_eo.message, 'NeuroElf - error', 'modal'));
            return;
        end

        % create config
        cc = struct('diaryfile', [tempname '.diary'], 'diarypos', 0, 'execonreturn', true);
        ch = struct('hFig', hFig, 'tags', hFig.TagStruct);

        % turn on diary
        diary(cc.diaryfile);

        % get workspace
        ws = evalin('base', 'whos');
        ws = {ws.name};
        ws = ws(:);
        ch.tags.LB_NEConsole_workspc.String = ws;
        ch.tags.LB_NEConsole_workspc.Value = [];
        if isempty(ws)
            ch.tags.LB_NEConsole_workspc.ListboxTop = 1;
        else
            ch.tags.LB_NEConsole_workspc.ListboxTop = max(1, numel(ws) - 10);
        end

        % get history
        chist = splittocell(asciiread([neuroelf_path('config') filesep 'history.txt']), char(10));
        chist = chist(:);
        if ~isempty(chist)
            chist(cellfun('isempty', chist)) = [];
        end
        cc.hist = chist;
        ch.tags.LB_NEConsole_history.String = chist;
        ch.tags.LB_NEConsole_history.Value = [];
        if isempty(chist)
            ch.tags.LB_NEConsole_history.ListboxTop = 1;
        else
            ch.tags.LB_NEConsole_history.ListboxTop = max(1, numel(chist) - 10);
        end

        % set command to empty string
        ch.tags.ED_NEConsole_command.String = '';

        % store config
        ne_gcfg.fcfg.Console = cc;
        ne_gcfg.h.Console = ch;

        % make visible
        figure(hFig.MLHandle);
        drawnow;
        hFig.Visible = 'on';

    % execute command
    case 'command'

        % get command
        cmd = ch.tags.ED_NEConsole_command.String;
        ch.tags.ED_NEConsole_command.String = '';

        % execute 
        try
            evalin('base', cmd);
            ne_gcfg.fcfg.Console.hist{end+1} = cmd;
            ch.tags.LB_NEConsole_workspc.String = ne_gcfg.fcfg.Console.hist;

            % update workspace
            cws = ch.tags.LB_NEConsole_workspc.String(ch.tags.LB_NEConsole_workspc.Value);
            ws = evalin('base', 'whos');
            ws = {ws.name};
            ws = ws(:);
            ch.tags.LB_NEConsole_workspc.String = ws;
            if ~isempty(cws) && ~isempty(ws)
                cwi = multimatch(cws, ws);
                cwi = cwi(cwi > 0);
            else
                cwi = [];
            end
            ch.tags.LB_NEConsole_workspc.Value = cwi;
            if ~isempty(cwi)
                ch.tags.LB_NEConsole_workspc.ListboxTop = max(1, max(cwi) - 10);
            else
                ch.tags.LB_NEConsole_workspc.ListboxTop = 1;
            end
        catch ne_eo;
            ch.tags.ED_NEConsole_output.String = ...
                [ch.tags.ED_NEConsole_output.String char(10) 'Error: ' ne_eo.message];
        end

        % re-read diary
        diary('off');
        dfp = fopen(cc.diaryfile, 'r');
        if dfp < 1
            diary('on');
            return;
        end
        fseek(dfp, cc.diarypos, -1);
        newdata = fread(dfp, [1, Inf], 'uint8=>char');
        ne_gcfg.fcfg.Console.diarypos = ftell(dfp);
        fclose(dfp);
        diary('on');
        cstring = ch.tags.ED_NEConsole_output.String;
        if iscell(cstring)
            cstring = gluetostring(cstring, char(10));
        elseif size(cstring, 1) > 1
            cstring = gluetostring(deblank(cellstr(cstring)), char(10));
        end
        if cstring(end) == char(10)
            cstring(end) = [];
        end
        cstring = [cstring char(10) '>> ' cmd char(10) newdata];
        ncstring = numel(splittocell(cstring), char(10));
        ch.tags.ED_NEConsole_output.String = cstring;
        ch.tags.ED_NEConsole_output.ListboxTop = max(1, ncstring - 10);

    % unknown action
    otherwise
        fprintf('Unknown console action: %s.\n', action);
end

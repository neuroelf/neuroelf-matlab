% FUNCTION ne_savetextoutput: saves the text output as a text file
function ne_savetextoutput(varargin)
%
% FORMAT:       ne_savetextoutput(src, evt [, filename])
%
% Input fields:
%
%       src         handle of UI object issuing the call (0 for GUI)
%       evt         event data (currently unused)
%       filename    optional filename (otherwise requested)

% Version:  v0.9c
% Build:    11051414
% Date:     May-14 2011, 2:03 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
ch = ne_gcfg.h;

% output filename
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
    isempty(regexpi(varargin{3}(:)', '\.(csv|txt)$'))
    tlist = {'*.txt', 'Tab-delimited text (*.txt)'; ...
         '*.csv', 'Comma-separated-values (*.csv)'};
    [savefile, savepath, saveidx] = uiputfile(tlist, ...
        'Save text output as...');
    if isequal(savefile, 0) || ...
        isequal(savepath, 0) || ...
        isempty(savefile)
        return;
    end
    if isempty(savepath)
        savepath = pwd;
    end
    [nullpath, savefile, saveext] = fileparts(savefile);
    if isempty(saveext) || ...
       ~any(strcmpi(saveext, {'.csv', '.txt'}))
        saveext = {'.csv', '.txt'};
        saveext = saveext{saveidx};
    end
    savefile = [savepath, '/', savefile, saveext];
else
    savefile = varargin{4}(:)';
end

% try to save text
try
    text = ch.ClusterTable.String;
    text(text == '|') = char(9);
    text(text == char(13)) = char(10);
    if size(text, 1) > 1
        text = cellstr(text);
    else
        text = splittocellc(text, char(10), true);
    end
    text = ddeblank(text);
    if ~isempty(text) && ...
       ~isempty(regexpi(text{1}, 'clustertable'))
        text = regexprep(text, ...
            '^([\+\-]?\d+)\s+([\+\-]?\d+)\s+([\+\-]?\d+)', ...
            ['$1' char(9) '$2' char(9) '$3']);
        text = strrep(text, [' L' char(9)], [char(9) 'L']);
        for lc = 1:numel(text)
            if ~isempty(regexpi(text{lc}, '^x\s+y\s+z'))
                text{lc} = gluetostringc(splittocellc( ...
                    text{lc}, char([9, 32]), true, true), char(9));
                break;
            end
        end
    end
    text = strrep(text, char(9), '__TAB__');
    text = regexprep(text, ' +__TAB__', '__TAB__');
    text = regexprep(text, '__TAB__ +', '__TAB__');
    text = strrep(text, '__TAB__', char(9));
    text = gluetostringc(text, char(10), true);
    if strcmpi(savefile(end-2:end), 'csv')
        text(text == char(9)) = ';';
    end
    asciiwrite(savefile, text);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg(['Error writing text output file: ', ne_eo.message], ...
        'NeuroElf - error', 'modal'));
end

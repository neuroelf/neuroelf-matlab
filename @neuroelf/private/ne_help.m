function varargout = ne_help(varargin)
% ne_help  - print help of NeuroElf GUI sub-function
%
% FORMAT:       [helptext = ] ne_help(SRC, EVT, fname)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       fname       function name on which to retrieve help
%
% Output fields:
%
%       helptext    if requested, return helptext instead of printing it
%
% Example:
%
%       ne_help(0, 0, 'srf_tools')
%
%       will retrieve the help text (comment) in function ne_srf_tools.m
%
% Applies to:
%
%       extended help is currently available for the following functions
%
%       closesatwindow  - close satellite window (scripted)
%       draw            - draw into currently selected dataset (VMR/HDR)
%       openfile        - add a file or object to the GUI
%       screenshot      - create a screenshot file
%       srf_save        - save currently selected SurfVar
%       srf_tools       - tools for currently selected SurfVar

% Version:  v1.0
% Build:    16011513
% Date:     Jan-15 2016, 1:56 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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

% pre-set output
if nargout > 0
    varargout = cell(1, nargout);
end

% no input
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3})
    fname = 'help';
else
    fname = lower(varargin{3}(:)');
end

% retrieve help
selfpath = fileparts(which(mfilename));
fullname = findfiles(selfpath, ['ne_' fname '.m'], 'depth=1');
if isempty(fullname)
    fullname = findfiles(selfpath, [fname '.m'], 'depth=1');
end
if isempty(fullname)
    helptext = sprintf('%% Function %s not found.', fname);
else
    
    % read file
    helptext = asciiread(fullname{1});
end

% process text
helptext = splittocellc(helptext, char(10));
lc = 1;
while numel(helptext) >= lc && ...
   (isempty(helptext{lc}) || ...
    helptext{lc}(1) ~= '%')
    lc = lc + 1;
end
flc = lc;
while numel(helptext) >= lc && ...
   ~isempty(helptext{lc}) && ...
    helptext{lc}(1) == '%'
    lc = lc + 1;
end
helptext = helptext([flc, flc:lc]);
helptext{1} = '';

% remove '% '
helptext = regexprep(helptext, '^%\s?', '');

% back to single string
helptext = gluetostringc(helptext, char(10));

% replace for ne_XXXX functions
if ~isempty(strfind(helptext, ['ne_' fname]))
    helptext = strrep(helptext, ['ne_' fname '(0, 0'], ...
        ['neuroelf_gui(''' fname '''']);
end

% what to do
if nargout == 0
    disp(helptext);
else
    varargout{1} = helptext;
end

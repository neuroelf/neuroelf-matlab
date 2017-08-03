function glued = gluetostring(textincells, linedelim, endterm)
% gluetostring  - glues a cell array of chars to one delimited char array
%
% FORMAT:         glued = gluetostring(textincells ,[linedelim,endterm])
%
% Input Fields:
%    textincells  cell array of char arrays to glue
%    linedelim    optional char array with glueing array
%                 defaults to system dependend CR/LF or LF
%    endterm      if given, do not remove end terminator
%
% if no linedelim is given, lines are limited system dependend
%
% See also splittocell.
%
% Note: this is the function equivalent to gluetostringc, written in
%       Matlab code for situations where the MEX file is not yet
%       compiled. Once this is done, functions use the faster one.

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

% enough/good arguments ?
if nargin < 1 || ...
   (~iscell(textincells) && ...
    ~ischar(textincells))
    error( ...
        'neuroelf:BadArguments',...
        'Bad or too few arguments. Try ''help %s''.',...
        mfilename ...
    );
end

% make sure we have a cell array
if ischar(textincells)
    textincells = cellstr(textincells);
end

% what delimiter
if nargin < 2
    if ispc
        linedelim = char([13,10]);
    else
        linedelim = char(10);
    end
end

% how many items
ilen = numel(textincells);

% for few items ...
if ilen < 4096
    if isempty(textincells)
        glued = '';
        return;
    end
    glued = sprintf(['%s' ...
       strrep(strrep(linedelim, '%', '%%'), '\', '\\')], textincells{:});

% many items ...
else
    tlen = min(256, floor(sqrt(ilen)));
    tvec = 1:tlen:ilen;
    tcnt = cell(1, numel(tvec));
    tac  = 1;
    for tc = 1:tlen:ilen
        tcnt{tac} = ...
            gluetostring(textincells(tc:min(ilen, (tc+tlen-1))), linedelim);
        tac = tac + 1;
    end
    glued = gluetostring(tcnt, linedelim, 1);
end

% end terminator
if (nargin < 3 || ...
    isempty(endterm)) && ...
   ~isempty(glued)
    glued((end+1-numel(linedelim)):end) = [];
end

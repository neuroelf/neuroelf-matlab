function asciiwrite(filename, filecont, append, linedelim)
% asciiwrite  - writes a textfile from a char array to file
%
% FORMAT:       asciiwrite(filename, filecont [, append, linedelim])
%
% Input Fields:
%       filename    name to a file, preferably absolute path
%       filecont    char array to write
%       append      if true, content will be appended to file
%       linedelim   optional char array with line delimiter
%
% Note: if no linedelim is given, lines are limited system dependend
%
% See also asciiread, binwrite.

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

% enough arguments ?
if nargin < 2 || ...
   ~ischar(filename) || ...
    isempty(filename) || ...
   ~ischar(filecont)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or too few arguments. Try ''help %s''.', ...
        mfilename ...
    );
end

% options
if nargin < 3 || ...
    isempty(append) || ...
   ~append(1)
    omode = 'w';
else
    omode = 'a';
end
if nargin < 4 || ...
   ~ischar(linedelim)
    if ispc
        linedelim = char([13, 10]);
    else
        linedelim = char(10);
    end
end

% standard check on filename
if ispc
    filename = strrep(filename, '/', filesep);
else
    filename = strrep(filename, '\', filesep);
end

% get file pointer
ofp = fopen(filename, omode);
if ofp < 1
    error( ...
        'neuroelf:FileNotWritable', ...
        'Couldn''t write to file ''%s''.', ...
        filename ...
    );
end
frewind(ofp);

% some more sanity checks on content
if iscell(filecont)
    filecont = char(filecont);
end

% make a 2-D char array a cell -> glue with linedelim
if all(size(filecont) > 1)
    filecont = gluetostringc(cellstr(filecont), linedelim);
end

% write output
fwrite(ofp, filecont, 'uchar');
fclose(ofp);

function fileswapendian(filename, wordsize, iooff, iosize)
% fileswapendian  - swap endian type of (partial) file content
%
% FORMAT:       fileswapendian(filename, wordsize [, iooff [, iosize]])
%
% Input fields:
%
%       filename    filename
%       wordsize    endian wordsize, either of [2, 4, 8, 16]
%       iooff       1x1 value, default 0
%       iosize      1x1 value, default til end of file

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 2:11 PM EST
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

% argument check
if nargin < 2 || ...
   ~ischar(filename) || ...
    isempty(filename) || ...
    exist(filename(:)', 'file') ~= 2 || ...
   ~isa(wordsize, 'double') || ...
    numel(wordsize) ~= 1 || ...
   ~any([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64] == wordsize)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
filename = filename(:)';
filedir = dir(filename);
filesiz = filedir.bytes;
if nargin < 3
    iooff = 0;
end
if nargin < 4
    iosize = Inf;
end
if ~isa(iooff, 'double') || ...
    numel(iooff) ~= 1 || ...
    isinf(iooff) || ...
    isnan(iooff) || ...
    iooff < 0 || ...
    iooff ~= fix(iooff) || ...
   ~isa(iosize, 'double') || ...
    numel(iosize) ~= 1 || ...
    isnan(iosize) || ...
    iosize < wordsize || ...
    (iosize / wordsize) ~= fix(iosize / wordsize) || ...
    (~isinf(iosize) && ...
     (iooff + iosize) > filesiz)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument combination.' ...
    );
end
if isinf(iosize)
    iosize = filesiz - iooff;
end
ioend = iooff + iosize;

% open file for read/write access
try
    fid = fopen(filename, 'r+', 'ieee-le');
    if fid < 1
        error('BAD_FILE_ID');
    end
    fseek(fid, iooff, -1);
    if ftell(fid) ~= iooff
        error('SEEK_FAILED');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:FileReadError', ...
        'Error opening/seeking file.' ...
    );
end

% create a loop to allow working on large files
lrgbufs = 3 * 2 ^ 20;
lrgread = lrgbufs / wordsize;
for cpos = iooff:lrgbufs:ioend

    % seek to position
    fseek(fid, cpos, -1);

    % for larger chunks
    if (cpos + lrgbufs) <= ioend

        % read content
        bytecont = fread(fid, [wordsize, lrgread], '*uint8');

        % seek back
        fseek(fid, cpos, -1);

        % write with order swapped
        fwrite(fid, bytecont(end:-1:1, :), 'uint8');

    % for smaller chunks
    else

        % calc number of elements
        smlread = (ioend - cpos) / wordsize;

        % check number
        if smlread ~= fix(smlread)
            fclose(fid);
            error( ...
                'neuroelf:InternalError', ...
                'Error accessing file.' ...
            );
        end

        % then read/seek/write
        bytecont = fread(fid, [wordsize, smlread], '*uint8');
        fseek(fid, cpos, -1);
        fwrite(fid, bytecont(end:-1:1, :), 'uint8');
    end
end

% close file
fclose(fid);

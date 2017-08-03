function isle = fileguessendian(filename, wordsize, iooff, iosize)
% fileguessendian  - guess endian type of (partial) file content
%
% FORMAT:       isle = fileguessendian(filename [, ws [, iooff [, iosize]]])
%
% Input fields:
%
%       filename    filename
%       ws          endian wordsize, either of [2, 4, 8, 16]
%       iooff       1x1 value, default 0
%       iosize      1x1 value, default til end of file
%
% Output fields:
%
%       isle        true if content was guessed as little endian
%                   false if content was guessed as big endian
%                   NaN if content could not be guessed as either
%
% Note: - if wordsize is not given, the file will be scanned from 40%
%         of it's size onwards for the usual wordsizes 4 and 8

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
if nargin < 1 || ...
   ~ischar(filename) || ...
    isempty(filename) || ...
    exist(filename(:)', 'file') ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
filename = filename(:)';
filedir = dir(filename);
filesiz = filedir.bytes;

% detect completely automatic
isle = NaN;
if nargin < 3 || ...
   ~isa(wordsize, 'double') || ...
    numel(wordsize) ~= 1 || ...
   ~any([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64] == wordsize)

    % check iopos and iosize
    if nargin < 3
        iooff = floor(0.4 * filesiz);
    end
    if nargin < 4
        iosize = Inf;
    end

    % repeat this at different positions
    for k = 1:6
        islev = zeros(1, 13);
        if nargin > 2
            islev(1) = fileguessendian(filename, 2, iooff, iosize);
        else
            islev(1) = NaN;
        end
        for c = 1:4
            islev(c + 1) = fileguessendian(filename, 4, iooff + c - 1, iosize);
        end
        for c = 1:8
            islev(c + 5) = fileguessendian(filename, 8, iooff + c - 1, iosize);
        end
        isled = (~isnan(islev));
        isles = sum(isled);
        if isles > 0
            islew = sum(islev(isled));
            if islew >= (0.75 * isles)
                isle = true;
            elseif islew <= (0.25 * isles)
                isle = false;
            end
        end
        if ~isnan(isle) && ...
            nargin > 2
            return;
        end

        % if file large enough get new IOPos
        if filesiz > 65536
            iooff = floor(rand(1) * filesiz);

        % otherwise give up
        else
            return;
        end
    end
    return;
end

% wordsize is known
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

% open file for read access
try
    fid = fopen(filename, 'r', 'ieee-le');
    if fid < 1
        error('BAD_FILE_ID');
    end
    fseek(fid, iooff, -1);
    if ftell(fid) ~= iooff
        error('SEEK_FAILED');
    end

    % look for first byte deviating from zero
    ciooff = iooff + 32768;
    iszero = 16;
    while iszero > 0
        if ciooff <= filesiz
            ftr = fread(fid, [32768, 1], '*uint8');
        else
            ftr = fread(fid, [1, Inf], '*uint8');
            iszero = 0;
        end
        ftn = find(ftr > 0);
        if ~isempty(ftn)
            iszero = 0;
            iooff = (ftell(fid) - numel(ftr)) + ...
                wordsize * floor((ftn(1) - 1) / wordsize);
            if iooff < 0
                iooff = 0;
            end
        elseif iszero == 0
            fclose(fid);
            return;
        end
        iszero = iszero - 1;
    end
catch ne_eo;
    error( ...
        'neuroelf:FileReadError', ...
        'Error opening/seeking file: %s.', ...
        ne_eo.message ...
    );
end
if isinf(iosize)
    iosize = filesiz - iooff;
end
ioend = iooff + wordsize * floor(iosize / wordsize);

% create a loop to allow working on large files
guesses = 0;
ilrgread = 1536 / wordsize;
lrgread = ilrgread;
lrgth = ceil(lrgread / 5 + 1);
for cpos = iooff:1536:ioend

    % seek to position
    fseek(fid, cpos, -1);

    % for smaller chunks
    if (cpos + 1536) > ioend

        % re-calc number of elements
        lrgread = (ioend - cpos) / wordsize;
        lrgth = ceil(lrgread / 4 + 1);

        % check number
        if lrgread ~= fix(lrgread)
            fclose(fid);
            warning( ...
                'neuroelf:InternalError', ...
                'Error accessing file.' ...
            );
            return;
        end
    end

    % read content
    bytecont = fread(fid, [wordsize, lrgread], '*uint8');

    % check some larger wordsize things first
    if any([4, 8] == wordsize)

        % small floats start with hex 0x3f
        fldet = sum(bytecont == 63, 2);
        nulsum = sum(bytecont == 0, 2);

        % check for float 1's
        if wordsize == 4
            if all([fldet(4), sum(bytecont(3, :) == 128), nulsum(1)] > lrgth)
                isle = true;
                break;
            elseif all([fldet(1), sum(bytecont(2, :) == 128), nulsum(end)] > lrgth)
                isle = false;
                break;
            end
        else
            if all([fldet(8), sum(bytecont(7, :) == 255)] > lrgth)
                isle = true;
                break;
            elseif all([fldet(1), sum(bytecont(2, :) == 255)] > lrgth)
                isle = false;
                break;
            end
        end

        % calc difference in meaningful range
        lhdiff = sum(bytecont(1, :) > 59 & bytecont(1, :) < 73) - ...
            sum(bytecont(end, :) > 59 & bytecont(end, :) < 73);
        fldmed = double(median(bytecont, 2));
        fldmed(fldmed == 0) = NaN;
        fldsum = sum(bytecont == repmat(fldmed, [1, lrgread]), 2);

        % check difference, LE?
        if wordsize == 4
            if fldsum(4) > lrgth && ...
                lhdiff < -lrgth

                % make sure it's in the right direction
                if fldsum(3) > (lrgth / 2)
                    isle = true;
                    break;
                elseif fldsum(1) > (lrgth / 2)
                    isle = false;
                    break;
                end
            elseif fldsum(1) > lrgth && ...
                lhdiff > lrgth

                % make sure it's in the right direction
                if fldsum(2) > (lrgth / 2)
                    isle = false;
                    break;
                elseif fldsum(4) > (lrgth / 2)
                    isle = true;
                    break;
                end
            end
        else
            if fldsum(8) > lrgth && ...
                lhdiff < -lrgth

                % make sure it's in the right direction
                isle = true;
                if fldsum(1) > (lrgth / 2)
                    isle = false;
                end
                break;
            elseif fldsum(1) > lrgth && ...
                lhdiff > lrgth

                % make sure it's in the right direction
                isle = false;
                if fldsum(2) > (lrgth / 2)
                    isle = true;
                end
                break;
            end
        end

        % for larger integrals
        if nulsum(end) == lrgread && ...
            nulsum(end - 1) > (0.75 * lrgread) && ...
            nulsum(end - 2) > 0 && ...
            nulsum(1) < (0.1 * lrgread)
            isle = true;
            break;
        elseif nulsum(1) == lrgread && ...
            nulsum(2) > (0.75 * lrgread) && ...
            nulsum(3) > 0 && ...
            nulsum(end) < (0.1 * lrgread)
            isle = false;
            break;
        end
    end

    % look for intresting entries
    lowon = (bytecont(1, :) > 3 & bytecont(1, :) < 252);
    highon = (bytecont(end, :) > 3 & bytecont(end, :) < 252);

    % check simple difference
    lhdiff = sum(highon(~lowon)) - sum(lowon(~highon));
    if lhdiff < -lrgth
        isle = true;
        break;
    elseif lhdiff > lrgth
        isle = false;
        break;
    end

    % some last resort
    if any([4, 8] == wordsize)
        if fldsum(end) > (0.8 * lrgread + 1) && ...
            all(fldsum(1:end-1) < lrgth)
            isle = true;
            break;
        elseif fldsum(1) > (0.8 * lrgread + 1) && ...
            all(fldsum(2:end) < lrgth)
            isle = false;
            break;
        end
    end

    % increase and check number of guesses
    guesses = guesses + 1;
    if guesses > 3 || ...
        lrgread ~= ilrgread
        break;
    end
end

% close file
fclose(fid);

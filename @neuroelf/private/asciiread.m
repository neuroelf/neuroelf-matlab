function [textinfile, starts, ends, lnum] = asciiread(filename, linedelim, uc)
% asciiread  - reads a textfile into one char array
%
% FORMAT:       textinfile = asciiread(filename [, linedelim [, uc]])
%
% Input Fields:
%
%       filename    name to a file, preferably absolute path
%       linedelim   optional char array with line delimiter
%       uc          unicode (detect unicode files, default: true)
%
% Output will be delimited as is written in the file. Hence, an
% optional line delimiter can be given...
%
% See also asciiwrite

% Version:  v0.9c
% Build:    14012811
% Date:     Jan-28 2014, 11:18 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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
    isempty(filename)
    error( ...
        'neuroelf:BadArgument',...
        'Bad argument given. Try ''help %s''.',...
        mfilename ...
    );
end
if nargin < 3 || ...
   ~islogical(uc) || ...
    numel(uc) ~= 1
    uc = true;
end

% open file (check readability)
filename = filename(:)';
fp = fopen(filename, 'r');
if fp < 1
    error( ...
        'neuroelf:FileNotReadable',...
        'File not readable: %s',...
        filename ...
    );
end

% read file contents
textinfile = fread(fp, [1, Inf], 'uchar=>char');

% unicode (discarding non-ascii characters!)
if uc && ...
    numel(textinfile) > 2 && ...
    mod(numel(textinfile), 2) == 0 && ...
    isequal(textinfile(1:2), char([255, 254]))

    % re-read file
    fseek(fp, 0, -1);
    textinfile = fread(fp, [1, Inf], 'uint16=>double');
    textinfile(1) = [];
    textinfile = char(textinfile);
end

% close file
fclose(fp);

% detect linedelim if not given
if (nargin > 1 && ...
     ischar(linedelim)) || ...
    nargout > 2

    % get first 16kB (if that many)
    arp  = textinfile(1:min(16384, length(textinfile)));

    % find Mac, Windows and Linux <RETURN> codes
    LfCr = strfind(arp, char([10, 13]));
    CrLf = strfind(arp, char([13, 10]));
    Lf   = strfind(arp, char(10));
    Cr   = strfind(arp, char(13));

    % is there a combined code (not unique, happens with blank lines!)
    if ~isempty(LfCr)

        % is there ONLY ONE type of combined code
        if isempty(CrLf)

            % then we know it
            srcld = char([10, 13]);

        % otherwise
        else

            % decide according to the number of matches
            if length(LfCr) > length(CrLf)
                srcld = char([10, 13]);
            else
                srcld = char([13, 10]);
            end
        end

    % if ONLY THE OTHER combination is found
    elseif ~isempty(CrLf)
        srcld = char([13, 10]);

    % only Linux linefeeds
    elseif ~isempty(Lf)
        srcld = char(10);

    % only carriage returns (RARE!)
    elseif ~isempty(Cr)
        srcld = char(13);

    % no <RETURN> codes, delimiter must then be given externally!
    else
        error( ...
            'neuroelf:AutoDetectFailed', ...
            'Could not determine line delimiter in source.' ...
        );
    end
end

% replace line delimiter before continuing
if nargin > 1 && ...
    ischar(linedelim)
    textinfile = strrep(textinfile, srcld, linedelim);
    srcld = linedelim;
end

% given second output argument if needed
if nargout > 2
    lsrcld = length(srcld);
    ltext = length(textinfile);
    if lsrcld == 1
        ends = find(textinfile == srcld);
    else
        ends = strfind(textinfile, srcld);
    end
    if (ends(end) + lsrcld) < ltext
        ends(end + 1) = ltext + 1;
    end
    starts = [1 (ends(1:end-1) + lsrcld)];
    ends = ends - 1;
    lnum = length(starts);
end

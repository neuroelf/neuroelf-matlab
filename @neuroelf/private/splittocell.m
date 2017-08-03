function [linetocell, cellcount] = splittocell(line, delimiter, varargin)
%SPLITTOCELL  Split a delimited string into a cell array.
%   SPLIT = SPLITTOCELL(GLUEDSTRING) splits the GLUEDSTRING into its parts
%   assuming the delimiter to be the tab character, CHAR(9).
%
%   SPLIT = SPLITTOCELL(GLUEDSTRING, DELIM) uses DELIM instead as a
%   delimiter (or set of delimiters).
%
%   SPLIT = SPLITTOCELL(GLUEDSTRING, DELIM, MD) treats all consecutive
%   delimiters as one, i.e. no empty tokens will be returned in SPLIT
%   besides a possible first and last token, if MD is set to 1.
%
% FORMAT:         [out, cnt] = splittocell(string [,delims, md, ad, hq])
%
% Input fields:
%    string       1xN char array to split
%    delims       char array containing one or more delimiters
%                 if left empty -> char(9) == <TAB>
%    md           must be '1' (numeric) to be effective, if set
%                 multiple delimiters will be treated as one
%    ad           match any of the delimiters (for multi-char
%                 delimiters strings)
%    hq           must be '1' (numeric) to be effective, if set
%                 delimiters within quotes will be ignored, sets
%                 md and ad to false!
%
% Output fields:
%    out          cell array containing the tokens after split
%    cnt          number of tokens in result
%
% See also gluetostring.
%
% Note: this is the function equivalent to splittocellc, written in
%       Matlab code for situations where the MEX file is not yet
%       compiled. Once this is done, functions use the faster one.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:10 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
if nargin < 1 || nargin > 5
    error('neuroelf:splittocell:badNumberOfInputs', 'Invalid number of inputs.');
end

% initialize return values and varargin{3}
linetocell = cell(0);
cellcount  = 0;

% return if empty
if isempty(line)
    return;
end

% input must be a line (1xN char)
if ~ischar(line) || size(line, 2) ~= numel(line)
    error('neuroelf:splittocell:badArgument', 'Input must be a 1xN shaped char array!');
end

% initialize
multidelim = false;
anydelim   = false;
heedquotes = false;
cellcount = 1;

% parse additional arguments
if nargin < 2 || ~ischar(delimiter)
    delimiter = char(9);
else
    delimiter = delimiter(:)';
end
if nargin > 2 && numel(varargin{1}) == 1 && (isnumeric(varargin{1}) || islogical(varargin{1})) && varargin{1}
    multidelim = true;
end
if nargin > 3 && numel(varargin{2}) == 1 && (isnumeric(varargin{2}) || islogical(varargin{2})) && varargin{2}
    anydelim = true;
end
if nargin > 4 && numel(varargin{3}) == 1 && (isnumeric(varargin{3}) || islogical(varargin{3})) && varargin{3}
    heedquotes = true;
end

% safety check
if heedquotes
    multidelim = false;
    anydelim = false;
end

% set initial parameters
lline  = size(line, 2);
ldelim = size(delimiter, 2);

% standard approach (delim can be 1x1 or a longer token)
if ~anydelim

    % delimiter is shorter than the text
    if ldelim < lline

        % last string in text is delimiter
        if strcmp(line(end+1-ldelim:end), delimiter)

            % find each occurrence
            cpos = [(1 - ldelim), strfind(line, delimiter)];

        % otherwise
        else

            % do the same, but add one additional virtual one at the end
            cpos = [(1 - ldelim), strfind(line, delimiter), lline + 1];
        end

    % delimiter has the same length
    elseif ldelim == lline

        % the text is the delimiter
        if strcmp(line, delimiter)

            % return empty string
            linetocell = {''};

        % otherwise
        else

            % the delimiter doesn't occur, so return the text instead
            linetocell = {line};
        end
        return;

    % the delimiter is longer, so it cannot occur
    else

        % just return the text
        linetocell = {line};
        return;
    end

% any of the given delimiters (e.g. white spaces)
else

    % last char is a delimiter
    if ~any(delimiter == line(end))
        cpos = [0, lline + 1];
    else
        cpos = 0;
    end

    % get all delimiter positions
    for pchar = delimiter
        cpos = union(cpos, strfind(line, pchar));
    end

    % set ldelim to 1!
    ldelim = 1;
end

% number of delimiters
lcpos = length(cpos);

% any delimiter found at all ?
if lcpos < 2

    % raise error
    error('neuroelf:splittocell:internalError', 'Error working with pattern.');

% exactly one token
elseif lcpos == 2

    % return that token
    linetocell = {line(cpos(1) + ldelim:cpos(2) - 1)};
    return;
end

% concatenate in case of multidelims
ecpos = cpos(2:end) - 1;
if multidelim
    if ecpos(1) == 0
        mcpos = [1, find(diff(ecpos) <= ldelim) + 1];
    else
        mcpos = find(diff(ecpos) <= ldelim) + 1;
    end
    cpos(mcpos) = [];
    ecpos(mcpos) = [];
end
cpos = cpos + ldelim;
ncpos = numel(ecpos);
cellcount = ncpos;
linetocell = cell(1, ncpos);

% extract substrings
if ~heedquotes
    for dpos = 1:ncpos
        if cpos(dpos) <= ecpos(dpos)
            linetocell{dpos} = line(cpos(dpos):ecpos(dpos));
        else
            linetocell{dpos} = '';
        end
    end
else
    dpos = 1;
    tpos = 1;
    while dpos <= ncpos
        cellc = line(cpos(dpos):ecpos(dpos));
        hq = false;
        if ~isempty(cellc) && ...
            any(cellc(1) == '"''')
            if numel(cellc) > 1 && ...
                cellc(1) == cellc(end)
                if numel(cellc) > 2
                    cellc = cellc(2:end-1);
                else
                    cellc = '';
                end
            else
                for xpos = dpos+1:ncpos
                    xcellc = line(cpos(xpos):ecpos(xpos));
                    if ~isempty(xcellc) && ...
                        xcellc(end) == cellc(1)
                        hq = true;
                        break;
                    end
                end
            end
        elseif isempty(cellc)
            cellc = '';
        end
        if hq
            linetocell{tpos} = line(cpos(dpos)+1:ecpos(xpos)-1);
            dpos = xpos;
        else
            linetocell{tpos} = cellc;
        end
        dpos = dpos + 1;
        tpos = tpos + 1;
    end
    if tpos < dpos
        linetocell(tpos:end) = [];
    end
end

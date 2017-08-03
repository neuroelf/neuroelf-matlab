function [grepped] = grep(togrep, varargin)
% grep  - GNU utils like grep
%
% FORMAT:       grepped = grep(togrep, params, pattern, ...)
%
% Input fields:
%
%       togrep      either a filename, a text file contents, or a list
%                   (cell array) of strings
%       params      list of options (must be preceded by '-')
%                   multiple options can be mixed by adding them to the
%                   params argument, e.g. '-inr', available options are:
%                   'aN' - match N lines after match
%                   'bN' - match N lines before match
%                   'hN' - only consider first N lines of input
%                   'HN' - only return first N lines of output
%                   'i'  - ignore case
%                   'l'  - return line numbers instead of contents
%                   'n'  - enumerate (prepend line number to strings)
%                   'tN' - only consider last N lines of input
%                   'TN' - only return last N lines of output
%                   'v'  - inverse match
%                   'x'  - use regexp instead of strfind
%                   '[STRING]' - tag the output lines with STRING
%       pattern     one or more patterns given to match input with
%
% Output fields:
%
%       grepped     list of strings or line numbers with matched lines
%
% See also strfind, strmatch, regexp.

% Version:  v0.9a
% Build:    11043015
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
   ~ischar(varargin{1})
    error( ...
        'neuroelf:BadArgument',...
        'Bad or missing argument.' ...
    );
end

% get option content right
firstarg = varargin{1}(:)';
if ~isempty(firstarg) && ...
    firstarg(1) == '-'
    opts = firstarg(2:end);
    patterns = {varargin{2:end}};
else
    opts = '';
    patterns = varargin;
end
pats = numel(patterns);
for pat = 1:pats
    patterns{pat} = patterns{pat}(:)';
end

% char content
if ischar(togrep)
    ichar = 1;
    togrep = togrep(:)';

    % filename
    if ~isempty(togrep) && ...
        numel(togrep) < 256 && ...
        exist(togrep, 'file') == 2
        togrep = asciiread(togrep);
    end

    % split according to line delim
    if ~isempty(strfind(togrep, char([13, 10])))
        [togrep, lines] = splittocellc(togrep, char([13, 10]), false, false);
    else
        [togrep, lines] = splittocellc(togrep, char(10), false, false);
    end

% cell array already
elseif iscell(togrep)
    ichar = 0;
    lines = numel(togrep);

    % check line content
    for line = 1:lines

        % only allow chars
        if ~ischar(togrep{line})
            error( ...
                'neuroelf:BadArgument',...
                'Bad cell input for grep.' ...
            );
        else
            togrep{line} = togrep{line}(:)';
        end
    end

% other content
else
    error( ...
        'neuroelf:BadArgument',...
        'Bad togrep argument.' ...
    );
end

% check option: [tag]
tagname = '';
if any(opts == '[') && ...
    any(opts == ']')
    opos = find(opts == '[');
    opos2 = find(opts == ']');
    if opos(1) < opos2(1)
        tagname = opts(opos(1)+1:opos2(1)-1);
        opts = strrep(opts, ['[' tagname ']'],'');
        tagname = [tagname ', '];
    end
end

% check option: a(fter)
doafter = 0;
if any(opts == 'a')
    opos = find(opts == 'a');
    if opos(1) < numel(opts)
        opos = opts(opos(1)+1:end);
        oposa = find(opos < 48 | opos > 57);
        if isempty(oposa)
            doafter = str2double(opos);
        elseif oposa(1) > 1
            doafter = str2double(opos(1:oposa(1)-1));
        end
    end
end

% check option: b(efore)
dobefore = 0;
if any(opts == 'b')
    opos = find(opts == 'b');
    if opos(1) < numel(opts)
        opos = opts(opos(1)+1:end);
        oposa = find(opos < 48 | opos > 57);
        if isempty(oposa)
            dobefore = str2double(opos);
        elseif oposa(1) > 1
            dobefore = str2double(opos(1:oposa(1)-1));
        end
    end
end

% check option: h(ead before grep)
dohead = -1;
if any(opts == 'h')
    opos = find(opts == 'h');
    if opos(1) < numel(opts)
        opos = opts(opos(1)+1:end);
        oposa = find(opos < 48 | opos > 57);
        if isempty(oposa)
            dohead = str2double(opos);
        elseif oposa(1) > 1
            dohead = str2double(opos(1:oposa(1)-1));
        end
    end
end

% check option: H(ead after grep)
doheada = -1;
if any(opts == 'H')
    opos = find(opts == 'H');
    if opos(1) < numel(opts)
        opos  = opts(opos(1)+1:end);
        oposa = find(opos < 48 | opos > 57);
        if isempty(oposa)
            doheada = str2double(opos);
        elseif oposa(1) > 1
            doheada = str2double(opos(1:oposa(1)-1));
        end
    end
end

% check option: t(ail before grep)
dotail = -1;
if any(opts == 't')
    opos = find(opts == 't');
    if opos(1) < numel(opts)
        opos = opts(opos(1)+1:end);
        oposa = find(opos < 48 | opos > 57);
        if isempty(oposa)
            dotail = str2double(opos);
        elseif oposa(1) > 1
            dotail = str2double(opos(1:oposa(1)-1));
        end
    end
end

% check option: T(ail after grep)
dotaila = -1;
if any(opts == 'T')
    opos = find(opts == 'T');
    if opos(1) < numel(opts)
        opos = opts(opos(1)+1:end);
        oposa = find(opos < 48 | opos > 57);
        if isempty(oposa)
            dotaila = str2double(opos);
        elseif oposa(1) > 1
            dotaila = str2double(opos(1:oposa(1)-1));
        end
    end
end

% check other options: case insensity, lines, enumberation
% and between options, reverse logic, regexp pattern match
if any(opts == 'i')
    dolower = 1;
else
    dolower = 0;
end
if any(opts == 'l')
    dolines = 0;
else
    dolines = 1;
end
if any(opts == 'n')
    doenum = 1;
else
    doenum = 0;
end
if any(opts == 'r')
    doand = 0;
else
    doand = 1;
end
if any(opts == 'v')
    doreverse = 1;
else
    doreverse = 0;
end
if any(opts == 'x')
    doregexp = 1;
else
    doregexp = 0;
end

% calculate complex option value
opts = 4 * doregexp + 2 * dolower + doreverse;

% preset output
mlines = [];

% lower patterns?
if dolower
    for pat = 1:pats
        patterns{pat} = lower(patterns{pat});
    end
end

% heading and/or tailing before grepping
if dohead > -1
    togrep = togrep(1:min(lines,dohead));
    lines = numel(togrep);
end
if dotail > -1
    togrep = togrep((lines+1)-min(lines,dotail):lines);
    lines = numel(togrep);
end

% iterate over lines
for line = 1:lines

    % ANDing patterns
    if doand

        % default 1
        match = 1;

        % which complex option
        % then check pattern and say no match if any failes
        switch (opts), case {0}
            for pat = 1:pats
                if  isempty(strfind(togrep{line}, patterns{pat}))
                    match = 0;
                    break;
                end
            end
        case {1}
            for pat = 1:pats
                if ~isempty(strfind(togrep{line}, patterns{pat}))
                    match = 0;
                    break;
                end
            end
        case {2}
            for pat = 1:pats
                if  isempty(strfind(lower(togrep{line}), patterns{pat}))
                    match = 0;
                    break;
                end
            end
        case {3}
            for pat = 1:pats
                if ~isempty(strfind(lower(togrep{line}), patterns{pat}))
                    match = 0;
                    break;
                end
            end
        case {4}
            for pat = 1:pats
                if  isempty(regexp(togrep{line}, patterns{pat}, 'once'))
                    match = 0;
                    break;
                end
            end
        case {5}
            for pat = 1:pats
                if ~isempty(regexp(togrep{line}, patterns{pat}, 'once'))
                    match = 0;
                    break;
                end
            end
        case {6}
            for pat = 1:pats
                if  isempty(regexpi(togrep{line}, patterns{pat}, 'once'))
                    match = 0;
                    break;
                end
            end
        case {7}
            for pat = 1:pats
                if ~isempty(regexpi(togrep{line}, patterns{pat}, 'once'))
                    match = 0;
                    break;
                end
            end
        end

    % ORing patterns
    else

        % default 0
        match = 0;

        % which complex option
        % then check pattern and say match if any is true
        switch (opts), case {0}
            for pat = 1:pats
                if ~isempty(strfind(togrep{line}, patterns{pat}))
                    match=1;
                    break;
                end
            end
        case {1}
            for pat = 1:pats
                if  isempty(strfind(togrep{line}, patterns{pat}))
                    match=1;
                    break;
                end
            end
        case {2}
            for pat = 1:pats
                if ~isempty(strfind(lower(togrep{line}), patterns{pat}))
                    match=1;
                    break;
                end
            end
        case {3}
            for pat = 1:pats
                if  isempty(strfind(lower(togrep{line}), patterns{pat}))
                    match=1;
                    break;
                end
            end
        case {4}
            for pat = 1:pats
                if ~isempty(regexp(togrep{line}, patterns{pat}, 'once'))
                    match=1;
                    break;
                end
            end
        case {5}
            for pat = 1:pats
                if  isempty(regexp(togrep{line}, patterns{pat}, 'once'))
                    match=1;
                    break;
                end
            end
        case {6}
            for pat = 1:pats
                if ~isempty(regexpi(togrep{line}, patterns{pat}, 'once'))
                    match=1;
                    break;
                end
            end
        case {7}
            for pat = 1:pats
                if  isempty(regexpi(togrep{line}, patterns{pat}, 'once'))
                    match=1;
                    break;
                end
            end
        end
    end

    % match found
    if match
        mlines(end+1) = line;
        if dobefore > 0
            mlines = union(mlines, max(1,line - dobefore):(line - 1));
        end
        if doafter  > 0
            mlines = union(mlines, (line + 1):min(lines, line + doafter));
        end
    end
end

% line enumberation
if dolines
    if doenum
        for ml = mlines
            togrep{ml} = sprintf('%s%d: %s', tagname, ml, togrep{ml});
        end
    end
    grepped = togrep(mlines);
else
    grepped = mlines;
    return;
end

% head or tail after grepping
if doheada > -1
    lines = numel(grepped);
    grepped = grepped(1:min(lines,doheada));
end
if dotaila > -1
    lines = numel(grepped);
    grepped = grepped((lines + 1) - min(lines, dotaila):lines);
end

% result as char
if ichar
    grepped = gluetostringc(grepped);
end

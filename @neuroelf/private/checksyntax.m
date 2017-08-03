function [errormsg] = checksyntax(snippet)
% checksyntax  - checks a snippet of code for syntax validity
%
% FORMAT:       [errmsg] = checksyntax(snippet)
%       OR      [errmsg] = checksyntax(filename)
%
% Input fields:
%
%       snippet     CR/LF separated lines of code (or single line)
%       filename    filename of the M-file to test
%
% Output fields:
%
%       errmsg      if requested (nargout > 0), any error detected
%                   by the MATLAB parser is returned, otherwise
%                   printed to stdout
%
% See also eval, try.

% Version:  v1.1
% Build:    16011612
% Date:     Jan-16 2016, 12:22 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
if nargin < 1 || ...
   ~ischar(snippet) || ...
    numel(snippet) ~= size(snippet, 2)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument' ...
    );
end

% for filenames (of existing files)
if ~any(snippet == ';' | snippet == ',' | snippet == ' ') && ...
    exist(snippet, 'file') == 2

    % read and split to lines, remove leading and trailing spaces
    filecont = ddeblank(splittocellc(asciiread(snippet), char([13, 10]), true, true));

    % remove empty and comment lines
    filecont(cellfun('isempty', filecont)) = [];
    filecont(~cellfun('isempty', regexpi(filecont, '^\%'))) = [];

    % replace initial "function\w+" with special token
    filecont = regexprep(filecont, '^function\s+', ['FUNC' 'TION_DEF__ ']);

    % parse functions specifically
    snippets  = splittocellc(gluetostringc(filecont, char(10)), ...
        ['FUNC' 'TION_DEF__ ']);
    ierrormsg = '';
    hasfuncs  = '';

    % for regular functions, the first snippet will be empty
    if ~isempty(snippets) && ...
       isempty(snippets{1})
        snippets(1) = [];
        hasfuncs = 'function ';

    % batch M-file (no sub-functions)
    elseif ~isempty(snippets) && ...
        numel(snippets) == 1
        ierrormsg = checksyntax(snippets{1});
        snippets  = cell(0);

    % function defitions within a non-function M-file (class?)
    elseif numel(snippets) > 1

        % split first part (before functions)
        firstsnip = splittocellc(snippets{1}, char([13, 10]), true, true);
        for lc = 1:numel(firstsnip)
            iscomm = find(firstsnip{lc} == '%');
            if ~isempty(iscomm)
                firstsnip{lc} = firstsnip{lc}(1:iscomm(1)-1);
            end
            if ~isempty(deblank(firstsnip{lc}))
                ierrormsg = ...
                    'A function declaration may not appear in a script M-file.';
                snippets  = {''};
                break;
            end
        end
        snippets = snippets(2:end);

    % otherwise nothing to check
    else
        snippets = cell(0);
    end

    % check snippets
    for snc = 1:numel(snippets)
        ierrormsg = checksyntax([hasfuncs snippets{snc}]);
        if ~isempty(ierrormsg)
            sniplines = splittocellc(snippets{snc}, char([13, 10]), true, true);
            ierrormsg = [ierrormsg ' (in function ' sniplines{1} ')'];
            break;
        end
    end

% code directly given (also used for filewise calls)
else

    % with functions
    ffound = strfind(snippet, 'function');

    % function found and at beginning
    if ~isempty(ffound) && ...
        ffound(1) < 2

        % split to lines
        snippet = splittocellc(snippet, char([13, 10]), true, true);
        lc = 1;

        % go to next function
        while isempty(strfind(snippet{lc}, 'function'))
            lc = lc + 1;
        end

        % get function line and snippet
        fline   = snippet{lc};
        snippet = gluetostringc(snippet((lc+1):end), char(10));

        % find argument delimiters
        flinea  = find(fline == '(');
        flinee  = find(fline == ')');

        % check for valid function declaration
        if ~isempty(flinea) && ...
           ~isempty(flinee) && ...
            flinee(end) > flinea(end)

            % get arguments
            fargs = splittocellc(fline((flinea(end)+1):(flinee(end)-1)), ...
                [', ' char(9)], true, true);

            % check arguments
            for fac = numel(fargs):-1:1
                if isempty(fargs{fac})
                    fargs(fac) = [];
                end
            end
        else

            % otherwise say 0 arguments
            fargs = cell(0);
        end

        % check arguments syntax by setting all args to [];
        snippet = ['varargin=cell(0);varargout=cell(0);' sprintf('%s=[];',fargs{:}) char(10) ...
                   snippet char(10)  'clear varargin varargout' sprintf(' %s',fargs{:}) ';'];
    end

    % try the actual code by never true if
    try
        eval(['if 0==1, ' snippet char(10) char(10) ' end']);
        ierrormsg = '';
    catch ne_eo;
        ierrormsg = ne_eo.message;
    end
end

% only print error if raised
if nargout > 0
    errormsg = ierrormsg;
else
    disp(['Reported: ' ierrormsg]);
end


% internal function to remove comment
function cremoved = rcomment(cremoved)
    pct = find(cremoved ~= char(9) & cremoved ~= char(32));
    if isempty(pct) || ...
        cremoved(pct(1)) == '%'
        cremoved = '';
    else
        pcc = find(cremoved == '%');
        if isempty(pcc)
            return;
        end
        if ~any(cremoved == '''')
            cremoved = cremoved(1:pcc(1)-1);
        end
    end
% end of function cremoved

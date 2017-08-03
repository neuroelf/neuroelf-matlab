function parsed = tfgparse(varargin)
% tfgparse  - parses a xfigure TFG file into a struct
%
% FORMAT:       TFGstruct = tfgparse(filename [, options])
%         or
%               TFGstruct = tfgparse(TFGstruct, options)
%
% Input fields:
%
%       filename    name of the TFG file to parse
%       TFGstruct   if given, write to TFG file or evaluate contents
%       options     struct with optional fields
%        .check     runs checksyntax on all Callback* fields
%        .evaluate  evaluate VARIABLES and *,] fields
%
% The second format creates the file specified by filename from the
% contents of TFGstruct.

% Version:  v1.1
% Build:    16042016
% Date:     Apr-20 2016, 4:11 PM EST
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

% global methods
global ne_methods;

% persistent memory
persistent tfgparser;
if isempty(tfgparser)

    % initialize struct
    tfgparser = struct;

    % callback functions to check
    tfgparser.cbcheck = { ...
        'callback', ...
        'callbackclreq', ...
        'callbackdblclick' ...
    };

    % default options
    tfgparser.defopts = { ...
        'check',    'logical', {false, true}, false; ...
        'evaluate', 'logical', {false, true}, true   ...
    };
end

% argument check
if nargin < 1 || (~ischar(varargin{1}) && ~isstruct(varargin{1})) || isempty(varargin{1})
    error('neuroelf:general:badArgument', 'Invalid number/type of arguments.');
end
if nargin > 1 && isstruct(varargin{2})
    options = checkstruct(varargin{2}, tfgparser.defopts, true);
else
    options = checkstruct(struct, tfgparser.defopts);
end

% try to read TFG file
if ischar(varargin{1})
    filename = varargin{1}(:)';
    if exist(filename, 'file') ~= 2
        error('neuroelf:general:badArgument', 'TFG File not found.');
    end
    try
        [figtext, nlines] = splittocell(asciiread(filename), char([10, 13]), 1, 1);
    catch ne_eo;
        error('neuroelf:general:badArgument', ...
            'TFG file not readable: %s.', ne_eo.message);
    end
else
    filename = '';
end

% parsing file
if ~isempty(filename)

    % initialize output structure, just in case...
    parsed = struct( ...
        'POPTIONS',     struct('CommentFieldWidth', 16), ...
        'COMMENTS',     struct, ...
        'VARIABLES',    emptystruct({'VarName', 'VarContent'}), ...
        'FIGURE',       emptystruct({'Position'}), ...
        'UICONTROLS',   emptystruct({'Position', 'Tag', 'Type'}), ...
        'UIRESIZE',     emptystruct({'Tag', 'Reference', 'RelPosition'}), ...
        'MENU',         emptystruct({'Callback', 'Caption', 'Level', 'Tag'}), ...
        'CONTEXTMENUS', emptystruct({'IsCM', 'Callback', 'Caption', 'Level', 'Tag'}));
    figsbegins = zeros(1, nlines);
    figsends   = zeros(1, nlines);
    figsnames  = cell(1, nlines);

    % iterate through lines
    fc = 0;
    lc = 0;
    while lc < nlines
        lc = lc + 1;

        % section begin ?
        lcf = strfind(figtext{lc}, 'BEGIN_');
        if ~isempty(lcf)

            % get section name tag
            psecbegin = lc;
            psecnamep = lcf(1) + 6;
            psectname = figtext{psecbegin}(psecnamep:end);

            % valid name tag ?
            tr = find(psectname < 65 | (psectname > 90 & psectname < 97) | ...
                psectname > 122);
            if ~isempty(tr)
                psectname(tr(1):end) = [];
            end
            if isempty(psectname)
                continue;
            end

            % find section end
            lec = lc;
            while lec < nlines
                lec = lec + 1;
                lecf = strfind(figtext{lec}, ['END_' psectname]);
                if ~isempty(lecf)
                    break;
                end
            end

            % if we're beyond EOF, don't use it!
            if lec > nlines
                lc = nlines + 1;
                continue;
            end

            % accept section
            psecend = lec;
            fc = fc + 1;
            figsbegins(fc) = psecbegin;
            figsends(fc) = psecend;
            figsnames{fc} = psectname;
            lc = lec;
        end
    end
    figsbegins = figsbegins(1:fc);
    figsends = figsends(1:fc);
    figsnames = figsnames(1:fc);

    % set final end marker
    figsbegins(end+1) = nlines + 1;

    % iterate over sections
    readopts = struct('convert', 'deblank', 'headline', '');
    numsects = length(figsnames);
    for sc = 1:numsects
        lrb = figsbegins(sc);
        lre = figsends(sc);
        sectname = figsnames{sc};

        % handle comments specially
        if strcmpi(sectname, 'COMMENTS')
            for lc = (lrb+1):(lre-1)
                [comv, comc] = splittocell(figtext{lc}, ':');
                if comc > 2
                    comv{2} = gluetostringc(comv(2:end), ':');
                end
                if comc < 2
                    comv{2} = '';
                end
                while numel(comv{2} > 0) && comv{2}(1) == ' '
                    comv{2}(1) = [];
                end
                parsed.COMMENTS.(makelabel(comv{1})) = comv{2};
            end
            parsed.COMMENTS.ALLCOMMENTS = ...
                gluetostringc(figtext(lrb:lre), char(10), 1);

        % all other sections here
        else
            parsed.(upper(sectname)) = ...
                acsvread(figtext(lrb+1:lre-1), '|', readopts);
        end
    end

    % numbers are parsed definitely
    fn = fieldnames(parsed);
    for fc1 = 1:numel(fn)
        fp = parsed.(fn{fc1});
        sfn = fieldnames(fp);
        for sc = 1:numel(fp)
            for fc2 = 1:numel(sfn)
                if ~isempty(fp(sc).(sfn{fc2})) && fp(sc).(sfn{fc2})(1) == '$'
                    try
                        fp(sc).(sfn{fc2}) = eval(['[' fp(sc).(sfn{fc2})(2:end) ']']);
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        fp(sc).(sfn{fc2}) = [];
                    end
                end
            end
        end
        parsed.(fn{fc1}) = fp;
    end

% otherwise
else
    parsed = varargin{1};
end

% evaluate the contents of the parsed TFG
if options.evaluate
    if isfield(parsed, 'VARIABLES') && isstruct(parsed.VARIABLES) && ...
        isfield(parsed.VARIABLES, 'VarName') && isfield(parsed.VARIABLES, 'VarContent')
        tvars = parsed.VARIABLES;
        assignin('base', 'tfgtv', struct);
        assignin('base', 'tfgne', ne_methods);

        % iterate over loaded variables
        for tvar = 1:numel(tvars)

            % no special character sequence (leading '*')
            if ~ischar(tvars(tvar).VarContent) || isempty(tvars(tvar).VarContent) || ...
               ~any('*$]' == tvars(tvar).VarContent(1))
                try
                    evalin('base', ['tfgtv.' tvars(tvar).VarName '=' ...
                         tvars(tvar).VarContent]);
                catch ne_eo;
                    error('neuroelf:general:evaluationError', ...
                        'Error evaluating expression for %s (%s).', ...
                        tvars(tvar).VarName, ne_eo.message);
                end

            % evaluation
            elseif any(']$' == tvars(tvar).VarContent(1))
                evalin('base', ['tfgtv.' tvars(tvar).VarName '=' ...
                    tvars(tvar).VarContent(2:end) ';'], '');

            % substitution
            else
                evalin('base', ['tfgtv.' tvars(tvar).VarName '=' ...
                    'tfgtv.' tvars(tvar).VarContent(2:end) ';'], '');
            end
        end
    end
    fn = fieldnames(parsed);
    for fc1 = 1:numel(fn)
        if strcmpi(fn{fc1}, 'variables')
            continue;
        end
        fp = parsed.(fn{fc1});
        sfn = fieldnames(fp);
        for sc = 1:numel(fp)
            for fc2 = 1:numel(sfn)
                if ischar(fp(sc).(sfn{fc2})) && ~isempty(fp(sc).(sfn{fc2})) && ...
                    any('$]*' == fp(sc).(sfn{fc2})(1))
                    switch (fp(sc).(sfn{fc2})(1))
                        case {'$', ']'}
                            fp(sc).(sfn{fc2}) = evalin('base', ['[' ...
                                fp(sc).(sfn{fc2})(2:end) ']'], '[]');
                        case '*'
                            fp(sc).(sfn{fc2}) = evalin('base', ['[tfgtv.' ...
                                fp(sc).(sfn{fc2})(2:end) ']'], '[]');
                    end
                end
            end
        end
        parsed.(fn{fc1}) = fp;
    end
    try
        evalin('base', 'clear tfgne tfgtv;', '');
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

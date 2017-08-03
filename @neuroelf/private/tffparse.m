function tffspec = tffparse(tffcont)
% tffparse  - parse TFF text file format description
%
% FORMAT:       tffspec = tffparse(tffcont)
%
% Input fields:
%
%       tffcont     1xN char: TFF specification file name or content
%
% Output fields:
%
%       tffspec     1x1 struct with TFF specification
%
% See also tffdocu, tffio

% Version:  v1.0
% Build:    16011416
% Date:     Jan-14 2016, 4:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

tffversion = 'v0.9c';

% argument check
if ...
    nargin < 1 || ...
   ~ischar(tffcont) || ...
    isempty(tffcont)
    error( ...
        'xff:BadArgument', ...
        'Invalid TFF specification argument given.' ...
    );
end

% check for line feeds and try to read file
if numel(tffcont) ~= length(tffcont) || ...
   ~any(tffcont == char(10))
    if exist(tffcont, 'file') ~= 2
        error( ...
            'xff:FileNotFound', ...
            'Either invalid content or file not found.' ...
        );
    end

    % read tffspec
    try
        tffcont = asciiread(tffcont);
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% check content
if isempty(strfind(tffcont, 'TextFileFormat'))
    error( ...
        'xff:BadTFFSpec', ...
        'Invalid TFF specification file/content given.' ...
    );
end

% split into lines
tfflines = splittocellc(tffcont(:)', char([10, 13]), true, true);

% remove comments, deblank, and then discard empty lines
tfflines = regexprep(tfflines, '\#.*$', '');
tfflines = deblank(tfflines);
tfflines(cellfun('isempty', tfflines)) = [];

% initialize content
tffspec = struct;
tffspec.FFTYPE             = 'TFF';
tffspec.TFFVERSION         = tffversion;
tffspec.AfterReadCode      = '';
tffspec.ArrayFormat        = '%-12s';
tffspec.BeforeWriteCode    = '';
tffspec.BinaryIO           = false;
tffspec.CheckForNaNs       = false;
tffspec.CustomDelimiters   = {};
tffspec.DefaultProperty    = {};
tffspec.Description        = 'All files';
tffspec.Extensions         = {};
tffspec.FieldDelimCollapse = true;
tffspec.FieldDelimiters    = {char(9)};
tffspec.Filetype           = 'Custom TFF file';
tffspec.FilenameMatch      = {};
tffspec.LineDelimiters     = {char([13,10]), char(10)};
tffspec.ListOfFields       = emptystruct({});
tffspec.Loops              = struct;
tffspec.Magic              = emptystruct({});
tffspec.NewFileCode        = '';
tffspec.ParagraphArrays    = true;
tffspec.SkipEmptyLines     = true;
tffspec.ValidFileCode      = '';
tffspec.Variables          = struct;

% parse content
lc  = 1;
llc = length(tfflines);
while (lc <= llc)

    % regexpi line matches '<TOKEN>:<VALUE>' ?
    tffline = tfflines{lc};
    [tlmatcht{1:3}] = regexpi( ...
        tffline, '^([a-z][a-z_0-9]*)\:\s*(.*)\s*$');
    tlmatcht = tlmatcht{3};

    % continue with next line if no match
    lc = lc + 1;
    if isempty(tlmatcht)
        continue;
    end

    % handle those sections/keywords:
    tffkey = tffline(tlmatcht{1}(1, 1):tlmatcht{1}(1, 2));
    tffval = tffline(tlmatcht{1}(2, 1):tlmatcht{1}(2, 2));
    switch (lower(tffkey))

        % after read code snippet
        case {'afterreadcode'}

            % find "EndAfterReadCode"
            endline = llc + 1;
            for slc = lc:llc
                if ~isempty(regexpi(tfflines{slc}, '^endafterreadcode'))
                    endline = slc - 1;
                    break;
                end
            end

            % check endline
            if endline > llc
                error( ...
                    'xff:BadTFFSpec', ...
                    'Unclosed AfterReadCode section in TFF spec.' ...
                );
            end

            % generate and check code snippet
            codesnip = gluetostringc(tfflines(lc:endline), char(10), true);
            codesnip = regexprep(codesnip, '\.{3,}\s+', ' ');
            tcodesnip = codesnip;
            tcodesnip(tcodesnip == '$' | tcodesnip == '@') = [];
            try
                eval(['if 1==0,' tcodesnip ';end']);
            catch ne_eo;
                error( ...
                    'xff:BadTFFSpec', ...
                    'Syntax error detected in AfterReadCode: ''%s''.', ...
                    ne_eo.message ...
                );
            end

            % store code
            tffspec.AfterReadCode = ...
                tff_parsecode(codesnip, 'namevars', 'tffcont');

        % array element format (sprintf)
        case {'arrayformat'}

            % check format
            testformat = sprintf(tffval, sprintf('%.3f', pi));
            if isempty(regexp(testformat, '3\.14', 'once'))
                error( ...
                    'xff:BadTFFSpec', ...
                    'Bad ArrayFormat specified: ''%s''.', ...
                    tffval ...
                );
            end

            % set format
            tffspec.ArrayFormat = tffval;

        % before write code snippet
        case {'beforewritecode'}

            % find "EndBeforeWriteCode"
            endline = llc + 1;
            for slc = lc:llc
                if ~isempty(regexpi(tfflines{slc}, '^endbeforewritecode'))
                    endline = slc - 1;
                    break;
                end
            end

            % check endline
            if endline > llc
                error( ...
                    'xff:BadTFFSpec', ...
                    'Unclosed BeforeWriteCode section in TFF spec.' ...
                );
            end

            % generate and check code snippet
            codesnip  = gluetostringc(tfflines(lc:endline), char(10), true);
            tcodesnip = codesnip;
            tcodesnip(tcodesnip == '$' | tcodesnip == '@') = [];
            try
                eval(['if 1==0,' tcodesnip ';end']);
            catch ne_eo;
                error( ...
                    'xff:BadTFFSpec', ...
                    'Syntax error detected in BeforeWriteCode: ''%s''.', ...
                    ne_eo.message ...
                );
            end

            % store code
            tffspec.BeforeWriteCode = ...
                tff_parsecode(codesnip, 'namevars', 'tffcont');

        % treat file as binary
        case {'binaryio'}
            try
                tffspec.BinaryIO = logical(eval(tffval));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xff:BadTFFSpec', ...
                    'Invalid BinaryIO value given: ''%s''.', ...
                    tffval ...
                );
            end

        % treat file as binary
        case {'checkfornans'}
            try
                tffspec.CheckForNaNs = logical(eval(tffval));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xff:BadTFFSpec', ...
                    'Invalid CheckForNaNs value given: ''%s''.', ...
                    tffval ...
                );
            end

        % custom delimiters
        case {'customdelimiters'}
            try

                % evaluate list
                eval(['tffspec.CustomDelimiters=' tffval ';']);

                % check list syntax
                if isempty(tffspec.CustomDelimiters) || ...
                   ~iscell(tffspec.CustomDelimiters)
                    tffspec.CustomDelimiters = {};
                    warning( ...
                        'xff:BadTFFSpec', ...
                        'Empty or malformed CustomDelimiters: ''%s''.', ...
                        tffval ...
                    );
                end

                % check list content
                delimlst = tffspec.CustomDelimiters(:)';
                for clcc = 1:length(delimlst)
                    if ~ischar(delimlst{clcc})
                        delimlst{clcc} = char(delimlst{clcc}(:)');
                    end
                end

                % put back to struct
                tffspec.CustomDelimiters = delimlst;

            % on error
            catch ne_eo;

                % bail out
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xff:BadTFFSpec', ...
                    'Invalid CustomDelimiters list: ''%s''.', ...
                    tffval ...
                );
            end

        % default property (for access)
        case {'defaultproperty'}
            if ~isempty(tffval) && ...
                tffval(1) == '{' && ...
                tffval(end) == '}'
                try
                    tffval = eval(tffval);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    error( ...
                        'xff:BadTFFSpec', ...
                        'Invalid DefaultProperty value given: ''%s''.', ...
                        tffval ...
                    );
                end
            end
            if ~iscell(tffval)
                tffval = {tffval};
            end
            tffspec.DefaultProperty = tffval;

        % file open description
        case {'description'}
            tffspec.Description = splittocellc(tffval, ';,', true, true);

        % valid extensions list
        case {'extensions'}
            tffspec.Extensions = splittocellc(tffval, ';,. ', true, true);

        % collapse field delimiters
        case {'fielddelimcollapse'}
            try
                tffspec.FieldDelimCollapse = logical(eval(tffval));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xff:BadTFFSpec', ...
                    'Invalid FieldDelimCollapse value given: ''%s''.', ...
                    tffval ...
                );
            end

        % field delimiters
        case {'fielddelimiters'}
            try

                % evaluate list
                eval(['tffspec.FieldDelimiters=' tffval ';']);

                % check list syntax
                if isempty(tffspec.FieldDelimiters) || ...
                   ~iscell(tffspec.FieldDelimiters)
                    tffspec.FieldDelimiters = {char(9)};
                    warning( ...
                        'xff:BadTFFSpec', ...
                        'Empty or malformed FieldDelimiters: ''%s''.', ...
                        tffval ...
                    );
                end

                % check list content
                delimlst = tffspec.FieldDelimiters(:)';
                for clcc = 1:length(delimlst)
                    if ~ischar(delimlst{clcc})
                        delimlst{clcc} = char(delimlst{clcc}(:)');
                    end
                end

                % put back to struct
                tffspec.FieldDelimiters = delimlst;

            % on error
            catch ne_eo;

                % bail out
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xff:BadTFFSpec', ...
                    'Invalid FieldDelimiters list: ''%s''.', ...
                    tffval ...
                );
            end

        % valid extensions list
        case {'filenamematch'}
            tffspec.FilenameMatch = splittocellc(tffval, ';, ', true, true);

        % line delimiters
        case {'linedelimiters'}
            try

                % evaluate list
                eval(['tffspec.LineDelimiters=' tffval ';']);

                % check list syntax
                if isempty(tffspec.LineDelimiters) || ...
                   ~iscell(tffspec.LineDelimiters)
                    tffspec.LineDelimiters = {char([13, 10]),char(9)};
                    warning( ...
                        'xff:BadTFFSpec', ...
                        'Empty or malformed LineDelimiters: ''%s''.', ...
                        tffval ...
                    );
                end

                % check list content
                delimlst = tffspec.LineDelimiters(:)';
                for clcc = 1:length(delimlst)
                    if ~ischar(delimlst{clcc})
                        delimlst{clcc} = char(delimlst{clcc}(:)');
                    end
                end

                % put back to struct
                tffspec.LineDelimiters = delimlst;

            % on error
            catch ne_eo;

                % bail out
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xff:BadTFFSpec', ...
                    'Invalid LineDelimiters list: ''%s''.', ...
                    tffval ...
                );
            end

        % list of fields
        case {'listoffields'}

            % find "EndListOfFields" line
            endline = llc + 1;
            for slc = (lc + 1):llc
                if ~isempty(regexpi(tfflines{slc}, '^endlistoffields'))
                    endline = slc - 1;
                    break;
                end
            end
            blc = lc;
            lc = endline + 1;

            % invalid ListOfFields
            if lc > llc
                error( ...
                    'xff:BadTFFSpec', ...
                    'Unclosed or bad ListOfFields block.' ...
                );
            end

            % list/table field separator
            listsep = tffval;

            % get rule headers
            rhead = splittocellc(tfflines{blc}, listsep, false);

            % build header struct
            hstruct = struct;
            for hfc = 1:length(rhead)
                hfield = lower(makelabel(rhead{hfc}));
                hstruct.(hfield) = hfc;
                rhead{hfc} = hfield;
            end
            blc = blc + 1;

            % bail out if invalid header given
            if ...
               ~isfield(hstruct, 'type') || ...
               ~isfield(hstruct, 'cond') || ...
               ~isfield(hstruct, 'field') || ...
               ~isfield(hstruct, 'datatype') || ...
               ~isfield(hstruct, 'format') || ...
               ~isfield(hstruct, 'dim') || ...
               ~isfield(hstruct, 'default') || ...
               ~isfield(hstruct, 'varname')
                warning( ...
                    'xff:BadTFFSpec', ...
                    'ListOfFields with bad headers.' ...
                );
                continue;
            end

            % get list of header field names (in their order)
            hfields = fieldnames(hstruct);
            nfields = numel(hfields);

            % build empty field list struct
            frules = emptystruct(hfields);

            % check global fieldlist struct
            tfields = fieldnames(tffspec.ListOfFields);

            % no rules yet -> OK
            if isempty(tfields)
                tffspec.ListOfFields = frules;

            % otherwise compare header fields
            else

                % check number and content of fields
                if length(tfields) ~= length(hfields) || ...
                   ~all(strcmp(tfields(:), hfields(:)))

                    % give warning and continue with next block
                    error( ...
                        'xff:BadTFFSpec', ...
                        'ListOfFields blocks must match in their headers.' ...
                    );
                end
            end

            % build list of rules to consider
            actrules = frules;
            slc = blc;
            while slc <= endline

                % split to fields
                rcont = splittocellc(tfflines{slc}, listsep, false);

                % increase counter
                slc = slc + 1;

                % continuation?
                while ~isempty(rcont) && ...
                    numel(rcont{end}) > 2 && ...
                    strcmp(rcont{end}(end-2:end), '...') && ...
                    slc <= endline

                    % split the next line also
                    rcontnext = splittocellc(tfflines{slc}, listsep, false);

                    % increase counter
                    slc = slc + 1;

                    % and join contents
                    rcont{end} = [rcont{end}(1:end-3) ddeblank(rcontnext{1})];
                    rcont = [rcont, rcontnext(2:end)];

                end

                % reject too short arrays
                if numel(rcont) ~= numel(rhead)
                    if numel(rcont) > 1
                        warning( ...
                            'xff:TFFSpecError', ...
                            'Invalid number of fields in line.' ...
                        );
                    end
                    continue;
                end

                % deal into struct
                rstruct = cell2struct(deblank(rcont(:)), hfields, 1);

                % deal with empty lines
                rstrtype = lower(rstruct.type);
                if isempty(rstrtype)
                    continue;
                elseif ...
                   ~any(strcmp(rstrtype, {'skipn', 'wrtln'})) && ...
                    isempty(rstruct.varname)
                    warning( ...
                        'xff:AmbiguousTFFSpec', ...
                        'Non recognized TYPE token: ''%s''.', ...
                        rstruct.type ...
                    );
                    continue;
                end

                % put non-empty type/varname rules into actrules
                actrules(end + 1) = rstruct;

            end

            % init loop checking variables and flist array
            looplist = {};
            loopused = {};

            % parse rules (using for; if any rules fail, error out!)
            for slc = 1:numel(actrules)

                % get rstruct from actrules
                rstruct = actrules(slc);
                rstrtype = lower(rstruct.type);

                % copy varname into empty fields
                if isempty(rstruct.field)
                    rstruct.field = rstruct.varname;
                    actrules(slc).field = rstruct.field;
                end

                % check syntax of fields, type
                if ~any(strcmp(rstrtype, ...
                    {'bloop', 'eloop', 'expre', 'skipn', 'xloop', ...
                     'array', 'field', 'flist', 'wrtln'}))
                    error( ...
                        'xff:BadTFFSpec', ...
                        'Invalid ListOfFields.type token: ''%s''.', ...
                        rstruct.type ...
                    );
                end

                % ... cond
                if ~isempty(rstruct.cond)
                    try
                        eval([ ...
                            'if 1==0,if ' ...
                            strrep(strrep(rstruct.cond, '$', ''), '@', '') ...
                            ',disp(1);end,end']);
                    catch ne_eo;
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid ListOfFields.cond content: ''%s'' (%s).', ...
                            rstruct.cond, ne_eo.message ...
                        );
                    end
                    rstruct.cond = ...
                        tff_parsecode(rstruct.cond, 'namevars', 'tffcont');
                end

                % ... field
                if ...
                    any(strcmp(rstrtype, ...
                        {'array', 'field', 'flist'})) && ...
                    isempty(regexpi(rstruct.field, '^[a-z][a-z_0-9\-]*$'))
                    error( ...
                        'xff:BadTFFSpec', ...
                        'Invalid ListOfFields.field name: ''%s''.', ...
                        rstruct.field ...
                    );
                end

                % ... datatype / format
                if any(strcmp(rstrtype, {'array', 'field', 'flist'})) && ...
                    (isempty(regexpi(rstruct.datatype, '^[a-z][a-z_0-9]+$')) || ...
                     ~any(rstruct.format == '%'))
                    error( ...
                        'xff:BadTFFSpec', ...
                        ['Invalid ListOfFields.datatype/format tag: ' ...
                         '''%s''/''%s''.'], ...
                        rstruct.datatype, rstruct.format ...
                    );
                end

                % ... dim (loops, fields and flists)
                if any (strcmp(rstrtype, {'bloop', 'field', 'flist'}))

                    % only single number OR variable
                    if isempty(regexpi(rstruct.dim, ...
                        '^(\d+|[\$\@][a-z][a-z_0-9\.]*(\((\d+|[\$\@][a-z][a-z_0-9\.]*)\))?)$'))
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid LOOP/FIELD/FLIST.DIM given: ''%s''.', ...
                            rstruct.dim ...
                        );
                    end

                % ... dim (arrays)
                elseif strcmp(rstrtype, 'array')

                    % multiple numbers AND/OR variables
                    if isempty(regexpi(rstruct.dim, ...
                        ['^((\d+|[\$\@][a-z][a-z_0-9]*)(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?' ...
                         '(\.[a-z][a-z_0-9]*(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?)*\,\s*)*' ...
                         '((\d+|[\$\@][a-z][a-z_0-9]*)(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?)$']))
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid ARRAY.dim given: ''%s''.', ...
                            rstruct.dim ...
                        );
                    end

                % ... dim (skipn)
                elseif strcmp(rstrtype, 'skipn')

                    % only a 1-D numeric value accepted
                    if isempty(regexp(rstruct.dim, '^\-?\d+$', 'once'))
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid SKIPN.dim given: ''%s''.', ...
                            rstruct.dim ...
                        );
                    end
                end
                rstruct.dim = ...
                    tff_parsecode(rstruct.dim, 'namevars', 'tffcont');
                if all(rstruct.dim == ',' | rstruct.dim == ' ' | ...
                      (rstruct.dim >= '0' & rstruct.dim <= '9'))
                    try
                        rstruct.dim = eval(['[' rstruct.dim ']']);
                    catch ne_eo;
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid numeric dimension given: %s.', ...
                            ne_eo.message ...
                        );
                    end
                end

                % ... default (if non-empty)
                if ~isempty(rstruct.default)
                    try
                        eval(['if 1==0,' ...
                              strrep(strrep(rstruct.default, '$', ''), '@', '') ',end']);
                    catch ne_eo;
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid ListOfFields.default value: ''%s'' (%s).', ...
                            rstruct.default, ne_eo.message ...
                        );
                    end
                end

                % ... varname (loops, fields)
                if any(strcmp(rstrtype, ...
                    {'array', 'bloop', 'eloop', 'field', 'flist', 'xloop'}))

                    % fields can be complex
                    if any(strcmp(rstrtype, {'array', 'field', 'flist'}))
                        if isempty(regexpi(rstruct.varname, ...
                            ['^([a-z][a-z_0-9]*(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?\.)*' ...
                             '[a-z][a-z_0-9]*(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?$']))
                            error( ...
                                'xff:BadTFFSpec', ...
                                'Invalid ListOfFields.varname: ''%s''.', ...
                                rstruct.varname ...
                            );
                        end

                    % loop variables MUST be simple
                    else
                        if isempty(regexpi(rstruct.varname, ...
                            '^[a-z][a-z_0-9]*(\(\d+\))?$'))
                            error ( ...
                                'xff:BadTFFSpec', ...
                                'Invalid ListOfFields.varname: ''%s''.', ...
                                rstruct.varname ...
                            );
                        end

                    end

                % ... varname (expressions, WRTLN must not be check)
                elseif ~strcmp(rstrtype, 'wrtln')
                    try
                        eval(['if 1==0,if 1==0,' ...
                            strrep(strrep(strrep(rstruct.varname, ...
                            '@@', 'cvar'), '$', 'nvar.'), '@', 'cvar.') ...
                            ',end,end']);
                    catch ne_eo;
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid EXPRE in ListOfFields.varname: ''%s'' (%s).', ...
                            rstruct.varname, ne_eo.message ...
                        );
                    end

                    % then discard other information!
                    rstruct.field = [];
                    rstruct.datatype = [];
                    rstruct.format = [];
                    rstruct.default = [];
                end

                % look at loop details
                if strcmp(rstrtype, 'bloop')

                    % varname might not be part of current looplist
                    if any(strcmp(rstruct.varname, loopused))
                        error( ...
                            'xff:BadTFFSpec', ...
                            'LOOP variable name reused: ''%s''.', ...
                            rstruct.varname ...
                        );
                    end

                    % put varname at the end of looplist
                    lvname = rstruct.varname;
                    looplist{end + 1} = lvname;
                    loopused{end + 1} = lvname;
                    slooplist = looplist;

                    % scan for loop end
                    eloopfound = false;
                    subloops = 1;
                    for sslc = (slc+1):length(actrules)

                        % get another shortcut
                        tstruct = actrules(sslc);
                        tstrtype = lower(tstruct.type);

                        % subloops (BLOOP)
                        if strcmp(tstrtype, 'bloop')
                            subloops = subloops + 1;
                            slooplist{end + 1} = tstruct.varname;

                        % end-of-loop (ELOOP)
                        elseif strcmp(tstrtype, 'eloop')
                            subloops = subloops - 1;
                            try
                                if ~strcmp(slooplist{end}, tstruct.varname)
                                    error('ILLEGALLOOP');
                                end
                                slooplist(end) = [];
                            catch ne_eo;
                                neuroelf_lasterr(ne_eo);
                                error( ...
                                    'xff:BadTFFSpec', ...
                                    'Illegal loop nesting found.' ...
                                );
                            end

                        % checking loop names for xloop
                        elseif strcmp(tstrtype, 'xloop')

                            % must be found in slooplist !
                            if ~any(strcmp(slooplist, tstruct.varname))
                                error( ...
                                    'xff:BadTFFSpec', ...
                                    'Unknown XLOOP token: ''%s''.', ...
                                    tstruct.varname ...
                                );
                            end
                        end

                        % type must be ELOOP and varname match
                        if ...
                            strcmpi(actrules(sslc).type, 'eloop') && ...
                            strcmp(actrules(sslc).varname, lvname)
                            eloopfound = true;
                            elooprule = sslc;
                            break;
                        end
                    end

                    % if no end of loop found...
                    if ~eloopfound
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Missing closing tag for LOOP %s.', ...
                            lvname ...
                        );
                    end

                    % subloops MUST be zero now
                    if subloops ~= 0
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid LOOP nesting in TFF spec.' ...
                        );
                    end

                    % build special loop struct
                    lstruct = struct( ...
                        'loopvar',   lvname, ...
                        'firstrule', slc, ...
                        'lastrule', elooprule, ...
                        'cond', rstruct.cond, ...
                        'dim', rstruct.dim);

                    % split dim from varname
                    ldim = 1;
                    fname = rstruct.varname;
                    [ldimmatcht{1:3}] = regexpi( ...
                        fname, '^[a-z][a-z_0-9]*\((\d+)\)$');
                    ldimmatcht = ldimmatcht{3};
                    if ~isempty(ldimmatcht)
                        ldim = str2double( ...
                            fname(ldimmatcht{1}(1, 1):ldimmatcht{1}(1, 1)));
                        fname = regexprep(fname, '\(.*\)', '');
                    end

                    % update known variables
                    if isfield(tffspec.Variables, fname)
                        tffspec.Variables.(fname) = max( ...
                            tffspec.Variables.(fname), ldim);
                    else
                        tffspec.Variables.(fname) = ldim;
                    end

                    % put loop struct into Loops
                    tffspec.Loops.(fname)(ldim) = lstruct;

                % loops (end)
                elseif strcmp(rstrtype, 'eloop')

                    % check if name matches last pushed name
                    if isempty(looplist) || ...
                       ~strcmp(rstruct.varname, looplist{end})
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid LOOP (end) detected in rule %d.', ...
                            slc ...
                        );
                    end

                    % pop loop from looplist
                    looplist(end) = [];

                % loops (exit)
                elseif strcmp(rstrtype, 'xloop')

                    % check if name matches any current loops
                    if isempty(looplist) || ...
                       ~any(strcmp(rstruct.varname, looplist))
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid XLOOP token: ''%s''.', ...
                            rstruct.varname ...
                        );
                    end

                    % don't do anything else !!! (see at bloop)

                % skip N rules
                elseif strcmp(rstrtype, 'skipn')

                    % dim MUST be a simple number
                    if ~isa(rstruct.dim, 'double')
                        error( ...
                            'xff:BadTFFSpec', ...
                            'SKIPN requires a simple numeric dim.' ...
                        );
                    end

                end
                rstruct.varname = ...
                    tff_parsecode(rstruct.varname, 'namevars', 'tffcont');

                % put back into actrules
                actrules(slc) = rstruct;

            end

            % put actrules at the end of ListOfFields
            tffspec.ListOfFields(end + 1:end + numel(actrules)) = actrules;

        % list of magic tokens
        case {'magic'}

            % find "EndMagic" line
            endline = llc + 1;
            for slc = (lc + 2):llc
                if ~isempty(regexpi(tfflines{slc}, '^endmagic'))
                    endline = slc - 1;
                    break;
                end
            end
            blc = lc;
            lc = endline + 1;

            % invalid Magic
            if lc > llc
                error( ...
                    'xff:BadTFFSpec', ...
                    'Unclosed Magic block in specification.' ...
                );
            end

            % list/table field separator
            listsep = tffval;

            % get rule headers
            rhead = splittocellc(tfflines{blc}, listsep, false);

            % build header struct
            hstruct = struct;
            for hfc = 1:length(rhead)
                hfield = makelabel(rhead{hfc});
                hstruct.(hfield) = hfc;
                rhead{hfc} = hfield;
            end
            blc = blc + 1;

            % bail out if invalid header given
            if ~isfield(hstruct, 'name') || ...
               ~isfield(hstruct, 'range') || ...
               ~isfield(hstruct, 'type') || ...
               ~isfield(hstruct, 'magic')
                error( ...
                    'xff:BadTFFSpec', ...
                    'Magic with bad headers.' ...
                );
            end

            % get list of header field names (in their order)
            hfields = fieldnames(hstruct);
            nfields = length(hfields);

            % build empty field list struct
            fmagic = emptystruct(hfields);

            % check global fieldlist struct
            tfields = fieldnames(tffspec.Magic);

            % no rules yet -> OK
            if isempty(tfields)
                tffspec.Magic = fmagic;

            % otherwise compare header fields
            else

                % assume no mismatch
                headermismatch = false;

                % check number of fields
                if length(tfields) ~= length(hfields)
                    tfields = {};
                    headermismatch = true;
                end

                % only if still content in tfields, check names
                for tfc = 1:length(tfields)
                    if ~strcmp(tfields{tfc}, hfields{tfc})
                        headermismatch = true;
                        break;
                    end
                end

                % if mismatch give warning and continue with next block
                if headermismatch
                    warning( ...
                        'xff:BadTFFSpec', ...
                        'Magic blocks must match in their headers.' ...
                    );
                    continue;
                end
            end

            % build list of magics to consider
            actmagic = fmagic;
            for slc = blc:endline

                % split to fields
                rcont = splittocellc(tfflines{slc}, listsep, false);

                % reject too short arrays
                if length(rcont) < length(rhead)
                    continue;
                end

                % deal into struct
                rstruct = hstruct;
                for hfc = 1:length(hfields)
                    hfield = hfields{hfc};
                    rstruct.(hfield) = deblank(rcont{hstruct.(hfield)});
                end

                % deal with empty lines
                if isempty(rstruct.name) || ...
                    isempty(rstruct.range) || ...
                    isempty(rstruct.type) || ...
                    isempty(rstruct.magic)
                    continue;
                end

                % put non-empty type/varname rules into actrules
                actmagic(end + 1) = rstruct;

            end

            % parse magics (using for; if any magic fails, error out!)
            for slc = 1:length(actmagic)

                % get rstruct from actrules
                rstruct = actmagic(slc);
                rstrtype = lower(rstruct.type);

                % check syntax of fields, name
                if ~strcmp(rstruct.name, makelabel(rstruct.name))
                    error( ...
                        'xff:BadTFFSpec', ...
                        'Invalid Magic.name token: ''%s''.', ...
                        rstruct.name ...
                    );
                end

                % ..., type
                if ~any(strcmp(rstrtype, ...
                    {'hex', 'regexp', 'regexpi', 'strfind'}))
                    error( ...
                        'xff:BadTFFSpec', ...
                        'Invalid Magic.type token: ''%s''.', ...
                        rstruct.type ...
                    );
                end

                % ..., range
                try
                    rrange = eval(['[' rstruct.range ']']);
                    if numel(rrange) ~= 2
                        error('INVALIDRANGE');
                    end
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    error( ...
                        'xff:BadTFFSpec', ...
                        'Invalid Magic.range specification: ''%s''.', ...
                        rstruct.range ...
                    );
                end
                actmagic(slc).range = rrange;

                % ..., magic -> parse hex codes
                if strcmp(rstrtype, 'hex')

                    % split at comma, semicolon or spaces
                    rhexcodes = splittocellc(rstruct.magic, ',; ', true, true);

                    % and convert
                    rhexvals  = zeros(1, numel(rhexcodes));
                    for rhc = 1:numel(rhexcodes)
                        rhexvals(rhc) = hex2dec(rhexcodes{rhc});
                    end

                    % put back into array
                    actmagic(slc).magic = rhexvals(:)';

                % allow hexadecimal content in other fields
                else

                    % pack %C{xx} and 0x{xx} into hex vars
                    [hexpackt{1:3}] = regexpi(rstruct.magic,  ...
                        '(\%c\{|0x\{)([0-9a-f][0-9a-f])(\})', 'once');
                    hexpackt = hexpackt{3};
                    while ~isempty(hexpackt)
                        rstruct.magic = strrep(rstruct.magic, ...
                            rstruct.magic(hexpackt(1, 1):hexpackt(3, 2)), ...
                            char(hex2dec( ...
                            rstruct.magic(hexpackt(2, 1):hexpackt(2, 2)))));
                        [hexpackt{1:3}] = regexpi(rstruct.magic, ...
                            '(\%c\{|0x\{)([0-9a-f][0-9a-f])(\})', 'once');
                    end

                    % put back into actmagic
                    actmagic(slc).magic = rstruct.magic;

                end

            end

            % put actrules at the end of ListOfFields
            tffspec.Magic(end + 1:end + length(actmagic)) = actmagic;

        % new file code snippet
        case {'newfilecode'}

            % find "EndNewFileCode"
            endline = llc + 1;
            for slc = lc:llc
                if ~isempty(regexpi(tfflines{slc}, '^endnewfilecode'))
                    endline = slc - 1;
                    break;
                end
            end

            % check endline
            if endline > llc
                error( ...
                    'xff:BadTFFSpec', ...
                    'Unclosed NewFileCode section in TFF spec.' ...
                );
            end

            % generate and check code snippet
            codesnip  = gluetostringc(tfflines(lc:endline), char(10), true);
            tcodesnip = codesnip;
            tcodesnip(tcodesnip == '@') = [];
            try
                eval(['if 1==0,' tcodesnip ';end']);
            catch ne_eo;
                error( ...
                    'xff:BadTFFSpec', ...
                    'Syntax error detected in NewFileCode: ''%s''.', ...
                    ne_eo.message ...
                );
            end

            % store code
            tffspec.NewFileCode = ...
                tff_parsecode(codesnip, 'namevars', 'tffcont');

        % paragraph arrays
        case {'paragrapharrays'}
            try
                tffspec.ParagraphArrays = logical(eval(tffval));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xff:BadTFFSpec', ...
                    'Invalid ParagraphArrays value given: ''%s''.', ...
                    tffval ...
                );
            end

        % skip empty lines
        case {'skipemptylines'}
            try
                tffspec.SkipEmptyLines = logical(eval(tffval));
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'xff:BadTFFSpec', ...
                    'Invalid SkipEmptyLines value given: ''%s''.', ...
                    tffval ...
                );
            end

        % valid file code snippet
        case {'validfilecode'}

            % find "EndValidFileCode"
            endline = llc + 1;
            for slc = lc:llc
                if ~isempty(regexpi(tfflines{slc}, '^endvalidfilecode'))
                    endline = slc - 1;
                    break;
                end
            end

            % check endline
            if endline > llc
                error( ...
                    'xff:BadTFFSpec', ...
                    'Unclosed ValidFileCode section in TFF spec.' ...
                );
            end

            % generate and check code snippet
            codesnip  = gluetostringc(tfflines(lc:endline), char(10), true);
            tcodesnip = codesnip;
            tcodesnip(tcodesnip == '@') = [];
            try
                eval(['if 1==0,' tcodesnip ';end']);
            catch ne_eo;
                error( ...
                    'xff:BadTFFSpec', ...
                    'Syntax error detected in ValidFileCode: ''%s''.', ...
                    ne_eo.message ...
                );
            end

            % store code
            tffspec.NewFileCode = codesnip;

        % unrecognized token, give warning
        otherwise
            warning( ...
                'xff:InvalidToken', ...
                'Invalid token in TFF file: %s.', ...
                tffkey ...
            );
    end
end

% check longest field name
fll = 8; % minimum field name length in ANY case
for flc = 1:length(tffspec.ListOfFields)

    % only for FIELD or FLIST
    if any(strcmpi(tffspec.ListOfFields(flc).type, ...
           {'field', 'flist'}))
        fll = max(fll, length(tffspec.ListOfFields(flc).field));
    end
end
tffspec.MaxFieldNameLength = fll;


function parsedcode = tff_parsecode(unparsed, nv, cnt)

    % any occurance of $$/@@
    parsedcode = strrep(strrep(unparsed, '@@', cnt), '$$', nv);

    % regexprep
    parsedcode = regexprep(parsedcode, '\@([a-zA-Z])', [cnt '.$1']);
    parsedcode = regexprep(parsedcode, '\$([a-zA-Z])', [nv '.$1']);

% end of function parsedcode = tff_parsecode(unparsed, nv, cnt)

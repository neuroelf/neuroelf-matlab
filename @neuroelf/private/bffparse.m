function bffspec = bffparse(bffcont)
% bffparse  - parse BFF binary file format description
%
% FORMAT:       bffspec = bffparse(bffcont)
%
% Input fields:
%
%       bffcont     1xN char: BFF specification file name or content
%
% Output fields:
%
%       bffspec     1x1 struct with BFF specification
%
% See also bffdocu, bffio

% Version:  v1.0
% Build:    16011416
% Date:     Jan-14 2016, 4:51 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2013, 2014, 2016, Jochen Weber
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

bffversion = 'v0.9c';
default_tiosize = Inf;

% create persistent datatype sizes struct if needed
persistent bff_dtsize;
if isempty(bff_dtsize)
    bff_dtsize = struct( ...
        'double', 8, ...
        'int16',  2, ...
        'int32',  4, ...
        'int64',  8, ...
        'int8',   1, ...
        'single', 4, ...
        'uchar',  1, ...
        'uint16', 2, ...
        'uint32', 4, ...
        'uint64', 8, ...
        'uint8',  1  ...
    );
end

% argument check
if nargin < 1 || ...
   ~ischar(bffcont) || ...
    isempty(bffcont)
    error( ...
        'xff:BadArgument', ...
        'Invalid BFF specification argument given.' ...
    );
end

% check for line feeds and try to read file
if ~any(bffcont == char(10))
    if exist(bffcont, 'file') ~= 2
        error( ...
            'xff:FileNotFound', ...
            'Either invalid content or file not found.' ...
        );
    end

    % read bffspec
    try
        bffcont = asciiread(bffcont);
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% check content
if isempty(strfind(bffcont, 'BinaryFileFormat'))
    error( ...
        'xff:BadArgument', ...
        'Invalid BFF specification file/content given.' ...
    );
end

% split into lines
bfflines = splittocellc(bffcont(:)', char([10, 13]), true, true);

% remove comments, deblank, and then discard empty lines
bfflines = regexprep(bfflines, '\#.*$', '');
bfflines = deblank(bfflines);
bfflines(cellfun('isempty', bfflines)) = [];

% initialize content
bffspec = struct;
bffspec.FFTYPE          = 'BFF';
bffspec.BFFVERSION      = bffversion;
bffspec.AfterReadCode   = '';
bffspec.BeforeWriteCode = '';
bffspec.DefaultProperty = {};
bffspec.Description     = 'All files';
bffspec.EncodingSyntax  = 'native';
bffspec.Extensions      = {};
bffspec.FilenameMatch   = {};
bffspec.Filetype        = 'Custom BFF file';
bffspec.ListOfFields    = emptystruct({});
bffspec.Loops           = struct;
bffspec.Magic           = emptystruct({});
bffspec.NewFileCode     = '';
bffspec.TransIOSize     = default_tiosize;
bffspec.ValidFileCode   = '';
bffspec.Variables       = struct;

% parse content
lc  = 1;
llc = length(bfflines);
while (lc <= llc)

    % regexpi line matches '<TOKEN>:<VALUE>' ?
    bffline = bfflines{lc};
    [tlmatcht{1:3}] = regexpi( ...
        bffline, '^([a-z][a-z_0-9]*)\:\s*(.*)\s*$');
    tlmatcht = tlmatcht{3};

    % continue with next line if no match
    lc = lc + 1;
    if isempty(tlmatcht)
        continue;
    end

    % handle those sections/keywords:
    bffkey = bffline(tlmatcht{1}(1, 1):tlmatcht{1}(1, 2));
    bffval = bffline(tlmatcht{1}(2, 1):tlmatcht{1}(2, 2));
    switch (lower(bffkey))

        % after read code snippet
        case {'afterreadcode'}

            % find "EndAfterReadCode"
            endline = llc + 1;
            for slc = lc:llc
                if ~isempty(regexpi(bfflines{slc}, '^endafterreadcode'))
                    endline = slc - 1;
                    break;
                end
            end

            % check endline
            if endline > llc
                error( ...
                    'xff:BadBFFSpec', ...
                    'Unclosed AfterReadCode section in BFF spec.' ...
                );
            end

            % generate and check code snippet
            codesnip  = gluetostringc(bfflines(lc:endline), char(10), true);
            tcodesnip = codesnip;
            tcodesnip(tcodesnip == '$' | tcodesnip == '@') = [];
            try
                eval(['if 1==0,' tcodesnip ';end']);
            catch ne_eo;
                error( ...
                    'xff:BadBFFSpec', ...
                    'Syntax error detected in AfterReadCode: ''%s''.', ...
                    ne_eo.message ...
                );
            end

            % store code
            bffspec.AfterReadCode = ...
                bff_parsecode(codesnip, 'namevars', 'bffcont');

        % before write code snippet
        case {'beforewritecode'}

            % find "EndBeforeWriteCode"
            endline = llc + 1;
            for slc = lc:llc
                if ~isempty(regexpi(bfflines{slc}, '^endbeforewritecode'))
                    endline = slc - 1;
                    break;
                end
            end

            % check endline
            if endline > llc
                error( ...
                    'xff:BadBFFSpec', ...
                    'Unclosed BeforeWriteCode section in BFF spec.' ...
                );
            end

            % generate and check code snippet
            codesnip  = gluetostringc(bfflines(lc:endline), char(10), true);
            tcodesnip = codesnip;
            tcodesnip(tcodesnip == '$' | tcodesnip == '@') = [];
            try
                eval(['if 1==0,' tcodesnip ';end']);
            catch ne_eo;
                error( ...
                    'xff:BadBFFSpec', ...
                    'Syntax error detected in BeforeWriteCode: ''%s''.', ...
                    ne_eo.message ...
                );
            end

            % store code
            bffspec.BeforeWriteCode = ...
                bff_parsecode(codesnip, 'namevars', 'bffcont');

        % default property (for access)
        case {'defaultproperty'}
            if ~isempty(bffval) && ...
                bffval(1) == '{' && ...
                bffval(end) == '}'
                try
                    bffval = eval(bffval);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    error( ...
                        'xff:BadBFFSpec', ...
                        'Invalid DefaultProperty value given: ''%s''.', ...
                        bffval ...
                    );
                end
            end
            if ~iscell(bffval)
                bffval = {bffval};
            end
            bffspec.DefaultProperty = bffval;

        % file open description
        case {'description'}
            bffspec.Description = splittocellc(bffval, ';,', true, true);

        % encoding syntax
        case {'encodingsyntax'}
            bffspec.EncodingSyntax = bffval;

        % valid extensions list
        case {'extensions'}
            bffspec.Extensions = splittocellc(bffval, ';,. ', true, true);

        % valid extensions list
        case {'filenamematch'}
            bffspec.FilenameMatch = splittocellc(bffval, ';, ', true, true);

        % list of fields
        case {'listoffields'}

            % find "EndListOfFields" line
            endline = llc + 1;
            for slc = (lc + 1):llc
                if ~isempty(regexpi(bfflines{slc}, '^endlistoffields'))
                    endline = slc - 1;
                    break;
                end
            end
            blc = lc;
            lc = endline + 1;

            % invalid ListOfFields
            if lc > llc
                error( ...
                    'xff:BadBFFSpec', ...
                    'Unclosed or bad ListOfFields block.' ...
                );
            end

            % list/table field separator
            listsep = bffval;

            % get rule headers
            rhead = splittocellc(bfflines{blc}, listsep, false);

            % build header struct
            hstruct = struct;
            for hfc = 1:length(rhead)
                hfield = lower(makelabel(rhead{hfc}));
                hstruct.(hfield) = hfc;
                rhead{hfc} = hfield;
            end
            blc = blc + 1;

            % bail out if invalid header given
            if ~isfield(hstruct, 'type') || ...
               ~isfield(hstruct, 'cond') || ...
               ~isfield(hstruct, 'disktype') || ...
               ~isfield(hstruct, 'datatype') || ...
               ~isfield(hstruct, 'dim') || ...
               ~isfield(hstruct, 'default') || ...
               ~isfield(hstruct, 'varname')
                warning( ...
                    'xff:BadBFFSpec', ...
                    'ListOfFields with bad headers.' ...
                );
                continue;
            end

            % get list of header field names (in their order)
            hfields = fieldnames(hstruct);
            nfields = numel(hfields);

            % build empty field list struct
            frules = cell2struct(cell(1, 1, nfields), hfields, 3);
            frules.disksize = 0;
            frules(:) = [];

            % check global fieldlist struct
            tfields = fieldnames(bffspec.ListOfFields);

            % no rules yet -> OK
            if isempty(tfields)
                bffspec.ListOfFields = frules;

            % otherwise compare header fields
            else

                % check number of fields and names
                if numel(tfields) ~= numel(fieldnames(frules)) || ...
                   ~all(strcmp(tfields, hfields))

                    % if mismatch give warning and continue with next block
                    warning( ...
                        'xff:BadBFFSpec', ...
                        'ListOfFields blocks must match in their headers.' ...
                    );
                    continue;
                end
            end

            % build list of rules to consider
            actrules = frules;
            slc = blc;
            while slc <= endline

                % split to fields
                rcont = splittocellc(bfflines{slc}, listsep, false);

                % increase counter
                slc = slc + 1;

                % continuation?
                while ~isempty(rcont) && ...
                    numel(rcont{end}) > 2 && ...
                    strcmp(rcont{end}(end-2:end), '...') && ...
                    slc <= endline

                    % split the next line also
                    rcontnext = splittocellc(bfflines{slc}, listsep, false);

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
                            'xff:BFFSpecError', ...
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
                elseif ~strcmp(rstrtype, 'skipn') && ...
                    isempty(rstruct.varname)
                    warning( ...
                        'xff:AmbiguousBFFSpec', ...
                        'Empty variable only allowed for SKIPN toked.' ...
                    );
                    continue;
                end

                % put non-empty type/varname rules into actrules
                rstruct.disksize = 0;
                actrules(end + 1) = rstruct;
            end

            % init loop checking variables
            looplist = {};
            loopused = {};

            % parse rules (using for; if any rules fail, error out!)
            for slc = 1:numel(actrules)

                % get rstruct from actrules
                rstruct = actrules(slc);

                % check syntax of fields, type
                rstrtype = lower(rstruct.type);
                if ~any(strcmp(rstrtype, ...
                    {'bloop', 'eloop', 'expre', 'field', 'skipn', 'xloop'}))
                    error( ...
                        'xff:BadBFFSpec', ...
                        'Invalid ListOfFields.type token: ''%s''.', ...
                        rstruct.type ...
                    );
                end

                % ... cond
                if ~isempty(rstruct.cond)
                    try
                        eval([ ...
                            'if 1==0,if ' ...
                            mstrrep(rstruct.cond, ...
                                {'(\@\@|\$\$)', '(\@|\$)'}, ...
                                {'tvar', 'tvar.'}, 1) ...
                            ',disp(1);end,end']);
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid ListOfFields.cond content: ''%s''.', ...
                            rstruct.cond ...
                        );
                    end
                    rstruct.cond = ...
                        bff_parsecode(rstruct.cond, 'namevars', 'bffcont');
                end

                % ... disktype
                if strcmp(rstrtype, 'field')

                    % get shortcut
                    rtype = rstruct.disktype;

                    % check type syntax
                    if isempty(regexpi(rtype, '^[a-z][a-z_0-9]+$'))
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid ListOfFields.disktype tag: ''%s''.', ...
                            rtype ...
                        );
                    end

                    % check disktype
                    switch lower(rtype)

                        % disktype is built in
                        case { ...
                                 'char',  'int8',  'int16',  'int32', ...
                                'uchar', 'uint8', 'uint16', 'uint32', ...
                                'int64', 'uint64','single', 'double', ...
                              'cstring', 'clinea'}
                            rtype = lower(rtype);

                        % disktype string, cstring, char*
                        case {'char*', 'string'}
                            rtype = 'cstring';

                        % unsupported disktype
                        otherwise
                            error( ...
                                'xff:BadBFFSpec', ...
                                'Invalid disktype in FIELD(%d).disktype', ...
                                slc ...
                            );
                    end

                    % put back into actrules(slc)
                    rstruct.disktype = rtype;
                    if isfield(bff_dtsize, lower(rtype))
                        rstruct.disksize = bff_dtsize.(lower(rtype));
                    else
                        rstruct.disksize = 0;
                    end
                end

                % ... datatype
                if ...
                    strcmp(rstrtype, 'field') && ...
                    isempty(regexpi(rstruct.datatype, '^[a-z][a-z_0-9]+$'))
                    error( ...
                        'xff:BadBFFSpec', ...
                        'Invalid ListOfFields.datatype tag: ''%s''.', ...
                        rstruct.datatype ...
                    );
                end

                % ... dim (loops)
                if strcmp(rstrtype, 'bloop')

                    % only single number OR variable
                    if isempty(regexpi(rstruct.dim, ...
                        '^(\d+|[\$\@][a-z][a-z_0-9\.]*(\((\d+|[\$\@][a-z][a-z_0-9\.]*)\))?)$'))
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid LOOP.DIM given: ''%s''.', ...
                            rstruct.dim ...
                        );
                    end

                % ... dim (fields)
                elseif strcmp(rstrtype, 'field')

                    % multiple numbers AND/OR variables
                    if isempty(regexpi(rstruct.dim, ...
                        ['^((\d+|[\$\@][a-z][a-z_0-9]*)(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?' ...
                         '(\.[a-z][a-z_0-9]*(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?)*\,\s*)*' ...
                         '((\d+|[\$\@][a-z][a-z_0-9]*)(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?)$']))
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid ListOfFields.dim given: ''%s''.', ...
                            rstruct.dim ...
                        );
                    end

                % ... dim (skipX)
                elseif strcmp(rstrtype(1:4), 'skip')

                    % only a 1-D numeric value accepted
                    if ~all(rstruct.dim > 47 & rstruct.dim < 58)
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid SKIPx.dim given: ''%s''.', ...
                            rstruct.dim ...
                        );
                    end
                end
                rstruct.dim = ...
                    bff_parsecode(rstruct.dim, 'namevars', 'bffcont');
                if all(rstruct.dim == ',' | rstruct.dim == ' ' | ...
                      (rstruct.dim >= '0' & rstruct.dim <= '9'))
                    try
                        rstruct.dim = eval(['[' rstruct.dim ']']);
                    catch ne_eo;
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid numeric dimension given: %s.', ...
                            ne_eo.message ...
                        );
                    end
                    if numel(rstruct.dim) == 1
                        rstruct.dim = [1, rstruct.dim];
                    end
                end

                % ... default (if non-empty)
                if ~isempty(rstruct.default)
                    try
                        eval(['if 1==0,' ...
                              mstrrep(rstruct.default, ...
                                  {'(\@\@|\$\$)', '(\@|\$)'}, ...
                                  {'tvar', 'tvar.'}, 1) ...
                              ',end']);
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid ListOfFields.default value: ''%s''.', ...
                            rstruct.default ...
                        );
                    end
                    rstruct.default = ...
                        bff_parsecode(rstruct.default, 'namevars', 'bffcont');
                end

                % ... varname (loops, fields)
                if any(strcmp(rstrtype, ...
                    {'bloop', 'eloop', 'field', 'xloop'}))

                    % fields can be complex
                    if strcmp(rstrtype, 'field')
                        if isempty(regexpi(rstruct.varname, ...
                            ['^([a-z][a-z_0-9]*(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?\.)*' ...
                             '[a-z][a-z_0-9]*(\((\d+|[\$\@][a-z][a-z_0-9]*)\))?$']))
                            error( ...
                                'xff:BadBFFSpec', ...
                                'Invalid ListOfFields.varname: ''%s''.', ...
                                rstruct.varname ...
                            );
                        end

                    % loop variables MUST be simple
                    else
                        if isempty(regexpi(rstruct.varname, ...
                            '^[a-z][a-z_0-9]*(\(\d+\))?$'))
                            error ( ...
                                'xff:BadBFFSpec', ...
                                'Invalid ListOfFields.varname: ''%s''.', ...
                                rstruct.varname ...
                            );
                        end
                    end

                % ... varname (expressions)
                else
                    try
                        eval(['if 1==0,' ...
                            mstrrep(rstruct.varname, ...
                                {'(\@\@|\$\$)', '(\@|\$)'}, ...
                                {'tvar', 'tvar.'}, 1) ...
                            ',end']);
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid EXPRE in ListOfFields.varname: ''%s'' / ''%s''.', ...
                            rstruct.varname, ne_eo.message ...
                        );
                    end
                end

                % look at loop details
                if strcmp(rstrtype, 'bloop')

                    % varname might not be part of current looplist
                    if any(strcmp(rstruct.varname, loopused))
                        error( ...
                            'xff:BadBFFSpec', ...
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
                    for sslc = (slc+1):numel(actrules)

                        % get another shortcut
                        tstruct = actrules(sslc);

                        % subloops (BLOOP)
                        tstrtype = lower(tstruct.type);
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
                                    'xff:BadBFFSpec', ...
                                    'Illegal loop nesting found.' ...
                                );
                            end

                        % checking loop names for xloop
                        elseif strcmp(tstrtype, 'xloop')

                            % must be found in slooplist !
                            if ~any(strcmp(slooplist, tstruct.varname))
                                error( ...
                                    'xff:BadBFFSpec', ...
                                    'Unknown XLOOP token: ''%s''.', ...
                                    tstruct.varname ...
                                );
                            end
                        end

                        % type must be ELOOP and varname match
                        if numel(actrules(sslc).type) == 5 && ...
                            all(lower(actrules(sslc).type) == 'eloop') && ...
                            strcmp(actrules(sslc).varname, lvname)
                            eloopfound = true;
                            elooprule = sslc;
                            break;
                        end
                    end

                    % if no end of loop found...
                    if ~eloopfound
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Missing closing tag for LOOP %s.', ...
                            lvname ...
                        );
                    end

                    % subloops MUST be zero now
                    if subloops ~= 0
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid LOOP nesting in BFF spec.' ...
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
                    if ~isempty(ldimmatcht) && ...
                       ~isempty(ldimmatcht{1})
                        ldim = str2double( ...
                            fname(ldimmatcht{1}(1, 1):ldimmatcht{1}(1, 1)));
                        fname = regexprep(fname, '\(.*\)', '');
                    end

                    % update known variables
                    if isfield(bffspec.Variables, fname)
                        bffspec.Variables.(fname) = max( ...
                            bffspec.Variables.(fname), ldim);
                    else
                        bffspec.Variables.(fname) = ldim;
                    end

                    % put loop struct into Loops
                    bffspec.Loops.(fname)(ldim) = lstruct;

                % loops (end)
                elseif strcmp(rstrtype, 'eloop')

                    % check if name matches last pushed name
                    if isempty(looplist) || ...
                       ~strcmp(rstruct.varname, looplist{end})
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid LOOP (end) detected in rule %d.', ...
                            slc ...
                        );
                    end

                    % pop loop from looplist
                    looplist(end) = [];

                % loops (exit)
                elseif strcmp(rstrtype, 'xloop')

                    % check if name matches any current loops
                    if ...
                        isempty(looplist) || ...
                       ~any(strcmp(rstruct.varname, looplist))
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid XLOOP token: ''%s''.', ...
                            rstruct.varname ...
                        );
                    end

                    % don't do anything else !!! (see at bloop)

                % skip N rules
                elseif strcmp(rstrtype, 'skipn')

                    % dim MUST be a simple number
                    if ~all(rstruct.dim > 47 & rstruct.dim < 58)
                        error( ...
                            'xff:BadBFFSpec', ...
                            'SKIPN requires a simple numeric dim.' ...
                        );
                    end

                end
                rstruct.varname = ...
                    bff_parsecode(rstruct.varname, 'namevars', 'bffcont');

                % put back into actrules
                actrules(slc) = rstruct;

            end

            % put actrules at the end of ListOfFields
            bffspec.ListOfFields(end + 1:end + length(actrules)) = actrules;

        % list of magic tokens
        case {'magic'}

            % find "EndMagic" line
            endline = llc + 1;
            for slc = (lc + 2):llc
                if ~isempty(regexpi(bfflines{slc}, '^endmagic'))
                    endline = slc - 1;
                    break;
                end
            end
            blc = lc;
            lc = endline + 1;

            % invalid Magic
            if lc > llc
                error( ...
                    'xff:BadBFFSpec', ...
                    'Unclosed Magic block in specification.' ...
                );
            end

            % list/table field separator
            listsep = bffval;

            % get rule headers
            rhead = splittocellc(bfflines{blc}, listsep, false);

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
                    'xff:BadBFFSpec', ...
                    'Magic with bad headers.' ...
                );
            end

            % get list of header field names (in their order)
            hfields = fieldnames(hstruct);
            nfields = length(hfields);

            % build empty field list struct
            fmagic = emptystruct(hfields);

            % check global fieldlist struct
            tfields = fieldnames(bffspec.Magic);

            % no rules yet -> OK
            if isempty(tfields)
                bffspec.Magic = fmagic;

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
                        'xff:BadBFFSpec', ...
                        'Magic blocks must match in their headers.' ...
                    );
                    continue;
                end
            end

            % build list of magics to consider
            actmagic = fmagic;
            for slc = blc:endline

                % split to fields
                rcont = splittocellc(bfflines{slc}, listsep, false);

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

                % check syntax of fields, name
                if ~strcmp(rstruct.name, makelabel(rstruct.name))
                    error( ...
                        'xff:BadBFFSpec', ...
                        'Invalid Magic.name token: ''%s''.', ...
                        rstruct.name ...
                    );
                end

                % ..., type
                rstrtype = lower(rstruct.type);
                if ~any(strcmp(rstrtype, ...
                    {'hex', 'regexp', 'regexpi', 'strfind'}))
                    error( ...
                        'xff:BadBFFSpec', ...
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
                        'xff:BadBFFSpec', ...
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
                    rhexvals  = [];
                    for rhc = 1:length(rhexcodes)
                        rhexvals(end+1) = hex2dec(rhexcodes{rhc});
                    end

                    % put back into array
                    actmagic(slc).magic = rhexvals(:)';

                % allow hexadecimal content in other fields
                else

                    % pack %C{xx} and 0x{xx} into hex vars
                    [hexpackt{1:3}] = regexpi( ...
                        rstruct.magic,  ...
                        '(\%c\{|0x\{)([0-9a-f][0-9a-f])(\})', 'once');
                    while ~isempty(hexpackt{3})
                        hexpackt = hexpackt{3};
                        rstruct.magic = strrep(rstruct.magic, ...
                            rstruct.magic(hexpackt(1, 1):hexpackt(3, 2)), ...
                            char(hex2dec( ...
                            rstruct.magic(hexpackt(2, 1):hexpackt(2, 2)))));
                        [hexpackt{1:3}] = regexpi( ...
                            rstruct.magic, ...
                            '(\%c\{|0x\{)([0-9a-f][0-9a-f])(\})', 'once');
                    end

                    % put back into actmagic
                    actmagic(slc).magic = rstruct.magic;

                end

            end

            % put actrules at the end of ListOfFields
            bffspec.Magic(end + 1:end + length(actmagic)) = actmagic;

        % create new file snippet
        case {'newfilecode'}

            % find "EndNewFileCode"
            endline = llc + 1;
            for slc = lc:llc
                if ~isempty(regexpi(bfflines{slc}, '^endnewfilecode'))
                    endline = slc - 1;
                    break;
                end
            end

            % check endline
            if endline > llc
                error( ...
                    'xff:BadBFFSpec', ...
                    'Unclosed NewFileCode section in BFF spec.' ...
                );
            end

            % generate and check code snippet
            codesnip = gluetostringc(bfflines(lc:endline), char(10), true);
            codesnip = regexprep(codesnip, '\.{3,}\s+', ' ');
            tcodesnip = codesnip;
            tcodesnip(tcodesnip == '@') = [];
            try
                eval(['if 1==0,' tcodesnip ';end']);
            catch ne_eo;
                error( ...
                    'xff:BadBFFSpec', ...
                    'Syntax error detected in NewFileCode: ''%s''.', ...
                    ne_eo.message ...
                );
            end

            % store code
            bffspec.NewFileCode = ...
                bff_parsecode(codesnip, 'namevars', 'bffcont');

        % critical transio object size
        case {'transiosize'}
            try
                tiosize = splittocellc(bffval, ';, ', true, true);
                if ~isempty(tiosize)
                    bffspec.TransIOSize = str2double(tiosize{1});
                end
                if length(bffspec.TransIOSize) ~= 1 || ...
                    isnan(bffspec.TransIOSize) || ...
                    bffspec.TransIOSize < 1
                    bffspec.TransIOSize = default_tiosize;
                else
                    bffspec.TransIOSize = fix(bffspec.TransIOSize);
                end
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                bffspec.TransIOSize = default_tiosize;
            end

        % create valid file snippet
        case {'validfilecode'}

            % find "EndValidFileCode"
            endline = llc + 1;
            for slc = lc:llc
                if ~isempty(regexpi(bfflines{slc}, '^endvalidfilecode'))
                    endline = slc - 1;
                    break;
                end
            end

            % check endline
            if endline > llc
                error( ...
                    'xff:BadBFFSpec', ...
                    'Unclosed ValidFileCode section in BFF spec.' ...
                );
            end

            % generate and check code snippet
            codesnip = gluetostringc(bfflines(lc:endline), char(10), true);
            codesnip = regexprep(codesnip, '\.{3,}\s+', ' ');
            tcodesnip = codesnip;
            tcodesnip(tcodesnip == '@') = [];
            try
                eval(['if 1==0,' tcodesnip ';end']);
            catch ne_eo;
                error( ...
                    'xff:BadBFFSpec', ...
                    'Syntax error detected in ValidFileCode: ''%s''.', ...
                    ne_eo.message ...
                );
            end

            % store code (unparsed, must be done later!)
            bffspec.ValidFileCode = codesnip;

        % unrecognized token, give a warning
        otherwise
            warning( ...
                'xff:InvalidToken', ...
                'Invalid token in BFF file: %s.', ...
                bffkey ...
            );
    end
end

function parsedcode = bff_parsecode(unparsed, nv, cnt)

    % any occurance of $$/@@
    parsedcode = strrep(strrep(unparsed, '@@', cnt), '$$', nv);

    % regexprep
    parsedcode = regexprep(parsedcode, '\@([a-zA-Z])', [cnt '.$1']);
    parsedcode = regexprep(parsedcode, '\$([a-zA-Z])', [nv '.$1']);

% end of function parsedcode = bff_parsecode(unparsed, nv, cnt)

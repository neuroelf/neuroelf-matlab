function [varargout] = tffio(varargin)
% tffio  - read/write binary file with TFF spec (and content)
%
% FORMAT:       tffcont = tffio(filename, tffspec [, options])
%          OR   tffio(wfilename, tffspec, tffcont [, options])
%
% Input fields:
%
%       filename    filename of binary file to read
%       tffspec     either filename or content or spec struct of TFF
%       options     optional arguments, one of
%                    'force' , 'f', or '-f' (for returning partial content)
%                   'header' , 'h', or '-h' (for reading the header only)
%                   'verbose', 'v', or '-v' (for verbose reading/writing)
%
%       wfilename   filename of binary file to write
%       tffcont     file content (struct)
%
% Output fields:
%
%       tffcont     file contents struct (depending on TFF)
%
% See also tffdocu, tffparse

% Version:  v1.1
% Build:    16012320
% Date:     Jan-23 2016, 8:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% TFF version number
tffversion = 'v1.0';

% argument check
if nargin < 2 || ...
    ~ischar(varargin{1}) || ...
    isempty(varargin{1}) || ...
   ((~isstruct(varargin{2}) || ...
     numel(varargin{2}) ~= 1) && ...
    ~ischar(varargin{2})) || ...
    isempty(varargin{2})
    error( ...
        'xff:BadArgument', ...
        'Bad or missing argument for tffio.' ...
    );
end
filename = varargin{1}(:)';
if ischar(varargin{2})
    tffspec  = varargin{2}(:)';
else
    tffspec  = varargin{2};
end
tffcont  = [];

% default options
forcemode   = false;
headeronly  = false;
verbosemode = false;
writemode   = false;

% tff content structure given
if nargin > 2 && ...
    isstruct(varargin{3}) && ...
   ~isempty(varargin{3})
    tffcont = varargin{3};
    writemode = true;
end

% check additional arguments
if nargin > 2
    for argc = 3:nargin

        % what is given
        argv = varargin{argc};

        % character argument
        if ischar(argv)

            % known argument
            switch lower(argv)

                % force mode
                case {'f', 'force'}
                    forcemode = true;

                % force mode
                case {'h', 'header'}
                    headeronly = true;

                % verbose mode
                case {'v', 'verbose'}
                    verbosemode = true;
            end

            % check GNU-style option
            if numel(argv) > 1 && ...
                argv(1) == '-'

                % force mode
                if any(argv == 'f')
                    forcemode = true;
                end

                % header only
                if any(argv == 'h')
                    headeronly = true;
                end

                % verbose mode
                if any(argv == 'v')
                    verbosemode = true;
                end
            end

        % discard other arguments/options
        end
    end
end

% TFF is a valid SPEC struct
if isstruct(tffspec) && ...
   ~isempty(tffspec) && ...
    isfield(tffspec, 'ArrayFormat') && ...
    ischar(tffspec.ArrayFormat) && ...
    any(tffspec.ArrayFormat(:) == '%') && ...
    isfield(tffspec, 'BinaryIO') && ...
   ~isempty(tffspec.BinaryIO) && ...
    islogical(tffspec.BinaryIO) && ...
    isfield(tffspec, 'CheckForNaNs') && ...
   ~isempty(tffspec.CheckForNaNs) && ...
    islogical(tffspec.CheckForNaNs) && ...
    isfield(tffspec, 'CustomDelimiters') && ...
    iscell(tffspec.CustomDelimiters) && ...
    isfield(tffspec, 'Extensions') && ...
    iscell(tffspec.Extensions) && ...
    isfield(tffspec, 'FieldDelimCollapse') && ...
   ~isempty(tffspec.FieldDelimCollapse) && ...
    islogical(tffspec.FieldDelimCollapse) && ...
    isfield(tffspec, 'FieldDelimiters') && ...
    iscell(tffspec.FieldDelimiters) && ...
   ~isempty(tffspec.FieldDelimiters) && ...
    isfield(tffspec, 'LineDelimiters') && ...
    iscell(tffspec.LineDelimiters) && ...
   ~isempty(tffspec.LineDelimiters) && ...
    isfield(tffspec, 'ListOfFields') && ...
    isstruct(tffspec.ListOfFields) && ...
    isfield(tffspec, 'Loops') && ...
    isstruct(tffspec.Loops) && ...
    isfield(tffspec, 'Magic') && ...
    isstruct(tffspec.Magic) && ...
    isfield(tffspec, 'MaxFieldNameLength') && ...
    isnumeric(tffspec.MaxFieldNameLength) && ...
   ~isempty(tffspec.MaxFieldNameLength) && ...
    isfield(tffspec, 'ParagraphArrays') && ...
   ~isempty(tffspec.ParagraphArrays) && ...
    islogical(tffspec.ParagraphArrays) && ...
    isfield(tffspec, 'SkipEmptyLines') && ...
   ~isempty(tffspec.SkipEmptyLines) && ...
    islogical(tffspec.SkipEmptyLines) && ...
    isfield(tffspec, 'Variables') && ...
    isstruct(tffspec.Variables)

% any other struct is rejected
elseif isstruct(tffspec)
    error( ...
        'xff:BadArgument', ...
        'Invalid struct TFF specification.' ...
    );

% for character arrays
elseif ischar(tffspec)

    % valid file contents (0x0a/0x0d delimited)
    if any(tffspec == char(10) | tffspec == char(13)) || ...
        exist(tffspec, 'file') == 2

        % parse content
        try
            if writemode
                tffcont = tffio( ...
                    filename, ...
                    tffparse(tffspec), ...
                    varargin{:});
            else
                tffcont = tffio( ...
                    filename, ...
                    tffparse(tffspec), ...
                    varargin{:});
            end
            if ~writemode
                varargout{1} = tffcont;
            else
                if nargout > 0
                    varargout = cell(1, nargout);
                end
            end
            return;
        catch ne_eo;
            rethrow(ne_eo);
        end
    end

    % otherwise
    error( ...
        'xff:BadArgument', ...
        'Invalid tffspec argument given.' ...
    );

end

% get shortcuts for specs
arf = tffspec.ArrayFormat;
aio = ~tffspec.BinaryIO;
cfn = tffspec.CheckForNaNs;
fds = tffspec.FieldDelimiters;
ldc = double(tffspec.SkipEmptyLines(1));
lds = tffspec.LineDelimiters;
mfl = num2str(abs(floor(tffspec.MaxFieldNameLength(1))) + 1);
par = tffspec.ParagraphArrays(1);

% check specs
if isempty(regexpi(mfl, '^\d+$'))
    error( ...
        'xff:BadTFFSpec', ...
        'Invalid MaxFieldNameLength content.' ...
    );
end

% initialize rules variables
rule  = tffspec.ListOfFields;
rules = length(rule);

% get loops shortcut and prepare loopi struct, loopx
loops = tffspec.Loops;
loopf = fieldnames(loops);
if ~isempty(loopf)
    loopf = fieldnames(loops.(loopf{1}));
    loopi = emptystruct(loopf);
else
    loopi = struct;
end
loopx = {};

% get file extension
[fnamex{1:3}] = fileparts(filename);
fnshort = [fnamex{2} fnamex{3}];
fnamex = fnamex{3};
fnamex(fnamex == '.') = [];

% initialize output
if isempty(tffcont)
    tffcont = struct;
end

% initialize internal variables
namevars = struct;
namevars.EXTENSION  = fnamex;
namevars.FILENAME   = filename;
namevars.FORCEMODE  = forcemode;
namevars.HEADERONLY = headeronly && ~writemode;
namevars.NEXTLINE = 1;
namevars.TFFREAD    = ~writemode;
namevars.TFFVERSION = tffversion;
namevars.TFFWRITE   =  writemode;
namevars.VERBOSE    = verbosemode;

% reading TFF content
if ~writemode

    % read entire file into 1xN char array
    try
        if aio
            textcont = asciiread(filename);
        else
            textcont = binread(filename);
            textlen = numel(textcont);
            textpos = 0;
        end
    catch ne_eo;
        error( ...
            'xff:FileNotReadable', ...
            'File not readable: ''%s'' (%s).', ...
            filename, ne_eo.message ...
        );
    end

    % the rest only for ASCII-based IO
    if aio

        % make sure line delimiter is correct
        for ldsc = 1:numel(lds)
            if ldsc < numel(lds)
                textcont = strrep(textcont, lds{ldsc}, lds{end});
            else
                lds = lds{ldsc};
            end
        end

        % no line delimiter found
        if ~ischar(lds)
            error( ...
                'xff:BadTFFCont', ...
                'No specified line delimiter detected.' ...
            );
        end

        % split to lines
        linecont  = splittocellc(textcont, lds, logical(ldc));
        if ldc && ...
           ~isempty(linecont) && ...
            isempty(linecont{1})
            linecont(1) = [];
        end
        if ldc
            linecont = regexprep(linecont, '^\s+$', '');
            elines = false(1, numel(linecont));
            for linecount = 1:numel(elines)
                if isempty(linecont{linecount})
                    elines(linecount) = true;
                end
            end
            linecont(elines) = [];
        end
        linecont  = linecont(:);
        linecount = length(linecont);

        % get file "size"
        namevars.LINECOUNT = linecount;
    end

% writing TFF content
else

    % BeforeWriteCode ?
    if ~isempty(tffspec.BeforeWriteCode)
        try
            eval(tff_parsecode(tffspec.BeforeWriteCode(:)'));
        catch ne_eo;
            error( ...
                'xff:EvaluationError', ...
                'Error raised by BeforeWriteCode: ''%s''.', ...
                ne_eo.message ...
            );
        end
    end

    % create empty variables
    linecont  = {};
    linecount = 0;

    % get standard variables
    fds = fds{1};
    lds = lds{1};

end

% make sure to parse code only once
for rulec = 1:rules
    if ischar(rule(rulec).dim) && ...
        isempty(strfind(rule(rulec).dim, 'namevars')) && ...
        isempty(strfind(rule(rulec).dim, 'tffcont'))
        try
            eval(['rule(rulec).dim=[' rule(rulec).dim '];']);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            warning( ...
                'xff:tffio:DimConversionError', ...
                'Error converting dim statement ''%s'' in rule %d: %s.', ...
                rule(rulec).dim, rulec, ne_eo.message ...
            );
        end
    end
end

% grand loop
linec = 1;  % number of next line in file
rulec = 1;  % pointer to next rule to process
loopc = 0;  % counter of loops (in which loop are we)

% try used for forcemode (no indentation !!!)
try % forcemode

while rulec <= rules

    % get shortcut to current rule to process (and subfields)
    prule = rule(rulec);
    rtype = lower(prule.type);
    rdata = lower(prule.datatype);
    rcond = prule.cond;
    rdim  = prule.dim;
    rexpr = prule.varname;
    rform = prule.format;

    % try to resolve dim
    if isempty(rdim)
        rdim = 1;
    elseif ischar(rdim)
        try
            eval(['rdim=[' rdim '];']);
            if isempty(rdim) || ...
                size(rdim, 2) > 2
                error('BADDIM');
            end
        catch ne_eo;
            error( ...
                'xff:BadExpression', ...
                'Bad DIM expression: ''%s'' (%s).', ...
                prule.dim, ne_eo.message ...
            );
        end
    end

    % special case, flist & write is done by field!
    if writemode && ...
        strcmp(rtype, 'flist')
        rtype = 'field';
    end

    % what to do
    switch (rtype)

        % array
        case {'array'}

            % check cond
            if ~isempty(rcond)
                try
                    eval(['iofield=false;if ' rcond ',iofield=true;end']);
                catch ne_eo;
                    error( ...
                        'xff:BadExpression', ...
                        'Couldn''t evaluation COND expression: ''%s'' (%s).', ...
                        rcond, ne_eo.message ...
                    );
                end
                if ~iofield
                    rulec = rulec + 1;
                    continue;
                end
            end

            % check dim
            if length(rdim) ~= 2
                error( ...
                    'xff:BadTFFSpec', ...
                    'Dim for ARRAY must be 2-D.' ...
                );
            end
            rdimrows = rdim(1);
            rdimcols = rdim(2);

            % read access
            if ~writemode

                % for ASCII access
                if aio

                    % initialize array
                    rarray = zeros(rdim);

                    % skip paragraphing line?
                    if par && ...
                        ldc == 0
                        linec = linec + 1;
                    end

                    % skip empty column paragraph
                    if rdimcols == 0 && ...
                        ldc
                        rdimrows = 0;
                    end

                    % check number of lines/columns to come
                    if (linec + rdimrows - 1) > linecount
                        error( ...
                            'xff:BadTTFCont', ...
                            'Too few lines to process ARRAY for %s in %s.', ...
                            rexpr, fnshort ...
                        );
                    end

                % binary access
                else

                    % for anything but strings
                    if ~any(strcmpi(rdata, {'char', 'cstring', 'string'}))

                        % read rarray
                        [rarray, textnen] = ...
                            u8str2double(textcont, rdimrows, rdimcols, textpos);

                        % new position
                        textpos = textlen - textnen;

                        % change class?
                        if strcmpi(rdata, 'double')
                            eval(['tffcont.'  rexpr '=rarray; clear rarray;']);
                        else
                            eval(['tffcont.'  rexpr ...
                                '=' rdata '(rarray); clear rarray;']);
                        end

                    % strings, etc.
                    else
                    end

                    % don't do the rest
                    continue;
                end

                % depending on datatype
                switch (lower(rdata))

                    % built-in
                    case { ...
                         'char',  'int8',  'int16',  'int32', ...
                                 'uint8', 'uint16', 'uint32', ...
                        'single', 'double', 'logical'}

                        % set function tokens
                        aconvopen  = lower(rdata);

                    % specially deal with strings
                    case {'cstring', 'string'}

                        % get data to read
                        readdata = linecont(linec:linec + rdimrows - 1);

                        % create new array
                        rdcell = cell(rdim(1), rdimcols);

                        % treat correctly formatted lines differently
                        if ~isempty(rform) && ...
                            rform(1) ~= '%' && ...
                            rform(1) == rform(end)

                            % get formating encloser
                            nrf = rform(1);

                            % iterate over lines
                            for arc = 1:rdimrows

                                % match with \%([^\%]*)\%\s*
                                [sarra{1:3}] = regexp(readdata{arc}, ...
                                    ['\' nrf '([^\' nrf ']*)\' nrf '\s*']);
                                sarrt = sarra{3};

                                % no match found
                                if isempty(sarrt)

                                    % fill cells
                                    rdcell{arc, 1} = readdata{arc};
                                    for acc = 2:rdimcols
                                        rdcell{arc, acc} = '';
                                    end

                                % matches found
                                else

                                    % over cells
                                    for acc = 1:rdimcols

                                        % if match for cell
                                        if length(sarrt) >= acc

                                            % fill cee
                                            rdcell{arc, acc} = ...
                                                readdata{arc}( ...
                                                sarrt{acc}(1, 1):sarrt{acc}(1, 2));

                                        % no match for cell
                                        else

                                            % default
                                            rdcell{arc, acc} = '';
                                        end
                                    end
                                end
                            end

                        % no correct delimiter
                        else

                            % set rdcell to data
                            rdcell = readdata;

                        end

                        % simply copy string lines (discard second dim)
                        try
                            eval(['tffcont.'  rexpr '=rdcell; clear rdcell;']);
                        catch ne_eo;
                            error( ...
                                'xff:BadTFFCont', ...
                                'Could not copy string list into %s (%s).', ...
                                rexpr, ne_eo.message ...
                            );
                        end
                        linec = linec + rdimrows;
                        rulec = rulec + 1;
                        continue;

                    % special
                    otherwise
                        aconvopen  = ['array2' rdata];
                end

                % try conversion to double first
                try
                    if rdimrows > 0
                        %eval(['rarray(:, :) = [' ...
                        %    gluetostringc(linecont(linec:linec + rdimrows - 1), ';') ...
                        %    '];']);
                        rarray(:, :) = u8str2double(gluetostringc( ...
                            linecont(linec:linec + rdimrows - 1), ';'), ...
                            rdim(1), rdim(2));
                        if cfn && ...
                            any(isnan(rarray(:)))
                            warning( ...
                                'xff:BadTFFCont', ...
                                'Could not parse numeric ARRAY %s with u8str2double; retrying.', ...
                                rexpr ...
                            );
                            error('RETRY');
                        end
                    end

                % possibly, the lines are badly formatted (Brain Voyager QX v1.10/SDM)
                catch ne_eo;

                    % keep track
                    neuroelf_lasterr(ne_eo);

                    % check lines in length and then try to split
                    llen = numel(linecont{linec});
                    for linecc = (linec + 1):(linec + rdimrows - 1)
                        if numel(linecont{linecc}) ~= llen
                            llen = -1;
                            break;
                        end
                    end

                    % see if all lines had the same length
                    if llen > 0 && ...
                        mod(llen, rdimcols) == 0
                        llen = llen / rdimcols;
                    else
                        llen = 0;
                    end

                    % if so
                    if llen > 0

                        % add another space after each llen
                        for linecc = linec:(linec + rdimrows - 1)
                            linecn = reshape(linecont{linecc}, llen, rdimcols);
                            linecn(llen + 1, :) = ' ';
                            linecont{linecc} = linecn(:)';
                        end

                        % another conversion attempt, set llen = 0 if fails
                        try
                            eval(['rarray(:, :) = [', ...
                                gluetostringc(linecont(linec:linec + rdimrows - 1), ';') ...
                            '];']);
                        catch ne_eo;
                            neuroelf_lasterr(ne_eo);
                            llen = 0;
                        end

                        % conversion OK?
                        if llen == 0
                            warning( ...
                                'xff:BadTFFCont', ...
                                'Could not parse numeric ARRAY %s.', ...
                                rexpr ...
                            );
                        end
                    end
                end

                % try conversion now
                if ~strcmp(aconvopen, 'double')
                    try
                        eval(['rarray=' aconvopen '(rarray);']);
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        error( ...
                            'xff:ConversionFailed', ...
                            'Couldn''t convert ARRAY from double=>%s.', ...
                            aconvopen ...
                        );
                    end
                end

                % store array
                try
                    eval(['tffcont.'  rexpr '=rarray;']);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    error( ...
                        'xff:EvaluationError', ...
                        'Couldn''t store ARRAY into struct: ''%s''.', ...
                        rexpr ...
                    );
                end

                % increase linec
                linec = linec + rdimrows;

            % write access
            else

                % initialize array
                larray = cell(rdimrows, 1);
                carray = cell(1, rdimcols);

                % get data to write
                try
                    eval(['writedata=tffcont.' rexpr ';']);
                catch ne_eo;
                    error( ...
                        'xff:EvaluationError', ...
                        'Data field evaluation error: ''%s'' (%s).', ...
                        rexpr, ne_eo.message ...
                     );
                end

                % check dim
                if length(size(writedata)) ~= 2 || ...
                    any(size(writedata) ~= rdim)
                    error( ...
                        'xff:BadTFFSpec', ...
                        'ARRAY dimension mismatch.' ...
                    );
                end

                % paragraph array ?
                if par
                    linecont{end+1} = '';
                    linecount = linecount + 1;
                end

                % for full format notation
                if sum(rform == '%') == size(writedata, 2) && ...
                    size(writedata, 2) > 1 && ...
                   ~strcmpi(rdata, 'string')

                    % use clever sprintf and continue
                    linecont{end+1} = sprintf(rform, double(writedata'));
                    linecount = linecount + 1;
                    rulec = rulec + 1;
                    continue;
                end

                % depending on datatype
                switch (lower(rdata))

                    % built-in
                    case { ...
                         'char',  'int8',  'int16',  'int32', ...
                                 'uint8', 'uint16', 'uint32', ...
                        'single', 'double', 'logical'}

                        % set function tokens
                        aconvopen  = lower(rdata);

                    % handle strings specifically
                    case {'cstring', 'string'}

                        % check array type
                        if ~iscell(writedata)
                            error( ...
                                'xff:BadTFFCont', ...
                                'String ARRAY must be of type cell.' ...
                            );
                        end

                        % just one columns
                        if rdimcols == 1 && ...
                            strcmp(rform, '%s')

                            % put strings into file
                            try
                                larray(1:end) = writedata;
                            catch ne_eo;
                                error( ...
                                    'xff:AssignmentFailed', ...
                                    'Could not copy string array from %s: %s.', ...
                                    rexpr, ne_eo.message ...
                                );
                            end

                        % more columns
                        else

                            % iterate over rows
                            numcols = size(writedata, 2);
                            for arc = 1:rdimrows

                                % format line
                                try
                                    if numcols > 0
                                        larray{arc} = deblank(sprintf([rform ' '], ...
                                            writedata{arc, :}));
                                    else
                                        larray{arc} = '';
                                    end
                                catch ne_eo;
                                    error( ...
                                        'xff:SprintfError', ...
                                        'Invalid sprintf call: ''%s''->''%s''.', ...
                                        ['sprintf(''' rform ' '', ...)'], ...
                                        ne_eo.message ...
                                    );
                                end

                            end
                        end

                        % put lines into file content
                        linecont(end + 1:end + rdimrows) = larray;
                        linecount = linecount + rdimrows;
                        rulec = rulec + 1;
                        continue;

                    % special
                    otherwise
                        aconvopen  = [rdata '2array'];
                end

                % try conversion
                if ~strcmpi(rdata, 'double')
                    try
                        eval(['writedata=double(' aconvopen '(writedata));']);
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        error( ...
                            'xff:ConversionFailed', ...
                            'Couldn''t convert ARRAY from %s=>double.', ...
                            aconvopen ...
                        );
                    end
                end

                % combine data over rows
                try
                    for rowc = 1:rdimrows

                        % build array elements
                        for colc = 1:rdimcols
                            carray{colc} = sprintf( ...
                                rform, writedata(rowc, colc));
                        end

                        % build entry with sprintf
                        larray{rowc} = deblank(sprintf([arf fds], carray{:}));
                    end
                catch ne_eo;
                    error( ...
                        'xff:FormatError', ...
                        'Error formatting ARRAY Element: ''%s''.', ...
                        ne_eo.message ...
                    );
                end

                % put lines into file content
                linecont(end + 1:end + rdimrows) = larray;
                linecount = linecount + rdimrows;

            end

        % begin of a loop
        case {'bloop'}

            % check loop var syntax (and extract dim, if needed)
            [ldimmatcha{1:3}] = regexpi( ...
                rexpr, '^([a-z][a-z_0-9]*)(\(\d+\))?');
            ldimmatcht = ldimmatcha{3};
            if isempty(ldimmatcht)
                error( ...
                    'xff:BadBFFSpec', ...
                    'Invalid LOOP variable name given: ''%s''.', ...
                    rexpr ...
                );
            end
            lvname = rexpr(ldimmatcht{1}(1, 1):ldimmatcht{1}(1, 2));
            if size(ldimmatcht{1}, 1) > 1 && ...
                ldimmatcht{1}(2, 2) > ldimmatcht{1}(2, 1)
                try
                    eval(['ldim = [' rexpr(ldimmatcht{1}(2, 1):ldimmatcht{1}(2, 2)) '];']);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    ldim = 1;
                end
            else
                ldim = 1;
            end

            % check whether appropriate entry in Loops.(...) exists
            if ~isfield(loops, lvname) || ...
                length(loops.(lvname)) < ldim
                error( ...
                    'xff:BadBFFSpec', ...
                    'Invalid LOOP variable given (not found: ''%s''.', ...
                    rexpr ...
                );
            end

            % get specific loop info
            loopinfo = loops.(lvname)(ldim);

            % check whether the loop should be entered at all
            enterloop = true;
            if ~isempty(rcond)
                enterloop = false;
                try
                    eval(['if ' rcond ',enterloop = true;end']);
                catch ne_eo;
                    error( ...
                        'xff:BadExpression', ...
                        'Couldn''t evaluate COND expression: ''%s'' (%s).', ...
                        rcond, ne_eo.message ...
                    );
                end
            end

            % if we're entering the loop
            if enterloop

                % try to resolve loop dim
                try
                    if ischar(loopinfo.dim)
                        loopinfo.dim = eval(tff_parsecode(loopinfo.dim));
                    end
                catch ne_eo;
                    error( ...
                        'xff:BadExpression', ...
                        'Couldn''t evaluate LOOP.DIM expression: ''%s'' (%s).', ...
                        loopinfo.dim, ne_eo.message ...
                    );
                end

                % initialize loop counter
                try
                    eval(['namevars.' rexpr '=1;']);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    error( ...
                        'xff:BadExpression', ...
                        'Error initializing LOOP counter ''%s'' to 1.', ...
                        rexpr ...
                    );
                end

                % only truly enter loop if dim is > 0
                if loopinfo.dim > 0

                    % increase loop counter and keep track of loop variable name
                    loopc = loopc + 1;
                    loopi(loopc) = loopinfo;
                    loopx{loopc} = rexpr;

                % otherwise
                else

                    % skip loop
                    rulec = loopinfo.lastrule;
                end

            % we're not entering the loop
            else

                % set current rule pointer to end of loop
                rulec = loopinfo.lastrule;

            end

        % end of a loop
        case {'eloop'}

            % check name of LAST ENTERED loop
            if ~strcmp(rexpr, loopx{end})
                error( ...
                    'xff:BadBFFSpec', ...
                    'Invalid LOOP end token found: ''%s''.', ...
                    rexpr ...
                );
            end

            % increase loop counter and check if the loop is done
            try
                leaveloop = false;
                eval([...
                    'namevars.' rexpr '=namevars.'  rexpr '+1;' ...
                    'if namevars.' rexpr '>loopi(loopc).dim,leaveloop=true;end']);
            catch ne_eo;
                error( ...
                    'xff:BadExpression', ...
                    'Couldn''t increase/check LOOP counter: ''%s'' (%s).', ...
                    rexpr, ne_eo.message ...
                );
            end

            % if we're NOT leaving loop
            if ~leaveloop

                % set next rule to process to first loop line
                rulec = loopi(loopc).firstrule;

            % if we're leaving loop
            else

                % pop loopi/x(loopc)
                loopi(loopc) = [];
                loopx(loopc) = [];
                loopc = loopc - 1;

            end

        % expression
        case {'expre'}

            % if no given condition is given
            if isempty(rcond)

                % simply execute expression
                try
                    eval(rexpr);
                catch ne_eo;
                    error( ...
                        'xff:BadExpression', ...
                        'Couldn''t evaluate EXPRE: ''%s'' (%s).', ...
                        rexpr, ne_eo.message ...
                    );
                end

            % with a given condition
            else

                % execute upon condition
                try
                    eval(['if ' rcond ',' rexpr ';end']);
                catch ne_eo;
                    error( ...
                        'xff:BadExpression', ...
                        'Couldn''t evaluate EXPRE: ''%s'' (%s).', ...
                        ['if ' rcond ',' rexpr ';end'], ne_eo.message ...
                    );
                end
            end

        % field (read/write)
        case {'field'}

            % check cond
            if ~isempty(rcond)
                try
                    eval(['iofield=false;if ' rcond ',iofield=true;end']);
                catch ne_eo;
                    error( ...
                        'xff:BadExpression', ...
                        'Couldn''t evaluation COND expression: ''%s'' (%s).', ...
                        rcond, ne_eo.message ...
                    );
                end
                if ~iofield
                    rulec = rulec + 1;
                    continue;
                end
            end

            % reading
            if ~writemode

                % check number of lines/columns to come
                if linec > linecount
                    error( ...
                        'xff:BadTTFCont', ...
                        'No more lines left to read FIELD.' ...
                    );
                end

                % get deblanked field line
                fieldline = deblank(linecont{linec});

                % check for name
                [linema{1:3}] = regexpi(fieldline, ...
                    '^\s*([a-z][a-z_0-9]*)\:?\s*(.*)$');
                linemt = linema{3};
                if ~isempty(linemt)
                    fieldname = fieldline(linemt{1}(1, 1):linemt{1}(1, 2));
                    fieldval  = fieldline(linemt{1}(2, 1):linemt{1}(2, 2));
                else
                    fieldname = '';
                    fieldval  = fieldline;
                end

                % give (only) a warning if name does NOT match
                if isempty(strfind(prule.field, fieldname))
                    warning( ...
                        'xff:BadTFFSpec', ...
                        ['Possibly reading wrong field, ' ...
                         '''%s'' instead of ''%s''.'], ...
                        fieldname, rexpr ...
                    );
                end

                % depending on datatype
                switch (lower(rdata))

                    % built-in
                    case { ...
                         'char',  'int8',  'int16',  'int32', ...
                                 'uint8', 'uint16', 'uint32', ...
                        'single', 'double', 'logical'}

                        % first to double then to value
                        try
                            eval(['tffcont.'  rexpr '=' ...
                                   lower(rdata) '(str2num(fieldval));']);
                        catch ne_eo;
                            error( ...
                                'xff:EvaluationError', ...
                                'Could not store %s field in %s (%s).', ...
                                rdata, rexpr, ne_eo.message ...
                            );
                        end

                    % handle strings specifically
                    case {'cstring', 'string'}

                        % format with delimiter, remove it then
                        if ~isempty(rform) && ...
                            length(fieldval) > 1 && ...
                            rform(1) == rform(end) && ...
                            fieldval(1) == rform(1) && ...
                            fieldval(end) == rform(1)
                            fieldval = fieldval(2:end - 1);
                        end

                        % simply copy content
                        try
                            eval(['tffcont.'  rexpr '=fieldval;']);
                        catch ne_eo;
                            error( ...
                                'xff:EvaluationError', ...
                                'Could not store string field in %s (%s).', ...
                                rexpr, ne_eo.message ...
                            );
                        end

                    % special
                    otherwise

                        % try custom conversion (only once)
                        try
                            eval(['tffcont.'  rexpr ...
                                  '=string2' rdata '(fieldval);']);
                        catch ne_eo;
                            error( ...
                                'xff:ConversionError', ...
                                'Could not convert string=>%s field for %s (%s).', ...
                                rdata, rexpr, ne_eo.message ...
                            );
                        end
                end

                % increase line counter
                linec = linec + 1;

            % writing
            else

                % get data to write
                try
                    eval(['writedata=tffcont.' rexpr ';']);
                catch ne_eo;
                    error( ...
                        'xff:EvaluationError', ...
                        'Data field evaluation error: ''%s'' (%s).', ...
                        rexpr, ne_eo.message ...
                     );
                end

                % check dim
                if length(size(writedata)) ~= 2 || ...
                    (~ischar(writedata) && ...
                      length(writedata) ~= rdim(1))
                    error( ...
                        'xff:BadTFFSpec', ...
                        'FIELD dimension mismatch.' ...
                    );
                end

                % get line start
                linestart = sprintf(['%-' mfl 's'], [prule.field ':']);

                % depending on datatype
                switch (lower(rdata))

                    % built-in but hopefully unconverted
                    case {'double'}

                        % first to double then to value
                        try
                            if ~isa(writedata, 'double')
                                writedata=double(writedata);
                            end
                            linevalue = deblank(sprintf([rform ' '], writedata));
                        catch ne_eo;
                            error( ...
                                'xff:ConversionFailed', ...
                                'Couldn''t convert from %s=>double (%s).', ...
                                class(writedata), ne_eo.message ...
                            );
                        end

                    % built-in but definitely converted
                    case { ...
                         'char',  'int8',  'int16',  'int32', ...
                                 'uint8', 'uint16', 'uint32', ...
                        'single', 'logical'}

                        % first to double then to value
                        try
                            eval(['writedata=double(' ...
                                   lower(rdata) '(writedata));']);
                            linevalue = deblank(sprintf([rform ' '], writedata));
                        catch ne_eo;
                            error( ...
                                'xff:ConversionFailed', ...
                                'Couldn''t convert from %s=>double (%s).', ...
                                rdata, ne_eo.message ...
                            );
                        end

                    % specifically handle strings
                    case {'cstring', 'string'}

                        % format
                        if ~isempty(rform) && ...
                            any(rform == '%')
                            try
                                writedata = ...
                                    deblank(sprintf([rform ' '], writedata));
                            catch ne_eo;
                                error( ...
                                    'xff:SprintfFailed', ...
                                    'Error using sprintf formatting: %s.', ...
                                    ne_eo.message ...
                                );
                            end
                        end

                        % simply assign value
                        linevalue = writedata;

                    % otherwise use handler function
                    otherwise

                        % try conversion
                        try
                            eval(['linevalue=double(' ...
                                   rtype '2string(writedata));']);
                        catch ne_eo;
                            error( ...
                                'xff:ConversionFailed', ...
                                'Couldn''t convert %s from %s=>string (%s).', ...
                                rexpr, rdata, ne_eo.message ...
                            );
                        end

                end

                % put line into stream
                linecont{end + 1} = [linestart fds linevalue];
                linecount = linecount + 1;

            end

        % field list (ONLY READ, write by field)
        case {'flist'}

            % get list of fields to come
            fieldlist = {};
            fieldlookup = struct;
            while rulec <= rules
                if ~strcmpi(rule(rulec).type, 'flist')
                    rulec = rulec - 1;
                    break;
                end

                % get shortcut to rule
                sprule = rule(rulec);
                srcond = sprule.cond;
                rulec  = rulec + 1;

                % check condition
                if ~isempty(srcond)
                    try
                        eval(['flinclude=false;if ' ...
                               srcond ',flinclude=true;end']);
                        if ~flinclude
                            continue;
                        end
                    catch ne_eo;
                        error( ...
                            'xff:EvaluationError', ...
                            'Error evaluating COND on FLIST ''%s'' (%s).', ...
                            sprule.field, ne_eo.message ...
                        );
                    end
                end

                fieldlist{end + 1} = lower(sprule.field);
                fieldlookup.(makelabel(fieldlist{end})) = rulec - 1;
            end

            % don't do anything on empty fieldlist
            if isempty(fieldlist)
                continue;
            end

            % build pattern matching terms
            fieldmlist = ['^\s*(' lower(gluetostringc(fieldlist, '\:|')) ...
                '\:)\s*(.*)$'];

            % check number of lines/columns to come
            if linec > linecount
                error( ...
                    'xff:BadTTFCont', ...
                    'No more lines left to read FLIST.' ...
                );
            end

            % read as long as pattern matches
            fieldline = deblank(linecont{linec});
            [flinea{1:3}] = ...
                regexpi(fieldline, fieldmlist);
            flinet = flinea{3};
            while ~isempty(fieldmlist) && ...
               ~isempty(flinet) && ...
                linec <= linecount

                % get fieldname that matched (without ':')
                fieldname = fieldline(flinet{1}(1, 1):flinet{1}(1, 2) - 1);
                fieldval  = fieldline(flinet{1}(2, 1):flinet{1}(2, 2));

                % lookup rule
                try
                    srulec = fieldlookup.(makelabel(lower(fieldname)));
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    error( ...
                        'xff:LookupFailed', ...
                        'Field lookup from struct failed: ''%s''.', ...
                        fieldname ...
                    );
                end

                % remove fieldname from fieldlist and rebuild fieldmlist
                fieldlist(strcmpi(fieldname, fieldlist)) = [];
                fieldmlist = ['^(' lower(gluetostringc(fieldlist, '\:|')) ...
                    '\:)\s*(.*)$'];

                % get sub rule
                sprule = rule(srulec);
                srtype = sprule.type;
                srdata = lower(sprule.datatype);
                srcond = sprule.cond;
                srexpr = sprule.varname;
                srform = sprule.format;

                % double check type and field
                if ~strcmpi(srtype, 'flist') || ...
                   ~strcmpi(sprule.field, fieldname)
                    error( ...
                        'xff:BadTFFSpec', ...
                        'Field type or name don''t match FLIST spec.' ...
                    );
                end

                % depending on what data
                switch (lower(srdata))

                    % built-in
                    case { ...
                         'char',  'int8',  'int16',  'int32', ...
                                 'uint8', 'uint16', 'uint32', ...
                        'single', 'double', 'logical'}

                        % first to double then to value
                        try
                            eval(['tffcont.'  srexpr '=' ...
                                   lower(srdata) '(str2num(fieldval));']);
                        catch ne_eo;
                            error( ...
                                'xff:EvaluationError', ...
                                'Could not store %s field in %s: %s.', ...
                                srdata, srexpr, ne_eo.message ...
                            );
                        end

                    % handle strings specifically
                    case {'cstring', 'string'}

                        % format with delimiter, remove it then
                        if ~isempty(srform) && ...
                            length(fieldval) > 1 && ...
                            srform(1) == srform(end) && ...
                            fieldval(1) == srform(1) && ...
                            fieldval(end) == srform(1)
                            fieldval = fieldval(2:end - 1);
                        end

                        % simply copy content
                        try
                            eval(['tffcont.'  srexpr '=fieldval;']);
                        catch ne_eo;
                            error( ...
                                'xff:EvaluationError', ...
                                'Could not store string field in %s: %s.', ...
                                rexpr, ne_eo.message ...
                            );
                        end

                    % special
                    otherwise

                        % try custom conversion (only once)
                        try
                            eval(['tffcont.' srexpr ...
                                  '=string2' srdata '(fieldval);']);
                        catch ne_eo;
                            error( ...
                                'xff:ConversionError', ...
                                'Could not convert string=>%s field for %s (%s).', ...
                                srdata, srexpr, ne_eo.message ...
                            );
                        end
                end

                % increase line counter and check
                linec = linec + 1;
                if linec > linecount
                    break;
                end

                % get next content to check and perform check
                fieldline = deblank(linecont{linec});
                [flinea{1:3}] = ...
                    regexpi(fieldline, fieldmlist);
                flinet = flinea{3};
            end

        % skip N rules
        case {'skipn'}

            % conditional skip
            if ~isempty(rcond)
                performskip = false;
                try
                    eval(['if ' rcond ',performskip=true;end']);
                catch ne_eo;
                    error( ...
                        'xff:EvaluationError', ...
                        'Error evaluating SKIPN condition: ''%s'' (%s).', ...
                        rcond, ne_eo.message ...
                    );
                end

                % if not skip, continue with next rule
                if ~performskip
                    rulec = rulec + 1;
                    continue;
                end
            end

            % perform skip
            nskiptgt = rulec + prule.dim;
            nskiplst = nskiptgt - 1;

            % skip too high ?
            if nskiplst > rules
                error( ...
                    'xff:BadTFFSpec', ...
                    'Can''t skip %d rules at rule %d/%d.', ...
                    nskipdim, rulec, rules ...
                );
            end

            % correctly leave any loops on our way
            for trulec = rulec:nskiplst

                % get shortcut for rule
                tstruct = rule(trulec);

                % eloop
                if strcmpi(tstruct.type, 'eloop')

                    % matching with last entered loop
                    if ~strcmp(tstruct.varname, loopx{end})
                        error( ...
                            'xff:BadTFFSpec', ...
                            'Invalid LOOP (ELOOP) nesting found: ''%s''.', ...
                            tstruct.varname ...
                        );
                    end

                    % pop loop
                    loopi(end) = [];
                    loopx(end) = [];
                    loopc = loopc - 1;

                end
            end

            % set new rulec
            rulec = nskiptgt;

        % write a line
        case {'wrtln'}

            % check condition first
            if ~isempty(rcond)
                try
                    eval(['dowrtln=false;if ' rcond ',dowrtln=true;end']);
                catch ne_eo;
                    error( ...
                        'xff:EvaluationError', ...
                        'Error evaluating WRTLN condition: ''%s'' (%s).', ...
                        rcond, ne_eo.message ...
                    );
                end

                % if not skip, continue with next rule
                if ~dowrtln
                    rulec = rulec + 1;
                    continue;
                end
            end

            % reading
            if ~writemode

                % empty line & skip empty lines, do nothing
                if ldc && ...
                    isempty(rexpr)
                    rulec = rulec + 1;
                    continue;
                end

                % simply skip next line
                linec = linec + 1;

            % writing
            else

                % evaluate line?
                if ~isempty(rexpr) && ...
                    rexpr(1) == '%' && ...
                    rexpr(end) == '%'
                    try
                        rexpr = eval(rexpr(2:end-1));
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        rexpr = '';
                    end
                    if isempty(rexpr)
                        rexpr = strrep(ne_eo.message, char(10), ' ');
                    end
                end

                % write line
                linecont{end + 1} = rexpr;
                linecount = linecount + 1;
            end

        % exit of a loop
        case {'xloop'}

            % check name of ENTERED loops
            if ~any(strcmp(rexpr, loopx))
                error( ...
                    'xff:BadBFFSpec', ...
                    'Invalid XLOOP token found: ''%s''.', ...
                    rexpr ...
                );
            end

            % check cond, whether loop has to be left
            try
                if isempty(rcond)
                    leaveloop = true;
                else
                    eval([ ...
                        'leaveloop=false;if ' rcond ',leaveloop=true;end']);
                end
            catch ne_eo;
                error( ...
                    'xff:EvaluationFailed', ...
                    'Error evaluating XLOOP condition: ''%s'' (%s).', ...
                    rcond, ne_eo.message ...
                );
            end

            % if we're leaving loop
            if leaveloop

                % find loop to leave
                leaveloop = find(strcmp(rexpr, loopx));

                % set rulec to lastloop of this loop
                rulec = loopi(leaveloop).lastrule;

                % clear loopi and loopx, set loopc
                loopi(leaveloop:end) = [];
                loopx(leaveloop:end) = [];
                loopc = leaveloop - 1;

            end

        % otherwise bail out
        otherwise
            error( ...
                'xff:BadBFFSpec', ...
                'Invalid rule type ''%s''.', ...
                rtype ...
            );
    end

    % increase next rule counter
    rulec = rulec + 1;
end

% remaining input in read mode
if ~writemode && ...
    linec <= linecount
    tffcont.REMAININGLINES = linecont(linec:end);

% remaining content in write mode
elseif writemode && ...
    isfield(tffcont, 'REMAININGLINES') && ...
    iscell(tffcont.REMAININGLINES) && ...
   ~isempty(tffcont.REMAININGLINES) && ...
    ischar(tffcont.REMAININGLINES{1})
    linecont(end + 1:end + numel(tffcont.REMAININGLINES)) = ...
        tffcont.REMAININGLINES(:);
end

% forcemode error handling
catch ne_eo;

    % keep track of errors
    neuroelf_lasterr(ne_eo);

    if forcemode

        % only warning
        warning( ...
            'xff:ForcedOverride', ...
            'Forceing error override: ''%s''.', ...
            ne_eo.message ...
        );

    % otherwise
    else

        % simply don't allow this
        rethrow(ne_eo);

    end

end % of try forcemode

% writing file
if writemode
    filename = namevars.FILENAME;
    try
        binwrite(filename, gluetostringc(linecont, lds, true));
    catch ne_eo;
        error( ...
            'xff:ErrorWritingFile', ...
            'Could not write output file, %s: ''%s''.', ...
            filename, ne_eo.message ...
        );
    end
end

% print in verbosemode
if ~writemode && ...
    verbosemode
    disp(tffcont);
end

% give correct output
if ~writemode

    % make sure it has run time vars
    if ~isfield(tffcont, 'RunTimeVars') || ...
       ~isstruct(tffcont.RunTimeVars) || ...
        numel(tffcont.RunTimeVars) ~= 1
        tffcont.RunTimeVars = struct;
    end

    % and set in output
    varargout{1} = tffcont;
else
    if nargout > 0
        varargout = cell(1, nargout);
        varargout{1} = tffcont;
        if nargout > 1
            varargout{2} = filename;
        end
    end
end



% %%%% sub function



function parsed = tff_parsecode(parsed)

    % shortcut
    if isempty(parsed) || ...
       ~any(parsed == '$' | parsed == '@')
        return;
    end

    % find first occurance
    pattern = '(\$|\@)([a-z][a-z_0-9]*)';
    [pcm{1:3}] = regexpi(parsed, pattern);
    pcm = pcm{3};

    % loop until tokens are gone
    while ~isempty(pcm)

        % what kind of token
        switch (parsed(pcm{1}(1, 1)))

            % replace $var
            case {'$'}
                parsed = [ ...
                    parsed(1:pcm{1}(1, 1) - 1) 'namevars.' ...
                    parsed(pcm{1}(2, 1):end)];
            % @var
            case {'@'}
                parsed = [ ...
                    parsed(1:pcm{1}(1, 1) - 1) 'tffcont.' ...
                    parsed(pcm{1}(2, 1):end)];

        end

        % rematch
        [pcm{1:3}] = regexpi(parsed, pattern);
        pcm = pcm{3};
    end

    % any occurance of $$/@@
    parsed = strrep(strrep(parsed, '@@', 'tffcont'), '$$', 'namevars');

% end of function parsedcode = tff_parsecode(unparsed, nv, cnt)

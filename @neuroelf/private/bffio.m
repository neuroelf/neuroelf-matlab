function [varargout] = bffio(varargin)
% bffio  - read/write binary file with BFF spec (and content)
%
% FORMAT:       bffcont = bffio(filename, bffspec [, options ...])
%          OR   bffio(wfilename, bffspec, bffcont [, options ...])
%
% Input fields:
%
%       filename    filename of binary file to read
%       bffspec     either filename or content or spec struct of BFF
%       options     optional arguments, one of
%                    'force' , 'f', or '-f' (for returning partial content)
%                   'header' , 'h', or '-h' (for reading the header only)
%                   'verbose', 'v', or '-v' (for verbose reading/writing)
%
%       wfilename   filename of binary file to write
%       bffcont     file content (struct)
%
% Output fields:
%
%       bffcont     file contents struct (depending on BFF)
%
% See also bffdocu, bffparse
%
% Note: options can also be GNU-style combine ('-fh')

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

% BFF version number
bffversion = 'v1.0';

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
        'Bad or missing argument for bffio.' ...
    );
end
bffname = varargin{1}(:)';
if ischar(varargin{2})
    bffspec = varargin{2}(:)';
else
    bffspec = varargin{2};
end
bffcont = [];

% default options
forcemode   = false;
headeronly  = false;
verbosemode = false;
writemode   = false;

% bff content structure given
if nargin > 2 && ...
    isstruct(varargin{3}) && ...
   ~isempty(varargin{3})
    bffcont = varargin{3};
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

                % header only
                case {'h', 'header'}
                    headeronly = true;

                % verbose mode
                case {'v', 'verbose'}
                    verbosemode = true;

            end

            % also check option style
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

% BFF is a valid SPEC struct
if isstruct(bffspec) && ...
   ~isempty(bffspec) && ...
    isfield(bffspec, 'AfterReadCode') && ...
    ischar(bffspec.AfterReadCode) && ...
    isfield(bffspec, 'BeforeWriteCode') && ...
    ischar(bffspec.BeforeWriteCode) && ...
    isfield(bffspec, 'EncodingSyntax') && ...
    ischar(bffspec.EncodingSyntax) && ...
   ~isempty(bffspec.EncodingSyntax) && ...
    isfield(bffspec, 'Extensions') && ...
    iscell(bffspec.Extensions) && ...
    isfield(bffspec, 'FilenameMatch') && ...
    iscell(bffspec.FilenameMatch) && ...
    isfield(bffspec, 'ListOfFields') && ...
    isstruct(bffspec.ListOfFields) && ...
    isfield(bffspec, 'Loops') && ...
    isstruct(bffspec.Loops) && ...
    isfield(bffspec, 'Magic') && ...
    isstruct(bffspec.Magic) && ...
    isfield(bffspec, 'TransIOSize') && ...
    isa(bffspec.TransIOSize, 'double') && ...
    numel(bffspec.TransIOSize) == 1 && ...
    isfield(bffspec, 'Variables') && ...
    isstruct(bffspec.Variables)

% any other struct is rejected
elseif isstruct(bffspec)
    error( ...
        'xff:BadArgument', ...
        'Invalid struct BFF specification.' ...
    );

% for character arrays
elseif ischar(bffspec)

    % valid file contents (0x0a/0x0d delimited)
    if ...
        any(bffspec == char(10) | bffspec == char(13)) || ...
        exist(bffspec, 'file') == 2

        % parse content
        try
            if ~writemode
                varargout{1} = bffio( ...
                    bffname, ...
                    bffparse(bffspec), ...
                    varargin{:});
            else
                [varargout{:}] = bffio( ...
                    bffname, ...
                    bffparse(bffspec), ...
                    varargin{:});
            end
            return;
        catch ne_eo;
            rethrow(ne_eo);
        end
    end

    % otherwise
    error( ...
        'xff:BadArgument', ...
        'Invalid bffspec argument given.' ...
    );

end

% initialize rules variables
rule  = bffspec.ListOfFields;
rules = numel(rule);

% get critical transio size
tiosz = bffspec.TransIOSize;
tiole = bffspec.EncodingSyntax;

% get loops shortcut and prepare loopi struct, loopx
loops = bffspec.Loops;
loopf = fieldnames(loops);
if ~isempty(loopf)
    loopf = fieldnames(loops.(loopf{1}));
    loopi = emptystruct(loopf);
else
    loopi = struct;
end
loopx = {};

% get file extension
[fnamex{1:3}] = fileparts(bffname);
fnamex = fnamex{3};
fnamex(fnamex == '.') = [];

% initialize output
if isempty(bffcont)
    bffcont = struct;
end

% initialize internal variables
namevars = struct;
namevars.BFFVERSION = bffversion;
namevars.BFFREAD    = ~writemode;
namevars.BFFWRITE   =  writemode;
namevars.EXTENSION  = fnamex;
namevars.FILENAME   = bffname;
namevars.FILESIZE   = -1;
namevars.FORCEMODE  = forcemode;
namevars.HEADERONLY = headeronly && ~writemode;
namevars.VERBOSE    = verbosemode;

% reading BFF content
if ~writemode

    % create transio object template if size suggests it
    try
        if ~isinf(tiosz) && ...
            tiosz > 0
            tioobjt = struct(transio(bffname, tiole, 'uint32', 0, [1, 1]));
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        if exist(bffname, 'file') > 0
            warning( ...
                'xff:TransIOError', ...
                'Error using transio for file ''%s''.', ...
                bffname ...
            );
        end
        tiosz = -Inf;
    end

    % check file (opening in little endian syntax by default)
    fid = fopen(bffname, 'r', tiole);
    if fid < 1
        error( ...
            'xff:FileNotReadable', ...
            'File not readable: ''%s''.', ...
            bffname ...
        );
    end

    % get file size
    fseek(fid, 0, 1);
    flen = ftell(fid);
    fseek(fid, 0, -1);
    namevars.FILESIZE = flen;

% writing BFF content
else

    % BeforeWriteCode ?
    if ~isempty(bffspec.BeforeWriteCode)
        try
            eval(bffspec.BeforeWriteCode);
        catch ne_eo;
            error( ...
                'xff:EvaluationError', ...
                'Error raised by BeforeWriteCode: ''%s''.', ...
                ne_eo.message ...
            );
        end
    end

    % check open writable file (read and write now standard!)
    wbffname = [bffname '.tmp'];
    fid = fopen(wbffname, 'w+', tiole);
    if fid < 1
        error( ...
            'xff:FileNotWritable', ...
            'File not writable: ''%s''.', ...
            wbffname ...
        );
    end
end

% grand loop
loopc = 0;  % counter of loops (in which loop are we)

% try used for forcemode (on indentation !!!)
try % forcemode

% try to preparse rule dims
for rulec = 1:numel(rule)
    if ischar(rule(rulec).dim) && ...
        isempty(strfind(rule(rulec).dim, 'namevars')) && ...
        isempty(strfind(rule(rulec).dim, 'bffcont'))
        try
            eval(['rule(rulec).dim=[' rule(rulec).dim '];']);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            warning( ...
                'xff:bffio:DimConversionError', ...
                'Error converting dim statement ''%s'' in rule %d: %s.', ...
                rule(rulec).dim, rulec, ne_eo.message ...
            );
        end
    end
end

% loop until all rules are parsed
rulec = 1;
while rulec <= rules

    % get shortcut to current rule to process (and subfields)
    prule = rule(rulec);
    rdisk = lower(prule.disktype);
    rdata = lower(prule.datatype);
    rcond = prule.cond;
    rdim  = prule.dim;
    rexpr = prule.varname;

    % try to resolve dim
    if isempty(rdim)
        rdim = [1, 1];
    elseif ischar(rdim)
        try
            eval(['rdim=double([' rdim ']);']);
            if isempty(rdim)
                error('EMPTYDIM');
            elseif numel(rdim) == 1
                rdim = [1, rdim(1)];
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            error( ...
                'xff:BadExpression', ...
                'Couldn''t evaluate DIM expression: ''%s''.', ...
                prule.dim ...
            );
        end
    end

    % what to do
    switch (lower(prule.type))

        % begin of a loop
        case {'bloop'}

            % check loop var syntax (and extract dim, if needed)
            [ldimmatcht{1:3}] = regexpi( ...
                rexpr, '^([a-z][a-z_0-9]*)(\(\d+\))?');
            ldimmatcht = ldimmatcht{3};
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
                ldim = str2double(rexpr(ldimmatcht{1}(2, 1):ldimmatcht{1}(2, 2)));
            else
                ldim = 1;
            end

            % check whether appropriate entry in Loops.(...) exists
            if ~isfield(loops, lvname) || ...
                numel(loops.(lvname)) < ldim
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
                        loopinfo.dim = eval(loopinfo.dim);
                    else
                        loopinfo.dim = loopinfo.dim(end);
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
                        'Couldn''t evaluate EXPRE: ''%s''\nMessage: %s.', ...
                        ['if ' rcond ',' rexpr ';end'], ne_eo.message ...
                    );
                end
            end

        % field (read/write)
        case {'field'}

            % check condition
            readfield  = ~writemode;
            writefield =  writemode;
            if ~isempty(rcond)
                try
                    readfield  = false;
                    writefield = false;
                    eval(['if ' rcond ',' ...
                          'readfield=~writemode;' ...
                          'writefield=writemode;' ...
                          'end']);
                catch ne_eo;
                    error( ...
                        'xff:BadExpression', ...
                        'Couldn''t evaluation COND expression: ''%s'' (%s).', ...
                        rcond, ne_eo.message ...
                    );
                end
            end

            % if we're reading the data
            if readfield

                % treat transio correctly
                if (~isinf(tiosz) && ...
                    (prule.disksize * prod(rdim)) >= tiosz && ...
                     (strcmp(rdisk, rdata) || ...
                      tiosz == 1)) || ...
                   (strcmp(rdata, 'transio') && ...
                    tiosz > 0)
                    try
                        cfpos = ftell(fid);
                        rsize = prule.disksize * prod(rdim);
                        fseek(fid, rsize, 0);
                        if ftell(fid) ~= (cfpos + rsize)
                            error('SEEK_ERROR');
                        end
                        tio_obj = tioobjt;
                        tio_obj.DataType = rdisk;
                        tio_obj.TypeSize = prule.disksize;
                        tio_obj.IOOffset = cfpos;
                        obs = rdim;
                        while ~isempty(obs) && ...
                            obs(end) == 1
                            obs(end) = [];
                        end
                        if numel(obs) < 2
                            obs(numel(obs)+1:2) = 1;
                        end
                        tio_obj.DataDims = obs;
                        tio_obj = transio(0, 'makeobject', tio_obj);
                        eval(['bffcont.' rexpr '=tio_obj;']);
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        error( ...
                            'xff:TransIOFailed', ...
                            'Error creating transio object.' ...
                        );
                    end
                    rulec = rulec + 1;
                    continue;
                end

                % what disk-bound datatype
                switch (rdisk)

                    % built-in types
                    case { ...
                             'char',  'int8',  'int16',  'int32', ...
                            'uchar', 'uint8', 'uint16', 'uint32', ...
                            'int64', 'uint64','single', 'double'}

                        % do read verbosely ?
                        if verbosemode
                            readdata = verboseread(fid, rdim, rdisk, rexpr);
                        else
                            readdata = ...
                                fread(fid, [1, prod(rdim)], ['*' rdisk]);
                            readdata = reshape(readdata, rdim);
                        end

                    % C-style strings
                    case {'cstring', 'clinea'}

                        % string or line
                        if strcmp(rdisk, 'cstring')
                            ltchar = char(0);
                        else
                            ltchar = char(10);
                        end

                        % for multiple strings (dim > 1)
                        if prod(rdim) > 1

                            % generate cell structure for strings
                            readdata = cell(rdim);

                            % do read verbosely ?
                            if verbosemode

                                % read strings into structure
                                for strc = 1:prod(rdim)
                                    readdata{strc} = ...
                                        verbosereadcstring(fid, ...
                                        sprintf('%s{%d}', rexpr, strc), ltchar, flen);
                                end

                            % otherwise
                            else

                                % read strings into structure
                                for strc = 1:prod(rdim)
                                    readdata{strc} = freadcstring(fid, ltchar, flen);
                                end
                            end

                        % only one string
                        else

                            % do read verbosely ?
                            if verbosemode
                                readdata = verbosereadcstring(fid, rexpr, ltchar, flen);
                            else
                                readdata = freadcstring(fid, ltchar, flen);
                            end
                        end

                    % unknown datatype
                    otherwise
                        error( ...
                            'xff:BadBFFSpec', ...
                            'Invalid DISKTYPE specified: ''%s''.', ...
                            rdisk ...
                        );
                end

                % check if we need type conversion
                if ~strcmp(rdisk, rdata)

                    % what type
                    switch (rdata)

                        % built-in types
                        case { ...
                                 'int8',  'int16',  'int32',  'int64', ...
                                'uint8', 'uint16', 'uint32', 'uint64', ...
                                 'char', 'single', 'double', 'logical'}

                            % try to convert via internal function
                            try
                               % eval(['readdata=' rdata '(readdata);']);
                                readdata = typeswitch(readdata, rdata);
                            catch ne_eo;
                                error( ...
                                    'xff:ConversionFailed', ...
                                    'Built-in conversion failed: %s', ...
                                    ne_eo.message ...
                                );
                            end

                        % all other types
                        otherwise

                            % try to convert data
                            try
                                eval(['readdata=' ...
                                    rdisk '2' rdata '(readdata);']);
                            catch ne_eo;
                                error( ...
                                    'xff:EvaluationFailed', ...
                                    'Failed data conversion, %s=>%s: %s', ...
                                    rdisk, rdata, ne_eo.message ...
                                );
                            end
                    end
                end

                % put the data read into the structure
                try
                    eval(['bffcont.' rexpr '=readdata;']);
                catch ne_eo;
                    error( ...
                        'xff:EvaluationFailed', ...
                        'Error storing read data into struct: %s.', ...
                        ne_eo.message ...
                    );
                end

            % are we else writing the data
            elseif writefield

                % get data to write
                try
                    towritedata = [];
                    eval(['towritedata=bffcont.' rexpr ';']);
                catch ne_eo;
                    error( ...
                        'xff:EvaluationFailed', ...
                        'Error retrieving data from struct: ''%s''.', ...
                        ne_eo.message ...
                    );
                end

                % treat transio correctly
                if istransio(towritedata, true)
                    stio = struct(towritedata);
                    try
                        cfpos = ftell(fid);
                        tfid = 0;
                        twcls = stio.DataType;
                        twicl = ['*' twcls];
                        twsiz = stio.TypeSize;
                        twchk = 2097152 / twsiz;
                        if ischar(stio.FileName)
                            if stio.LittleND
                                tfid = fopen(stio.FileName, 'rb', 'ieee-le');
                            else
                                tfid = fopen(stio.FileName, 'rb', 'ieee-be');
                            end
                            fseek(tfid, stio.IOOffset, -1);
                            twnum = prod(stio.DataDims);
                            while twnum >= twchk
                                fwrite(fid, fread(tfid, [twchk, 1], twicl), twcls);
                                twnum = twnum - twchk;
                            end
                            if twnum > 0
                                fwrite(fid, fread(tfid, [twnum, 1], twicl), twcls);
                            end
                        elseif iscell(stio.FileName)
                            for tiofc = 1:numel(stio.FileName)
                                if stio.LittleND
                                    tfid = fopen(stio.FileName{tiofc}, 'rb', 'ieee-le');
                                else
                                    tfid = fopen(stio.FileName{tiofc}, 'rb', 'ieee-be');
                                end
                                fseek(tfid, stio.IOOffset(tiofc), -1);
                                twnum = prod(stio.DataDims(1:end-1));
                                while twnum >= twchk
                                    fwrite(fid, fread(tfid, [twchk, 1], twicl), twcls);
                                    twnum = twnum - twchk;
                                end
                                if twnum > 0
                                    fwrite(fid, fread(tfid, [twnum, 1], twicl), twcls);
                                end
                                if tfid > 0
                                    fclose(tfid);
                                    tfid = 0;
                                end
                            end
                        else
                            error( ...
                                'xff:TransIOFailed', ...
                                'Invalid filename in transio member %s.', ...
                                rexpr ...
                            );
                        end
                    catch ne_eo;
                        if tfid > 0
                            fclose(tfid);
                        end
                        error( ...
                            'xff:TransIOFailed', ...
                            'Error skipping transio object size in file: %s', ...
                            ne_eo.message ...
                        );
                    end
                    if tfid > 0
                        fclose(tfid);
                        tfid = 0;
                    end
                    stio.FileName = bffname;
                    stio.IOOffset = cfpos;
                    eval(['bffcont.' rexpr '=transio(0,''makeobject'',stio);']);
                    rulec = rulec + 1;
                    continue;
                end

                % assume built-in type
                builtintype = true;

                % what current datatype
                switch (rdata)

                    % built-in datatypes
                    case { ...
                             'int8',  'int16',  'int32',  'int64', ...
                            'uint8', 'uint16', 'uint32', 'uint64', ...
                             'char', 'single', 'double', 'logical', ...
                          'cstring', 'clinea'}

                        % do nothing here

                    % otherwise we need conversion!
                    otherwise
                        builtintype = false;
                end

                % if datatype isn't built-in
                if ~builtintype

                    % try to convert data
                    try
                        eval(['towritedata=' ...
                            rdata '2' rdisk '(towritedata);']);
                    catch ne_eo;
                        error( ...
                            'xff:EvaluationFailed', ...
                            'Failed data conversion, %s=>%s: %s', ...
                            rdata, rdisk, ne_eo.message ...
                        );
                    end

                % datatypes still must be converted
                elseif ~strcmp(rdisk, rdata)

                    % try internal conversion
                    try
                        eval(['towritedata=' rdisk '(towritedata);']);
                    catch ne_eo;
                        error( ...
                            'xff:EvaluationFailed', ...
                            'Failed data conversion, %s=>%s: %s', ...
                            rdata, rdisk, ne_eo.message ...
                        );
                    end
                end

                % write REAL built-in datatypes
                if ~any(strcmp(rdisk, {'cstring', 'clinea'}))

                    % verbose
                    if verbosemode
                        try
                            verbosewrite(fid, towritedata, rdisk, rexpr);
                        catch ne_eo;
                            error( ...
                                'xff:WriteFailed', ...
                                'Writing to file failed: ''%s''.', ...
                                ne_eo.message ...
                            );
                        end

                    % non-verbose
                    else
                        try
                            fwrite(fid, towritedata, rdisk);
                        catch ne_eo;
                            error( ...
                                'xff:WriteFailed', ...
                                'Writing to file failed: ''%s''.', ...
                                ne_eo.message ...
                            );
                        end
                    end

                % c-style strings
                else

                    % which term char
                    if strcmp(rdisk, 'cstring')
                        ltchar = 0;
                    else
                        ltchar = 10;
                    end

                    % one string
                    if ischar(towritedata)

                        % verbose
                        if verbosemode
                            try
                                verbosewritecstring(fid, towritedata, rexpr, ltchar);
                            catch ne_eo;
                                error( ...
                                    'xff:WriteFailed', ...
                                    'Writing to file failed: ''%s''.', ...
                                    ne_eo.message ...
                                );
                            end

                        % non-verbose
                        else
                            try
                                fwrite(fid, ...
                                    uint8(double(towritedata(:)')), ...
                                    'uint8');
                                fwrite(fid, ltchar, 'uint8');
                            catch ne_eo;
                                error( ...
                                    'xff:WriteFailed', ...
                                    'Writing to file failed: ''%s''.', ...
                                    ne_eo.message ...
                                );
                            end
                        end

                    % more strings
                    elseif iscell(towritedata)

                        % verbose
                        if verbosemode

                            % loop over strings
                            try
                                for strc = 1:numel(towritedata)
                                    verbosewritecstring( ...
                                        fid, towritedata{strc}, ...
                                        sprintf('%s{%d}', rexpr, strc), ltchar);
                                end
                            catch ne_eo;
                                error( ...
                                    'xff:WriteFailed', ...
                                    'Writing to file failed: ''%s''.', ...
                                    ne_eo.message ...
                                );
                            end

                        % non-verbose
                        else

                            % loop over strings
                            try
                                for strc = 1:numel(towritedata)
                                    fwrite(fid, ...
                                        uint8(double(towritedata{strc})), ...
                                        'uint8');
                                    fwrite(fid, ltchar, 'uint8');
                                end
                            catch ne_eo;
                                error( ...
                                    'xff:WriteFailed', ...
                                    'Writing to file failed: ''%s''.', ...
                                    ne_eo.message ...
                                );
                            end
                        end
                    else
                        error( ...
                            'xff:WriteFailed', ...
                            'Writing cstring(s) requires char or cell array.' ...
                        );
                    end
                end

            % end of 'if readfield, elseif writefield'
            end

        % skip N rules
        case {'skipn'}

            % extra check dim field
            if ~isa(prule.dim, 'double') || ...
                numel(prule.dim) ~= 1
                error( ...
                    'xff:BadBFFSpec', ...
                    'Invalid SKIPN dim (N) given.' ...
                );
            end

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

            % get dim
            nskiplst = min(rulec + prule.dim - 1, rules);

            % correctly leave any loops on our way
            for trulec = rulec:nskiplst

                % get shortcut for rule
                tstruct = rule(trulec);

                % eloop
                switch (lower(tstruct.type)), case {'eloop'}

                    % matching with last entered loop
                    if ~strcmp(tstruct.varname, loopx{end})
                        error( ...
                            'xff:BadBFFSpec', ...
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
            rulec = nskiplst + 1;

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
                lower(prule.type) ...
            );
    end

    % increase next rule counter
    rulec = rulec + 1;
end

% forcemode error handler
catch ne_eo;
    neuroelf_lasterr(ne_eo);

    % for writemode, close file, try to delete .tmp, and give error
    if writemode
        fclose(fid);
        ine_eo = ne_eo;
        try
            delete(wbffname);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
        rethrow(ine_eo);
    end

    % forcemode requested
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

end % try forcemode

% AfterReadCode ?
if ~isempty(bffspec.AfterReadCode)
    try
        eval(bffspec.AfterReadCode);
    catch ne_eo;
        error( ...
            'xff:EvaluationError', ...
            'Error raised by AfterReadCode: ''%s''.', ...
            ne_eo.message ...
        );
    end
end

% fclose on file
fclose(fid);

% print in verbosemode
if ~writemode && ...
    verbosemode
    disp(bffcont);
end

% give correct output
if ~writemode

    % make sure bffcont has runtime variables
    if ~isfield(bffcont, 'RunTimeVars') || ...
       ~isstruct(bffcont.RunTimeVars) || ...
        numel(bffcont.RunTimeVars) ~= 1
        bffcont.RunTimeVars = struct;
    end

    % and set in output
    varargout{1} = bffcont;
else

    % rename file
    try
        succ = renamefile(wbffname, bffname);
        if succ ~= 0
            error('MOVE_ERROR');
        end
    catch ne_eo;
        error( ...
            'xff:FileNotMovable', ...
            'Error moving file ''%s'' to ''%s'': %s.', ...
            wbffname, bffname, ne_eo.message ...
        );
    end

    if nargout > 0
        varargout{1} = bffcont;
        if nargout > 1
            varargout(2:nargout) = cell(1, nargout - 1);
        end
    end
end



% %%% subfunctions



function cstring = freadcstring(fid, tchar, flen)

    % read up to 256 bytes
    fpos = ftell(fid);
    rsz = flen - fpos;
    cstring = fread(fid, [1, min(256, rsz)], 'uint8=>char');

    % find occurrence of terminator
    cstrcmp = findfirst(cstring == tchar);

    % end mark is found
    if ~isempty(cstrcmp)

        % cut string
        if cstrcmp == 1
            cstring = '';
        else
            cstring = cstring(1:cstrcmp-1);
        end

        % seek back and return
        fseek(fid, fpos + cstrcmp, -1);
        return;

    % no more data
    elseif rsz <= 256

        % return anyway
        return;
    end

    % keep reading
    rsz = rsz - 256;
    while isempty(cstrcmp)

        % read next (up to) 256 bytes
        cbyte = fread(fid, [1, min(256, rsz)], 'uint8=>char');

        % search for terminator
        cstrcmp = findfirst(cbyte == tchar);

        % found
        if ~isempty(cstrcmp)

            % paste
            cstring = [cstring, cbyte(1:cstrcmp-1)];

            % seek back and return
            fseek(fid, fpos + numel(cstring) + 1, -1);
            return;
        end

        % simply paste and move on
        cstring = [cstring, cbyte];
        rsz = rsz - 256;
    end

% end of function cstring = freadcstring(fid)


function readdata = verboseread(fid, rdim, rtype, target)

    % make output
    fpos = ftell(fid);

    % read (at max) 16 bytes
    vout = fread(fid, [1, 16], 'uint8=>double');
    fseek(fid, fpos, -1);

    % make output
    disp(sprintf('0x%08X:%-48s -> %s (%s)', fpos, ...
         sprintf(' %02X', vout), target, rtype));

    % read data
    readdata = reshape(fread(fid, [1, prod(rdim)], ['*' rtype]), rdim);

    % print some data ?
    if ~isempty(readdata)
        disp(sprintf(' ->%s ...', ...
             sprintf(' %d', double(readdata(1:min(8, numel(readdata)))))));
    end

% end of function readdata = verboseread(fid, rdim, type)


function readdata = verbosereadcstring(fid, target, tchar, flen)

    % make output
    fpos = ftell(fid);

    % read (at max) 16 bytes
    vout = fread(fid, [1, 16], 'uint8=>double');
    fseek(fid, fpos, -1);

    % make output
    disp(sprintf('0x%08X:%-48s -> %s (%s)', fpos, ...
         sprintf(' %02X', vout), target, 'cstring'));

    % read data
    readdata = freadcstring(fid, tchar, flen);

    % print some data
    disp(sprintf(' -> %s...', readdata(1:min(60, numel(readdata)))));

% end of function readdata = verbosereadcstring(fid, rdim, target)


function verbosewrite(fid, wdata, wtype, target)

    % make output
    fpos = ftell(fid);

    % write data
    fwrite(fid, wdata, wtype);

    % rewind to previous pos
    if fseek(fid, fpos, -1) < 0
        warning( ...
            'xff:FileIOError', ...
            'Error seeking to %d in file ID %d.', ...
            fpos, fid ...
        );
        return;
    end

    % read at max 16 bytes
    bdata = fread(fid, [1, 16], 'uint8=>double');

    % and read at max 8 elements of type datatype
    fseek(fid, fpos, -1);
    tdata = fread(fid, [1, min(numel(wdata), 8)], [wtype '=>double']);

    % go back to end of file
    fseek(fid, 0, 1);

    % make output
    disp(sprintf('0x%08X:%-48s -> %s (%s)', fpos, ...
         sprintf(' %02X', bdata), target, wtype));

    % print some data ?
    if ~isempty(tdata)
        disp(sprintf(' ->%s ...', ...
             sprintf(' %d', tdata(1:min(8, numel(tdata))))));
    end

% end of function verbosewrite(fid, wdata, wtype, target)


function verbosewritecstring(fid, wstring, target, tchar)

    % use verbosewrite(...)
    verbosewrite(fid, uint8([double(wstring(:)'), tchar]), target);

% end of function verbosewritecstring(fid, wstring, target)


function data = typeswitch(data, type)
    switch (lower(type))
        case {'double'}
            data = double(data);
        case {'single'}
            data = single(data);
        case {'int32'}
            data = int32(data);
        case {'uint32'}
            data = uint32(data);
        case {'int16'}
            data = int16(data);
        case {'uint16'}
            data = uint16(data);
        case {'char'}
            data = char(data);
        case {'int8'}
            data = int8(data);
        case {'uint8'}
            data = uint8(data);
        case {'logical'}
            % first to double -> Octave 3.4.0 bug
            data = logical(double(data));
        case {'int64'}
            data = int64(data);
        case {'uint64'}
            data = uint64(data);
    end
% end of function typeswitch(data, type)

function [varargout] = xini(varargin)
% xini (class)
%
% creates an object for ini-file handling
%
% Input: for meaningful object construction, give filename!
%
% ___ Constructor methods: ___
% ----------------------------
%
% xini          returns a handle to an xini object
%       FORMAT: IniFileObject = xini(Filename [,convert]);
%       FORMAT: NewObject = xini('convert' | 'new');
%       FORMAT: FactoryObject = xini;
%
% ReleaseAllFiles       issue release for all open objects
% ReloadAllFiles        issue reload for all open objects
% ResetClass            resets the internal variable, useful for debugging
% SetIniString          create an IniFile object from an inistring
%
% ___ NOTES: ___
% --------------
%
%  *)   with conversion active, xini uses any2ascii, another part
%       of this toolbox to write non-character variables to disk;
%       when not active, only valid one-line char arrays are allowed!
%
%  **)  storage for the ini file settings is NOT part of the object, but
%       is requested and maintained by the singleton object itself, thus
%       only little memory is needed for argument handling, even with
%       bigger ini file, plus multiple functions can share one ini file!
%
% Using: asciiread, gluetostring, splittocell.

% Version:  v1.1b
% Build:    25011123
% Date:     Jan-11 2025, 11:20 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010 - 2016, 2025, Jochen Weber
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

% persistent storage for xini object and factory variables
persistent xinic;
persistent xinif;

% xinic                 internal object with subfields
%    .caller            who opened this file initially
%    .children          internal list of file-IDs that depend on this file
%    .convert           is the content converted (evaled/hxdouble'd)
%    .filename          holds the filename (or a temporary name)
%    .icontent          content of the file as struct representation
%    .parent            either empty or internal file-ID of parent file
%    .postfile          string of characters after terminating sequence
%    .prefile           string of characters before first section
%    .readwrite         single character option for write protection,
%                       whereas 'a' stands for all, 'c' for caller,
%                       'm' for memory-only, and 'n' for no write access
%
% xinif                 factory object with subfields
%    .fterm             ini-file terminator
%    .is_init           initialization complete flag
%    .lineout           configurable line break character sequence
%    .matrix            ID lookup table
%    .raction           read-only actions
%    .uaction           list of actions that don't need a valid ID
%    .waction           write-file actions

% check for class initialization -> do if necessary
if isempty(xinif) || ...
   ~xinif.is_init

    % configuration settings
    elf = neuroelf;
    xinif.asciiread = elf.asciiread;
    xinif.fterm = '[/EndOfXINI]';
    xinif.gluetostring = elf.gluetostring;
    xinif.is_init = false;
    xinif.lineout = char(10);
    xinif.matrix = [];

    % array of read-only actions for protection scheme
    xinif.raction = struct( ...
        'display',            true, ...
        'getcaller',          true, ...
        'getchildren',        true, ...
        'getcomplete',        true, ...
        'getfilename',        true, ...
        'getfoot',            true, ...
        'gethead',            true, ...
        'getid',              true, ...
        'getinisection',      true, ...
        'getinisetting',      true, ...
        'getinistring',       true, ...
        'getparents',         true, ...
        'getprotection',      true, ...
        'getsections',        true, ...
        'getsectionsettings', true, ...
        'issection',          true, ...
        'issetting',          true, ...
        'isvalid',            true  ...
    );

    % array of available actions for calls with object ID = 0
    xinif.uaction = struct( ...
        'display',         true, ...
        'isvalid',         true, ...
        'loadinifile',     true, ...
        'newinifile',      true, ...
        'parseinistring',  true, ...
        'releaseallfiles', true, ...
        'reloadallfiles',  true, ...
        'resetclass',      true  ...
    );

    % array of file-write actions for protection scheme
    xinif.waction = struct( ...
        'saveinifile',   true, ...
        'saveinifileas', true  ...
    );

    % initialize object structure
    xinic = struct(    ...
        'caller'   , '',     ...
        'children' , [],     ...
        'convert'  , 0,      ...
        'filename' , '',     ...
        'icontent' , struct, ...
        'parent'   , [],     ...
        'postfile' , '',     ...
        'prefile'  , '',     ...
        'readwrite', 'a',     ...
        'sectvals',  false   ... 
    );

    % set initialization successful to true
    xinif.is_init = true;

    % and get factory handle
    xinif.xhfactory = class(struct('L', 0), 'xini');
end
hf = xinif.xhfactory;

% if no arguments given return factory
if nargin < 1

    % initialize basic return object and return
    varargout{1} = hf;
    return;
end

% input argument check
if nargin < 3 && ...
   ~isa(varargin{1}, 'xini')

    % first is string
    if ischar(varargin{1}) && ...
       ~isempty(varargin{1})

        % but no filename...
        if exist(varargin{1}(:)','file') ~= 2

            % if it is neither 'convert', 'exact', or 'new'
            if ~any(strcmpi(varargin{1}(:)', {'convert', 'exact', 'new'}))
                error( ...
                    'xini:FileNotFound', ...
                    'The specified file wasn''t found.' ...
                );

            % assume NewIniFile constructor meaning
            else
                varargout{1} = ...
                    xini(hf, 'newinifile', lower(varargin{1}(:)'));
                return;
            end
        end

        % perform filename call
        varargout{1} = xini(hf, 'loadinifile', varargin{:});
        return;

    % if we get a numeric argument...
    elseif isa(varargin{1}, 'double') && ...
        numel(varargin{1}) == 1

        % is an object in list
        if any(xinif.matrix == varargin{1})
            varargout{1} = class(struct('L', varargin{1}), 'xini');
            return;

        % give an error otherwise
        else
            error( ...
                'xini:InvalidFileID', ...
                'The given ID wasn''t found in the lookup matrix.' ...
            );
        end
    end

% otherwise bark out if first argument in NOT an xini object
elseif ~isa(varargin{1}, 'xini') || ...
    numel(varargin{1}) ~= 1 || ...
    nargin < 2 || ...
   ~ischar(varargin{2})
    error( ...
        'xini:CallingConvention', ...
        'Illegal call to xini.' ...
    );
end
hIniFile = varargin{1};
action = lower(varargin{2}(:)');
varargout{1} = hIniFile;

% get ID and look up the matrix position of file...
cfid = hIniFile.L;
if cfid ~= 0
    cfid = find(xinif.matrix == cfid);
    if isempty(cfid)
        if ~strcmp(action, 'isvalid')
            error( ...
                'xini:InvalidFileID', ...
                'File lookup error. Possible programming error?' ...
            );
        else
            varargout{1} = false;
            return;
        end
    end
end

% disallow all actions that do need a "vadid" ID
if ~isfield(xinif.uaction, action) && ...
    cfid == 0
    error( ...
        'xini:InvalidAction', ...
        'The requested action needs a valid object to work.' ...
    );
end

% get current object
if cfid > 0
    cfile = xinic(cfid);

    % check protection first
    switch cfile.readwrite

        % only caller
        case {'c'}

            % test caller / action pair
            ncaller = extcaller(1);
            ocaller = cfile.caller;
            if ~isfield(xinif.raction, action) && ...
               ~strcmp(ncaller, ocaller)
                warning( ...
                    'xini:ProtectionViolation', ...
                    'The object has ''caller-only'' write access.' ...
                );
                return;
            end

        % file protection only
        case {'m'}

            % test for writefile action strings
            if isfield(xinif.waction, action)
                warning( ...
                    'xini:ProtectionViolation', ...
                    'The object has ''memory-only'' write access.' ...
                );
                return;
            end

        % no write access at all
        case {'n'}

            % test for read-only action string
            if ~isfield(xinif.raction, action)
                warning( ...
                    'xini:ProtectionViolation', ...
                    'The object has write access disabled.' ...
                );
                return;
            end
    end
end

%  switch over action
switch (action)

    % display (also used for writing!)
    case {'display'}

        % if sections can have values, use standard display on struct
        if cfile.sectvals
            display(cfile.icontent);
            varargout{1} = '';
            return;
        end

        % initialize output variable
        lineout = xinif.lineout;
        dispout = '';
        filepost = '';
        indent = '   ';
        equalchar = ' = ';

        % if we have a valid object ID (other than the factory)
        if cfid > 0

            % get conversion state
            tconvert = cfile.convert;
            issave = false;

            % if we're going to write the file to disk reset some vars!
            if nargin > 2 && ...
                ischar(varargin{3}) && ...
                strcmpi(varargin{3}(:)', 'saveinifile')
                dispout = cfile.prefile(:)';
                filepost = cfile.postfile(:)';
                indent = '';
                equalchar = '=';
                lineout = char([13, 10]);
                issave = true;

                % make sure to add correct footer if is invalid or empty
                if ~ischar(filepost) || ...
                    numel(filepost) < 1
                    filepost = [xinif.fterm lineout];
                end

                % do we have a correct footer now?
                hasend = strfind(filepost,xinif.fterm);
                if isempty(hasend)
                    filepost = [xinif.fterm lineout filepost];

                % remove any leading chars before good terminator
                elseif hasend(1) ~= 1
                    filepost(1:(hasend(1)-1)) = [];
                end
            end

            % correctly terminate header and footer
            if numel(dispout) > 0 && ...
                dispout(end) ~= char(10) && ...
                dispout(end) ~= char(13)
                dispout = [dispout lineout];
            end
            if numel(filepost) > 0 && ...
                filepost(end) ~= char(10) && ...
                filepost(end) ~= char(13)
                filepost = [filepost lineout];
            end

            % iterate over sections
            cfs = fieldnames(cfile.icontent);
            for sectcount = 1:numel(cfs)

                % initialize section output
                sdispout = '';
                extralarge = {};

                % get section name and content, and add name to dispout
                sectname = cfs{sectcount};
                csection = cfile.icontent.(cfs{sectcount});
                sdispout = [sdispout '[' sectname ']' lineout];

                % re-set non-structs
                if ~isstruct(csection)
                    dispout = [dispout '[' sectname ':] = ' any2ascii(csection) lineout lineout];
                    continue;
                end

                % iterate over settings in section
                namsettings = fieldnames(csection);
                for namecount = 1:numel(namsettings)

                    % get value for setting
                    fvalue = csection.(namsettings{namecount});

                    % do we have non character, invalid characters or not 1xN
                    if tconvert < 1 && ...
                       (~ischar(fvalue) || ...
                        any(fvalue < 32 | fvalue > 127) || ...
                        size(fvalue, 1) > 1 || ...
                        ndims(fvalue) > 2)
                        if issave
                            error( ...
                                'xini:BadSettingValue', ...
                                'Bad value for %s.%s (no conversion).', ...
                                sectname, ...
                                namsettings{namecount} ...
                            );
                        else
                            warning( ...
                                'xini:BadSettingValue', ...
                                'Bad value for %s.%s (no conversion).', ...
                                sectname, ...
                                namsettings{namecount} ...
                            );
                            continue;
                        end
                    end

                    % for non-struct value, try conversion
                    if ~isstruct(fvalue)
                        try
                            % for special character arrays ([...], {...})
                            if ischar(fvalue) && ...
                                tconvert && ...
                               ~isempty(fvalue) && ...
                               ((fvalue(1) == '[' && ...
                                 fvalue(end) == ']') || ...
                                (fvalue(1) == '{' && ...
                                 fvalue(end) == '}'))
                                fvalue = any2ascii(fvalue);

                            % for exact double conversion
                            elseif tconvert == 2 && ...
                                isa(fvalue, 'double') && ...
                                isreal(fvalue) && ...
                                any(floor(fvalue(:)) ~= fvalue(:))
                                fvalue = any2ascii(fvalue, 'exact');

                            % otherwise, for non chars
                            elseif ~ischar(fvalue)
                                fvalue = any2ascii(fvalue);
                            end

                        % if conversion fails, log error and continue
                        catch ne_eo;
                            if issave
                                error( ...
                                    'xini:ConversionFailed', ...
                                    'Bad value for %s.%s (%s).', ...
                                    sectname, ...
                                    namsettings{namecount}, ...
                                    ne_eo.message ...
                                );
                            else
                                warning( ...
                                    'xini:ConversionFailed', ...
                                    'Bad value for %s.%s (%s).', ...
                                    sectname, ...
                                    namsettings{namecount}, ...
                                    ne_eo.message ...
                                );
                                continue;
                            end
                        end

                    % in case of a struct, insert into special array
                    else

                        extralarge{end + 1} = namsettings{namecount};

                        % for display, omit the cell2struct line
                        if ~issave
                            continue;
                        end

                        % print cell2struct line for other structs
                        fvalfields = fieldnames(fvalue);
                        fvalue = sprintf('[cell2struct(cell(%.0f,%.0f,%.0f),%s,3)]', ...
                                     size(fvalue, 1), ...
                                     size(fvalue, 2), ...
                                     numel(fvalfields), ...
                                     any2ascii(fvalfields));
                    end

                    % put Fieldname = <Value> line into soutput
                    sdispout = [sdispout indent ...
                                namsettings{namecount} equalchar fvalue lineout];
                end

                % do we have any structs in section?
                if ~isempty(extralarge)

                    % iterate over structs
                    for namecount = 1:length(extralarge)

                        % try conversion and log error if failure
                        try
                            sdispout = [sdispout lineout ...
                                     i_dispxl( ...
                                         [sectname '.' extralarge{namecount}], ...
                                         csection.(extralarge{namecount}),     ...
                                         indent,    ...
                                         lineout,   ...
                                         equalchar, ...
                                         tconvert   ...
                                     )];
                        catch ne_eo;
                            if issave
                                error( ...
                                    'xini:ConversionFailed', ...
                                    'Bad struct for %s.%s (%s).', ...
                                    sectname, ...
                                    extralarge{namecount}, ...
                                    ne_eo.message ...
                                );
                            else
                                warning( ...
                                    'xini:ConversionFailed', ...
                                    'Bad struct for %s.%s (%s).', ...
                                    sectname, ...
                                    extralarge{namecount}, ...
                                    ne_eo.message ...
                                );
                                continue;
                            end
                        end
                    end
                end

                % insert section output into global out
                dispout = [dispout sdispout lineout];
            end

        % for base class we show a list of files with IDs
        else

            % print nice header
            dispout = ['Currently active ini files (' ...
                num2str(length(xinif.matrix)) '):' lineout lineout];
            dispout = [dispout '-------------------------------------------------------------------------' lineout];
            dispout = [dispout 'OBJECT_FID (PARENT_FID)  C  P  Caller                          Filename  ' lineout];
            dispout = [dispout '-------------------------------------------------------------------------' lineout];

            % iterate over all loaded IniFiles
            for fiter = 1:numel(xinif.matrix)
                fifid = xinif.matrix(fiter);

                % get conversion and protection states, parent ID, and caller
                if xinic(fiter).convert
                    isconv = 'y';
                else
                    isconv = 'n';
                end
                protstat = xinic(fiter).readwrite;
                if isempty(xinic(fiter).parent)
                    p_id = '(    --    )';
                else
                    p_id = sprintf('(%10.0f)', xinic(fiter).parent(1));
                end
                fcaller = xinic(fiter).caller;
                if length(fcaller) > 30
                    fcaller = ['...' fcaller(end-26:end)];
                end

                % show ID, matrix position, and settings
                % file contents can then be obtained via xini(ID)
                dispout = [dispout ...
                    sprintf('%10.0f %12s  %s  %s  %-30s  %s', ...
                    fifid, p_id, isconv, protstat, fcaller, ...
                    xinic(fiter).filename) lineout];
            end
        end

        % replace varargout{1} and display if should be displayed!
        varargout{1} = [dispout filepost];
        if nargout == 0
            disp(dispout);
        end


    % get name of caller
    case {'getcaller'}
        varargout{1} = cfile.caller;


    % get list of children IDs
    case {'getchildren'}
        varargout{1} = cfile.children;


    % give complete IniFile contents as struct
    case {'getcomplete'}
        varargout{1} = cfile.icontent;


    % give filename of object
    case {'getfilename'}

        % set file name
        varargout{1} = cfile.filename;


    % get file footer
    case {'getfoot'}

        % get footer
        varargout{1} = strrep(cfile.postfile, xinif.fterm, '');

        % remove extra linefeeds...
        while numel(varargout{1}) > 0 && ...
           (varargout{1}(1) == char(10) || ...
            varargout{1}(1) == char(13))
            varargout{1}(1) = [];
        end


    % get file header
    case {'gethead'}
        varargout{1} = cfile.prefile;


    % give id of object
    case {'getid'}
        varargout{1} = hIniFile.L;


    % give filename of object
    case {'getinistring'}

        % set file name
        varargout{1} = xini(hIniFile, 'display', 'saveinifile');


    % get list of parent IDs
    case {'getparents'}

        % initialize parent ID list
        loparents = [];
        tfid = cfid;

        % iterate until no more parent found
        while ~isempty(xinic(tfid).parent)
            loparents(end + 1) = xinic(tfid).parent;
            tfid = find(xinif.matrix == loparents(end));

            % bail out if parent disappeared
            if isempty(tfid)
                error( ...
                    'xini:InvalidFileID', ...
                    'Parent object disappeared, memory glitch?' ...
                );
            end
        end

        % return list
        varargout{1} = loparents;


    % give id of object
    case {'getprotection'}
        varargout{1} = cfile.readwrite;


    % give one section of IniFile as struct
    case {'getinisection'}

        % check for section field in call
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3})
            error( ...
                'xini:CallingConvention', ...
                'GetIniSection needs a CHAR section name in call.' ...
            );
        end

        % if section doesn't exist
        if ~isfield(cfile.icontent, varargin{3}(:)')
            error( ...
                'xini:InvalidSection', ...
                'Invalid section in call: %s in file %s.', ...
                varargin{3}(:)', ...
                cfile.filename ...
            );
        end
        varargout{1} = cfile.icontent.(varargin{3}(:)');


    % return sections names
    case {'getsections'}
        varargout{1} = fieldnames(cfile.icontent);


    % return sections names
    case {'getsectvals'}
        varargout{1} = cfile.sectvals;


    % return section's setting names
    case {'getsectionsettings'}

        % check for section field in call
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3})
            error( ...
                'xini:CallingConvention', ...
                'IniSectionSettings needs a CHAR section name in call.' ...
            );
        end

        % if section doesn't exist, return empty list
        if ~isfield(cfile.icontent, varargin{3}(:)')
            error( ...
                'xini:InvalidSection', ...
                'The given section (%s) is not part of this file.', ...
                varargin{3}(:)' ...
            );
        end

        % get section and if it's struct get fieldnames
        csection = cfile.icontent.(varargin{3}(:)');
        if isstruct(csection)
            varargout{1} = fieldnames(csection);
        else
            varargout{1} = {'Value'};
        end


    % give one setting in section of IniFile
    case {'getinisetting'}

        % check for section and setting fields in call
        if nargin < 4 || ...
          ~ischar(varargin{3}) || ...
          ~ischar(varargin{4}) || ...
           isempty(varargin{3}) || ...
           isempty(varargin{4})
            error( ...
                'xini:CallingConvention', ...
                'GetIniSetting needs two CHAR arguments in call.' ...
            );
        end

        % if section doesn't exist and no parent is available
        if ~isfield(cfile.icontent, varargin{3}(:)') && ...
            isempty(cfile.parent)
            error( ...
                'xini:InvalidSection', ...
                'Invalid section in call: %s in file %s.', ...
                varargin{3}(:)', ...
                xinic(cfid).filename ...
            );
        end

        % try to retrieve section
        if isfield(cfile.icontent, varargin{3}(:)')
            csection = cfile.icontent.(varargin{3}(:)');
            if ~isstruct(csection)
                varargout{1} = csection;
                return;
            end
        else
            csection = struct;
        end

        % if we didn't succeed, check parents
        if ~isfield(csection, varargin{4}(:)')
            snotfound = true;
            tfid = cfid;

            % follow parents until found or no more parent
            while snotfound && ...
               ~isempty(cfile.parent)

                % lookup parent ID
                tfid = find(xinif.matrix == xinic(tfid).parent(1));

                % if not found, bail out
                if isempty(tfid)
                    error( ...
                        'xini:InvalidFileId', ...
                        'Parent object disappeared, memory glitch?' ...
                    );
                end

                % try to retrieve section from parent
                if isfield(xinic(tfid).icontent, varargin{3}(:)')
                    tsection = xinic(tfid).icontent.(varargin{3}(:)');
                    if ~isstruct(tsection)
                        if ~iscell(tsection)
                            tsection = struct('Value', tsection);
                        else
                            tsection = struct('Value', {tsection});
                        end
                    end
                else
                    tsection = struct;
                end

                % if setting exists, leave loop
                if isfield(tsection, varargin{4}(:)')
                    sfound = tsection.(varargin{4}(:)');
                    snotfound = false;
                end
            end

            % setting still not found?
            if snotfound
                error( ...
                    'xini:InvalidSetting', ...
                    'Invalid setting in call: %s.%s in file %s.', ...
                    varargin{3}(:)', ...
                    varargin{4}(:)', ...
                    cfile.filename ...
                );
            end

            % found -> return
            varargout{1} = sfound;

        % direct hit
        else
            varargout{1} = csection.(varargin{4}(:)');
        end


    % is given argument an existing section ?
    case {'issection'}

        % check for section field in call
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3})
            error( ...
                'xini:CallingConvention', ...
                'IsSection needs a CHAR section name in call.' ...
            );
        end

        % get info...
        varargout{1} = isfield(cfile.icontent, varargin{3}(:)');


    % are given arguments an existing section/setting pair ?
    case {'issetting'}

        % check for section and setting fields in call
        if nargin < 4 || ...
           ~ischar(varargin{3}) || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{3}) || ...
            isempty(varargin{4})
            error( ...
                'xini:CallingConvention', ...
                'IsSetting needs two CHAR arguments in call.' ...
            );
        end

        % init output variable
        varargout{1} = false;

        % only further checks if section exists
        if isfield(cfile.icontent, varargin{3}(:)')

            % get test section
            tsection = cfile.icontent.(varargin{3}(:)');
            if ~isstruct(tsection)
                varargout{1} = true;
                return;
            end

            % reset output variable
            varargout{1} = isfield(tsection, varargin{4}(:)');

            % if exists or no parent should be inspected then return
            if varargout{1}
                return;
            end
        end

        % if we have a parent and should inspect them
        if ~isempty(cfile.parent) && ...
           (nargin < 5 || ...
            ~islogical(varargin{5}) || ...
            isempty(varargin{5}) || ...
            ~varargin{5}(1))
            sfound = false;
            tfid = cfid;

            % iterate until found or more parent
            while ~sfound && ...
               ~isempty(xinic(tfid).parent)

                % lookup parent
                tfid = find(xinif.matrix == xinic(tfid).parent(1));

                % bail out if gone
                if isempty(tfid)
                    error( ...
                        'xini:InvalidFileId', ...
                        'Parent object disappeared, memory glitch?' ...
                    );
                end

                % check for field
                if isfield(xinic(tfid).icontent, varargin{3}(:)')
                    tsection = xinic(tfid).icontent.(varargin{3}(:)');
                    if ~isstruct(tsection)
                        if ~iscell(tsection)
                            tsection = struct('Value', tsection);
                        else
                            tsection = struct('Value', {tsection});
                        end
                    end
                else
                    tsection = struct;
                end

                if isfield(tsection, varargin{4}(:)')
                    sfound = true;
                end
            end
            varargout{1} = sfound;
        end


    % get validity flag
    case {'isvalid'}
        varargout{1} = true;


    % read IniFile from disk
    case {'loadinifile'}

        % do we have a valid filename
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3}) || ...
            exist(varargin{3}(:)', 'file') ~= 2
            error( ...
                'xini:CallingConvention', ...
                'LoadIniFile needs a valid filename in call.' ...
            );
        end

        % only works for the factory object
        if cfid ~= 0
            error( ...
                'xini:CallingConvention', ...
                'LoadIniFile only works for the factory object.' ...
            );
        end

        % get filename, and make sure to have absolute path
        [isabs{1:2}] = isabsolute(varargin{3}(:)');
        filename = isabs{2};

        % try to read file
        try
            filecont = xinif.asciiread(filename);
        catch ne_eo;
            error( ...
                'xini:FileNotReadable', ...
                'The specified file is not readable: %s.', ...
                ne_eo.message ...
            );
        end

        % conversion activated?
        iconvert = 0;
        if nargin > 3
            sconvert = varargin{4};
            if isnumeric(sconvert) && ...
                numel(sconvert) == 1
                if sconvert == 1
                    iconvert = 1;
                elseif sconvert == 2
                    iconvert = 2;
                end
            else
                iconvert = 1;
            end
        end

        % check for current files in matrix (re-open)
        for ofc = 1:length(xinif.matrix)

            % do we already have this file opened?
            if strcmp(filename, xinic(ofc).filename) || ...
               (ispc && ...
                strcmpi(filename, xinic(ofc).filename))

                % then initialize output object
                varargout{1} = ...
                    class(struct('L', xinif.matrix(ofc)), 'xini');

                % conversion on top? -> reconvert now
                if iconvert > 0 && ...
                    xinic(ofc).convert == 0

                    % get current content
                    ccontent = xinic(ofc).icontent;
                    sectnames = fieldnames(ccontent);

                    % try to convert any strings with [...] or {...} content
                    try
                        % iterate over sections
                        for csect = 1:length(sectnames)

                            % get setting names
                            tsection = ccontent.(sectnames{csect});
                            if ~isstruct(tsection)
                                if ~iscell(tsection)
                                    tsection = struct('Value', tsection);
                                else
                                    tsection = struct('Value', {tsection});
                                end
                            end
                            ssettings = fieldnames(tsection);

                            % iterate over settings
                            for cset = 1:length(ssettings)

                                % get value
                                svalue = tsection.(ssettings{cset});

                                % if convertable
                                if ischar(svalue) && ...
                                   ~isempty(svalue) && ...
                                   ((svalue(1) == '[' && ...
                                     svalue(end) == ']') || ...
                                    (svalue(1) == '{' && ...
                                     svalue(end) == '}'))
                                    try
                                        tsection.(ssettings{cset}) = ...
                                            eval(strrep(svalue, ...
                                                 '$SELF.', 'ccontent.'));
                                    catch ne_eo;
                                        neuroelf_lasterr(ne_eo);
                                        error('CONVERSION_FAILED');
                                    end
                                end
                            end
                            ccontent.(sectnames{csect}) = tsection;
                        end

                    % conversion failed, bail out
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        error( ...
                            'xini:ConversionFailed', ...
                            'Couldn''t convert already open file %s.', ...
                            xinic(ofc).filename ...
                        );
                    end

                    % now set new conversion mode and content
                    xinic(ofc).convert = iconvert;
                    xinic(ofc).icontent = ccontent;
                end

                % do return
                return;
            end
        end

        % file not already open -> get unique object fid
        cfid = floor(2 ^ 31 * rand(1) + 1);
        while ~isempty(xinif.matrix) && ...
           any(xinif.matrix == cfid)
            cfid = floor(2 ^ 31 * rand(1) + 1);
        end

        % read file contents and cut trailing lines
        [icontent, prefile, postfile, newerrors] = i_parse( ...
            filecont, xinif.fterm, iconvert);

        % errors during parse
        if ~isempty(newerrors)
            error( ...
                'xini:IniParseError', ...
                'Parsing errors occurred:\n\n%s', ...
                xinif.gluetostring(newerrors, char(10)) ...
            );
        end

        % add new file to matrix
        xinif.matrix(end + 1) = cfid;
        ifid = length(xinif.matrix);
        xinic(ifid).caller = extcaller(1);
        xinic(ifid).children = [];
        xinic(ifid).convert = iconvert;
        xinic(ifid).filename = filename;
        xinic(ifid).icontent = icontent;
        xinic(ifid).parent = [];
        xinic(ifid).postfile = postfile;
        xinic(ifid).prefile = prefile;
        xinic(ifid).readwrite = 'a';

        % set output object
        varargout{1} = class(struct('L', cfid),'xini');


    % create IniFile object from scratch
    case {'newinifile'}

        % only works for factory object
        if cfid ~= 0
            error( ...
                'xini:CallingConvention', ...
                'NewIniFile only works for the factory object.' ...
            );
        end

        % get unique object fid
        cfid = floor(2 ^ 31 * rand(1) + 1);
        while ~isempty(xinif.matrix) && ...
            any(xinif.matrix == cfid)
            cfid = floor(2 ^ 31 * rand(1) + 1);
        end

        % update matrix and init object fields
        xinif.matrix(end + 1) = cfid;
        ifid = numel(xinif.matrix);
        xinic(ifid).caller = extcaller(1);
        xinic(ifid).children = [];
        xinic(ifid).convert = 0;
        xinic(ifid).filename = '';
        xinic(ifid).icontent = struct;
        xinic(ifid).parent = [];
        xinic(ifid).postfile = [xinif.fterm xinif.lineout];
        xinic(ifid).prefile = '';
        xinic(ifid).readwrite = 'a';

        % set output object values
        varargout{1} = class(struct('L', cfid), 'xini');

        % conversion?
        if nargin > 2
            sconvert = varargin{3}(:)';

            % numeric convert -> [0..2]
            if isnumeric(sconvert) && ...
               ~numel(sconvert) == 1 && ...
               ~isnan(sconvert) && ...
               ~isinf(sconvert)
                xinic(ifid).convert = min(2, max(0, floor(sconvert)));

            % for 'exact' set to 2
            elseif ischar(sconvert) && ...
                strcmpi(sconvert, 'exact')
                xinic(ifid).convert = 2;

            % else to 1
            else
                xinic(ifid).convert = 1;
            end
        end


    % parse IniFile from string
    case {'parseinistring'}

        % do we have a valid ini-string?
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3})
            error( ...
                'xini:CallingConvention', ...
                'ParseIniString needs a non-empty inistring to parse.' ...
            );
        end

        % only works for the factory object
        if cfid ~= 0
            error( ...
                'xini:CallingConvention', ...
                'NewIniFile only works for the factory object.' ...
            );
        end

        % conversion activated?
        iconvert = 0;
        if nargin > 3
            sconvert = varargin{4}(:)';
            if isnumeric(sconvert) && ...
               numel(sconvert) == 1
                if sconvert(1) == 1
                    iconvert = 1;
                elseif sconvert(1) == 2
                    iconvert = 2;
                end
            elseif ischar(sconvert) && ...
                strcmpi(sconvert, 'exact')
                iconvert = 2;
            else
                iconvert = 1;
            end
        end

        % get unique object fid
        cfid = floor(2 ^ 31 * rand(1) + 1);
        while ~isempty(xinif.matrix) && ...
            any(xinif.matrix == cfid)
            cfid = floor(2 ^ 31 * rand(1) + 1);
        end

        % parse file contents
        [icontent, prefile, postfile, newerrors] = i_parse( ...
            varargin{2}(:)', xinif.fterm, iconvert);

        % errors during parse
        if ~isempty(newerrors)
            error( ...
                'xini:IniParseError', ...
                'Parsing errors occurred:\n\n%s', ...
                xinif.gluetostring(newerrors, char(10)) ...
            );
        end

        % add new file to matrix
        xinif.matrix(end + 1) = cfid;
        ifid = length(xinif.matrix);
        xinic(ifid).caller = extcaller(1);
        xinic(ifid).children = [];
        xinic(ifid).convert = iconvert;
        xinic(ifid).filename = '';
        xinic(ifid).icontent = icontent;
        xinic(ifid).parent = [];
        xinic(ifid).postfile = postfile;
        xinic(ifid).prefile = prefile;
        xinic(ifid).readwrite = 'a';

        % init output object
        varargout{1} = class(struct('L', cfid), 'xini');


    % user request to release IniFile
    case {'release'}

        % test for right caller or forced mode
        ocaller = cfile.caller;
        ncaller = extcaller(1);
        if ~strcmp(ocaller, ncaller) && ...
            (nargin < 3 || ...
            ~ischar(varargin{3}) || ...
             isempty(varargin{3}) || ...
             lower(varargin{3}(1)) ~= 'f')
            warning( ...
                'xini:CantReleaseObject', ...
                'ReleaseIniFile not forced or issued by caller.' ...
            );
            return;
        end

        % we're going to release this file, so make good output variable
        varargout{1} = xinif.xhfactory;

        % inform parent
        pfid = cfile.parent;
        if ~isempty(pfid)

            % lookup parent object in matrix
            pfmpos = find(xinif.matrix == pfid(1));

            % bail out if parent no longer available
            if isempty(pfmpos)
                error( ...
                    'xini:InvalidFileID', ...
                    'Parent object disappeared, memory glitch?' ...
                );
            end

            % remove self from parent's children list
            ifid = hIniFile.L;
            xinic(pfmpos).children = ...
                xinic(pfmpos).children(xinic(pfmpos).children ~= ifid);
        end

        % inform children
        chid = cfile.children;
        if ~isempty(chid)

            % iterate over children
            for schid = chid

                % lookup child object in matrix
                cfmpos = find(xinif.matrix == schid);

                % bail out if child no longer available
                if isempty(cfmpos)
                    error( ...
                        'xini:InvalidFileID', ...
                        'Child object disappeared, memory glitch?' ...
                    );
                end

                % unset parent object in child
                xinic(cfmpos).parent = [];
            end
        end

        % finally remove object itself
        xinic(cfid) = [];
        xinif.matrix(cfid) = [];


    % user request to release all files from memory
    case {'releaseallfiles'}

        % only works for the factory object
        if cfid ~= 0
            error( ...
                'xini:CallingConvention', ...
                'ReleaseAllFiles only works for the factory object.' ...
            );
        end

        % iterate over all open files in matrix
        xinif.matrix = [];
        xinic(:) = [];
        varargout{1} = xinif.xhfactory;


    % reload all files in memory from disk
    case {'reloadallfiles'}

        % only works for the factory object
        if cfid ~= 0
            error( ...
                'xini:CallingConvention', ...
                'ReloadAllFiles only works for the factory object.' ...
            );
        end

        % init output object
        varargout{1} = xinif.xhfactory;

        % iterate over all open files in matrix
        for rfid = xinif.matrix
            rfob.L = rfid;
            xini(class(rfob, 'xini'), 'reloadinifile');
        end


    % reload ini file contents
    case {'reloadinifile'}

        % get filename and convert state
        cfilename = cfile.filename;
        cconvert = cfile.convert;

        % check for file status
        if isempty(cfilename)
            return;
        end

        % mark filename, so we can reload it!
        xinic(cfid).filename = '!BEING_RELOADED!';

        % try to load file as new
        try
            if cconvert > 1
                reloaded = xini(cfilename, 'exact');
            elseif cconvert > 0
                reloaded = xini(cfilename, 'convert');
            else
                reloaded = xini(cfilename);
            end
        catch ne_eo;
            xinic(cfid).filename = cfilename;
            warning( ...
                'xini:ReloadFailed', ...
                'Reload of file %s failed (%s).', ...
                cfilename, ...
                ne_eo.message ...
            );
            return;
        end

        % restore filename
        xinic(cfid).filename = cfilename;

        % locate new object
        rfid = find(xinif.matrix == reloaded.L(1));

        % bail out if not found
        if isempty(rfid)
            error( ...
                'xini:InvalidFileID', ...
                'Re-read file %s disappeared from memory. Glitch?', ...
                 cfilename ...
            );
        end

        % set new file content
        xinic(cfid).icontent = xinic(rfid(1)).icontent;
        xinic(cfid).postfile = xinic(rfid(1)).postfile;
        xinic(cfid).prefile = xinic(rfid(1)).prefile;

        % clear re-loaded in matrix
        xinic(rfid(1)) = [];
        xinif.matrix(rfid(1)) = [];


    % delete contents (!) of IniFile in memory
    case {'removeinicomplete'}

        % clear file contents
        xinic(cfid).icontent = [];


    % remove one section from memory
    case {'removeinisection'}

        % check for section field in call
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3})
            error( ...
                'xini:CallingConvention', ...
                'RemoveIniSection needs a CHAR section name in call.' ...
            );
        end

        % test if section exists
        if ~isstruct(cfile.icontent) || ...
           ~isfield(cfile.icontent, varargin{3}(:)')
            return;
        end

        % remove section and get new names
        xinic(cfid).icontent = ...
            rmfield(xinic(cfid).icontent, varargin{3}(:)');


    % remove one setting from an IniFile section
    case {'removeinisetting'}

        % check for section and setting fields in call
        if nargin < 4 || ...
           ~ischar(varargin{3}) || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{3}) || ...
            isempty(varargin{4})
            error( ...
                'xini:CallingConvention', ...
                'RemoveIniSetting needs two CHAR arguments name in call.' ...
            );
        end

        % test if section and setting exist
        if ~isstruct(cfile.icontent) || ...
           ~isfield(cfile.icontent, varargin{3}(:)')
            return;
        end
        tsection = cfile.icontent.(varargin{3}(:)');
        if ~isstruct(tsection)
            xinic(cfid).icontent = struct;
            return;
        end
        if ~isfield(tsection, varargin{4}(:)')
            return;
        end

        % remove setting
        xinic(cfid).icontent.(varargin{3}(:)') = ...
            rmfield(tsection, varargin{4}(:)');


    % set file footer
    case {'resetclass'}

        % only works for the factory object
        if cfid ~= 0
            error( ...
                'xini:CallingConvention', ...
                'NewIniFile only works for the xini factory object.' ...
            );
        end

        % clear persistent memory and re-run init code
        xinif = [];
        xinic = [];
        xini;

        % return new factory as output
        varargout{1} = xinif.xhfactory;


    % write IniFile back to disk
    case {'save'}

        % test whether we have a good filename
        if isempty(cfile.filename)
            error( ...
                'xini:MemoryOnlyID', ...
                'Object is a new file, use SaveAs(FILENAME).' ...
            );
        end

        % test if conversion is requested
        if cfile.convert == 0
            tconvert = 0;
            if nargin > 2 && ...
               ~isempty(varargin{3})
                tconvert = varargin{3}(:)';
                if isnumeric(tconvert) && ...
                    numel(tconvert) == 1 && ...
                   ~isnan(tconvert) && ...
                   ~isinf(tconvert)
                    tconvert = min(2, max(0, floor(tconvert)));
                elseif ischar(tconvert)
                    if strcmpi(tconvert, 'exact')
                        tconvert = 2;
                    else
                        tconvert = 1;
                    end
                else
                    tconvert = 1;
                end
            end
            if tconvert > 0
                xinic(cfid).convert = tconvert;
            end
        end

        % open file
        rfid = fopen(xinic(cfid).filename, 'w');
        if rfid < 1
            error( ...
                'xini:FileNotWritable', ...
                'Couldn''t write ini-file %s.', ...
                xinic(cfid).filename ...
            );
        end

        % write file contents
        frewind(rfid);
        fwrite(rfid, xini(hIniFile, 'display', 'saveinifile'), 'uchar');
        fclose(rfid);


    % write IniFile under a new name
    case {'saveas'}

        % test for asfilename field
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3})
            error( ...
                'xini:CallingConvention', ...
                'SaveAs needs a non-empty CHAR filename.' ...
            );
        end

        % try to write file
        asfile = varargin{3}(:)';
        rfid = fopen(asfile, 'w');
        if rfid < 1
            error( ...
                'xini:FileNotWritable', ...
                'Couldn''t write xini as %s.', ...
                asfile ...
            );
        end
        fclose(rfid);

        % set new filename and write file
        xinic(cfid).filename = asfile;
        xini(hIniFile, 'save', varargin{4:end});


    % set file footer
    case {'setfoot'}

        % test for foot
        if nargin < 3 || ...
           ~ischar(varargin{3})
            error( ...
                'xini:CallingConvention', ...
                'SetFileFoot needs a CHAR foot option.' ...
            );
        end

        % set footer
        xinic(cfid).postfile = varargin{3}(:)';


    % set file header
    case {'sethead'}

        % test for head
        if nargin < 3 || ...
           ~ischar(varargin{3})
            error( ...
                'xini:CallingConvention', ...
                'SetFileHead needs a CHAR head option.' ...
            );
        end

        % set header
        newhead = varargin{3}(:)';
        newhok = find(newhead == '[');
        while ~isempty(newhok)
            newhead = newhead(newhok(1) + 1:end);
            newhok = find(newhead == ']');
            if ~isempty(newhok)
                newhead = newhead(newhok(1) + 1:end);
            end
            newhok = find(newhead == '[');
        end
        xinic(cfid).prefile = newhead;


    % set a complete IniFile structure
    case {'setcomplete'}

        % test complete field
        if nargin < 3 || ...
           ~isstruct(varargin{3}) || ...
            numel(varargin{3}) ~= 1
            error( ...
                'xini:CallingConvention', ...
                'SetIniComplete needs a 1x1 STRUCT argument.' ...
            );
        end

        % get shortcut
        cmpl = varargin{3};
        
        % shortcut
        if cfile.sectvals
            xinic(cfid).icontent = cmpl;
            return;
        end

        % iterate over sections
        sects = fieldnames(cmpl);
        for sc = 1:numel(sects)
            csection = cmpl.(sects{sc});
            if ~isstruct(csection)
                if ~iscell(csection)
                    csection = struct('Value', csection);
                else
                    csection = struct('Value', {csection});
                end
                cmpl.(sects{sc}) = csection;
                continue;
            elseif numel(csection) ~= 1
                error( ...
                    'xini:CallingConvention', ...
                    'Fields in complete struct must be 1x1 structs.' ...
                );
            end

            % if we need to check content, iterate over fields
            if cfile.convert == 0
                sfields = fieldnames(csection);
                for fc = 1:length(sfields)
                    ftest = csection.(sfields{fc});

                    % reject every non-char, illegal chars and bad dim'ed chars
                    if ~ischar(ftest) || ...
                        any(ftest(:) < 32 | ftest(:) > 127) || ...
                       (~isempty(ftest) && ...
                        numel(ftest) ~= size(ftest, 2))
                        error( ...
                            'xini:InvalidSetting', ...
                            'Illegal value for %s.%s (no conversion).', ...
                            sects{sc}, ...
                            sfields{fc} ...
                        );
                    end
                end
            end
        end

        % set new contents
        xinic(cfid).icontent = cmpl;


    % set an entire section
    case {'setinisection'}

        % test for section and inisection fields
        if nargin < 4 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3}) || ...
           ~isrealvarname(varargin{3}(:)') || ...
           (~cfile.sectvals && ...
            (~isstruct(varargin{4}) || ...
            numel(varargin{4}) ~= 1))
            error( ...
                'xini:CallingConvention', ...
                'SetIniSection needs CHAR section and a 1x1 struct.' ...
            );
        end

        % if no conversion iterate over fields
        csection = varargin{4};
        if cfile.convert == 0
            sfields = fieldnames(csection);
            for fc = 1:length(sfields)
                ftest = csection.(sfields{fc});

                % reject every non-char, illegal chars and bad dim'ed chars
                if ~ischar(ftest) || ...
                    any(ftest(:) < 32 | ftest(:) > 127) || ...
                   (~isempty(ftest) && ...
                    numel(ftest) ~= size(ftest, 2))
                    error( ...
                        'xini:InvalidSetting', ...
                        'Illegal value for %s.%s (no conversion).', ...
                        varargin{3}(:)', ...
                        sfields{fc} ...
                    );
                end
            end
        end

        % set section in matrix
        xinic(cfid).icontent.(varargin{3}(:)') = csection;


    % add/change an IniFile setting
    case {'setinisetting'}

        % test for section, setting, and value fields
        if nargin < 5 || ...
           ~ischar(varargin{3}) || ...
            isempty(varargin{3}) || ...
           ~isrealvarname(varargin{3}(:)') || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{4}) || ...
           ~isrealvarname(varargin{4}(:)')
            error( ...
                'xini:CallingConvention', ...
                'SetIniSetting needs a CHAR section and setting option.' ...
            );
        end

        % if no conversion, reject invalid content
        fval = varargin{5};
        if cfile.convert == 0 && ...
           (~ischar(fval) || ...
             any(fval(:) < 32 | fval(:) > 127) || ...
             (~isempty(fval) && ...
              numel(fval) ~= size(fval, 2)))
            error( ...
                'xini:InvalidSetting', ...
                'Illegal value for %s.%s (no conversion).', ...
                varargin{3}(:)', ...
                varargin{4}(:)' ...
            );
        end

        % add new section first?
        if ~isfield(cfile.icontent, varargin{3}(:)')
            xinic(cfid).icontent.(varargin{3}(:)') = struct;
        end

        % set value
        xinic(cfid).icontent.(varargin{3}(:)').(varargin{4}(:)') = fval;


    % set IniFile's parent object id
    case {'setparent'}

        % test for parent field
        if nargin < 3 || ...
            numel(varargin{3}) > 1 || ...
           ~any(strcmpi(class(varargin{3}), {'double', 'xini'}))
            error( ...
                'xini:CallingConvention', ...
                'SetParent needs a valid object as parent option.' ...
            );
        end

        % get parent object's ID and locate matrix position
        try
            % try to get requested FileID
            pobj = varargin{3};
            if isa(pobj, 'xini') && ...
                numel(pobj) == 1
                pobj = pobj.L;
            elseif isempty(pobj) == 0
                pobj = [];
            end

            % locate ID in matrix
            if ~isempty(pobj)
                pobj = pobj(1);
                pobjid = find(xinif.matrix == pobj);

                % handle disappeared parent
                if numel(pobjid) ~= 1
                    error('INVADID_FILE_ID');
                end
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            error( ...
                'xini:InvalidFileID', ...
                'The requested file ID doesn''t exist. Glitch?' ...
            );
        end

        % parent object already set
        if ~isempty(cfile.parent)
            ofid = find(xinif.matrix == cfile.parent(1));

            % unset parent object reference
            xinic(cfid).parent = [];

            % valid File-ID
            if ~isempty(ofid)

                % clear parent's child reference to this
                ifid = hIniFile.L;
                xinic(ofid).children(xinic(ofid).children == ifid) = [];

            % invalid File-ID
            else
                warning( ...
                    'xini:InvalidFileID', ...
                    'The former parent object is gone. Glitch?' ...
                );
            end
        end

        % no more parent...
        if isempty(pobj)
            return;
        end

        % self being parent is illegal
        if pobjid == cfid
            error( ...
                'xini:InvalidParentID', ...
                'The object can''t be it''s own parent.' ...
            );
        end

        % try to detect loop
        pparents = xini(class(struct('L', pobj), 'xini'), ...
            'getparents');
        if any(pparents == hIniFile.L)
            error( ...
                'xini:InvalidParentID', ...
                'Parent loop detected.' ...
            );
        end

        % set parent references
        xinic(cfid).parent = pobj;
        xinic(pobjid).children(end + 1) = hIniFile.L;


    % set protection scheme
    case {'setprotection'}

        % check for protection field in call
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
            numel(varargin{3}) ~= 1 || ...
           ~any('acmn' == lower(varargin{3}))
            error( ...
                'xini:CallingConvention', ...
                'SetProtection needs a valid 1x1 CHAR protection tag.' ...
            );
        end

        % re-set protection scheme
        xinic(cfid).readwrite = lower(varargin{3});


    % set section as values flag
    case {'setsectvals'}

        % check for protection field in call
        if nargin < 3 || ...
           ~islogical(varargin{3}) || ...
            numel(varargin{3}) ~= 1
            error( ...
                'xini:CallingConvention', ...
                'SetSectVals needs a valid 1x1 LOGICAL flag.' ...
            );
        end

        % re-set protection scheme
        xinic(cfid).sectvals = varargin{3};


    % otherwise log invalid action
    otherwise

        % error or not
        error( ...
            'xini:InvalidMethod', ...
            'Method %s unknown.', ...
            action(:)' ...
        );
end


% % % internal functions


function xldispout = i_dispxl(xlname, xlcont, xlind, xllo, xleq, xlconv)

    % test what size the struct has (indexing needed)
    xlsize = size(xlcont);
    if all(xlsize == 1)

        % start with section name and line feed
        xldispout = ['[' xlname ']' xllo];

        % iterate over fieldnames
        xlxl = {};
        xlfields = fieldnames(xlcont);
        for numf = 1:length(xlfields)

            % get and check content
            xlvalue = xlcont.(xlfields{numf});

            % simple string
            if ischar(xlvalue) && ...
               ~any(xlvalue < 32 | xlvalue > 127) && ...
                ndims(xlvalue) < 3 && ...
                size(xlvalue, 1) < 2

                % if conversion activated and [...] or {...}
                if xlconv && ...
                   ~isempty(xlvalue) && ...
                   ((xlvalue(1) == '[' && ...
                     xlvalue(end) == ']') || ...
                    (xlvalue(1) == '{' && ...
                     xlvalue(end) == '}'))
                    xlvalue = any2ascii(xlvalue);
                end

                % add to output
                xldispout = [xldispout xlind xlfields{numf} xleq xlvalue xllo];

            % complex variable (conv is always at least 1!)
            else

                % add structs to list
                if isstruct(xlvalue)
                    xlxl{end + 1} = xlfields{numf};

                    % print cell2struct line
                    xlvalfields = fieldnames(xlvalue);
                    xldispout = [xldispout xlind xlfields{numf} xleq ...
                        sprintf('[cell2struct(cell(%.0f,0,0),%s,1)]', ...
                        length(xlvalfields), any2ascii(xlvalfields)) xllo];
                    continue;
                end

                % if 'exact' conversion, treat non-empty, non-integer doubles
                if xlconv == 2 && ...
                   ~isempty(xlvalue) && ...
                    isa(xlvalue, 'double') && ...
                    isreal(xlvalue) && ...
                    any(floor(xlvalue(:)) ~= xlvalue(:))
                    xldispout = [xldispout xlind ...
                        xlfields{numf} xleq any2ascii(xlvalue,'exact') xllo];

                % otherwise simply any2ascii
                else
                    xldispout = [xldispout xlind ...
                        xlfields{numf} xleq any2ascii(xlvalue) xllo];
                end
            end
        end

        % add substructs to dispout
        for numf = 1:length(xlxl)
            xldispout = [xldispout xllo ...
                i_dispxl([xlname '.' xlxl{numf}], ...
                         xlcont.(xlxl{numf}), xlind, xllo, xleq, xlconv)];
        end

    % multi-dim structs
    else

        % initialize output
        xldispout = '';

        % get number of dims (N-D), number of elements and start with (1,1,...)
        xlmsize = size(xlsize, 2);
        xlpsize = prod(xlsize);
        xlisize = ones(1,xlmsize);

        % iterate over struct elements
        for xlpcount = 1:xlpsize

            % format indexing string
            xlpindex = strrep(['(' sprintf('%.0f,',xlisize) ')'],',)',')');

            % add substruct to dispout
            xldispout = [xldispout xllo ...
                i_dispxl([xlname xlpindex], ...
                         xlcont(xlpcount), xlind, xllo, xleq, xlconv)];

            % get new index info
            for uc = xlmsize:-1:1

                % increase index
                xlisize(uc) = xlisize(uc) + 1;

                % if index less or equal max index, end this loop
                if xlisize(uc) <= xlsize(uc)
                    break;
                end

                % set index back to one and go to next
                xlisize(uc)=1;
            end
        end

        % if we have output, remove final line feed
        if ~isempty(xldispout)
            xldispout(1:length(xllo)) = [];
        end
    end
% end of function i_dispxl


function [xlp, pre, post, errs] = i_parse(fcont, et, icnv)

    % persistent function handles
    persistent splittocell;
    if isempty(splittocell)
        nelf = neuroelf;
        splittocell = nelf.splittocell;
    end

    % initialize vars
    xlp = struct;
    pre = '';
    post = '';
    errs = {};

    % empty altogether
    if isempty(fcont)
        return;
    end

    % any post-file content ?
    endofini = strfind(fcont, et);
    if ~isempty(endofini)

        % no content at all
        if endofini(1) == 1
            post = fcont;
            return;
        end

        % only accept end terms starting on a new line
        while ~isempty(endofini) && ...
            fcont(endofini(1) - 1) ~= char(10) && ...
            fcont(endofini(1) - 1) ~= char(13)
            endofini(1) = [];
        end

        % if we have a valid terminator
        if ~isempty(endofini)

            % get footer
            post = fcont(endofini(1):end);

            % remove footer from file
            fcont(endofini(1):end) = [];
        end
    end

    % any pre-file content ?
    begofini = strfind(fcont, '[');

    % first section must start on a new line
    while ~isempty(begofini) && ...
        begofini(1) > 1 && ...
        fcont(begofini(1) - 1) ~= char(10) && ...
        fcont(begofini(1) - 1) ~= char(13)
        begofini(1) = [];
    end

    % no sections at all
    if isempty(begofini)
        pre = fcont;
        return;

    % first section not at beginning of file
    elseif begofini(1) > 1

        % set file header
        pre = fcont(1:begofini(1) - 1);

        % remove header from file
        fcont(1:begofini(1) - 1) = [];
    end

    % split lines to cell array (removing empty lines)
    [fcntcell, nolines] = splittocell(fcont, char([13, 10]), 1, 1);

    % but also make a char array for strmatch & find
    fcntchar = char(fcntcell);

    % on what lines do we have sections?
    sectline = find(fcntchar(:, 1) == '[');

    % number of sections
    sectlines = size(sectline,1);

    % terminate sections correctly
    sectline(sectlines + 1, 1) = nolines + 1;

    % iterate over sections (lines)
    for linecount = 1:sectlines

        % get section name line and name range
        nextsect = fcntcell{sectline(linecount)};
        linelength = numel(nextsect) - 1;

        % and extract section name and test it
        nextsect = nextsect(2:linelength);
        validsect = false;
        obrackets = find(nextsect == '(');
        if ~isempty(nextsect)
            try
                eval(['validsects = struct; validsects.' nextsect ' = 1;']);
                if isempty(obrackets) || ...
                   any(nextsect(1:obrackets(1) - 1) == '.')
                    validsect = true;
                end
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end

        % in case of a valid section
        if validsect

            % section start, end
            sectfrom = sectline(linecount) + 1;
            sectto = sectline(linecount+1) - 1;

            % initialize section structure
            sectstruct = struct;

            % section with any contents?
            if sectto >= sectfrom

                % iterate over section lines
                for insectline = sectfrom:sectto

                    % split into name and value
                    eqfound = find(fcntcell{insectline} == '=');
                    if isempty(eqfound)
                        setname = 'Value';
                        setvalue = fcntcell{insectline};
                    else
                        setname = fcntcell{insectline}(1:eqfound(1) - 1);
                        setvalue = fcntcell{insectline}(eqfound(1) + 1:end);
                    end

                    % check for validity of name
                    rsetname = deblank(setname);
                    if ~isvarname(rsetname)
                        errs{end + 1} = ...
                            [sprintf( ...
                             ['xini::i_parse -> ' ...
                             'bad line %d in file %s, section %s:'], ...
                             insectline, fname, nextsect) ...
                             char(10) ' -> ' ...
                             fcntcell{insectline}];
                        continue;
                    end

                    % check for equal sign on line
                    if strcmp(setname, 'Value') && ...
                        isempty(setvalue)
                        errs{end+1} = ...
                            [sprintf( ...
                             ['xini::i_parse -> ' ...
                             'bad line %d in file %s, section %s:'], ...
                             insectline, fname, nextsect) ...
                             char(10) ' -> ' ...
                             fcntcell{insectline}];
                        setvalue = '[]';
                    end

                    % value with at least 2 chars and conversion activated?
                    if numel(setvalue) > 1 && ...
                        icnv > 0

                        % valid bracketting?
                        if ((setvalue(1) == '[' && ...
                             setvalue(end) == ']') || ...
                            (setvalue(1) == '{' && ...
                             setvalue(end) == '}'))

                            % try conversion (replace $SELF.(...) with xlp.(...))
                            try
                                setvalue = ...
                                    eval(strrep(setvalue, '$SELF.', 'xlp.'));

                            % if fails, bark out or log error
                            catch ne_eo;
                                errs{end+1} = ...
                                    ['xini::i_parse -> ' ...
                                     'conversion failed: ' ...
                                     '[' fname '].' ...
                                     nextsect '.' rsetname ' = ' ...
                                     setvalue ' (' ne_eo.message ')'];

                                % remove further brackets for future use
                                while ((setvalue(1) == '[' && ...
                                    setvalue(end) == ']') || ...
                                   (setvalue(1) == '{' && ...
                                    setvalue(end) == '}'))
                                    setvalue = setvalue(2:(end - 1));
                                end
                            end
                        end
                    end

                    % set correct "empty" values
                    if isempty(setvalue)
                        if ischar(setvalue)
                            setvalue = '';
                        elseif isnumeric(setvalue)
                            setvalue = [];
                        end
                    end

                    % now put value into section structure
                    sectstruct.(rsetname) = setvalue;
                end
            end

            % and finally try to insert section structure into IniFile content
            try
                eval(['xlp.' nextsect ' = sectstruct;']);
            catch ne_eo;
                errs{end+1} = ...
                    ['xini::i_parse -> ' ...
                     'couldn''t add ' nextsect ' to the content (' ne_eo.message ')'];
            end
        end
    end
% end of function i_parse



% extcaller  - from where was a call issued
function fromwhere = extcaller(varargin)

    % persistent versioning information
    persistent ec_mlv;
    if isempty(ec_mlv)
        ec_mlv = version;
        if ec_mlv(1) > '6'
            ec_mlv = 'file';
        else
            ec_mlv = 'name';
        end
    end

    % how to get stack info
    if ec_mlv(1) == 'n'
        pathstack = dbstack;
    else
        pathstack = dbstack('-completenames');
        while ~isempty(pathstack) && strcmpi(pathstack(1).file, 'xini.m')
            pathstack(1) = [];
        end
    end

    % set default caller name
    fromwhere = 'CONSOLE_OR_GUI';
    if numel(pathstack) < 1
        return;
    end

    % input argument handling
    if nargin < 1
        func = lower(pathstack(1).(ec_mlv));

    elseif ischar(varargin{1})
        func = lower(varargin{1});

    elseif isnumeric(varargin{1}) && ...
       ~isempty(varargin{1}) && ...
       ~isnan(varargin{1}(1))
        func = lower(pathstack(1).(ec_mlv));
        if varargin{1}(1) == 1
            func = fileparts(func);
        end

    else
        error( ...
            'neuroelf:BadArgument', ...
            'Unsupported input argument.' ...
        );
    end

    % try to find matching caller
    for dbsi = 1:numel(pathstack)
        if isempty(strfind(lower(pathstack(dbsi).(ec_mlv)), func))
            break;
        end
    end
    if dbsi < numel(pathstack) || ...
        isempty(strfind(lower(pathstack(dbsi).(ec_mlv)), func))
        fromwhere = pathstack(dbsi).(ec_mlv);
    end
% end of function extcaller

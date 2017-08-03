classdef xff < handle
%XFF  Create and handle xff (eXtended File Formats) objects.
%
%   OBJ = XFF without any arguments gives access to the class factory,
%   which keeps track of all instantiated objects. It is important to
%   understand that this is somewhat similar to graphics handle objects,
%   which, once created, are not automatically deleted at the end of a
%   function (when the variable containing the handle is cleared from
%   a workspace/memory).
%
%   All objects, other than the factory, must be deleted to free the memory
%   they occupy. To do so, use OBJ.ClearObject; or clearxffobjects({OBJ});
%
%   To access any loaded (or created) object, the object identifier (either
%   a numeric ID, the hash ID, or its filename) can be used, similar to
%   VisualBasic syntax in Excel, when addressing WorkSheet objects (see
%   more details below).
%
%   OBJ = XFF(FILENAME) loads the data in file FILENAME (see supported
%   formats below), and makes them accessible as a handle-based object via
%   variable OBJ. Typically the file format is determined by the extension.
%   If the extension doesn't identify which format to use, a set of Magic
%   tokens is checked against the file. If no format can be identified, an
%   error is produced. The same applies if the file is not found or not
%   readable.
%
%   OBJ = XFF(FILENAME, FORMAT) forces the file format specified in FORMAT
%   to be used to interpret the data. So for instance to load a DICOM file,
%   use OBJ = XFF(FILENAME, 'dcm'); regardless of the filename. FORMAT must
%   be a char array (string) of either 1x3 or 1x4 size. For a list of all
%   available formats, use OBJ = XFF; disp(OBJ.Extensions);
%
%   OBJ = XFF(FILENAMES) creates an 1xF XFF array if FILENAMES is a cell
%   array of filenames. If (and only if) the list of FILENAMES refers to
%   a set of NIFTI (HDR) files with the same spatial layout (size, etc.),
%   the constructor attempts to load all data into a single 1x1 XFF object.
%
%   OBJ = XFF('*.FMT') shows a file selector for files of type FMT. If a
%   file is selected, this file will be loaded. If no file is selected,
%   XFF returns the empty array (0x0 of type XFF). If the given format does
%   not exist, a generic file selector ('*.*' files) is displayed instead.
%
%   OBJ = XFF('*.FMT', NUMBER) allows to pick multiple files, which will
%   be returned as 1xF XFF array. It is not necessary that all files are
%   of the same type, so long as all files can be loaded correctly, given
%   that XFF is a pseudo-abstract class.
%
%   OBJ = XFF('*.FMT', SELECTORTITLE) in addition sets the text for the
%   selector to inform people what file exactly needs to be loaded. In
%   addition, SELECTORTITLE can also be set to a 1x1 double (integer)
%   number representing the number of files to select instead.
%
%   OBJ = XFF('*.*') shows a generic file selector for the most typically
%   loaded files and then proceeds the same way.
%
%   OBJ = XFF(FILENAME, 'h') only loads the header of an object (only
%   applies to binary objects specified in a BFF file format).
%
%   OBJ = XFF(FILENAME, 't') forces large fields to be loaded as @transio
%   objects (if they aren't already; also only works for binary files). See
%   the help of @transio for more information. Importantly, content can be
%   accessed "on demand" and doesn't require as much free memory at once
%   for large files (such as GLM files with many or high resolution maps).
%
%   OBJ = XFF(FILENAME, 'T') forces every field to be read as @transio,
%   even if the content is only 1x1 in size (which allows to change header
%   values without loading the content of large files).
%
%   OBJ = XFF('new:FMT') creates a blank (empty) object of type FMT.
%
%   The contents of OBJ is accessible via OBJ.Field (struct subsref
%   notation). And contents can be changed via the same notation (subsasgn)
%   using OBJ.Field = NEWVALUE; and OBJ.Field(SUBSASGN_EXP) = NEWVALUES;
%
%   The same notation is ALSO used for methods, which at first may be
%   somewhat unusual for MATLAB code, but is relatively common in other
%   languages, e.g. VMR = XFF(VMRFILE); VMR.SmoothData3D([3, 3, 3]); which
%   alters the content of the VMR file (in memory). Methods can also
%   return values. If the returned value of a method is a 1x1 struct (e.g.
%   for the BoundingBox method), sub-fields can directly be accessed
%   transparently: VTC = XFF(VTCFILE); BBOX = VTC.BoundingBox.BBox;
%
%   To write an existing object, which has been previously read, back to
%   disk, use OBJ.Save; if the location is read-only, an error will be
%   raised.
%
%   If an object is "new" (it doesn't have an associated filename), use
%   OBJ.SaveAs(NEWFILENAME); instead.
%
%   To revert any changes made to an object that was loaded from disk, use
%   OBJ.ReloadFromDisk;
%
%   Next to values that are stored in the file format, users can add any
%   additional content in a field RunTimeVars, which is a 1x1 struct, using
%   OBJ.RunTimeVars.SomeField = 'additional content';
%   The contents of this field will be saved as a file next to the original
%   file with the same name but extension .rtv instead. If multiple files
%   share the same name but different extensions (e.g. somefile.glm and
%   somefile.vmp), the content of *both* files will be stored in the same
%   auxiliary file called somefile.rtv.
%   To save this additional file use OBJ.SaveRunTimeVars; and if the
%   RunTimeVars struct contains a 1x1 logical field AutoSave, this call
%   will be executed whenever the object itself is saved (both when
%   OBJ.Save; and OBJ.SaveAs(FILENAME); are called).
%
%   General methods (which are available for most* if not all objects)
%
%   * BoundingBox      - details for objects with spatial information
%   * Browse           - adds supported objects to the NeuroElf GUI
%     ClearObject      - remove the object from memory (necessary!)
%     CopyObject       - create a copy of an object (another handle!)
%     DeleteHandle     - remove object-specific information (internal use)
%     FilenameOnDisk   - returns the filename of the object (or '')
%     Filetype         - returns the format (3- or 4-letter extension)
%   * GetVolume(INDEX) - returns the data associated with a volume index
%     Handles          - inquire object-specific information (internal use)
%     Help             - prints a list of object-specific methods
%     Help(METHOD)     - prints the help of an object-specific method
%   * Layout           - 1x15 double for spatial compatibility checks
%     LoadTransIOData  - load all transio-coded fields into memory
%     ReloadFromDisk   - re-loads the content from disk (must be loaded!)
%     Save             - writes a temp file and renames it to the old one
%     SaveAs(NEWNAME)  - writes a temp file and renames it to NEWNAME
%     SetHandle        - set object-specific information (internal use)
%   * UnBrowse         - remove object from NeuroElf's GUI
%
%   In addition to these general methods, type-specific methods are also
%   available. A list of these methods can be printed to the console with
%   OBJ.Help, and method-specific help is available via OBJ.Help(METHOD)
%
%   Special methods are also available for the ROOT object:
%
%     ClearObjects(TYPE)  - clear all objects specified by TYPE
%     Config              - retrieve and set configuration (see help)
%     Document(SPEC)      - find document by name or ID
%     Documents           - retrieve list of documents (char names/IDs)
%     Help(SPEC)          - print help for method or object
%     TransIOSize(TYPE,S) - retrieve (and set) TransIOSize setting
%     UpdateState(TYPE,S) - retrieve (and set) UpdateState setting

% Version:  v1.1
% Build:    16061411
% Date:     Jun-14 2016, 11:53 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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

% hidden properties
properties (Access = 'private', Hidden = true)
    C = struct; % object actual content
    F = '';     % object filename (if set)
    H = struct; % handles (i.e. special fields)
    L = '------------------------'; % internal ID
    S = [];     % file format specification (extension, etc.)
end


% methods (no indentation for readability)
methods (Access = 'public', Hidden = true)


% constructor, bulk of work for OBJECT = xff(FILENAME); calls
function oo = xff(varargin)

    % get access to neuroelf methods
    global ne_methods xffsngl;

    % check for Class initialization
    if isempty(xffsngl) || ~isstruct(xffsngl) || ~xffsngl.INIT

        % run initialization, which will call xff() to create ROOT
        xffinit(sprintf('v%s;%d', neuroelf_version, neuroelf_build));
    end

    % default constructor
    if nargin == 0

        % on first call, this is not yet set
        if isempty(xffsngl.OBJS)

            % fill ROOT object (content)
            oo.C = struct('Extensions', xffsngl.EXT, 'Formats', xffsngl.FF, ...
                'Magic', xffsngl.FF.magic, 'Methods', xffsngl.FM);
            oo.F = '<ROOT>';
            oo.H = struct('xff', 0, 'CleanUp', {{}}, 'ShownInGUI', false, 'SourceObject', -1);
            oo.L = 'XXXXXXXXXXXXXXXXXXXXXXXX';
            oo.S = struct('DefaultProperty', {{'Document'}}, 'Extensions', {{'ROOT'}}, 'FFTYPE', 'ROOT');

        % otherwise (subsequent calls)
        else

            % return existing object (first element)
            oo = xffsngl.OBJS{1, 4};
        end
        return;
    end

    % shortcuts
    conf = xffsngl.CONF;
    ext = xffsngl.EXT;
    exn = fieldnames(ext);
    ffnd = ne_methods.findfirst;
    objs = xffsngl.OBJS;

    % first argument is a 0, perform a requested action
    if nargin > 1 && isequal(varargin{1}, 0) && ischar(varargin{2})

        % which action
        switch (lower(varargin{2}(:)'))

            % remove one object from the array
            case 'clearobj'

                % input argument (object ID) OK?
                if nargin < 3 || ~ischar(varargin{3}) || size(varargin{3}, 2) ~= 24 || ...
                    any(varargin{3}(:) < '0' | upper(varargin{3}(:)) > 'F')
                    error('neuroelf:xff:invalidLookup', 'Invalid object lookup');
                end

                % find object(s)
                for oc = 1:size(varargin{3}, 1)
                    olup = ffnd(strcmpi(objs(:, 3), varargin{3}(oc, :)));

                    % delete found object
                    if ~isempty(olup)
                        delete(objs{olup, 4});
                    end
                end

                % return ROOT (constructor compliance)
                oo = objs{1, 4};

            % get object from lookup
            case 'object'

                % input argument (object ID/document name) OK?
                if nargin < 3 || ~ischar(varargin{3}) || isempty(varargin{3})
                    error('neuroelf:xff:lookupError', 'Invalid object lookup.');
                end

                % file ID
                if isequal(size(varargin{3}), [1, 24]) && ...
                    all(varargin{3} >= '0' & upper(varargin{3}) <= 'F') && ...
                    any(strcmpi(objs(:, 3), varargin{3}))
                    oo = objs{strcmpi(objs(:, 3), varargin{3}), 4};

                % try filename
                else
                    lup = strrep(lower(varargin{3}(:)'), '\', '/');
                    lup = ffnd(strcmpi(objs(:, 1), lup));
                    if isempty(lup)
                        error('neuroelf:xff:lookupError', 'Object not found.');
                    end
                    oo = objs{lup, 4};
                end

            % create object from struct
            case 'makeobject'

                % any 3rd argument?
                if nargin < 3
                    error('neuroelf:xff:missingInput', 'Missing input argument.');
                end

                % generate object
                try
                    oo = xffmakeobject(varargin{3});
                catch xfferror
                    rethrow(xfferror);
                end

            % list of objects (copy of persistent array without factory)
            case 'objects'
                oo = cat(1, objs{:, 4});
                oo(1) = [];

            % clear storage completely
            case 'clearallobjects'

                % use delete (in case there are GZip files, etc.)
                delete(cat(1, objs{2:end, 4}));
                oo = objs{1, 4};

            % bail out on study commands
            otherwise
                error('neuroelf:xff:badAction', 'Invalid action given.');
        end
        return;
    end

    % single char argument (or char + double)
    if nargin > 0 && nargin < 3 && ischar(varargin{1}) && ~isempty(varargin{1}) && ...
       (nargin < 2 || (isa(varargin{2}, 'double') && numel(varargin{2}) == 1) || ...
        (ischar(varargin{2}) && ~isempty(varargin{2})))

        % linearize filename
        filename = varargin{1}(:)';

        % unique ID (24 hex)
        if nargin == 1 && numel(filename) == 24 && ...
            all(filename >= '0' & upper(filename) <= 'F') && size(objs, 1) > 1

            % try to locate ID
            xids = objs(:, 3);
            xidi = ffnd(strcmpi(xids, filename));

            % found?
            if ~isempty(xidi)

                % return object
                oo = objs{xidi, 4};
                return;
            end
        end

        % get filename components
        [fnparts{1:3}] = fileparts(filename);

        % support gzip-ed data
        isgziped = '';
        if conf.loadgziped && strcmpi(fnparts{3}, '.gz') && ...
            ~any(fnparts{2} == '*' | fnparts{2} == ':')
            isgziped = fnparts{3};
            filename(end-2:end) = [];
            [fnparts{1:3}] = fileparts(filename);
        end

        % check for "*.???"
        if ~isempty(fnparts{2}) && any(fnparts{2} == '*')

            % force use of extension for ?ff function
            extf = fnparts{3}(fnparts{3} ~= '.');

            % allow special case for img
            if any(strcmpi(extf, {'hdr', 'img', 'nii'})) && ...
               ~isempty(fnparts{1}) && numel(fnparts{2}) > 1

                % try to find images
                try

                    % requires full path
                    [imgfp, imgfn, imgfe] = fileparts(filename);
                    if isempty(imgfp)
                        error('neuroelf:xff:badArgument', ...
                            'Multi HDR/IMG lookup requires full path.');
                    end

                    % look up files
                    imgfs = dir([imgfp, filesep, imgfn, imgfe]);

                    % found
                    if ~isempty(imgfs)

                        % try to open
                        oo = xff(regexprep({imgfs.name}, '^(.)', [imgfp '/$1']));
                        if iscell(oo) && numel(oo) == 1
                            oo = oo{1};
                        end

                        % and return
                        return;
                    end

                % somthing went wrong...
                catch xfferror
                    rethrow(xfferror);
                end
            end

            % how many files
            if nargin < 2 || ~isa(varargin{2}, 'double') || isnan(varargin{2}) || varargin{2} < 1
                numreq = 1;
            else
                numreq = floor(varargin{2});
            end

            % title given
            if nargin > 1 && ischar(varargin{end}) && numel(varargin{end}) > 1
                dtitle = varargin(end);
            else
                dtitle = cell(0, 1);
            end

            % request filename
            filename = xffrequestfile(numreq, filename, ext, xffsngl.FF, dtitle{:});

            % no filename(s) returned
            if isempty(filename)

                % return empty outputs
                oo(:) = [];
                return;

            % for single filename
            elseif ischar(filename)

                % cell-ify
                filename = {filename};
            end
            nfiles = numel(filename);

            % still allow "reuse-same" syntax
            if nfiles == 1 && ~isempty(filename{1}) && nargin > 1 && ...
                isa(varargin{2}, 'double') && numel(varargin{2}) == 1 && varargin{2} == 0
                if ispc
                    namelup = ffnd(strcmpi(filename{1}, objs(:, 1)));
                else
                    namelup = ffnd(strcmp(filename{1}, objs(:, 1)));
                end
                if ~isempty(namelup)
                    oo = objs{namelup, 4};
                    return;
                end
            end

            % create necessary (blank) output object(s)
            oo = oo(1, ones(1, nfiles));
            ol = false(1, nfiles);

            % try (re-) loading objects
            for oc = 1:nfiles
                try

                    % re-load loaded objects?
                    if conf.reloadsame

                        % then pass on as single file (with correct ext)
                        oo(oc) = xff(filename{oc}, extf);
                        ol(oc) = true;

                    % otherwise
                    else

                        % see if already loaded?
                        if ispc
                            namelup = ffnd(strcmpi(filename{oc}, objs(:, 1)));
                        else
                            namelup = ffnd(strcmp(filename{oc}, objs(:, 1)));
                        end

                        % found?
                        if ~isempty(namelup)
                            oo(oc) = objs{namelup, 4};

                        % if not found
                        else

                            % pass on to load
                            oo(oc) = xff(filename{oc}, extf);
                            ol(oc) = true;
                        end
                    end

                % something went wrong with loading
                catch xfferror
                    clearxffobjects(oo(ol));
                    rethrow(xfferror);
                end
            end
            return;

        % patch extension only for OBJECT = xff(TYPE); calls
        elseif numel(filename) < 5 && ...
           (exist(filename, 'file') ~= 2 || exist([filename '.m'], 'file') == 2)

            % prepend 'new:' part
            filename = ['new:' filename];
        end

        % checking for 'new:???' (for valid types only)
        if numel(filename) > 5 && strcmpi(filename(1:4), 'new:') && ...
           ~any(filename == '\' | filename == '/' | filename == '.') && ...
            isfield(ext, lower(filename(5:end)))

            % get spec
            fftype = lower(filename(5:end));
            ffspec = ext.(fftype);
            xfft = ffspec{1}(end-2:end);
            spec = xffsngl.FF.(xfft)(ffspec{2});

            % no "new" code is available
            if isempty(spec.NewFileCode)
                error('neuroelf:xff:incompleteSpec', ...
                    'For %s type files, no NewFileCode is available.', fftype);
            end

            % create new object's content
            try
                oo.C = xffnewcont(fftype);
                nlup = oo.C.RunTimeVars.xffID;
                oo.F = '';
                oo.H = struct('xff', nlup, 'CleanUp', {{}}, 'ShownInGUI', false, 'SourceObject', []);
                oo.L = nlup;
                oo.S = spec;

                % for types with spatial content
                if any(strcmp(fftype, {'ava', 'cmp', 'dmr', 'fmr', 'glm', ...
                        'hdr', 'head', 'nlf', 'srf', 'vmp', 'vmr', 'vtc'})) && ...
                   ~isfield(oo.C.RunTimeVars, 'TrfPlus')

                    % add TrfPlus field
                    oo.C.RunTimeVars.TrfPlus = eye(4);
                end

                % add to global storage
                xffsngl.OBJS(end+1, :) = {oo.F, fftype, oo.L, oo};

            % something went wrong (with new:????)
            catch xfferror
                error('xff:EvaluationError', ...
                    'Couldn''t evaluate NewFileCode snippet for type %s: %s.', ...
                    fftype, xfferror.message);
            end
            return;
        end

        % make absolute!
        [isabs{1:2}] = isabsolute(filename);
        filename = isabs{2};

        % get file name parts
        [fx{1:3}] = fileparts(filename);
        fx = regexprep(fx{3}, ',\d+$', '');
        if ~isempty(fx) && fx(1) == '.'
            fx(1) = [];
        end

        % override extension
        if ~any(strcmpi(fx, exn)) && nargin > 1 && ischar(varargin{2}) && ...
            numel(varargin{2}) <= 5 && any(strcmpi(varargin{2}(:)', exn))
            fx = lower(varargin{2}(:)');
        end

        % any extension based match
        exf = strcmpi(fx, exn);
        if any(exf)

            % get match and matching extension field name
            exf = ffnd(exf);
            exn = exn{exf(1)};

            % look up in extensions
            if isfield(ext, exn)
                ffspec = ext.(exn);
                xfft = ffspec{1}(end-2:end);
                ff = xffsngl.FF.(xfft)(ffspec{2});

            % or error !
            else
                error('neuroelf:xff:badExtension', ...
                    'Extension given, but not in BFF/TFF list: %s.', exn);
            end

            % match found
            exf = true;
        else

            % no match found
            exf = false;
        end

        % not yet identified, try magic
        if ~exf

            % not supported with gzipped files
            if ~isempty(isgziped)
                error('neuroelf:xff:noMagicWithGzip', ...
                    'Magic token detection not supported with gziped files.');
            end
            maf = false;
            try
                detmag = xffdetectmagic(filename, xffsngl.MAG);
            catch xfferror
                warning(xfferror.identifier, xfferror.message);
                detmag = '';
            end
            if ~isempty(detmag)

                % either bff
                if isfield(ext, detmag)
                    ffspec = ext.(detmag);
                    xfft = ffspec{1}(end-2:end);
                    ff = xffsngl.FF.(xfft)(ffspec{2});

                % or error !
                else
                    error('neuroelf:xff:badExtension', ...
                        'Magic found, but type not in BFF/TFF list: %s.', detmag);
                end
                maf = true;
            end
        end

        % file type could not be identified
        if ~exf && ~maf
            error('neuroelf:xff:badFileContent', ...
                'Unknown file type. Cannot read ''%s''.', filename);
        end

        % get filetype (fix extension)
        fft = lower(ff.Extensions{1});

        % reload same ?
        if nargin < 2 || ischar(varargin{2})
            rls = conf.reloadsame;
        else
            rls = (varargin{2} ~= 0);
        end

        % try to find among loaded objects
        if ~rls

            % look up
            if ispc
                namelup = ffnd(strcmpi(filename, objs(:, 1)));
            else
                namelup = ffnd(strcmp(filename, objs(:, 1)));
            end

            % if found
            if ~isempty(namelup)

                % set and return
                oo = objs{namelup, 4};
                return;
            end
        end

        % special treatments (header only, transio/forced, reuse, etc.)
        fullf = true;
        if nargin > 1 && ischar(varargin{2}) && numel(varargin{2}) == 1

            % header only
            if varargin{2} == 'h'
                fullf = false;
                fullfa = '-fh';

            % transio (for fields > minsize)
            elseif varargin{2} == 't' && ...
                strcmpi(ff.FFTYPE, 'bff')
                ff.TransIOSize = 1024;

            % transio FOR ALL FIELDS (incl. 1x1!)
            elseif varargin{2} == 'T' && ...
                strcmpi(ff.FFTYPE, 'bff')
                ff.TransIOSize = 1;

            % verbose loading
            elseif varargin{2} == 'v'
                fullf = false;
                fullfa = '-v';

            % content only (don't register object in factory)
            elseif varargin{2} == 'c'
                fullf = false;
                fullfa = '';
            elseif varargin{2} == 'r'
                if ispc
                    namelup = ffnd(strcmpi(filename, objs(:, 1)));
                else
                    namelup = ffnd(strcmp(filename, objs(:, 1)));
                end
                if ~isempty(namelup)
                    oo = objs{namelup, 4};
                    return;
                end
            end
        end

        % uncompress if needed
        if ~isempty(isgziped)
            ofname = filename;
            try
                tdir = conf.settings.GZip.TempDir;
                gunzip([filename isgziped], tdir);
                filename = [tempname(tdir) fnparts{3}];
                if movefile([tdir '/' fnparts{2} fnparts{3}], filename) ~= 1
                    error('neuroelf:xff:moveFileError', ...
                        'Error renaming xff temporary file.');
                end
            catch xfferror
                rethrow(xfferror);
            end
        end

        % read file
        try

            % generate a 24-char hex string
            nlup = hxdouble(randn(1, 2));
            nlup = nlup([4:15, 20:31]);

            % re-generate until it's a unique one
            while any(strcmpi(objs(:, 3), nlup))
                nlup = hxdouble(randn(1, 2));
                nlup = nlup([4:15, 20:31]);
            end

            switch lower(ff.FFTYPE)

                % for binary files
                case 'bff'

                    % allow sub-volumes
                    subvolm = regexpi(filename, ',\s*\d+\s*$');
                    if ~isempty(subvolm)
                        if ~any(strcmp(fft, {'hdr', 'head'}))
                            error('neuroelf:xff:badArgument', ...
                                'Sub-volume selection only for HDR/NII/HEAD files.' ...
                            );
                        end
                        subvol = str2double(filename(subvolm+1:end));
                        filename = deblank(filename(1:subvolm-1));
                    else
                        subvol = [];
                    end

                    % full
                    if fullf
                        oo.C = ne_methods.bffio(filename, ff);
                        oo.C.RunTimeVars.xffID = nlup;
                        oo.F = filename;
                        oo.H = struct('xff', nlup, 'CleanUp', {{}}, ...
                            'ShownInGUI', false, 'SourceObject', []);
                        oo.L = nlup;
                        oo.S = ff;

                    % or header
                    else
                        oo.C = ne_methods.bffio(filename, ff, fullfa);
                        oo.C.RunTimeVars.xffID = 'HHHHHHHHHHHHHHHHHHHHHHHH';
                        oo.L = 'HHHHHHHHHHHHHHHHHHHHHHHH';
                        oo.S = struct('Extensions', {ff.Extensions});
                        return;
                    end

                    % sub-volume access
                    if ~isempty(subvol)

                        % depends on extension
                        try
                            switch (lower(ff.Extensions{1}))

                                % HDR (NII)
                                case 'hdr'

                                    % patch Map and Mat44 fields, if exists
                                    if isfield(oo.C.RunTimeVars, 'Map') && ...
                                        numel(oo.C.RunTimeVars.Map) == size(oo.C.VoxelData, 4)
                                        oo.C.RunTimeVars.Map = oo.C.RunTimeVars.Map(subvol);
                                    end
                                    if isfield(oo.C.RunTimeVars, 'Mat44') && ...
                                        size(oo.C.RunTimeVars.Mat44, 3) == size(oo.C.VoxelData, 4)
                                        oo.C.RunTimeVars.Mat44 = oo.C.RunTimeVars.Mat44(:, :, subvol);
                                    end

                                    % as well as ImgDim(5) and content
                                    oo.C.ImgDim.Dim(5) = 1;
                                    oo.C.VoxelData = oo.C.VoxelData(:, :, :, subvol, :);

                                % HEAD (BRIK)
                                case 'head'
                                    oo.C.NrOfVolumes = 1;
                                    oo.C.Brick = oo.C.Brick(subvol);
                            end

                        % something went wrong
                        catch xfferror
                            error('xff:BadArgument', ...
                                'Couldn''t access subvolume %d in file %s: %s.', ...
                                subvol, filename, xfferror.message);
                        end

                        % add to filename
                        oo.F = sprintf('%s,%d', filename, subvol);
                    end

                % for text-based files
                case 'tff'

                    % full file
                    if fullf

                        % read content
                        oo.C = ne_methods.tffio(filename, ff);

                        % new ID
                        oo.C.RunTimeVars.xffID = nlup;

                        % store filename
                        oo.F = filename;

                        % and default handles, ID, and file format
                        oo.H = struct('xff', nlup, 'CleanUp', {{}}, 'ShownInGUI', false, 'SourceObject', []);
                        oo.L = nlup;
                        oo.S = ff;

                    % or header (object without filename/factory tracking)
                    else
                        oo.C = ne_methods.tffio(filename, ff, fullfa);
                        oo.C.RunTimeVars.xffID = 'HHHHHHHHHHHHHHHHHHHHHHHH';
                        oo.L = 'HHHHHHHHHHHHHHHHHHHHHHHH';
                        oo.S = struct('Extensions', {ff.Extensions});
                        return;
                    end

                    % for later (RTV compatibility)
                    subvol = [];

                % this should never happen, but who knows
                otherwise
                    error('neuroelf:xff:badFFTYPE', ...
                        'FF type %s not supported yet.', ff.FFTYPE);
            end

            % to read RunTimeVars, check GZIP-ed status
            if isempty(isgziped)
                [filenp, filenn, filene] = fileparts(filename);
            else
                [filenp, filenn, filene] = fileparts(ofname);
            end

            % try
            try

                % open the RTV file (faster than try, load, catch)
                filenm = fopen([filenp '/' filenn '.rtv']);

                % exists and readable
                if filenm > 0

                    % close again
                    fclose(filenm);

                    % then load as MAT file
                    filenr = load([filenp '/' filenn '.rtv'], '-mat');

                    % correct content
                    if isfield(filenr, 'RunTimeVars') && isstruct(filenr.RunTimeVars) && ...
                        numel(filenr.RunTimeVars) == 1 && isfield(filenr.RunTimeVars, filene(2:end))

                        % get file-type specific content only
                        rtv = filenr.RunTimeVars.(filene(2:end));
                        rtvf = fieldnames(rtv);

                        % set in object RunTimeVars (override defaults)
                        for rtvfc = 1:numel(rtvf)
                            oo.C.RunTimeVars.(rtvf{rtvfc}) = rtv.(rtvf{rtvfc});
                        end

                        % no sub-vol?
                        if isempty(subvol) && ~any(strcmpi(objs(:, 3), oo.C.RunTimeVars.xffID))

                            % keep original ID
                            oo.L = oo.C.RunTimeVars.xffID;

                        % with sub-vol
                        else

                            % reverse (set new ID)
                            oo.C.RunTimeVars.xffID = oo.L;
                        end

                        % in addition, check fields from RunTimeVars
                        rtvf = fieldnames(oo.C.RunTimeVars);
                        for rtvfc = 1:numel(rtvf)

                            % if CONTENT has a struct field with same name & size
                            if isfield(oo.C, rtvf{rtvfc}) && isstruct(oo.C.RunTimeVars.(rtvf{rtvfc})) && ...
                                isstruct(oo.C.(rtvf{rtvfc})) && ...
                                numel(oo.C.RunTimeVars.(rtvf{rtvfc})) == numel(oo.C.(rtvf{rtvfc}))

                                % store into those as sub-RunTimeVars
                                for rtvfcc = 1:numel(oo.C.(rtvf{rtvfc}))
                                    oo.C.(rtvf{rtvfc})(rtvfcc).RunTimeVars = ...
                                        oo.C.RunTimeVars.(rtvf{rtvfc})(rtvfcc);
                                end

                                % and remove from the content
                                oo.C.RunTimeVars = rmfield(oo.C.RunTimeVars, rtvf{rtvfc});
                            end
                        end
                    end
                end

            % something went wrong for RunTimeVars
            catch xfferror
                fprintf('Error loading RunTimeVars: %s.\n', xfferror.message);
            end

            % validation code
            if ~isempty(ff.ValidFileCode)
                try

                    % assign to variable (to protect object)
                    ffcont = oo.C;

                    % execute code
                    eval(regexprep(ff.ValidFileCode, '\@([a-zA-Z])', 'ffcont.$1'));

                    % if no error occurred, re-assign back
                    oo.C = ffcont;

                % with error
                catch xfferror

                    % simply print warning
                    fprintf('Error executing ValidFileCode for %s: %s.\n', ...
                        filename, xfferror.message);
                end
            end

            % add field for per-object transformation where useful
            if any(strcmp(fft, {'ava', 'cmp', 'dmr', 'fmr', 'glm', 'hdr', ...
                    'head', 'nlf', 'srf', 'vmp', 'vmr', 'vtc'})) && ...
               ~isfield(oo.C.RunTimeVars, 'TrfPlus')
                oo.C.RunTimeVars.TrfPlus = eye(4);
            end

            % complete object specs
            if ispc
                oo.F = strrep(oo.F, '\', '/');
            end

            % tally all (non-header/content-only) objects (prevent clear)
            xffsngl.OBJS(end+1, :) = {oo.F, ff.Extensions{1}, oo.L, oo};

            % came from gziped file
            if ~isempty(isgziped)
                oo.H.GZIPext = isgziped;
                oo.H.GZIPfile = ofname;
            end

            % bring up in GUI
            if conf.loadingui && ...
                any(strcmpi(ff.Extensions{1}, {'dmr', 'fmr', 'glm', 'vmp', 'vmr'}))
                aft_Browse(oo);
            end

            % now return
            return;

        % general error in this branch
        catch xfferror
            error('xff:XFFioFailed', 'Error calling *ffio(...): ''%s''.', ...
                xfferror.message);
        end

    % elseif ... other input argument combinations
    elseif nargin == 2 && ischar(varargin{1}) && ...
       ~isempty(varargin{1}) && numel(varargin{2}) == 1 && ...
        xffisobject(varargin{2}, true)

        % try writing file
        try
            oo = varargin{2};

            % don't allow volume marker
            if ~isempty(regexpi(varargin{1}(:)', ',\s*\d+\s*$'))
                error('neuroelf:xff:badArgument', 'Saving of sub-volumes not permitted.');
            end

            % based on type
            switch lower(ostr.S.FFTYPE)

                % for binary files
                case 'bff'

                    % reassign after writing
                    oo.C = ne_methods.bffio(varargin{1}, oo.S, oo.C);

                % for text-based files
                case 'tff'

                    % also re-assign filename!
                    [oo.C, oo.F] = ne_methods.tffio(varargin{1}, oo.S, oo.C);

                % should never get to this...
                otherwise
                    error('neuroelf:xff:badFFTYPE', ...
                        'FF type %s not supported yet.', ostr.S.FFTYPE);
            end

        % file-write error
        catch xfferror
            error('neuroelf:xff:xFFioFailed', ...
                'Error calling ?ffio(...): ''%s''.', xfferror.message);
        end

    % just a lookup value
    elseif nargin == 1 && isa(varargin{1}, 'double') && numel(varargin{1}) == 1 && ...
       ~isinf(varargin{1}) && ~isnan(varargin{1}) && varargin{1} >= 1

        % return object
        if varargin{1} <= (size(objs, 1) - 1)

            % looked up
            oo = objs{round(varargin{1}) + 1, 4};

        % if not found
        else
            error('neuroelf:xff:badLookup', 'Invalid lookup for handle syntax.');
        end

    % list of filenames
    elseif nargin > 0 && iscell(varargin{1}) && ~isempty(varargin{1})

        % try to load objects
        ol = false(1, numel(varargin{1}));
        oo = oo(1, ones(1, numel(ol)));
        for cc = 1:numel(ol)

            % regular filename
            if ischar(varargin{1}{cc}) && ~isempty(varargin{1}{cc})
                try

                    % also allow .img instead of .hdr
                    filename = regexprep(varargin{1}{cc}, ...
                        '\.img(\s*\,\s*\d+\s*)?$', '.hdr$1', 'preservecase');

                    % with additional arguments
                    if nargin > 1
                        oo(cc) = xff(filename, varargin{2:end});

                    % or without
                    else
                        oo(cc) = xff(filename);
                    end
                    ol(cc) = true;

                    % special case for first of several hdr files
                    if cc == 1 && numel(varargin{1}) > 1 && ...
                        xffisobject(oo(1), true, 'hdr') && ...
                       ~any(cellfun('isempty', regexpi(varargin{1}(:), ...
                            '\.(hdr|img|nii)(\s*\,\s*\d+\s*)?$')))

                        % get content
                        hdrc = oo(1).C;

                        % anything but 3D transio voxel data
                        if isempty(hdrc.VoxelData) || ~istransio(hdrc.VoxelData) || ndims(hdrc.VoxelData) ~= 3

                            % then go on...
                            continue;
                        end

                        % get transio (as struct), endian type
                        hdrt = struct(hdrc.VoxelData);
                        if hdrt.LittleND
                            endtype = 'ieee-le';
                        else
                            endtype = 'ieee-be';
                        end

                        % get values
                        filenames = regexprep(varargin{1}(:), ...
                            '\.hdr(\s*\,\s*\d+\s*)?$', '.img$1', 'preservecase')';

                        % get offset and size
                        hdrbo = hdrc.ImgDim.VoxOffset;
                        hdro = hdrbo + zeros(1, numel(filenames));
                        hdro(1) = hdrt.IOOffset;
                        hdrsz = [hdrt.DataDims, numel(hdro)];
                        hdrby = prod(hdrt.DataDims) * hdrt.TypeSize;

                        % remove volume identifiers
                        hasvol = cellfun('isempty', regexpi(filenames, '\.(img|nii)$'));
                        hasvol(1) = false;
                        if any(hasvol)
                            for scc = 2:numel(hdro)
                                if hasvol(scc)
                                    hdro(scc) = hdrbo + (str2double(regexprep(filenames{scc}, ...
                                        '^.*\,\s*(\d+)\s*$', '$1')) - 1) * hdrby;
                                end
                            end
                        end

                        % update object
                        hdrc.ImgDim.Dim = [numel(hdrsz), hdrsz, ones(1, 7 - numel(hdrsz))];
                        hdrc.ImgDim.PixSpacing(numel(hdrsz) + 1) = 1;

                        % and try to resolve all filenames as transio
                        try
                            hdrc.VoxelData = transio(filenames, ...
                                endtype, hdrt.DataType, hdro, hdrsz);

                        % if error occurred, print error, and return
                        catch xfferror
                            fprintf('Error fast-mapping NIFTI file: %s.\n', xfferror.message);
                            continue;
                        end

                        % otherwise, prepare for further reading
                        hdrc.VoxelDataCT = cell(1, numel(hdro));
                        hdrc.RunTimeVars.Map = ...
                            hdrc.RunTimeVars.Map(ones(1, numel(hdro)));

                        % try to read specific header fields into uint8(:, :)
                        filenames = regexprep(filenames, '\s*\,\s*\d+\s*$', '');
                        filenames = regexprep(filenames, '\.img$', '.hdr', 'preservecase');
                        hfid = 0;
                        hcnt = uint8(0);
                        hcnt(numel(hdro), 256) = 0;

                        % for each header file
                        for scc = 1:numel(hdro)
                            try

                                % open file
                                hfid = fopen(filenames{scc}, 'r', hdrc.Endian);

                                % if failed
                                if hfid < 1

                                    % stop trying
                                    break;
                                end

                                % read first 256 bytes
                                hcnt(scc, :) = fread(hfid, [256, 1], '*uint8')';

                                % close and continue
                                fclose(hfid);
                                hfid = 0;

                            % if an error occurs
                            catch xfferror

                                % close last file
                                if hfid > 0
                                    fclose(hfid);
                                end

                                % print error
                                fprintf('Error fast-mapping NIFTI files: %s.\n', xfferror.message);

                                % but don't bother exiting
                                break;
                            end
                        end

                        % check crucial settings (must match)
                        if any(any(diff(hcnt(:, [41:48, 71:74, 77:120, 253:256])) ~= 0))

                            % continue without data concatenation
                            continue;
                        end

                        % all checks worked? then set descriptions
                        hdesc = deblank(cellstr(char(hcnt(:, 149:228))));
                        [hdrc.RunTimeVars.Map(:).Name] = deal(hdesc{:});

                        % and parse Mat44 information if necessary
                        cfr = hdrc.ImgDim.PixSpacing(2:4);
                        dimf = hdrc.ImgDim.Dim(2:4);
                        dimh = 0.5 + 0.5 * dimf;
                        mat44 = repmat(eye(4), [1, 1, prod(hdrc.ImgDim.Dim(5:8))]);
                        mat4t = 1;
                        mat4f = regexprep(filename, '\.(hdr|nii)$', '.mat', 'preservecase');
                        q = ne_methods.emptystruct({'QSFormCode', ...
                            'QuaternionBCDXYZ', 'AffineTransXYZ'}, [1, 1]);

                        % check for Mat44 files
                        for scc = 1:numel(hdro)

                            % mat file available
                            try
                                hfid = fopen(mat4f{scc}, 'r');
                                if hfid > 0
                                    fclose(hfid);

                                    % load file
                                    hmatcnt = load(mat4f{scc});

                                    % and copy respective content
                                    if isfield(hmatcnt, 'mat')
                                        mfilec = hmatcnt.mat;
                                    elseif isfield(hmatcnt, 'M')
                                        mfilec = hmatcnt.M;
                                    else
                                        mfilec = [];
                                    end
                                    if ~isempty(mfilec)
                                        mat44(:, :, mat4t:mat4t+size(mfilec, 3)-1) = mfilec;
                                        mat4t = mat4t + size(mfilec, 3);
                                        continue;
                                    end
                                end

                            % error during loading (ignore for now)
                            catch xfferror
                                neuroelf_lasterr(xfferror);
                            end

                            % read in NIFTI part of header
                            if hdrc.NIIFileType > 0
                                hfid = fopen(filenames{scc}, 'r', hdrc.Endian);
                                fseek(hfid, 252, -1);
                                q.QSFormCode = fread(hfid, [2, 1], 'int16=>double');
                                q.QuaternionBCDXYZ = fread(hfid, [6, 1], 'single=>double');
                                q.AffineTransXYZ = fread(hfid, [12, 1], 'single=>double');
                                fclose(hfid);

                                % first check SFormCode
                                if q.QSFormCode(2) > 0

                                    % use AffineTransX/Y/Z
                                    mfilec = [reshape(q.AffineTransXYZ, 4, 3)'; ...
                                        0, 0, 0, 1];
                                    mfilec(1:3, 4) = mfilec(1:3, 4) - mfilec(1:3, 1:3) * [1; 1; 1];

                                % next check QFormCode
                                elseif q.QSFormCode(1) > 0

                                    % use that information instead
                                    mfilec = ne_methods.spmtrf(q.QuaternionBCDXYZ(4:6)') * ...
                                        ne_methods.spmq2m(q.QuaternionBCDXYZ(1:3)') * ...
                                        ne_methods.spmtrf([0,0,0], [0,0,0], hdrc.ImgDim.PixSpacing(2:4));
                                    mfilec(1:3, 4) = mfilec(1:3, 4) - mfilec(1:3, 1:3) * [1; 1; 1];

                                % old-school stuff
                                else
                                    mfilec = [ ...
                                        cfr(1), 0 ,   0   , -cfr(1) * dimh(1); ...
                                          0 , cfr(2), 0   , -cfr(2) * dimh(2); ...
                                          0 ,   0 , cfr(3), -cfr(3) * dimh(3); ...
                                          0 ,   0 ,   0   ,    1];

                                    % support default flip
                                    if (hdrc.ImgDim.PixSpacing(1) < 0 || ...
                                        (hdrc.ImgDim.PixSpacing(1) == 0 && ...
                                         hdrc.DataHist.Orientation == 0 && ...
                                         conf.type.hdr.assumeflipped))

                                        % perform x-flip to get to TAL
                                        mfilec(1, :) = -mfilec(1, :);
                                    end
                                end

                            % really old-school stuff
                            else

                                % an originator is given (SPM2)
                                if ~all(hdrc.DataHist.OriginSPM(1:3) == 0)

                                    % use this information
                                    dho = hdrc.DataHist.OriginSPM(1:3);
                                    if all(hdrc.DataHist.OriginSPM(1:3) < 0)
                                        dho = dimf + dho;
                                    end

                                    % then create a transformation matrix
                                    mfilec = [ ...
                                        cfr(1), 0   , 0   , -cfr(1) * dho(1); ...
                                          0 , cfr(2), 0   , -cfr(2) * dho(2); ...
                                          0 ,   0 , cfr(3), -cfr(3) * dho(3); ...
                                          0 ,   0   , 0   ,    1];

                                % for older images
                                else

                                    % depending on orientation flag
                                    switch (double(hcnt(scc, 253)))
                                        case 1
                                            mfilec = [ ...
                                                cfr(1), 0 ,   0   , -cfr(1) * dimh(1); ...
                                                  0 ,   0 , cfr(3), -cfr(3) * dimh(3); ...
                                                  0 , cfr(2), 0   , -cfr(2) * dimh(2); ...
                                                  0 ,   0 ,   0   ,    1];
                                        case 2
                                            mfilec = [ ...
                                                  0 , cfr(2), 0   , -cfr(2) * dimh(2); ...
                                                  0 ,   0 , cfr(3), -cfr(3) * dimh(3); ...
                                                cfr(1), 0 ,   0   , -cfr(1) * dimh(1); ...
                                                  0 ,   0 ,   0   ,    1];
                                        case 3
                                            mfilec = [ ...
                                                cfr(1), 0 ,   0   , -cfr(1) * dimh(1); ...
                                                  0 ,-cfr(2), 0   ,  cfr(2) * dimh(2); ...
                                                  0 ,   0 , cfr(3), -cfr(3) * dimh(3); ...
                                                  0 ,   0 ,   0   ,    1];
                                        case 4
                                            mfilec = [ ...
                                                cfr(1), 0 ,   0   , -cfr(1) * dimh(1); ...
                                                  0 ,   0 ,-cfr(3),  cfr(3) * dimh(3); ...
                                                  0 , cfr(2), 0   , -cfr(2) * dimh(2); ...
                                                  0 ,   0 ,   0   ,    1];
                                        case 5
                                            mfilec = [ ...
                                                  0 , cfr(2), 0   , -cfr(2) * dimh(2); ...
                                                  0 ,   0 ,-cfr(3),  cfr(3) * dimh(3); ...
                                                cfr(1), 0 ,   0   , -cfr(1) * dimh(1); ...
                                                  0 ,   0 ,   0   ,    1];
                                        otherwise
                                            mfilec = [ ...
                                                cfr(1), 0 ,   0   , -cfr(1) * dimh(1); ...
                                                  0 , cfr(2), 0   , -cfr(2) * dimh(2); ...
                                                  0 ,   0 , cfr(3), -cfr(3) * dimh(3); ...
                                                  0 ,   0 ,   0   ,    1];
                                    end
                                end

                                % support default flip
                                if (hdrc.ImgDim.PixSpacing(1) < 0 || ...
                                    (hdrc.ImgDim.PixSpacing(1) == 0 && ...
                                     hdrc.DataHist.Orientation == 0 && ...
                                     conf.type.hdr.assumeflipped))

                                    % perform x-flip to get to TAL
                                    mfilec(1, :) = -mfilec(1, :);
                                end
                            end

                            % add to Mat44
                            mat44(:, :, mat4t:mat4t+size(mfilec, 3)-1) = mfilec;
                            mat4t = mat4t + size(mfilec, 3);
                        end

                        % and finally store in RTV
                        hdrc.RunTimeVars.Mat44 = mat44;

                        % set back
                        oo(1).C = hdrc;
                        oo = oo(1);
                        ol = true;

                        % and leave loop
                        break;
                    end

                % deal with general errors (any filename/object)
                catch xfferror
                    delete(oo);
                    rethrow(xfferror);
                end
            end
        end
        oo(~ol) = [];

    % else
    else
        error('neuroelf:xff:badArgument', 'Bad argument combination or file not writable.');

    % end of argument test if
    end
end

% destructor (delete)
function delete(xo)

    % input checks
    if ~isa(xo, 'xff') || numel(xo) ~= 1
        error('neuroelf:xff:badArgument', 'Bad or missing argument.');
    end

    % don't do anything for ROOT
    if xo.L(1) == 'X' || xo.L(1) == '-'
        return;
    end

    % call clear code and remove from factory
    xffclear(xo.L);
end

% for compatibility
function xo = bless(xo, varargin)
    warning('neuroelf:xff:deprecated', 'Use of BLESS(OBJECT) is deprecated.');
end

% fieldnames
function names = fieldnames(xo)

    % nothing for empty arrays
    if isempty(xo)
        names = cell(0, 1);
        return;
    end

    % names for first item
    names = fieldnames(xo(1).C);

    % and intersect with any other item
    for oc = 2:numel(xo)
        names = intersect(names, fieldnames(xo(oc).C));
    end
end

% methods, properties
function names = methods(xo)

    % keep track of last accessed object (for full help)
    global xffsngl;

    % default is nothing
    names = cell(0, 1);

    % for empty array
    if isempty(xo)
        return;
    end

    % same as last?
    if ~isequal(xo, xffsngl.LAST)
        pril = false;
    else
        pril = true;
    end

    % filetype
    stype = lower(xo(1).S.Extensions{1});

    % add methods if any
    if isfield(xffsngl.FM, stype)
        tm = xffsngl.FM.(stype);
        mf = fieldnames(tm);

        % iterate
        names = mf;
        for mc = 1:numel(mf)
            sm = tm.(mf{mc}){1};
            up = find(sm == '_');
            if pril
                names{mc} = [sm(up(1)+1:end) tm.(mf{mc}){2}];
            else
                names{mc} = sm(up(1)+1:end);
            end
        end
    end

    % add AFT methods
    tm = xffsngl.FM.aft;
    mf = fieldnames(tm);

    % iterate
    nnames = numel(names);
    names{nnames + numel(mf), 1} = '';
    for mc = 1:numel(mf)
        sm = tm.(mf{mc}){1};
        if pril
            names{nnames + mc} = [sm(5:end) tm.(mf{mc}){2}];
        else
            names{nnames + mc} = sm(5:end);
        end
    end

    % add further
    for oc = 2:numel(xo)
        names = intersect(names, methods(xo(oc)));
    end
    xffsngl.LAST = xo;
end

% properties (tab-press, etc.) is combination
function names = properties(xo)
    names = [fieldnames(xo); methods(xo)];
end

% isfield
function isf = isfield(xo, F)

    % if S is not 1x1 or F is a cell array error out
    if numel(xo) ~= 1 && iscell(F)
        error('neuroelf:xff:invalidSyntax', ...
            'fieldnames cannot be called with multiple object and fields.');
    end

    % sanity check
    if isempty(F) || (~ischar(F) && ~iscell(F))
        error('neuroelf:xff:badArgument', ...
            'F must be either a single fieldname or a list of fieldnames.');
    end
    if ischar(F)
        F = {F(:)'};
    end
    F = F(:)';
    for fc = 1:numel(F)
        if ~ischar(F{fc}) || isempty(F{fc})
            error('neuroelf:xff:badArgument', 'Invalid fieldname.');
        end
        F{fc} = F{fc}(:)';
    end

    % single object
    if numel(xo) == 1

        % create output
        isf = false(size(F));

        % return original names
        names = fieldnames(xo.C);

        % check for fieldnames
        for fc = 1:numel(F)
            isf(fc) = any(strcmp(names, F{fc}));
        end

    % multiple objects
    else

        % create output
        isf = false(1, numel(xo));

        % iterate across objects
        for oc = 1:numel(xo)

            % test
            isf(oc) = any(strcmp(fieldnames(xo(oc).C), F{1}));
        end
    end
end

% filename and filetype
function n = filename(xo)
    try
        if numel(xo) == 1
            n = xo.F;
        else
            n = cell(size(xo));
            for fc = 1:numel(n)
                n{fc} = xo(fc).F;
            end
        end
    catch xfferror
        rethrow(xfferror);
    end
end
function t = filetype(xo)
    try
        if numel(xo) == 1
            t = xo.S.Extensions{1};
        else
            t = cell(size(xo));
            for fc = 1:numel(t)
                t{fc} = xo(fc).S.Extensions{1};
            end
        end
    catch xfferror
        rethrow(xfferror);
    end
end

% get field or entire (super) contents
function fc = get(xo, F)
    if nargin ~= 2 || (~ischar(F) && (~iscell(F) || numel(xo) > 1))
        error('neuroelf:xff:badArgument', 'Invalid input argument.');
    end

    % pass on to default subsref (of struct C)
    fc = xo.C.(F(:)');
end
function c = getcont(xo)
    if numel(xo) == 1
        c = xo.C;
    else
        c = {xo.C};
    end
end
function c = getscont(xo)
    if numel(xo) == 1
        c = struct('C', xo.C, 'F', xo.F, 'H', xo.H, 'L', xo.L, 'S', xo.S);
    else
        c = struct('C', {xo.C}, 'F', {xo.F}, 'H', {xo.H}, 'L', {xo.L}, 'S', {xo.S});
    end
end

% get handles
function h = handles(xo)
    if numel(xo) == 1
        h = xo.H;
    else
        h = {xo.H};
    end
end

% help
function help(xo, varargin)
    if sfile.L(1) == 'X'
        disp(root_Help(xo));
    else
        disp(aft_Help(xo, varargin{:}));
    end
end

% object ID
function id = objectid(xo)
    if numel(xo) == 1
        id = xo.L;
    else
        id = {xo.L};
    end
end

% saveobj
function c = saveobj(xo)

    % factory
    global xffsngl;

    % argument check
    if nargin ~= 1 || numel(xo) ~= 1 || ~xffisobject(xo, true)
        error('neuroelf:xff:BadArgument', 'Invalid object or mismatching fieldnames.');
    end

    % depending on behavior
    b = lower(xffsngl.CONF.settings.Behavior.SaveObj);
    if isempty(xo.F)
        b = 'data';
    end
    switch (lower(b(1)))

        % saving as data
        case 'd'

            % generate data struct
            c = struct('C', struct('BEH', 'data', 'DAT', xo.C, ...
                'EXT', xo.S.Extensions{1}, 'TYP', xo.S.FFTYPE), ...
                'F', '', 'H', struct, 'L', '', 'S', struct);

        % saving as filename
        otherwise

            % generate filename base struct
            c = struct('C', struct('BEH', 'filename', 'DAT', xo.F, ...
                'EXT', xo.S.Extensions{1}, 'TYP', xo.S.FFTYPE), ...
                'F', '', 'H', struct, 'L', '', 'S', struct);
    end
end

% set, setcont
function xo = set(xo, varargin)

    % requires second and valid argument
    if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true) || ...
       ~ischar(varargin{1}) || isempty(varargin{1}) || ~isfield(xo.C, varargin{1}(:)')
        error('neuroelf:xff:badArgument', 'Invalid set-assignment.');
    end

    % try-make all assignments
    try
        bc = xo.C;
        for ac = 1:2:(nargin-1)
            bc.(varargin{ac}(:)') = varargin{ac+1};
        end
        xo.C = bc;
    catch xfferror
        error('neuroelf:xff:assignmentFailed', 'Assignment failed: %s.', xfferror.message);
    end
end
function xo = setcont(xo, bc)

    % requires second and valid argument
    if nargin ~= 2 || numel(xo) ~= 1 || numel(bc) ~= 1 || ~xffisobject(xo, true) || ...
       ~isstruct(bc) || xo.L(1) == 'X'
        error('neuroelf:xff:badArgument', 'Invalid object for assignment.');
    end

    % allow REMAININGCONTENT field
    xocf = fieldnames(xo.C);
    if any(strcmp(xocf, 'REMAININGCONTENT'))
        xocf(strcmp(xocf, 'REMAININGCONTENT')) = [];
    end
    if isfield(bc, 'REMAININGCONTENT')
        bcr = bc.REMAININGCONTENT;
        bc = rmfield(bc, 'REMAININGCONTENT');
    else
        bcr = [];
    end

    % get content first
    if numel(xocf) ~= numel(fieldnames(bc)) || ~all(strcmp(xocf, fieldnames(bc)))
        error('neuroelf:xff:badArgument', 'Mismatching fieldnames, content not set.');
    end

    % assignment (but keep original ID!)
    try
        xo.C = bc;
        if ~isempty(bcr)
            xo.C.REMAININGCONTENT = bcr;
        end
        xo.C.RunTimeVars.xffID = xo.L;
    catch xfferror
        rethrow(xfferror);
    end
end

% end of public, visible methods
end


% static methods
methods (Access = 'public', Hidden = true, Static = true)

% loadobj overload
function xo = loadobj(c)

    % make sure factory is loaded
    global xffsngl;
    if isempty(xffsngl)
        xff();
    end

    % argument check
    if nargin ~= 1 || numel(c) ~= 1 || ~isstruct(c) || ...
       ~isfield(c, 'C') || numel(c.C) ~= 1 || ~isstruct(c.C) || ...
       ~isfield(c.C, 'BEH') || ~isfield(c.C, 'DAT') || ...
       ~isfield(c.C, 'EXT') || ~isfield(c.C, 'TYP')
        error('neuroelf:xff:badArgument', 'Invalidly stored object.');
    end

    % root object
    if strcmpi(c.C.TYP, 'root')

        % simply return current root object
        xo = xffsngl.OBJS{1, 4};
        return;
    end

    % depending on behavior
    b = c.C.BEH;
    switch (b)

        % loading from data
        case 'data'

            % create new object
            xo = xff(['new:' c.C.EXT]);

            % and set content
            xo.C = c.C.DAT;
            xo.C.RunTimeVars.xffID = xo.L;

        % load from file
        case 'filename'

            % try loading
            try
                xo = xff(c.C.DAT, c.C.EXT);
            catch xfferror
                warning('neuroelf:xff:loadError', ...
                    'Error loading saved object from ''%s'' (using default): %s.', ...
                    c.C.DAT, xfferror.message);
                xo = xff(['new:' hstr.EXT]);
            end
    end
end


% end of static methods
end

% end of classdef
end

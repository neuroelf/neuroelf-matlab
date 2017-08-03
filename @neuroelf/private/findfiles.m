function [varargout] = findfiles(varargin)
%FINDFILES  locate files and directories, comparable to Linux/UNIX find.
%   FILES = FINDFILES(PATTERN) searches for files matching PATTERN in the
%   present working directory (and, by default, all sub-directories). The
%   pattern may contain the asterisk (*) for {0,} repetitions of any
%   character, as well as the question mark (?) for exactly one occurrence
%   of any character. With this syntax, the filenames returned in FILES are
%   absolute pathnames, including the PWD (see 'relative=#' option below).
%
%   FILES = FINDFILES(STARTFOLDER, PATTERN) searches in the given
%   start folder instead; using this syntax, both the pattern as well as
%   the start folder can also be a set of strings in a cell array, such as
%   FILES = FINDFILES(PWD, {'*.m', '*.mex*'}) -- or
%   FILES = FINDFILES({'/path/to/folder1', '/path/to/folder2'}, '*.m')
%
%   FILES = FINDFILES(STARTFOLDER, PATTERN, ...) allows to further specify
%   which files (or directories) are to be located with the following
%   'OPTION=VALUE' strings as additional inputs ('-' options combinable!)
%
%   'depth=X'    only locate files at exactly depth X from the start folder
%                whereas the start folder is considered depth=1; as an
%                alternative shortcut, pass in '-dX'
%   'dirs=1'     locate directories instead of files (alternative '-D')
%   'filesize=X' locate only files of a specific size (alternative '-sX')
%   'maxage=X'   restrict to files that have been changed at most X seconds
%                before; e.g. to search for files modified within the last
%                hour use 'maxage=3600' (alternative '-AX')
%   'maxdepth=X' only locate files up to depth X from start folder
%   'maxsize=X'  only locate files up to filesize X (alternative '-MX')
%   'minage=X'   restrict to files that have been changed at least X
%                seconds before (alternative '-aX')
%   'mindepth=X' only locate files with at least depth X from start folder
%   'minsize=X'  only locate files of filesize X or greater (alt. -'mX')
%   'oneperdir'  locate only the first matching file per folder
%   'relative=#' replace the start folder with an arbitrary string
%
%   FILES = FINDFILES(PWD, '*.m', '-d2m4096') would locate files with an
%   extension of '.m' in folder depth 2 (one below the PWD) and with at
%   least 4096 bytes filesize.
%
%   Alternative to a set of strings, options can also be passed into this
%   function as a struct:
%
%   FILES = FINDFILES(STARTFOLDER, PATTERN, OPTIONS) whereas each optional
%   field must be the name of one of the (long) options, such as in
%   FILES = FINDFILES(PWD, '*.m', STRUCT('mindepth', 3))
%
%   [FILES, NUMBER] = FINDFILES(...) also returns the number of files.
%
%   [FILES, NUMBER, SIZES] = FINDFILES(...) also returns the sizes for
%   each of the found files. For folders, it returns the number of entries
%   per folder.
%
%   [FILES, NUMBER, SIZES, DATES] = FINDFILES(...) also returns a date (as
%   a string) for each of the found files (or directories).
%
%   Note: the minage/maxage feature only fully works when the system
%   returns English-style month in calls to dir. i.e. under Linux, set the
%   LANG environmental setting to 'en_US' before starting up MATLAB.
%
%   See also: DIR.

% Version:  v1.1
% Build:    16060814
% Date:     Jun-08 2016, 2:49 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% number of inputs
if nargin == 0
    varargin{1} = '*';
end

% enough arguments ?
if ((~ischar(varargin{1}) || isempty(varargin{1})) && ...
    (~iscell(varargin{1}) || isempty(varargin{1}) || ...
     ~ischar(varargin{1}{1}) || isempty(varargin{1}{1})))
    error('neuroelf:findfiles:invalidInputType', 'Invalid input type: %s.', class(varargin{1}));
end

% shortcut variable
fsep = filesep;

% for single argument
if nargin < 2

    % check for full path
    [p, f, x] = fileparts(varargin{1}(:)');

    % path specification
    if isempty(p) || strcmp(p, '.')
        varargin{1} = pwd;
    else
        varargin{1} = p;
    end

    % extend varargin
    varargin{2} = [f, x];

% allow local path with options
elseif ischar(varargin{1}) && any(varargin{1} == '*') && ...
    ischar(varargin{2}) && any(varargin{2} == '=')

    % pass on
    [varargout{1:max(1, nargout)}] = findfiles(pwd, varargin{:});
    
    % \ filesep
    if fsep == '\' && ~isempty(varargout)
        varargout{1} = strrep(varargout{1}, '/', fsep);
    end
    return;
end

% sanity checks: startfolder
startfolder = varargin{1};
if ischar(startfolder) && ~isempty(startfolder)

    % does startfolder contain an asterisk?
    if any(startfolder == '?' | startfolder == '*')

        % then put startfolder into a cell array to treat this case
        [varargout{1:max(1, nargout)}] = findfiles({startfolder}, varargin{2:end});
        
        % \ filesep
        if fsep == '\' && ~isempty(varargout)
            varargout{1} = strrep(varargout{1}, '/', fsep);
        end
        return;
    end

% startfolder is a list of folders
elseif iscell(startfolder) && ~isempty(startfolder)

    % generate expanded list
    nstartfolder = cell(0, 1);

    % iterate over items
    for nelem = 1:numel(startfolder)

        % must be non-empty char
        if ~ischar(startfolder{nelem}) || isempty(startfolder{nelem})
            error('neuroelf:findfiles:badArgument', 'Bad startfolder argument.');
        end

        % expand pattern?
        if ~any(startfolder{nelem} == '?' | startfolder{nelem} == '*')

            % put into final list
            nstartfolder{end+1, 1} = startfolder{nelem};

        % or look up matching folders
        else

            % split along possible path separators
            [pparts, cparts] = splitbysep(startfolder{nelem}, '/\');

            % remove empty parts (double filesep) 
            if any('/\:' == startfolder{nelem}(1))
                pparts = [{''}, pparts];
                cparts = cparts + 1;
            end

            % look for first pattern
            for cpart = 1:cparts
                if any(pparts{cpart} == '?' | pparts{cpart} == '*')
                    break;
                end
            end

            % if the pattern occurs in first part, look in current dir
            if cpart == 1
                pfolders = findfiles('.', pparts{1}, struct('dirs', 1, 'depth', 1));

            % otherwise glue first parts together and look up matches
            else
                spart = gluebysep(pparts(1:(cpart-1)), fsep);
                if exist(spart, 'dir') > 0
                    pfolders = findfiles(spart, pparts{cpart}, ...
                        struct('dirs', 1, 'depth', 1));
                else
                    pfolders = cell(0, 1);
                end
            end

            % the pattern was not in last part
            if cpart < cparts

                % put remaining parts back on
                for ppart = 1:numel(pfolders)
                    pfolders{ppart} = [pfolders{ppart}, fsep, ...
                        gluebysep(pparts((cpart+1):end), fsep)];
                end
            end

            % put results at end of list
            nstartfolder = [nstartfolder; pfolders];
        end
    end

    % we start with no files found
    varargout{1} = cell(0, 1);
    varargout{2} = 0;
    if nargout > 2
        varargout{3} = zeros(0, 1);
        if nargout > 3
            varargout{4} = cell(0, 1);
        end
    end

    % allow different file separators
    if fsep == '\'
        nstartfolder = strrep(nstartfolder, '\', '/');
    end

    % for each folder in list
    for nelem = 1:numel(nstartfolder)

        % if there (iteratively) remains a pattern char, redo this
        if any(nstartfolder{nelem} == '?' | nstartfolder{nelem} == '*')
            [nff, nfn, nfs, nfd] = findfiles(nstartfolder(nelem), varargin{2:end});
            varargout{1} = [varargout{1}; nff];
            varargout{2} = varargout{2} + nfn;
            if nargout > 2
                varargout{3} = [varargout{3}; nfs];
                if nargout > 3
                    varargout{4} = [varargout{4}; nfd];
                end
            end

        % otherwise get files in this folder and put at end of array
        elseif exist(nstartfolder{nelem}, 'dir') == 7
            [nff, nfn, nfs, nfd] = findfiles(nstartfolder{nelem}, varargin{2:end});
            varargout{1} = [varargout{1}; nff];
            varargout{2} = varargout{2} + nfn;
            if nargout > 2
                varargout{3} = [varargout{3}; nfs];
                if nargout > 3
                    varargout{4} = [varargout{4}; nfd];
                end
            end
        end
    end

    % \ filesep
    if fsep == '\' && ~isempty(varargout)
        varargout{1} = strrep(varargout{1}, '/', fsep);
    end
    return;

% illegal first argument
else
    error('neuroelf:findfiles:badArgument', 'Bad startfolder argument.');
end

% we're now going for the single folder case and startfolder is OK
varargout{1} = cell(0, 1);
if nargout > 1
    varargout{2} = 0;
    if nargout > 2
        varargout{3} = zeros(0, 1);
        if nargout > 3
            varargout{4} = cell(0, 1);
        end
    end
end

% startfolder exists?
if exist(startfolder, 'dir') ~= 7
    return;
end

% append missing fsep if needed
if ~any(startfolder(end) == '\/')
    startfolder = [startfolder '/'];
end

% default is cell, otherwise
patterns = varargin{2};
if ~iscell(patterns)

    % if is character, put into cell
    if ischar(patterns)
        patterns = {patterns};

    % otherwise bail out
    else
        error('neuroelf:findfiles:badArgument', 'Patterns must be char or cell.');
    end
end

% check each pattern
for count = 1:numel(patterns)

    % only accept chars
    if ~ischar(patterns{count}) || isempty(patterns{count})
        patterns{count} = '*';
    end
end
startfolder = strrep(startfolder, '\', '/');
patterns = strrep(patterns, '\', '/');
patterns = unique(patterns(:));

% option argument parsing, default options
if nargin < 3
    opt.dirs     = 0;
    opt.filesize = -1;
    opt.maxage   = -1;
    opt.maxdepth = 0;
    opt.maxsize  = -1;
    opt.minage   = -1;
    opt.mindepth = 0;
    opt.minsize  = -1;
    opt.oneperdir = 0;
    opt.relative = 0;
    opt.rfolder  = startfolder;

% parse options
else
    opt = varargin{3};

    % non-struct options
    if ~isstruct(opt)

        % yet start with default struct
        opt = struct;
        opt.dirs     = 0;
        opt.filesize = -1;
        opt.maxage   = -1;
        opt.maxdepth = 0;
        opt.maxsize  = -1;
        opt.minage   = -1;
        opt.mindepth = 0;
        opt.minsize  = -1;
        opt.oneperdir = 0;
        opt.relative = 0;
        opt.rfolder  = startfolder;

        % parse all arguments
        for acount = 3:nargin

            % char option
            if ischar(varargin{acount})

                % special case: -a#A#d#Dors#
                % (minage, maxage, depth, dirs, oneperdir, relative=, size)
                if ~isempty(varargin{acount}) && varargin{acount}(1) == '-'
                    optarg = varargin{acount}(2:end);
                    while ~isempty(optarg)
                        switch (optarg(1))
                            case {'a'}
                                if numel(optarg) > 1 && optarg(2) >= '0' && optarg(2) <= '9'
                                    optpos = find(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos(1));
                                        optarg(1:optpos(1)) = [];
                                    end
                                    opt.minage = str2double(oval);
                                else
                                    warning('neuroelf:BadOption', ...
                                        'Option -a (minage) requires numeric input.');
                                    optarg(1) = [];
                                end
                            case {'A'}
                                if numel(optarg) > 1 && optarg(2) >= '0' && optarg(2) <= '9'
                                    optpos = find(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos(1));
                                        optarg(1:optpos(1)) = [];
                                    end
                                    opt.maxage = str2double(oval);
                                else
                                    warning('neuroelf:BadOption', ...
                                        'Option -A (maxage) requires numeric input.');
                                    optarg(1) = [];
                                end
                            case {'d'}
                                if numel(optarg) > 1 && optarg(2) >= '0' && optarg(2) <= '9'
                                    optpos = find(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos(1));
                                        optarg(1:optpos(1)) = [];
                                    end
                                    opt.maxdepth = str2double(oval);
                                    opt.mindepth = opt.maxdepth;
                                else
                                    warning('neuroelf:BadOption', ...
                                        'Option -d (depth) requires numeric input.');
                                    optarg(1) = [];
                                end
                            case {'D'}
                                opt.dirs = 1;
                                optarg(1) = [];
                            case {'m'}
                                if numel(optarg) > 1 && optarg(2) >= '0' && optarg(2) <= '9'
                                    optpos = find(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos(1));
                                        optarg(1:optpos(1)) = [];
                                    end
                                    opt.minsize = str2double(oval);
                                else
                                    warning('neuroelf:BadOption', ...
                                        'Option -m (minsize) requires numeric input.');
                                    optarg(1) = [];
                                end
                            case {'M'}
                                if numel(optarg) > 1 && optarg(2) >= '0' && optarg(2) <= '9'
                                    optpos = find(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos(1));
                                        optarg(1:optpos(1)) = [];
                                    end
                                    opt.maxsize = str2double(oval);
                                else
                                    warning('neuroelf:BadOption', ...
                                        'Option -M (maxsize) requires numeric input.');
                                    optarg(1) = [];
                                end
                            case {'o'}
                                opt.oneperdir = 1;
                                optarg(1) = [];
                            case {'r'}
                                opt.relative = 1;
                                opt.rfolder = '';
                                optarg(1) = [];
                            case {'s'}
                                if numel(optarg) > 1 && optarg(2) >= '0' && optarg(2) <= '9'
                                    optpos = find(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos(1));
                                        optarg(1:optpos(1)) = [];
                                    end
                                    opt.filesize = str2double(oval);
                                else
                                    warning('neuroelf:BadOption', ...
                                        'Option -s (filesize) requires numeric input.');
                                    optarg(1) = [];
                                end
                            otherwise
                                warning('neuroelf:BadOption', ...
                                    'Unknown findfiles option: %s.', optarg(1));
                                optarg(1) =[];
                        end
                    end
                    continue;
                end

                % get argument name and option value
                argnv = find(varargin{acount} == '=');
                if isempty(argnv)
                    oname = varargin{acount};
                    oval = '';
                else
                    oname = varargin{acount}(1:argnv(1)-1);
                    oval = varargin{acount}(argnv(1)+1:end);
                end

                % only accept known arguments
                switch lower(oname)

                    % option: depth (min and max)
                    case {'depth'}
                        if str2double(oval) >= 0
                            opt.maxdepth = str2double(oval);
                            opt.mindepth = str2double(oval);
                        else
                            opt.maxdepth = 0;
                            opt.mindepth = 0;
                        end

                    % option: dirs, set lookup type
                    case {'dirs'}
                        oval = str2double(oval);
                        if oval == 0
                            opt.dirs = 0;
                        else
                            opt.dirs = 1;
                        end

                    % option: filesize
                    case {'filesize'}
                        if str2double(oval) >= 0
                            opt.filesize = fix(str2double(oval));
                        else
                            opt.filesize = -1;
                        end

                    % option: maxage
                    case {'maxage'}
                        if str2double(oval) >= 0
                            opt.maxage = fix(str2double(oval));
                        else
                            opt.maxage = -1;
                        end

                    % option: maxdepth
                    case {'maxdepth'}
                        if str2double(oval) >= 0
                            opt.maxdepth = fix(str2double(oval));
                        else
                            opt.maxdepth = 0;
                        end

                    % option: maxsize
                    case {'maxsize'}
                        if str2double(oval) >= 0
                            opt.maxsize = fix(str2double(oval));
                        else
                            opt.maxsize = -1;
                        end

                    % option: minage
                    case {'minage'}
                        if str2double(oval) >= 0
                            opt.minage = fix(str2double(oval));
                        else
                            opt.minage = -1;
                        end

                    % option: mindepth
                    case {'mindepth'}
                        if str2double(oval) >= 0
                            opt.mindepth = fix(str2double(oval));
                        else
                            opt.mindepth = 0;
                        end

                    % option: minsize
                    case {'minsize'}
                        if str2double(oval) >= 0
                            opt.minsize = fix(str2double(oval));
                        else
                            opt.minsize = -1;
                        end

                    % option: oneperdir
                    case {'oneperdir'}
                        oval = str2double(oval);
                        if oval == 0
                            opt.oneperdir = 0;
                        else
                            opt.oneperdir = 1;
                        end

                    % option: relative
                    case 'relative'
                        noval = str2double(oval);
                        if ~isnan(noval)
                            opt.relative = noval;
                            if noval < 1
                                opt.rfolder = startfolder;
                            else
                                opt.rfolder = ['.' fsep];
                            end
                        else
                            opt.relative = 1;
                            opt.rfolder = oval;
                        end
                end
            end
        end

    % struct option argument
    else

        % make sure options are present
        if ~isfield(opt, 'dirs')
            opt.dirs = 0;
        end
        if ~isfield(opt, 'filesize')
            opt.filesize = -1;
        end
        if ~isfield(opt, 'maxage')
            opt.maxage = -1;
        end
        if ~isfield(opt, 'maxdepth')
            if isfield(opt, 'depth')
                opt.maxdepth = opt.depth;
            else
                opt.maxdepth = 0;
            end
        end
        if ~isfield(opt, 'maxsize')
            opt.maxsize = -1;
        end
        if ~isfield(opt, 'minage')
            opt.minage = -1;
        end
        if ~isfield(opt, 'mindepth')
            if isfield(opt, 'depth')
                opt.mindepth = opt.depth;
            else
                opt.mindepth = 0;
            end
        end
        if ~isfield(opt, 'minsize')
            opt.minsize = -1;
        end
        if ~isfield(opt, 'oneperdir')
            opt.oneperdir = 0;
        end
        if  isfield(opt, 'rfolder')
            opt = rmfield(opt, 'rfolder');
        end
        if ~isfield(opt, 'relative')
            opt.relative = 0;
            opt.rfolder = startfolder;
        else
            if ischar(opt.relative)
                opt.rfolder = opt.relative;
                opt.relative = 1;
            else
                if double(opt.relative) >= 1
                    opt.rfolder = ['.' fsep];
                    opt.relative = 1;
                else
                    opt.rfolder = startfolder;
                    opt.relative = 0;
                end
            end
        end
    end
end

% check option types
if opt.dirs ~= 0
    opt.dirs = 1;
end
if ~isa(opt.filesize, 'double')
    opt.filesize = -1;
end
if ~isa(opt.maxage, 'double')
    opt.maxage = -1;
end
if ~isa(opt.maxdepth, 'double')
    opt.maxdepth = 0;
end
if ~isa(opt.maxsize, 'double')
    opt.maxsize = -1;
end
if ~isa(opt.minage, 'double')
    opt.minage = -1;
end
if ~isa(opt.mindepth, 'double')
    opt.mindepth = 0;
end
if ~isa(opt.minsize, 'double')
    opt.minsize = -1;
end
if opt.oneperdir ~= 0
    opt.oneperdir = 1;
end
if opt.relative ~= 0
    opt.relative = 1;
else
    opt.rfolder = startfolder;
end

% calculate age here
opt.maxage = opt.maxage / 86400;
if opt.maxage < 0
    opt.maxage = -1;
end
opt.minage = opt.minage / 86400;
if opt.minage < 0
    opt.minage = -1;
end

% make call for files
if opt.dirs == 0
    [varargout{1:max(1, nargout)}] = findsubfiles(startfolder, patterns, 1, ...
        opt.filesize, opt.mindepth, opt.maxdepth, opt.minage, opt.maxage, ...
        opt.minsize, opt.maxsize, opt.oneperdir, opt.rfolder);

% make call for dirs
else
    [varargout{1:max(1, nargout)}] = findsubdirs(startfolder, patterns, 1, ...
        opt.mindepth, opt.maxdepth, opt.minage, opt.maxage, opt.minsize, ...
        opt.maxsize, opt.oneperdir, opt.rfolder);
end

% \ filesep
if fsep == '\' && ~isempty(varargout)
    varargout{1} = strrep(varargout{1}, '/', fsep);
end

% %%%%internal functions%%%%

% findsubfiles
function [found, nfound, sizes, dates] = findsubfiles(path, patterns, adepth, fsize, sdepth, mdepth, mnage, mxage, mnsize, mxsize, operdir, relative)

% start with zero files found
found = cell(0, 1);
nfound = 0;
sizes = zeros(0, 1);
dates = cell(0, 1);

% first, recursively handle all subfolders, if depth is still valid
if mdepth == 0 || adepth < mdepth

    % get list of files and folders, and size of list
    ilist = dir(path);
    ilist(~cat(1, ilist.isdir)) = [];

    % check items
    for count = 1:numel(ilist)

        % don't heed . and ..
        if strcmp(ilist(count).name, '.') || strcmp(ilist(count).name, '..')
            continue;
        end

        % find files in subdirs
        if nargout < 3
            filestoadd = findsubfiles([path ilist(count).name '/'], ...
                patterns, adepth + 1, fsize, sdepth, mdepth, mnage, mxage, ...
                mnsize, mxsize, operdir, [relative ilist(count).name '/']);
            sfound = numel(filestoadd);
        else
            [filestoadd, sfound, nfs, nfd] = findsubfiles([path ilist(count).name '/'], ...
                patterns, adepth + 1, fsize, sdepth, mdepth, mnage, mxage, ...
                mnsize, mxsize, operdir, [relative ilist(count).name '/']);
        end

        % if files found
        if sfound > 0
            nfoundfrm = nfound + 1;
            nfoundnew = nfound + sfound;
            found(nfoundfrm:nfoundnew, 1) = filestoadd(:);
            if nargout > 2
                sizes(nfoundfrm:nfoundnew, 1) = nfs;
                if nargout > 3
                    dates(nfoundfrm:nfoundnew, 1) = nfd;
                end
            end
            nfound = nfoundnew;
        end
    end
end

% how many files so far
nfoundb = nfound + 1;

% then, if depth is valid, add files to the output
if sdepth == 0 || sdepth <= adepth

    % only get time if needed
    if any([mnage, mxage] >= 0)
        rnow = now;
    end

    % number of patterns
    spatt = numel(patterns);
    for pcount = 1:spatt

        % no "*" pattern
        if ~any(patterns{pcount} == '*') && ~any(patterns{pcount} == '?')
            ilist = dir([path patterns{pcount} '*']);
            if isempty(ilist)
                continue;
            end
            ilist(cat(1, ilist.isdir)) = [];
            
            % patch bytes for links to -1, date to now
            if nargout > 2 || any([fsize, mnsize, mxsize, mnage, mxage] >= 0)
                iliste = cellfun('isempty', {ilist.bytes}) | cellfun('isempty', {ilist.date});
                if any(iliste)
                    [ilist(iliste).bytes] = deal(-1);
                    [ilist(iliste).date] = deal(datestr(now));
                end
            end

            % conditions
            if fsize >= 0
                ilist(cat(1, ilist.bytes ~= fsize)) = [];
            end
            if mnsize > 0
                ilist(cat(1, ilist.bytes) < mnsize) = [];
            end
            if mxsize > 0
                ilist(cat(1, ilist.bytes) > mxsize) = [];
            end
            if mnage > 0
                ilist((rnow - datenum(cat(1, ilist.date))) < mnage) = [];
            end
            if mxage > 0
                ilist((rnow - datenum(cat(1, ilist.date))) > mxage) = [];
            end
            
            % match found?
            ilistn = {ilist(:).name};
            ilistf = find(strcmp(ilistn, patterns{pcount}));
            if ~isempty(ilistf)
                nfound = nfound + 1;
                found{nfound, 1} = [relative patterns{pcount}];
                if nargout > 2
                    sizes(nfound, 1) = ilist(ilistf(1)).bytes;
                    if nargout > 3
                        dates{nfound, 1} = ilist(ilistf(1)).date;
                    end
                end
                if operdir == 1
                    return;
                end
            end
            continue;

        % find matching entries with ?
        elseif any(patterns{pcount} == '?')
            ilist = dir([path regexprep(strrep(patterns{pcount}, '?', '*'), '\*\*+', '*')]);
            ilistn = {ilist(:).name};
            ilist(cellfun('isempty', regexp(ilistn, ['^' strrep(strrep(strrep( ...
                patterns{pcount}, '.', '\.'), '?', '.'), '*', '.*') '$']))) = [];

        % and without ?
        else
            ilist = dir([path patterns{pcount}]);
        end
        if isempty(ilist)
            continue;
        end
        
        % remove dirs
        ilist(cat(1, ilist.isdir)) = [];
        
        % patch bytes for links to -1, date to now
        if nargout > 2 || any([fsize, mnsize, mxsize, mnage, mxage] >= 0)
            iliste = cellfun('isempty', {ilist.bytes}) | cellfun('isempty', {ilist.date});
            if any(iliste)
                [ilist(iliste).bytes] = deal(-1);
                [ilist(iliste).date] = deal(datestr(now));
            end
        end

        % remove by age and size
        if fsize >= 0
            ilist(cat(1, ilist.bytes) ~= fsize) = [];
        end
        if mnage > 0
            ilist((rnow - datenum(cat(1, ilist.date))) < mnage) = [];
        end
        if mxage > 0
            ilist((rnow - datenum(cat(1, ilist.date))) > mxage) = [];
        end
        if mnsize > 0
            ilist(cat(1, ilist.bytes) < mnsize) = [];
        end
        if mxsize > 0
            ilist(cat(1, ilist.bytes) > mxsize) = [];
        end

        % if only one per dir
        if operdir == 1
            if ~isempty(ilist)
                nfound = nfound + 1;
                found{nfound, 1} = [relative ilist(1).name];
                if nargout > 2
                    sizes(nfound, 1) = ilist(1).bytes;
                    if nargout > 3
                        dates{nfound, 1} = ilist(1).date;
                    end
                end
                return;
            end

        % otherwise
        else
            
            % patch names
            ilistn = regexprep({ilist.name}, '^(.*)$', [relative '$1']);

            % remove already found
            if pcount > 1 && nfound >= nfoundb
                [ilistn, ilistni] = setdiff(ilistn(:), found(nfoundb:nfound));
                if nargout > 2
                    ilist = ilist(ilistni);
                end
            end
            
            % add
            nfoundnew = nfound + numel(ilistn);
            found(nfound+1:nfoundnew, 1) = ilistn(:);
            if nargout > 2
                sizes(nfound+1:nfoundnew, 1) = cat(1, ilist.bytes);
                if nargout > 3
                    dates(nfound+1:nfoundnew, 1) = reshape({ilist.date}, numel(ilist), 1);
                end
            end
            nfound = nfoundnew;
        end
    end
end


% findsubdirs
function [found, nfound, sizes, dates] = findsubdirs(path, patterns, adepth, sdepth, mdepth, mnage, mxage, mnsize, mxsize, operdir, relative)

% start with zero dirs found
found = cell(0, 1);
nfound = 0;
sizes = zeros(0, 1);
dates = cell(0, 1);

% first, recursively handle all subfolders, if depth is still valid
if mdepth == 0 || adepth < mdepth

    % get list of files and folders, and size of list
    ilist = dir(path);
    ilist(~cat(1, ilist.isdir)) = [];

    % check items
    for count = 1:numel(ilist)

        % don't heed . and ..
        if strcmp(ilist(count).name, '.') || strcmp(ilist(count).name, '..')
            continue;
        end

        % find dirs in subdirs
        if nargout < 3
            filestoadd = findsubdirs([path ilist(count).name '/'], ...
                patterns, adepth + 1, sdepth, mdepth, mnage, mxage, ...
                mnsize, mxsize, operdir, [relative ilist(count).name '/']);
            sfound = numel(filestoadd);
        else
            [filestoadd, sfound, nfs, nfd] = findsubdirs([path ilist(count).name '/'], ...
                patterns, adepth + 1, sdepth, mdepth, mnage, mxage, ...
                mnsize, mxsize, operdir, [relative ilist(count).name '/']);
        end

        % if files found
        if sfound > 0
            nfoundfrm = nfound + 1;
            nfoundnew = nfound + sfound;
            found(nfoundfrm:nfoundnew, 1) = filestoadd(:);
            if nargout > 2
                sizes(nfoundfrm:nfoundnew, 1) = nfs;
                if nargout > 3
                    dates(nfoundfrm:nfoundnew, 1) = nfd;
                end
            end
            nfound = nfoundnew;
        end
    end
end

% how many files so far
nfoundb = nfound + 1;

% then, if depth is valid, add folders to the output
if sdepth == 0 || sdepth <= adepth

    % only get time if needed
    if any([mnage, mxage] >= 0)
        rnow = now;
    end

    % number of patterns
    spatt = numel(patterns);
    for pcount = 1:spatt

        % no "*" or "?" pattern
        if ~any(patterns{pcount} == '*') && ~any(patterns{pcount} == '?')
            if exist([path patterns{pcount}], 'dir') == 7
                nfound = nfound + 1;
                found{nfound, 1} = [relative patterns{pcount}];
                if nargout > 2
                    sizes(nfound, 1) = max(0, numel(dir([path patterns{pcount}])) - 2);
                    if nargout > 3
                        ldates = dir([path patterns{pcount} '*']);
                        ldates(~strcmp({ldates.name}, patterns{pcount})) = [];
                        if ~isempty(ldates)
                            dates{nfound, 1} = ldates(1).date;
                        else
                            dates{nfound, 1} = datestr(rnow);
                        end
                    end
                end
                if operdir == 1
                    return;
                end
            end
            continue;

        % "?" pattern/s
        elseif any(patterns{pcount} == '?')
            ilist = dir([path regexprep(strrep(patterns{pcount}, '?', '*'), '\*\*+', '*')]);
            ilistn = {ilist(:).name};
            ilist(cellfun('isempty', regexp(ilistn, ['^' strrep(strrep(strrep( ...
                patterns{pcount}, '.', '\.'), '?', '.'), '*', '.*') '$']))) = [];

        % "*" pattern/s
        else
            ilist = dir([path patterns{pcount}]);
        end
        
        % only directories
        if isempty(ilist)
            continue;
        end
        ilist(~cat(1, ilist.isdir)) = [];
        while ~isempty(ilist) && any(strcmpi(ilist(1).name, {'.', '..'}))
            ilist(1) = [];
        end
        if isempty(ilist)
            continue;
        end

        % patch bytes for links to -1, date to now
        if any([mnsize, mxsize, mnage, mxage] >= 0)
            iliste = cellfun('isempty', {ilist.bytes}) | cellfun('isempty', {ilist.date});
            if any(iliste)
                [ilist(iliste).bytes] = deal(-1);
                [ilist(iliste).date] = deal(datestr(now));
            end
        end

        % remove by age and size
        if mnage > 0
            ilist((rnow - datenum(cat(1, ilist.date))) < mnage) = [];
        end
        if mxage > 0
            ilist((rnow - datenum(cat(1, ilist.date))) > mxage) = [];
        end
        if mnsize > 0 || mxsize > 0
            ilists = zeros(numel(ilist), 1);
            for count = 1:numel(ilists)
                ilists(count) = numel(dir([path ilist(count).name]));
            end
            ilists = max(0, ilists - 2);
            if mnsize > 0
                ilist(ilists < mnsize) = [];
                ilists(ilists < mnsize) = [];
            end
            if mxsize > 0
                ilist(ilists > mxsize) = [];
            end
        end

        % if only one per dir
        if operdir == 1
            if ~isempty(ilist)

                % get next entry
                nfound = nfound + 1;
                found{nfound, 1} = [relative ilist(1).name];
                if nargout > 2
                    sizes(nfound, 1) = max(0, numel(dir([path ilist(1).name])) - 2);
                    if nargout > 3
                        dates{nfound, 1} = ilist(1).date;
                    end
                end
                return;
            end

        % otherwise check all
        else
            
            % patch names
            ilistn = regexprep({ilist.name}, '^(.*)$', [relative '$1']);

            % remove already found
            if pcount > 1 && nfound >= nfoundb
                [ilistn, ilistni] = setdiff(ilistn(:), found(nfoundb:nfound));
                if nargout > 2
                    ilist = ilist(ilistni);
                end
            end
            
            % add
            nfoundnew = nfound + numel(ilistn);
            found(nfound+1:nfoundnew, 1) = ilistn(:);
            if nargout > 2
                ilists = zeros(numel(ilist), 1);
                for count = 1:numel(ilists)
                    ilists(count) = numel(dir([path ilist(count).name]));
                end
                sizes(nfound+1:nfoundnew, 1) = max(0, ilists - 2);
                if nargout > 3
                    dates(nfound+1:nfoundnew, 1) = reshape({ilist.date}, numel(ilist), 1);
                end
            end
            nfound = nfoundnew;
        end
    end
end

% function to split along fileseps
function [sp, nsp] = splitbysep(p, s)
p = p(:)';
ps = (p == s(1));
for sc = 2:numel(s)
    ps = ps | (p == s(sc));
end
pss = sum(ps);
if pss == 0
    sp = {p};
    nsp = 1;
    return;
end
ps = [0, reshape(find(ps), 1, pss), numel(p) + 1];
pss = pss + 1;
sp = cell(1, pss);
for sc = 1:pss
    sp{sc} = p(ps(sc)+1:ps(sc+1)-1);
end
sp(cellfun('isempty', sp)) = [];
nsp = numel(sp);

% function to glue by filesep
function p = gluebysep(sp, s)
sp = regexprep(sp, '(.)$', ['$1' s]);
sp(cellfun('isempty', sp)) = {s};
p = cat(2, sp{:});

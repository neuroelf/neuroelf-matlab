function [pfnd, prel] = findfiledir(opath, spath, isafile, useff)
% findfiledir  - find a named file or directory
%
% FORMAT:       [pfnd, prel] = findfiledir(opath, spath, isafile, useff)
%
% Input fields:
%
%       opath       original path
%       spath       either string or 1xP cell array to search in
%       isafile     is it a file (true) or directory (false)
%       useff       make use of findfiles or not
%
% Output fields:
%
%       pfnd        either the found file/dir path or empty
%       prel        relative path from the matching spath

% Version:  v0.9a
% Build:    10051716
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

% argument check
if nargin < 4 || ...
   ~ischar(opath) || ...
    isempty(opath) || ...
    isempty(spath) || ...
   (~ischar(spath) && ...
    (~iscell(spath) || ...
     isempty(spath{1}))) || ...
   (~isnumeric(isafile) && ...
    ~islogical(isafile)) || ...
    isempty(isafile) || ...
   (~isnumeric(useff) && ...
    ~islogical(useff)) || ...
    isempty(useff)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
isafile = isafile(1);
if isafile
    isafile = 2;
    srcfile = 'file';
else
    isafile = 7;
    srcfile = 'dir';
end
pfnd = '';
prel = '';

% shortcut if original file is found
opath = strrep(opath(:)', '\', '/');
if exist(opath, srcfile) == isafile
    pfnd = opath;
    return;
end

% check OS
if ~ispc && ...
    opath(2) == ':'
    opath(1:2) = [];
end
if ispc && ...
    opath(2) ~= ':' && ...
    opath(1) == '/'
    [opabs{1:2}] = isabsolute(opath);
    opath = opabs{2};
end

% check whether path is absolute and get particles
opc = splittocell(opath, '/');
[opf{1:3}] = fileparts(opath);
opn = numel(opc);
opm = opn - 1;

% create good alternative paths
if ~iscell(spath)
    spath = {spath(:)'};
end
spath = spath(:)';
for sc = numel(spath):-1:1
    if isempty(spath{sc}) || ...
        exist(spath{sc}(:)', 'dir') ~= 7
        spath(sc) = [];
        continue;
    end
    spath{sc} = strrep(spath{sc}(:)', '\', '/');
end

% possibly add parents of original path
apath = cell(0, 1);
for c = opm:-1:1
    apath{end+1} = gluetostring(opc(1:c), '/');
end
spath = [spath(:)', apath(:)'];

% check spath (should never run into this error...)
if isempty(spath)
    error( ...
        'neuroelf:BadArgument', ...
        'Given spath is invalid, all non existing directories.' ...
    );
end

% get usefindfiles flag
if useff(1)
    useff = true;
    ffopts = struct;
    ffopts.oneperdir = 1;
    if isafile == 2
        ffopts.dirs = 0;
    else
        ffopts.dirs = 1;
    end
else
    useff = false;
end

% iterate over search paths
for sc = 1:numel(spath)

    % check if file/dir exists there
    tpath = [spath{sc} '/' opc{end}];
    if exist(tpath, srcfile) == isafile
        pfnd = tpath;
        break;
    end

    % check with added particles
    for pc = opm:-1:2
        tpath = [spath{sc} '/' gluetostring(opc(pc:end), '/')];
        if exist(tpath, srcfile) == isafile
            pfnd = tpath;
            prel = gluetostring(opc(pc:end - 1), '/');
            break;
        end
    end

    % only if findfiles is activated
    if useff

        % start findfiles
        ff = findfiles(spath{sc}, ['*' opf{2} '*' opf{3}], ffopts);
        if ~isempty(ff)
            pfnd = ff{1};
            prel = fileparts(pfnd);
            break;
        end
    end
end

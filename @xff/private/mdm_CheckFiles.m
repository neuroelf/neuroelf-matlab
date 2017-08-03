function xo = mdm_CheckFiles(xo, opts)
% MDM::CheckFiles  - check an MDM file's referential consistency
%
% FORMAT:       mdm.CheckFiles([options])
%
% Input fields:
%
%       options     optional 1x1 struct with fields
%        .altpath   alternative path for searching
%        .autofind  flag, try to locate matching files, default: true
%        .silent    flag, don't ask for user interaction, default: false
%
% Output fields:
%
%       mdm         checked object
%
% Using: findfiledir, findfiles.

% Version:  v1.1
% Build:    16053116
% Date:     May-31 2016, 4:25 PM EST
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

% neuroelf library
global ne_methods;
findfiledir = ne_methods.findfiledir;
findfiles   = ne_methods.findfiles;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if isempty(bc.XTC_RTC) || ~any(strcmpi(bc.TypeOfFunctionalData(:)', {'fmr', 'mtc', 'vtc'}))
    warning('neuroelf:xff:invalidObject', 'Unsupported MDM content (empty XTC_RTC).');
    return;
end
fdt = lower(bc.TypeOfFunctionalData(:)');

% check options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'altpath') || isempty(opts.altpath) || (~ischar(opts.altpath) && ...
    (~iscell(opts.altpath) || isempty(opts.altpath)))
    opts.altpath = {};
end
if ~iscell(opts.altpath)
    opts.altpath = {opts.altpath(:)'};
end
if ~isfield(opts, 'autofind') || ~islogical(opts.autofind) || isempty(opts.autofind)
    opts.autofind = true;
else
    opts.autofind = opts.autofind(1);
end
if ~isfield(opts, 'silent') || ~islogical(opts.silent) || isempty(opts.silent)
    opts.silent = false;
else
    opts.silent = opts.silent(1);
end

% get "main" directories for XTC/SDM/SSM files
rfiles = strrep(bc.XTC_RTC, '\', '/');
numstudy = size(rfiles, 1);
if fdt(1) == 'm'
    ssmpath = fileparts(rfiles{1});
else
    ssmpath = '';
end
xtcpath = fileparts(rfiles{1, end - 1});
if isempty(xtcpath)
    xtcpath = '.';
end
rtcpath = fileparts(rfiles{1, end});
if isempty(rtcpath)
    rtcpath = '.';
end
for fc = 2:numstudy
    if ~isempty(strfind(fileparts(rfiles{fc, end - 1}), xtcpath)) || ...
       (ispc && ~isempty(strfind(lower(fileparts(rfiles{fc, end - 1})), lower(xtcpath))))
        continue;
    end
    xtcpath = fileparts(xtcpath);
    if isempty(xtcpath) || strcmp(xtcpath, '.') || ...
       (ispc && numel(xtcpath) == 3 && xtcpath(2) == ':')
        xtcpath = '';
        break;
    end
end
for fc = 2:numstudy
    if ~isempty(strfind(fileparts(rfiles{fc, end}), rtcpath)) || ...
       (ispc && ~isempty(strfind(lower(fileparts(rfiles{fc, end})), lower(rtcpath))))
        continue;
    end
    rtcpath = fileparts(rtcpath);
    if isempty(rtcpath) || strcmp(rtcpath, '.') || ...
       (ispc && numel(rtcpath) == 3 && rtcpath(2) == ':')
        rtcpath = '';
        break;
    end
end
if fdt(1) == 'm'
    for fc = 2:numstudy
        if ~isempty(strfind(fileparts(rfiles{fc}), ssmpath)) || ...
           (ispc && ~isempty(strfind(lower(fileparts(rfiles{fc})), lower(ssmpath))))
            continue;
        end
        ssmpath = fileparts(ssmpath);
        if isempty(ssmpath) || strcmp(ssmpath, '.') || ...
           (ispc && numel(ssmpath) == 3 && ssmpath(2) == ':')
            ssmpath = '';
            break;
        end
    end
end

% alternative paths
altpath = opts.altpath(:)';
if ~isempty(fileparts(xo.F))
    altpath{end+1} = fileparts(xo.F);
end
altpath{end+1} = strrep(pwd, '\', '/');
for c = numel(altpath):-1:1
    if exist(altpath{c}, 'dir') ~= 7
        altpath(c) = [];
    end
end
% aps = numel(altpath);

% check SSM paths
if fdt(1) == 'm' && (exist(ssmpath, 'dir') ~= 7 || exist(rfiles{1}, 'file') ~= 2)

    % check whether subdir is found elsewhere
    [npath, nrel] = findfiledir(ssmpath, altpath, false, false);
    if isempty(npath)
        [npath, nrel] = findfiledir(rfiles{1}, altpath, true, opts.autofind);
        if ~isempty(npath)
            npath = fileparts(npath);
        end
    else
        [fp{1:3}] = fileparts(ssmpath);
        nrel = [nrel '/' fp{2} fp{3}];
    end
    if isempty(npath) && ~opts.autofind
        error('neuroelf:xff:missingFlag', ...
            'Feature autofind disabled. SSM reference path not found.');
    end
    for c = 1:numstudy
        if ~ispc
            rfiles{c} = strrep(rfiles{c}, ssmpath, npath);
        else
            rfiles{c} = strrep(lower(rfiles{c}), lower(ssmpath), npath);
        end
    end
end

% do the same for XTCs and RTCs
if ~isempty(xtcpath) && (exist(xtcpath, 'dir') ~= 7 || ...
    (~any(rfiles{1, end - 1} == '*') && exist(rfiles{1, end - 1}, 'file') ~= 2))

    % check whether subdir is found elsewhere
    [npath, nrel] = findfiledir(xtcpath, altpath, false, false);
    if isempty(npath)
        [npath, nrel] = findfiledir(rfiles{1, end - 1}, altpath, true, opts.autofind);
        if ~isempty(npath)
            npath = fileparts(npath);
        end
    else
        [fp{1:3}] = fileparts(xtcpath);
        nrel = [nrel '/' fp{2} fp{3}];
    end
    if isempty(npath) && ~opts.autofind
        error('neuroelf:xff:missingFlag', ...
            'Feature autofind disabled. XTC reference path not found.');
    end
    for c = 1:numstudy
        if ~ispc
            rfiles{c, end - 1} = strrep(rfiles{c, end - 1}, xtcpath, npath);
        else
            rfiles{c, end - 1} = strrep(lower(rfiles{c, end - 1}), lower(xtcpath), npath);
        end
    end
elseif exist(xtcpath, 'dir') ~= 7 && any(rfiles{1, end - 1} == '*')
    error('neuroelf:xff:badCombination', ...
        'File patterns (containing *) require correct path names.');
end
if ~isempty(rtcpath) && (exist(rtcpath, 'dir') ~= 7 || exist(rfiles{1, end}, 'file') ~= 2)

    % check whether subdir is found elsewhere
    [npath, nrel] = findfiledir(rtcpath, altpath, false, false);
    if isempty(npath)
        [npath, nrel] = findfiledir(rfiles{1, end}, altpath, true, opts.autofind);
        if ~isempty(npath)
            npath = fileparts(npath);
        end
    else
        [fp{1:3}] = fileparts(rtcpath);
        nrel = [nrel '/' fp{2} fp{3}];
    end
    if isempty(npath) && ~opts.autofind
        error('neuroelf:xff:missingFlag', ...
            'Feature autofind disabled. SDM reference path not found.');
    end
    for c = 1:numstudy
        if ~ispc
            rfiles{c, end} = strrep(rfiles{c, end}, rtcpath, npath);
        else
            rfiles{c, end} = strrep(lower(rfiles{c, end}), lower(rtcpath), npath);
        end
    end
end

% check if files exist
rexist = false(size(rfiles));
for fc = 1:numstudy
    for cc = 1:size(rexist, 2)
        if ~any(rfiles{fc, cc} == '*')
            rexist(fc, cc) = (exist(rfiles{fc, cc}, 'file') == 2);
        else
            [imgfp, imgfn, imgfe] = fileparts(rfiles{fc, cc});
            if isempty(imgfp)
                error('neuroelf:xff:badCombination', ...
                    'File patterns (containing *) require full paths.');
            end
            rexist(fc, cc) = ~isempty(findfiles(imgfp, [imgfn, imgfe], 'depth=1'));
        end
    end
end
if ~all(rexist(:))
    error('neuroelf:xff:fileNotFound', '%d or %d required file(s) not found.', ...
        sum(~rexist(:)), numel(rexist));
end

% store back
xo.C = bc;
xo.H.FilesChecked = true;

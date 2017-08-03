function samelengthfilenames(ext, folder, opts)
% samelengthfilenames  - ensure that all file names have the same length
%
% FORMAT:       samelengthfilenames(ext [, folder, [, opts]])
%
% Input fields:
%
%       ext         file extension to restrict (e.g. '.dcm')
%       folder      folder name to (begin) search in
%       opts        optional settings
%        .expchar   expansion character (default: '0')
%        .subdirs   logical flag, include sub-directories (default: true)
%
% No output fields.


% Version:  v1.0
% Build:    15111810
% Date:     Nov-18 2015, 10:15 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, Jochen Weber
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
if nargin < 1 || ...
   ~ischar(ext)
    ext = '.*';
end
if nargin < 2 || ...
   ~ischar(folder) || ...
    isempty(folder)
    folder = pwd;
elseif exist(folder(:)', 'dir') ~= 7
    error( ...
        'neuroelf:BadArgument', ...
        'Folder does not exist.' ...
    );
end
folder = folder(:)';
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'expchar') || ...
   ~ischar(opts.expchar) || ...
    numel(opts.expchar) ~= 1 || ...
    any(double(opts.expchar) == [0:47, 58:60, 62, 63, 91:94, 96, 123:127])
    opts.expchar = '0';
end
if ~isfield(opts, 'subdirs') || ...
   ~islogical(opts.subdirs) || ...
    numel(opts.subdirs) ~= 1
    opts.subdirs = true;
end

% filesep
fsep = filesep;

% sub-dirs
if opts.subdirs
    
    % get folders
    sfolders = dir(folder);
    sfolders = sfolders(cat(1, sfolders.isdir));
    sfolders(strcmp({sfolders.name}, '.')) = [];
    sfolders(strcmp({sfolders.name}, '..')) = [];
    if ~isempty(sfolders)
        for fc = 1:numel(sfolders)
            samelengthfilenames(ext, [folder fsep sfolders(fc).name], opts);
        end
    end
end

% get list of files
files = dir([folder fsep '*' ext]);

% ensure no sub-dirs
files(cat(1, files.isdir)) = [];

% no work left
if isempty(files)
    return;
end

% get names
names = {files.name};

% get lengths
nlens = cellfun('prodofsize', names);
mlens = max(nlens);

% no work left
if all(nlens == mlens)
    return;
end

% find first discrepancy
cnames = char(names);
cdiff = find(any(diff(double(cnames), 1, 1) ~= 0, 1));
cdiff = cdiff(1);
fdiff = cdiff - 1;

% nothing to do for longer names
lnames = (nlens == mlens);
names(lnames) = [];

% how much to expand
nlens = mlens - nlens(~lnames);

% expansion lists
elist = cell(1, max(nlens));
for fc = 1:numel(elist)
    elist{fc} = repmat(opts.expchar, 1, fc);
end

% loop
for fc = 1:numel(names)
    
    % expand filename and rename
    movefile( ...
        [folder fsep names{fc}], ...
        [folder fsep names{fc}(1:fdiff) elist{nlens(fc)} names{fc}(cdiff:end)]);
end

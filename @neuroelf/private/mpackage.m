function mpackage(varargin)
% mpackage  - MatLab Installation Package Generator
%
% FORMAT:       mpackage(destination, source [,source ...] [,opt])
%
% next to the destination filename and the source folders, you
% can optionally specify some options to enhance the installation
%
% Input field 'opt' 1x1 struct array with subfields...
%
%         .addpath  cell array naming the folders to add to the path
%         .banner   string specifying a banner text to show initially
%         .crypt    cell array with files to be crypted
%         .csyntax  if true check syntax of *.m files
%         .destlin  string specifying the default Linux target path
%  -and-  .destwin  string specifying the default Windows target path
%         .destref  string naming a dot-m-file to use as a reference
%  -and-  .destrup  number of dir levels to skip upwards on ref-file
%         .dontask  assume '-i' as parameter for executable M-file
%         .exclude  cell array with filenames to exclude
%         .finish   cell array with strings to eval after install
%         .ierrors  ignore errors and continue (default: false)
%         .include  cell array with filenames to include all the same
%         .maxage   option passed to findfiles when locating files
%         .minage   option passed to findfiles when locating files
%         .notwith  cell array with functions conflicting the package
%         .verbose  verbose packing (show info messages, default: false)
%         .postclr  do 'clear classes' after installation is complete
%         .release  minimum MATLAB major release number for install
%         .remove   list of files to remove from installation dir
%                   (useful for update packages, unimplemented)
%         .reqcond  cell array with conditions being evaled, whereas
%                   all must hold true for depackaging to procede
%         .require  cell array with filenames of prerequisites
%                   (M-file names which must be available)
%         .savev7   if true, use compressed saving on v7.x (default: true)
%         .strrep   Rx2 cell array with string replacements (compression)
%         .update   is package an update package
%
% once a package is created, you'll have an M-file and a MAT-file,
% whereas the M-file contains the install information for the MAT-file.
%
% with the current implementation, empty folders will NOT be packed!

% Version:  v1.0
% Build:    14092317
% Date:     Sep-23 2014, 5:28 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, Jochen Weber
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

% the dot-mat-file will hold a variable pkgcont, which is a struct.
% fields of pkgcont are
%
%         .autoi    automatically install (without -i switch)
%         .banner   greeting message
%         .date     date of creation (returned by now();)
%         .packs    cell array with sub-packages and file numbers
%         .pfiles   mXn cell array naming the files to pack
%         .pfcont   mXn cell array containing the files
%         .popt     copy of the opt struct on creation time
%         .reqrel   require Matlab release to run

% enough arguments ?
if nargin < 2 || ...
   ~ischar(varargin{1}) || ...
   (nargin < 3 && ...
    isstruct(varargin{end}))
    error( ...
        'neuroelf:BadArgument',...
        'Too few or bad arguments.' ...
    );
end

% get library functions
sclib = splittocell(asciiread( ...
    [fileparts(mfilename('fullpath')) filesep 'mpackage.m']), ...
    char([13, 10]), 1, 1);
nsclin = numel(sclib);
for lc = 1:nsclin
    if ~isempty(strfind(sclib{lc}, ['depackage ' 'functions']))
        break;
    end
end
if (lc + 5) > nsclin
    error( ...
        'neuroelf:InternalError', ...
        'Error locating depackaging functions.' ...
    );
end
sclib(1:lc) = [];
sclib = gluetostring(sclib, char(10));

% with or without options
if isstruct(varargin{end}) && ...
    numel(varargin{end}) == 1
    topt = varargin{end};
    folders = varargin(2:(end-1));
else
    topt = struct;
    folders = varargin(2:end);
end

% check given folders
for n = 1:numel(folders)
    if ~ischar(folders{n}) || ...
        exist(folders{n}, 'dir') ~= 7
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid folder selection' ...
        );
    end
end

% check individual options
opt = struct;
if isfield(topt, 'addpath') && ...
    iscell(topt.addpath) && ...
   ~isempty(topt.addpath) && ...
    ischar(topt.addpath{1})
    opt.addpath = topt.addpath;
else
    opt.addpath = {};
end
if isfield(topt, 'banner') && ...
    ischar(topt.banner)
    opt.banner = topt.banner(:)';
else
    opt.banner = '';
end
if isfield(topt, 'crypt') && ...
    iscell(topt.crypt) && ...
   ~isempty(topt.crypt) && ...
    ischar(topt.crypt{1}) && ...
   ~isempty(topt.crypt{1})
    opt.crypt = topt.crypt(:);
else
    opt.crypt = {};
end
if ~isfield(topt, 'csyntax') || ...
   ~islogical(topt.csyntax) || ...
    numel(topt.csyntax) ~= 1
    opt.csyntax = false;
else
    opt.csyntax = topt.csyntax;
end
if isfield(topt, 'destlin') && ...
    isfield(topt, 'destwin') && ...
    ischar(topt.destlin) && ...
    ischar(topt.destwin) && ...
   ~isempty(topt.destlin) && ...
   ~isempty(top.destwin)
    opt.destlin = topt.destlin;
    opt.destwin = topt.destwin;
end
if isfield(topt, 'destref') && ...
    isfield(topt, 'destrup') && ...
    ischar(topt.destref) && ...
   ~isempty(topt.destref) && ...
    isnumeric(topt.destrup) && ...
    numel(topt.destrup) == 1 && ...
   ~isinf(topt.destrup) && ...
   ~isnan(topt.destrup) && ...
    topt.destrup >= 0
    desttgt = which(topt.destref);
    if ~isempty(desttgt) && ...
        numel(desttgt) > 2 && ...
        strcmp('.m', desttgt((end-1):end))
        opt.destref = topt.destref;
        opt.destrup = round(topt.destrup);
    end
end
if ~isfield(topt, 'dontask') || ...
   ~islogical(topt.dontask) || ...
    numel(topt.dontask) ~= 1
    opt.dontask = false;
else
    opt.dontask = topt.dontask;
end
if isfield(topt, 'exclude') && ...
    iscell(topt.exclude) && ...
   ~isempty(topt.exclude) && ...
    ischar(topt.exclude{1})
    opt.exclude = topt.exclude;
else
    opt.exclude = {};
end
if isfield(topt, 'finish') && ...
    iscell(topt.finish) && ...
   ~isempty(topt.finish) && ...
    ischar(topt.finish{1})
    opt.finish = topt.finish;
else
    opt.finish = {};
end
if isfield(topt, 'ierrors') && ...
    islogical(topt.ierrors) && ...
    numel(topt.ierrors) == 1
    opt.ierrors = topt.ierrors;
else
    opt.ierrors = false;
end
if isfield(topt, 'include') && ...
    iscell(topt.include) && ...
   ~isempty(topt.include) && ...
   (ischar(topt.include{1}) || ...
    iscell(topt.include{1}))
    opt.include = topt.include;
else
    opt.include = {};
end
if isfield(topt, 'maxage') && ...
    isnumeric(topt.maxage) && ...
    numel(topt.maxage) == 1 && ...
    topt.maxage > 0
    opt.maxage = topt.maxage;
else
    opt.maxage = [];
end
if isfield(topt, 'minage') && ...
    isnumeric(topt.minage) && ...
    numel(topt.minage) == 1 && ...
    topt.minage >= 0
    opt.minage = topt.minage(1);
else
    opt.minage = [];
end
if isfield(topt, 'notwith') && ...
    iscell(topt.notwith) && ...
   ~isempty(topt.notwith) && ...
    ischar(topt.notwith{1})
    opt.notwith = topt.notwith;
else
    opt.notwith = {};
end
if ~isfield(topt, 'postclr') || ...
   ~islogical(topt.postclr) || ...
    numel(topt.postclr) ~= 1
    opt.postclr = true;
else
    opt.postclr = topt.postclr;
end
if isfield(topt, 'release') && ...
    isnumeric(topt.release) && ...
    numel(topt.release) == 1 && ...
    topt.release >= 0
    opt.release = topt.release;
else
    opt.release = 0;
end
if isfield(topt, 'remove') && ...
    iscell(topt.remove) && ...
   ~isempty(topt.remove) && ...
    ischar(topt.remove{1})
    opt.remove = topt.remove;
else
    opt.remove = {};
end
if isfield(topt, 'reqcond') && ...
    iscell(topt.reqcond) && ...
   ~isempty(topt.reqcond) && ...
    ischar(topt.reqcond{1})
    opt.reqcond = topt.reqcond;
else
    opt.reqcond = {};
end
if isfield(topt, 'require') && ...
    iscell(topt.require) && ...
   ~isempty(topt.require) && ...
    ischar(topt.require{1})
    opt.require = topt.require;
else
    opt.require = cell(0);
end
if ~isfield(topt, 'savev7') || ...
   ~islogical(topt.savev7) || ...
    numel(topt.savev7) ~= 1
    opt.savev7 = true;
else
    opt.savev7 = topt.savev7;
end
if isfield(topt, 'strrep') && ...
    iscell(topt.strrep) && ...
   ~isempty(topt.strrep) && ...
    size(topt.strrep, 2) == 2 && ...
    all(all(cellfun(@ischar, topt.strrep)))
    opt.strrep = topt.strrep;
else
    opt.strrep = cell(0, 2);
end
opt.ufilesep = filesep;
if ~isfield(topt, 'update') || ...
   ~islogical(topt.update) || ...
    isempty(topt.update)
    opt.update = true;
else
    opt.update = topt.update;
end
if isfield(topt, 'verbose') && ...
    islogical(topt.verbose) && ...
    numel(topt.verbose) == 1
    opt.verbose = topt.verbose;
else
    opt.verbose = false;
end

% encrypting requested files
cryptst = uint8([]);
if ~isempty(opt.crypt)

    % request password online
    passwd = input('Encryption password (leave empty to keep files plain): ', 's');
    if ~isempty(passwd)
        passwc = input('Re-enter password: ', 's');
        if ~isequal(passwd, passwc)
            return;
        end
        cryptst = mp_encrypt(uint8('CRYPTSTTOKEN'), passwd);
    else
        opt.crypt = {};
    end
end

% start filling output structure
pkgcont.autoi = opt.dontask;
pkgcont.banner = opt.banner;
pkgcont.crypted = {};
pkgcont.cryptst = cryptst;
pkgcont.date = now;
pkgcont.packs = {};
pkgcont.popt = opt;
pkgcont.reqrel = opt.release;
pkgcont.popt.strrep = opt.strrep(:, [2, 1]);

% keep track of errors
serrors = 0;
for n = 1:numel(folders)

    % check folders
    folder = folders{n};
    if folder(end) == filesep
        folder(end) = [];
    end

    % get single filename
    [qfold, package] = fileparts(folder);
    pkgcont.packs{n} = package;
    tfiles = mp_scanfolder(folder, ...
        opt.exclude, opt.include, opt.maxage, opt.minage, opt.verbose);
    tgt = 1;
    pkgcont.pfcont{n} = cell(0, 1);
    pkgcont.pfiles{n} = cell(0, 1);
    for m = 1:numel(tfiles)
        try
            pkgcont.pfiles{n}{tgt} = tfiles{m};
            pkgcont.pfcont{n}{tgt} = mp_readfile([folder tfiles{m}(2:end)]);
            for mc = 1:size(opt.strrep, 1)
                tcont = char(pkgcont.pfcont{n}{tgt}(:)');
                if ~isempty(strfind(tcont, opt.strrep{mc, 1}))
                    pkgcont.pfcont{n}{tgt} = uint8(strrep(tcont, ...
                        opt.strrep{mc, 1}, opt.strrep{mc, 2}))';
                end
            end
            for mc = 1:numel(opt.crypt)
                if ~isempty(regexpi(tfiles{m}, opt.crypt{mc}))
                    pkgcont.crypted{end+1} = tfiles{m};
                    disp(['Encrypting ' tfiles{m} '...']);
                    pkgcont.pfcont{n}{tgt} = mp_encrypt( ...
                        pkgcont.pfcont{n}{tgt}, passwd);
                    break;
                end
            end
            if opt.csyntax && ...
                numel(tfiles{m}) > 2 && ...
                strcmpi(tfiles{m}((end-1):end), '.m')
                csyntax = checksyntax([folder tfiles{m}(2:end)]);
                if ~isempty(csyntax)
                    disp(['syntax error in read file ' tfiles{m} ': ' csyntax]);
                    serrors = serrors + 1;
                elseif opt.verbose
                    disp(['read file: ' tfiles{m} ' -> syntax ok.']);
                end
            elseif opt.verbose
                disp(['read file: ' tfiles{m} ' (no M-file)']);
            end
            tgt = tgt + 1;
        catch ne_eo;
            warning('mpackage:warn', ne_eo.message);
        end
    end
end

% show syntax errors
if serrors > 0 && ...
   ~opt.ierrors
    disp('Some files contained syntax errors! Aborting...');
    return;
end

% put into package list
for n = 1:numel(pkgcont.packs)
    opt.packs{n}{1} = pkgcont.packs{n};
    opt.packs{n}{2} = numel(pkgcont.pfiles{n});
end
pkgcont.packs = opt.packs;

% get destination filename right
dest = regexprep(varargin{1}, '\.[mM][aA][tT]$', '');
[destsh{1:2}] = fileparts(dest);
destsh = destsh{2};

% save correctly
if mainver < 7 || ...
    opt.savev7

    % save normal/compressed
    save(dest, 'pkgcont');
else

    % save uncompressed
    save(dest, 'pkgcont', '-v6');
end

% write corresponding *.m file
odf = fopen([dest '.m'], 'w');
if odf < 1
    error( ...
        'neuroelf:FileWriteError', ...
        'Error writing package M-file ''%s''.', ...
        [dest '.m'] ...
    );
end

% start with function header, name, and the function call
fprintf(odf, ...
    ['func' 'tion ' destsh '(varargin)' char(10) ...
     '    mp_depackage(mp_readpkg(mfilename, varargin{:}), varargin{:});' ...
      char(10) char(10)]);

% write library as char and close file
fwrite(odf, sclib, 'char');
fclose(odf);


% library functions %


function myfiles = mp_scanfolder(foldername, excludes, includes, maxage, minage, verbose)

% basic checks
if nargin < 6 || ...
   ~ischar(foldername) || ...
   ~iscell(excludes)
    error('mpackage:mp_scanfolder', 'Calling convention misfit.');
end

% incorporate file age args
ffargs = {'*', 'relative=1'};
if ~isempty(maxage)
    ffargs{end+1} = sprintf('maxage=%.0f', maxage);
end
if ~isempty(minage)
    ffargs{end+1} = sprintf('minage=%.0f', minage);
end

% locate files
myfiles = findfiles(foldername, ffargs{:});

% remove excludes
for m = 1:numel(excludes)
    if ischar(excludes{m}) && ...
       ~isempty(excludes{m})
        if any(excludes{m} == '*')
            for n = numel(myfiles):-1:1
                if ~isempty(regexpi(myfiles{n}, excludes{m}))
                    if verbose
                        disp(sprintf('excluding %s', myfiles{n}));
                    end
                    myfiles(n) = [];
                end
            end
        else
            for n = numel(myfiles):-1:1
                if ~isempty(strfind(myfiles{n}, excludes{m}))
                    if verbose
                        disp(sprintf('excluding %s', myfiles{n}));
                    end
                    myfiles(n) = [];
                end
            end
        end
    end
end

if foldername(end) ~= filesep
    foldername(end+1) = filesep;
end

% forced includes
for m = 1:numel(includes)
    if ischar(includes{m}) && ...
       ~isempty(includes{m})

        % is file
        if exist([foldername includes{m}], 'file') == 2
            myfiles{end+1} = ['.' filesep includes{n}];

        % or folder
        elseif exist([foldername includes{m}], 'dir') == 7
            myifiles = findfiles( ...
                [foldername includes{m}], '*', ...
                ['relative=.' filesep includes{m} filesep]);
            myfiles(end+1:end+numel(myifiles)) = myifiles;
        end

    % for sub-includes as cell with folder / file pattern
    elseif iscell(includes{m}) && ...
        numel(includes{m}) > 1 && ...
        ischar(includes{m}{1}) && ...
        ischar(includes{m}{2})

        % check again and use findfiles
        if ~isempty(includes{m}{1})
            myifiles = findfiles( ...
                [foldername includes{m}{1}], includes{m}{2}, ...
                ['relative=.' filesep includes{m}{1} filesep]);
        else
            myifiles = findfiles(foldername(1:(end-1)), includes{m}{2}, ...
                'relative=');
        end
        myfiles(end+numel(myifiles)) = myifiles;
    end
end

function filecont = mp_readfile(filename)

if nargin < 1 || ...
   ~ischar(filename) || ...
   ~any([2,3,4,6] == exist(filename, 'file'))
    error('mpackage:readerror', 'Calling convention misfit.');
end

ifp = fopen(filename, 'r');
if ifp < 1
    error('mpackage:readerror', 'File not readable.');
end

filecont = fread(ifp, Inf, '*uint8');
fclose(ifp);

% this is simply a copy of textcrypt (which is forced also on binary data)
function t = mp_encrypt(t, k)
[ki, kn, ss] = mp_iki(k);
t = t(:)';
t(end+1:end+16) = uint8([103 56 80 98 116 123 84 46 47 57 112 56 109 55 120 65]);
tr = (t > 32 & t < 127);
r = double(t(tr)) - 32;
cp = 1;
nr = numel(r);
while cp <= nr
    ss = min(ss, nr - cp + 1);
    rs = r(cp:cp + ss - 1);
    rs = mod(ki(kn, rs) + ki(kn, 1:ss), 94) + 1;
    kr = mod(rs(end) + 13 * [rs(1):94, 1:rs(1)-1], 94) + 1;
    ki(kn, :) = ki(kn, kr);
    r(cp:cp + ss - 1) = rs;
    cp = cp + ss;
    knn = mod(sum(rs), 8) + 1;
    ss = mod(sum(rs), 32) + 63;
    ki(knn, :) = ki(knn, ki(kn, :));
    kn = knn;
end
t(tr) = uint8(r + 32);
t = t(:);


% depackage functions %

function mlr = mp_mlrelease

mlr = version('-release');
hasdots = find(mlr == '.');
if ~isempty(hasdots)
    mlr(hasdots:end) = [];
end
mlr = str2double(mlr);
if isempty(mlr)
    try
        mlr = ver('matlab');
        if ~isempty(mlr) && ...
            isstruct(mlr)
            mlr = mlr(1).Release;
        else
            mlr = ver('octave');
            mlr = ['r' mlr(1).Version(mlr(1).Version ~= '.')];
        end
        rpos = find(lower(mlr) == 'r');
        if isempty(rpos)
            error('Invalid release string');
        end
        mlr = mlr(rpos+1:end);
        rpos = find(mlr < 48 | mlr > 57);
        if ~isempty(rpos)
            mlr(rpos:end)=[];
        end
        mlr = str2double(mlr);
        if isempty(mlr)
            error('Invalid release string');
        end
    catch ne_eo;
        warning('mpackage:warn', ['Error detecting Matlab release number: ' ne_eo.message]);
        mlr = 1;
    end
end

function action = mp_printbanner(banner, pdate, packs, autoi, reqrel, varargin)

action = '';
disp(' ');
disp('MatLab Installation Package (created with NeuroElf v1.0)');
disp(['Package creation date: ' datestr(pdate)]);
disp('Package contents:');
disp('---------------------------------');
for nr = 1:numel(packs)
    disp(sprintf('* %s (%d files)', packs{nr}{1}, packs{nr}{2}));
end
disp('---------------------------------');
if ~isempty(banner)
    disp(banner);
    disp('---------------------------------');
end
disp(' ');

if nargin < 6 || ...
   ~ischar(varargin{1}) || ...
  (~strcmpi(varargin{1}, '-i') && ~strcmpi(varargin{1}, '-f'))
    if autoi == 1
        action = '-i';
    else
        [fn{1:2}] = fileparts(mfilename);
        fn = fn{2};
        disp(['Usage: ' fn ' <-i|-f>']);
        disp(['       ' fn ' -i : install with conditional checks']);
        disp(['       ' fn ' -f : install with forced mode (be careful!)']);
        disp(' ');
        return;
    end
else
    action = varargin{1};
end

if reqrel > 0
    disp(['Required MATLAB release: ' num2str(reqrel) '.']);
    if mp_mlrelease < reqrel && ...
       ~strcmpi(action, '-f')
        disp([' - MATLAB release insufficient for this installer.' char(10) ...
              '   Use forced installation to override!']);
        action = '';
        return;
    end
end

function pchars = mp_printline(line)
pchars = fprintf(1, strrep(strrep(line, '%', '%%'), '\', '\\'));

function mp_clearline(pchars)
fprintf(1, char(ones(1, pchars) * sprintf('\b')));

function pkgcont = mp_readpkg(filename, varargin)

if nargin < 1 || ...
   ~ischar(filename) || ...
    numel(filename) < 2
    error('mpackage:readpkg', 'Calling convention misfit.');
end

if strcmp('.m', filename((end-1):end))
    filename((end-1):end) = [];
end
filename = which(filename);
filename((end-1):end) = [];
if exist([filename '.mat'], 'file') ~= 2
    error('mpackage:readpkg', 'Requested file doesn''t exist.');
end

pkgtemp = load(filename);
if ~isfield(pkgtemp, 'pkgcont')
    error('mpackage:readpkg', 'Requested file has invalid content.');
end
pkgcont = pkgtemp.pkgcont;

action = mp_printbanner(pkgcont.banner, pkgcont.date, pkgcont.packs, ...
    pkgcont.autoi,pkgcont.reqrel,varargin{:});
if isempty(action)
    pkgcont = struct('noaction', 1);
    return;
end
pkgcont.action = action;
pkgcont.instpass = '';
for ac = 1:(nargin-2)
    if ischar(varargin{ac}) && ...
        strcmpi(varargin{ac}, '-p') && ...
        ischar(varargin{ac+1})
        pkgcont.instpass = varargin{ac+1};
    end
end

function [o_status, o_message] = mp_mkadir(dirname)

if nargin < 1 || ...
   ~ischar(dirname)
    error('mpackage:mkadir', 'Calling convention misfit.');
end

if dirname(end) == filesep
    dirname(end) = [];
end

[parent, name, ext] = fileparts(dirname);
if exist(parent, 'dir') ~= 7
    [status, message] = mp_mkadir(parent);
    if status < 1
        error('mpackage:mkadir', 'Couldn''t create necessary parent directory.');
    end
end

if ~isempty(ext)
    name = [name ext];
end
status = 0;
maxtries = 5;
while status < 1 && ...
    maxtries > 0
    maxtries = maxtries - 1;
    [status, message] = mkdir(parent, name);
end

if status < 1
    error('mpackage:mkadir', 'Couldn''t create target directory.');
end
if nargout > 0
    o_status = status;
end
if nargout > 1
    o_message = message;
end

function mp_writefile(filename, filecont)

if nargin < 2 || ...
   ~ischar(filename)
    error('mpackage:writefile', 'Calling convention misfit.');
end

tpath = fileparts(filename);
if exist(tpath, 'dir') ~= 7
    mp_mkadir(tpath);
end
ofp = fopen(filename, 'w');
if ofp < 1
    error('mpackage:writefile', 'File not writable.');
end
fwrite(ofp, filecont, 'uint8');
fclose(ofp);

function funcavail = mp_funcavail(funcname)

if nargin < 1 || ...
   ~ischar(funcname) || ...
    isempty(funcname)
    error('mpackage:funcavail', 'Calling convention misfit.');
end
h = which(funcname);
if isempty(h)
    funcavail = 0;
elseif strcmp('built-in', h)
    funcavail = 1;
else
    funcavail = 2;
end

function reffolder = mp_reffolder(tfilename, refup)

if nargin < 2 || ...
   ~ischar(tfilename) || ...
   ~isnumeric(refup)
    error('mpackage:reffolder', 'Calling convention misfit.');
end
h = which(tfilename);
if isempty(h)
    if ispc
        reffolder = 'C:\';
    else
        reffolder = '/';
    end
    return;
elseif strcmp('built-in', h)
    reffolder = fileparts(fileparts(which('mean')));
    return;
end

reffolder = fileparts(h);
for n = 1:refup
    reffolder = fileparts(reffolder);
end

function [ki, kn, ss] = mp_iki(k)
persistent tcki;
if isempty(tcki)
    tcki = double([ ...
        '"i@o><K8e''/-HyLX$URE4c|3b2j}x*v&wd%6`T^7q!thn.N#aWg]?mIz)=fk;9uplFQr1O([sYA~0SJM\_+G5B{PDC,ZV:'; ...
        'keK<h#La97d^$5xyf|VAI[2%B8;EjTX(H*bP3/Sc04N`pzg,)!n~l1v@oQqsGwJ{_''u\?.YtFW}=-ZrR6+]:DCiU>m"M&O'; ...
        '-SKw6@yIAb>#D{Qme79GsVUOc:(F&\XEapCZ0<)]Win''hHT.fBJdM$|^j"*?%Lo`P8;uNv,2t~/r1R5}Y[=gz_+k43lxq!'; ...
        '/X?%hg[sR<lK\up_"&|$yjtVA3>b`@)qO{MQ=F#G9Zc!*~Lf5,]^;z-EH10C}oiTr+m(kUN''n2xBdJW6D4e7vYPI:8Saw.'; ...
        'q@\W{wOCx.J%b5~!92zghG;6j]S"3=&`fI''LE/U}KsA|<dNM#DcnR0tV-H^T17iaom,YF_B4ylv:$Z?u+e(P>rQp*8[kX)'; ...
        ':{hiTG#Wr9myMuN!s)<?(OgKVp4kB|DS\l1jd&[Lt]^3-6''aAo_nxU"f*~X;5YP.bz8%eqZ+=@J7Ev,$RIwC>F`Q/c20H}'; ...
        '?OA,p^-`sXNz#C2.&{ME+]LvlDk|mZ%JU@86_/\*tYW)i=F[BwaPS7Toj}0KRu1x$e4gI>(;drQ"''39nbGyf:h~5HqV!c<'; ...
        '47W\[neDx;dB5MV{&qR1y-cJi]6HX`*$%ATwm^P~a<K+.u9}Ut/FEg''0rYl@N,:p>"Lz?#Zb3jf(ov)2|I8h=O!kGsC_QS'])' - 32;
end
ki = tcki;
k = mod(double(k(:)') - 33, 94) + 1;
kl = numel(k);
if kl < 16
    k = repmat(k, 1, ceil(16 / kl));
    k = k(1:16);
    kl = 16;
end
while kl > 48
    if mod(kl, 2) ~= 0
        k(end + 1) = k(1);
    end
    k = mod(k(1:2:end) + k(1:2:end) .* k(2:2:end) - k(end:-2:2), 94) + 1;
end

kp = true(1, 94);
kn = mod(sum(k), 8) + 1;
ss = mod(sum(k), 32) + 63;

for sc = 1:8
    kp(:) = true;
    k = mod(k + ki(k + 94 * mod(k * 5 + (1:kl), 8)), 94) + 1;
    kp(k) = false;
    ki = [ki(~kp, :); ki(kp, :)];
    ks = sum(k .* k);
    [ko, koi] = sort(ki(mod(ks:ks+7, 87) + 1, mod(ks, 8) + 1));
    ki = ki(:, koi);
end
ki = ki';

function t = mp_decrypt(t, k)
[ki, kn, ss] = mp_iki(k);
t = t(:)';
tr = (t > 32 & t < 127);
r = double(t(tr)) - 32;
cp = 1;
nr = numel(r);
while cp <= nr
    ss = min(ss, nr - cp + 1);
    rs = r(cp:cp + ss - 1);
    rc = mod(rs - (2 + ki(kn, 1:ss)), 94) + 1;
    knn = mod(sum(rs), 8) + 1;
    ssn = mod(sum(rs), 32) + 63;
    [dv, dc] = sort(ki(kn, :));
    kr = mod(rs(end) + 13 * [rs(1):94, 1:rs(1)-1], 94) + 1;
    ki(kn, :) = ki(kn, kr);
    rs = dc(rc);
    ki(knn, :) = ki(knn, ki(kn, :));
    kn = knn;
    r(cp:cp + ss - 1) = rs;
    cp = cp + ss;
    ss = ssn;
end
t(tr) = uint8(r + 32);
if ~isequal(t(end-15:end), ...
    uint8([103 56 80 98 116 123 84 46 47 57 112 56 109 55 120 65]))
    error('INCORRECT_PASSWORD');
end
t(end-15:end) = [];
t = t(:);

function mp_depackage(varargin)
if nargin < 1 || ...
   ~isstruct(varargin{1})
    error('mpackage:depackage', 'Calling convention misfit.');
end
pkgcontin = varargin{1};
if isfield(pkgcontin, 'noaction')
    return;
end

disp(' ');
disp('NOTE: !!  Please do not press any keys during file extraction.    !!');
disp('      !!  Due to some problem arising in MATLAB''s mkdir command,  !!');
disp('      !!  this can lead to an aborted installation!               !!');
disp(' ');

action = pkgcontin.action;
if ~strcmpi(action, '-f')

    disp(' ');
    disp('Performing initial checks before installation...');
    disp(' ');

    nr = numel(pkgcontin.popt.require);
    if nr > 0
        disp([' - required functions for this installer are: ' char(10), ...
              '( ' sprintf('%s ',pkgcontin.popt.require{:}) ')']);
        for n = 1:nr
            if mp_funcavail(pkgcontin.popt.require{n}) < 2
                disp([' -> installation cannot procede -- please install ''' ...
                    pkgcontin.popt.require{n} ''' functionality']);
                return;
            else
                disp([' -> ' pkgcontin.popt.require{n} ' found']);
            end
        end
        disp(' ');
    end

    nr = numel(pkgcontin.popt.reqcond);
    if nr > 0
        disp(' - checking required conditions for this installer...');
        for n = 1:nr
            fprintf(1, '%-64s -> ', ...
                sprintf(' -> condition (%s)',pkgcontin.popt.reqcond{n}));
            eval(['if ' pkgcontin.popt.reqcond{n} ...
                ', ct = true; else, ct = false; end'], 'ct = false;');
            if ct
                fprintf(1, 'passed\n');
            else
                fprintf(1, 'not passed!\n');
                return;
            end
        end
        disp(' ');
    end

    nr = numel(pkgcontin.popt.notwith);
    if nr > 0
        disp([' - conflicting files with this installer are: ' ...
             sprintf('%s ',pkgcontin.popt.notwith{:})]);
        for n = 1:nr
            if mp_funcavail(pkgcontin.popt.notwith{n}) > 0
                disp([' -> installation cannot procede -- ' ...
                      'conflicting M-file found: ' pkgcontin.popt.notwith{n} ...
                      ' (' which(pkgcontin.popt.notwith{n}) ')']);
                return;
            else
                disp([' -> ' pkgcontin.popt.notwith{n} ' not found (OK)']);
            end
        end
        disp(' ');
    end

else
    disp(['WARNING: this is a forced installation, ' ...
          'make sure you know what you''re doing!']);
    disp(' ');
end

tffound = false;
tfolder = '';
if ~tffound && ...
    isfield(pkgcontin.popt, 'destref') && ...
    isfield(pkgcontin.popt, 'destrup')
    disp(' - trying installation folder detection...');
    tfolder = mp_reffolder(pkgcontin.popt.destref, pkgcontin.popt.destrup);
    if numel(tfolder) > 3
        if ~pkgcontin.popt.update
            reply = input(['Do you want to install these packages into ' ...
                strrep(tfolder, '\', '\\') ' ? (y/n) '], 's');
            if strcmpi(reply, 'y')
                tffound = true;
            end
        else
            reply = input(['You''re about to update the installation in ' ...
                strrep(tfolder, '\', '\\') '.' char(10) ...
                'Do you wish to continue ? (y/n) '], 's');
            if ~strcmpi(reply, 'y')
                disp('Update aborted.');
                return;
            end
            tffound = true;
        end
    end
    disp(' ');
end
if ~tffound && ...
    isfield(pkgcontin.popt, 'destwin') && ...
    isfield(pkgcontin.popt, 'destlin')
    disp(' - using preconfigured target folder...');
    if ispc
        tfolder = pkgcontin.popt.destwin;
    else
        tfolder = pkgcontin.popt.destlin;
    end
    reply = input(['Do you want to install these packages into ' ...
        strrep(tfolder, '\', '\\') ' ? (y/n) '], 's');
    if strcmpi(reply, 'y')
        tffound = true;
    end
    disp(' ');
end

while ~tffound || ...
    numel(tfolder) < 4
    tfolder = input( ...
        'Please specify an installation location (leave blank to abort): ', 's');
    if isempty(tfolder)
        disp('aborted.');
        return;
    end
    if exist(tfolder, 'dir') ~= 7
        reply = input(['Folder ''' strrep(tfolder, '\', '\\') ...
            ''' doesn''t exist. Create it? (y/n) '], 's');
        if strcmpi(reply, 'y')
            try
                mp_mkadir(tfolder);
            catch ne_eo;
                warning('mpackage:warn', ['mkadir failed: ' ne_eo.message]);
            end
            if exist(tfolder, 'dir') ~= 7
                tfolder='';
            else
                tffound = true;
            end
        else
            tfolder = '';
        end
    else
        tffound = true;
    end
    disp(' ');
end

ufilesep = pkgcontin.popt.ufilesep;
if ~strcmp(ufilesep, filesep)
    rfilesep = true;
    tfilesep = filesep;
else
    rfilesep = false;
end

strreps = pkgcontin.popt.strrep;
xtfolders = cell(1, numel(pkgcontin.packs));
passwd = '';
if ~isempty(pkgcontin.cryptst) && ...
   ~isempty(pkgcontin.crypted) && ...
   ~isempty(pkgcontin.instpass)
    try
        decrypt = mp_decrypt(pkgcontin.cryptst, pkgcontin.instpass);
        if ~all(decrypt(:)' == uint8('CRYPTSTTOKEN'))
            disp('Wrong password.');
        else
            passwd = pkgcontin.instpass;
        end
    catch ne_eo;
        warning('mpackage:warn', ['Decryption error: ' ne_eo.message]);
    end
end
for n = 1:numel(xtfolders);
    xtfolder = [tfolder filesep pkgcontin.packs{n}{1} filesep];
    if exist(xtfolder, 'dir') == 7 && ...
       ~pkgcontin.popt.update
        reply = input([' -> the installation folder already exists! ' ...
            'do you wish to continue (y/n) '], 's');
        if ~strcmpi(reply, 'y')
            disp('installation aborted.');
            return;
        end
    end
    disp([' --> installing ' pkgcontin.packs{n}{1} ' into ' xtfolder]);
    fprintf(1, '     progress: ');
    pchars = 0;
    numfiles = length(pkgcontin.pfiles{n});
    for f = 1:numfiles
        xtfile = pkgcontin.pfiles{n}{f};
        xtcont = pkgcontin.pfcont{n}{f};
        if any(strcmp(xtfile, pkgcontin.crypted))
            if isempty(passwd)
                continue;
            end
            try
                xtcont = mp_decrypt(xtcont, passwd);
            catch ne_eo;
                passwd = '';
                if pchars > 0
                    mp_clearline(pchars);
                end
                warning('mpackage:warn', ['Yet incorrect password, skipping all further encrypted files! (' ne_eo.message ')']);
                continue;
            end
        end
        for mc = 1:size(strreps, 1)
            tcont = char(xtcont(:)');
            if ~isempty(strfind(tcont, strreps{mc, 1}))
                xtcont = uint8(strrep(tcont, strreps{mc, 1}, strreps{mc, 2}));
            end
        end
        if rfilesep
            xtfile = strrep(pkgcontin.pfiles{n}{f}, ufilesep, tfilesep);
        end
        if xtfile(1) == '.'
            xtfile(1) = [];
        end
        if xtfile(1) == filesep
            xtfile(1) = [];
        end
        mp_writefile([xtfolder xtfile], xtcont);
        if pchars > 0
            mp_clearline(pchars);
        end
        pchars = mp_printline(sprintf(' %3.0f%%, extracted: %-80s', ...
            (100*f/numfiles), strrep(xtfile, '\', '\\')));
    end
    disp(' ');
    xtfolders{n} = xtfolder(1:(end-1));
end

nr = numel(pkgcontin.popt.addpath);
addedfolders = cell(0, 1);
if nr > 0
    disp(' - temporarily altering path...');
    for n = nr:-1:1
        tx = -1;
        for t = 1:numel(pkgcontin.packs)
            if strcmpi(pkgcontin.packs{t}{1}, pkgcontin.popt.addpath{n})
                tx = t;
                break;
            end
        end
        if tx > 0
            disp([' -> prepending installation dir of ' ...
                pkgcontin.packs{tx}{1} ' (' xtfolders{tx} ') to path...']);
            addpath(xtfolders{tx}, '-begin');
            addedfolders{end+1} = xtfolders{tx};
        end
    end
    rehash path;
    try
        rehash toolboxcache;
    catch ne_eo;
        warning('mpackage:warn', ['Error rehashing toolbox cache: ' ne_eo.message]);
    end
    try
        savepath;
    catch ne_eo;
        warning('mpackage:warn', ['Error saving path: ' ne_eo.message]);
    end
end

if ~isempty(addedfolders)
    hline = char(45 * ones(1, 72));
    disp([char(10) hline]);
    disp('-  NOTE: Please ensure these folder(s) are in MatLab''s path permanently !! -');
    disp(hline);
    for n = 1:numel(addedfolders)
        disp([' *  ' addedfolders{n}]);
    end
    disp(hline);
    disp(' ');
end

nr = numel(pkgcontin.popt.finish);
if nr > 0
    disp(' - running finishing functions for this installer...');
    for n = 1:nr
        if iscell(pkgcontin.popt.finish{n})
            if length(pkgcontin.popt.finish{n}) > 1 && ...
                ischar(pkgcontin.popt.finish{n}{2}) && ...
               ~isempty(pkgcontin.popt.finish{n}{2})
                disp([' -> ' pkgcontin.popt.finish{n}{2}]);
            end
            try
                eval(pkgcontin.popt.finish{n}{1});
            catch ne_eo;
                disp(['    failure: ' ne_eo.message]);
            end
        else
            try
                eval(pkgcontin.popt.finish{n});
            catch ne_eo;
                disp(['    failure: ' ne_eo.message]);
            end
        end
        if exist('pkgcontin', 'var') ~= 1
            break;
        end
    end
    disp(' ');
end

if exist('addedfolders', 'var') == 1 && ...
   ~isempty(addedfolders)
    hline = char(45 * ones(1, 72));
    disp([char(10) hline]);
    disp('- AGAIN: Please add these folder(s) to your MatLab path permanently !! -');
    disp(hline);
    for n = 1:numel(addedfolders)
        disp([' *  ' addedfolders{n}]);
    end
    disp(hline);
    disp(' ');
end

disp('Installation completed!');

if exist('pkgcontin', 'var') ~= 1 || ...
    pkgcontin.popt.postclr
    try
        try
            delete(findobj('type', 'figure'));
        catch ne_eo;
            warning('mpackage:warn', ['Error removing figures: ' ne_eo.message]);
        end
        try
            fclose('all');
        catch ne_eo;
            warning('mpackage:warn', ['Error closing open files: ' ne_eo.message]);
        end
        clear persistent;
        clear global;
        clear functions;
        clear classes;
    catch ne_eo;
        warning('mpackage:warn', ['Error clearing workspaces: ' ne_eo.message]);
    end
    clear all;
    try
        clear classes;
    catch ne_eo;
        warning('mpackage:warn', ['Error clearing classes: ' ne_eo.message]);
    end
end

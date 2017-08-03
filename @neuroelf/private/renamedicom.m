function p = renamedicom(dcmfile, pattern, action)
% renamedicom  - rename a dicom file according to a pattern
%
% FORMAT:       [p = ] renamedicom(dcmfile, pattern [, action])
%       Or             renamedicom(dcmpack, 'target/%', 'unpack');
%
% Input fields:
%
%       dcmfile     dicom filename(s)
%       pattern     string containing keys, default:
%                   'Anon-%04k00200010-%04k00200011-%05k00200013.dcm'
%                   alternatively, parts of the original filename can
%                   be used with
%                   %fX:Y for filename(X:Y) (relative filename!) and
%                   %pX:Y for folder(X:Y) (parent foldername)
%       action      action to take, either of 'copy', 'pack', {'rename'},
%                   or 'unpack'
%
%       dcmpack     packed DICOM files (from previous call using 'pack')
%
% Output fields:
%
%       p           packed DICOM files (for action := 'pack')

% Version:  v0.9d
% Build:    14061709
% Date:     Jun-17 2014, 9:49 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2013, 2014, Jochen Weber
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

% global UI variable
global ne_ui;

% preset output
p = {[], [], []};

% UI call
nnargin = nargin;
if nnargin < 1 || ...
   (isstruct(dcmfile) && ...
    numel(dcmfile) == 1)

    % load figure
    try
        hFig = xfigure([neuroelf_path('tfg') '/renamedicom.tfg']);
    catch ne_eo;
        error( ...
            'neuroelf:xfigureError', ...
            'Error creating UI for renamedicom: %s.', ...
            ne_eo.message ...
        );
    end

    % get tags
    hTag = hFig.TagStruct;

    % set into global variable
    ne_ui.renamedicom = struct( ...
        'hFig', hFig, ...
        'hTag', hTag, ...
        'actn', 'rename', ...
        'files', {{}}, ...
        'info', [], ...
        'pbar', [], ...
        'patt', 'Anon-%04k00200010-%04k00200011-%05k00200013.dcm');

    % possibly copy pbar
    if nargin > 0 && ...
        isfield(dcmfile, 'pbar') && ...
        any(strcmpi(class(dcmfile.pbar), {'xfigure', 'xprogress'}))
        ne_ui.renamedicom.pbar = dcmfile.pbar;
    end

    % initialize controls
    hTag.LB_renamedicom_files.Value = [];
    hTag.LB_renamedicom_files.String = {};
    hTag.LB_renamedicom_fields.Value = [];
    hTag.LB_renamedicom_fields.String = {};
    hTag.LB_renamedicom_fieldc.Value = [];
    hTag.LB_renamedicom_fieldc.String = {};
    hTag.ED_renamedicom_pattern.String = ne_ui.renamedicom.patt;
    hTag.LB_renamedicom_pattest.String = 'Anon-0001-0001-00001.dcm';

    % set callbacks
    hTag.ED_renamedicom_folder.Callback = @rd_enablesearch;
    hTag.BT_renamedicom_folder.Callback = @rd_browse;
    hTag.BT_renamedicom_search.Callback = @rd_search;
    hTag.BT_renamedicom_fminus.Callback = @rd_delfiles;
    hTag.BT_renamedicom_getinfo.Callback = @rd_getinfo;
    hTag.LB_renamedicom_fields.Callback = @rd_fields2c;
    hTag.LB_renamedicom_fieldc.Callback = @rd_fieldc2s;
    hTag.BT_renamedicom_testpat.Callback = @rd_testpat;
    hTag.BT_renamedicom_copy.Callback = {@rd_rename, 'copy'};
    hTag.BT_renamedicom_pack.Callback = {@rd_rename, 'pack'};
    hTag.BT_renamedicom_rename.Callback = @rd_rename;
    hTag.BT_renamedicom_cancel.Callback = @rd_closeui;

    % set visible and modal
    hFig.HandleVisibility = 'callback';
    hFig.Visible = 'on';
    hFig.WindowStyle = 'modal';

    % wait for dialog to be closed
    uiwait(hFig.MLHandle);

    % no work to do?
    if isempty(ne_ui.renamedicom.files)
        ne_ui.renamedicom = [];
        return;
    end

    % setup arguments for further call
    rd = ne_ui.renamedicom;
    action = rd.actn;
    dcmfile = rd.files;
    pattern = rd.patt;
    opts = struct('pbar', rd.pbar);

    % override nargin
    nnargin = 3;

% regular call
else
    opts = struct('pbar', []);
end

% argument check
if nnargin < 2 || ...
   ((~ischar(dcmfile) || ...
     isempty(dcmfile) || ...
     exist(dcmfile, 'file') ~= 2) && ...
    (~iscell(dcmfile) || ...
     isempty(dcmfile))) || ...
   ~ischar(pattern) || ...
    isempty(pattern) || ...
   ~any(pattern(:)' == '%')
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing dicom filename.' ...
    );
end
if nnargin < 3 || ...
   ~ischar(action) || ...
    isempty(action) || ...
   ~any(strcmpi(action(:)', {'copy', 'pack', 'rename', 'unpack'}))
    action = 'rename';
else
    action = lower(action(:)');
end

% for cell arrays
if iscell(dcmfile)

    % unpacking
    if action(1) == 'u'

        % check packed data
        if size(dcmfile, 2) ~= 3 || ...
           ~all(cellfun(@ischar, dcmfile(:, 1))) || ...
           ~all(cellfun(@ischar, dcmfile(:, 2))) || ...
           ~all(strcmp(cellfun(@class, dcmfile(:, 3), 'UniformOutput', false), 'uint8'))
            error( ...
                'neuroelf:BadArgument', ...
                'Bad or missing renamedicom package.' ...
            );
        end

        % progress bar stuff
        if ~isempty(opts.pbar)
            opts.pbar.Visible = 'on';
            verb = 'Unpacking';
            opts.pbar.Progress(0, sprintf('%s %d DICOM files...', verb, size(dcmfile, 1)));
            lpr = now;
            lpd = 1 / 86400;
        end

        % write
        for fc = 1:size(dcmfile, 1)
            if ~isempty(opts.pbar) && ...
               (now - lpr) >= lpd
                try
                    lpr = now;
                    opts.pbar.Progress((fc - 1) / numel(dcmfile), ...
                        sprintf('%s (%d/%d) DICOM file...', verb, fc, size(dcmfile, 1)));
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    opts.pbar = [];
                end
            end
            try
                binwrite(strrep(pattern, '%', dcmfile{fc, 2}), dcmfile{fc, 3});
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end

    % renameing/packing/copying
    else

        % remove no-char and illegal things
        badidx = false(size(dcmfile));
        for fc = 1:numel(dcmfile)
            badidx = (~ischar(dcmfile{fc}) || isempty(dcmfile{fc}));
        end
        dcmfile(badidx) = [];

        % progress bar stuff
        if ~isempty(opts.pbar)
            opts.pbar.Visible = 'on';
            switch (action(1))
                case {'c'}
                    verb = 'Copying';
                case {'p'}
                    verb = 'Packing';
                case {'r'}
                    verb = 'Renaming';
            end
            opts.pbar.Progress(0, sprintf('%s %d DICOM files...', verb, numel(dcmfile)));
            lpr = now;
            lpd = 1 / 86400;
        end

        % rename single files
        p = cell(numel(dcmfile), 3);
        for fc = 1:numel(dcmfile)
            if ischar(dcmfile{fc})
                if ~isempty(opts.pbar) && ...
                   (now - lpr) >= lpd
                    try
                        lpr = now;
                        opts.pbar.Progress((fc - 1) / numel(dcmfile), ...
                            sprintf('%s (%d/%d) DICOM file...', verb, fc, numel(dcmfile)));
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        opts.pbar = [];
                    end
                end
                try
                    p(fc, :) = renamedicom(dcmfile{fc}, pattern, action);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
        end
        p(cellfun('isempty', p(:, 3)), :) = [];

    end

    % return
    return;
end

% try to read dicom file
try
    dcmobj = xff(dcmfile, 'h');
    if length(dcmobj.DataKeys) < 10
        error('BAD_DICOM_FILE');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:BadArgument', ...
        'Probably no DICOM file or not supported yet.' ...
    );
end

% get new filename
[newfile, warn] = parsepattern(pattern(:)', dcmobj.DataKeyLookup, dcmobj.Data, dcmfile);
if ~isempty(warn)
    warning( ...
        'neuroelf:BadArgument', ...
        'Bad pattern: %s', ...
        warn ...
    );
    return;
end

% copy/rename if not the same
if action(1) == 'p'
    p = {dcmfile, newfile, binread(dcmfile)};
elseif ~strcmpi(dcmfile, newfile)
    if action(1) == 'c'
        cpfile(dcmfile, newfile);
    else
        movefile(dcmfile, newfile);
    end
end



% internal functions
function [fpat, warn] = parsepattern(fpat, keyl, data, fname)
warn = '';
% parse pattern
[rbeg, rend] = regexp(fpat, '\%(\d*[kK]\d{8}|[fFpP]\-?\d+\:\-?\d+)');
if isempty(rend)
    warn = 'Invalid pattern input, no valid variable found.';
    return;
end

% parse filename
[ffpath, fname] = fileparts(fname);
[fnull, fpath] = fileparts(ffpath);

% keep parsing until complete
while ~isempty(rend)
    ppat = fpat(rbeg(1):rend(1));
    if any(ppat == 'k' | ppat == 'K')
        keyv = sprintf('k_%s_%s', ppat(end-7:end-4), ppat(end-3:end));
        if numel(ppat) > 10
            nummod = ['%' ppat(2:end-9) 'g'];
        else
            nummod = '%g';
        end
        if numel(nummod) > 2 && ...
            nummod(2) ~= '0'
            nummod = ['%0' nummod(2:end)];
        end
        if ~isfield(keyl, keyv)
            warn = ['The requested key wasn''t found: ' ppat(end-8:end) '.'];
            fpat = [fpat(1:rbeg(1)-1), 'NF' keyv, fpat(rend(1)+1:end)];
            return;
        end
        value = data(keyl.(keyv)).Value;
        if ischar(value) && ...
            all((value >= '0' & value <= '9') | value == '.') && ...
            sum(value == '.') <= 1
            value = str2double(value);
        end
        if isnumeric(value)
            if numel(value) == 1
                fpat = sprintf(['%s', nummod, '%s'], fpat(1:rbeg(1)-1), ...
                    value, fpat(rend(1)+1:end));
            elseif numel(value) == 2
                fpat = sprintf(['%s', nummod, nummod, '%s'], fpat(1:rbeg(1)-1), ...
                    value(1), value(2), fpat(rend(1)+1:end));
            elseif numel(value) == 3
                fpat = sprintf(['%s', nummod, nummod, nummod, '%s'], fpat(1:rbeg(1)-1), ...
                    value(1), value(2), value(3), fpat(rend(1)+1:end));
            else
                warn = ['Unsupported numeric value size in key: ' ppat(end-8:end) '.'];
                fpat = [fpat(1:rbeg(1)-1), 'IV' keyv, fpat(rend(1)+1:end)];
                return;
            end
        elseif ischar(value)
            fpat = [fpat(1:rbeg(1)-1), value(:)', fpat(rend(1)+1:end)];
        else
            warn = ['Unsupported value class in key: ' ppat(end-8:end) '.'];
            fpat = [fpat(1:rbeg(1)-1), 'IV' keyv, fpat(rend(1)+1:end)];
            return;
        end
    else
        colpos = find(ppat == ':');
        vfrom = str2double(ppat(3:colpos-1));
        vto = str2double(ppat(colpos+1:end));
        if any(ppat == 'f' | ppat == 'F')
            value = fname;
        elseif any(ppat == 'p')
            value = fpath;
        else
            value = ffpath;
        end
        if sign(vfrom) <= 0
            vfrom = numel(value) + vfrom;
        end
        if sign(vto) <= 0
            vto = numel(value) + vto;
        end
        if vfrom > 0 && ...
            vfrom <= numel(value) && ...
            vto >= vfrom && ...
            vto <= numel(value)
            value = value(vfrom:vto);
        end
        fpat = [fpat(1:rbeg(1)-1), value(:)', fpat(rend(1)+1:end)];
    end
    [rbeg, rend] = regexp(fpat, '\%(\d*[kK]\d{8}|[fFpP]\-?\d+\:\-?\d+)');
end
if any(fpat == '%')
    warn = 'Other pattern specifiers remain (unsupported).';
end



% UI functions
function rd_closeui(varargin)
global ne_ui;
ne_ui.renamedicom.hFig.Delete;

function rd_enablesearch(varargin)
global ne_ui;
hTag = ne_ui.renamedicom.hTag;
fstring = hTag.ED_renamedicom_folder.String;
if any(fstring == '*')
    [fpath, ffind] = fileparts(fstring);
    fexist = ~isempty(findfiles(fpath, ffind, 'dirs=1', 'depth=1'));
else
    fexist = (exist(fstring, 'dir') > 0);
end
if fexist
    hTag.BT_renamedicom_search.Enable = 'on';
else
    hTag.BT_renamedicom_search.Enable = 'off';
end

function rd_browse(varargin)
global ne_ui;

% browse for directory
dd = uigetdir(pwd, 'Please select the folder containing the DICOM files...');
if isequal(dd, 0)
    return;
end

% put into edit field
hTag = ne_ui.renamedicom.hTag;
hTag.ED_renamedicom_folder.String = dd;
hTag.BT_renamedicom_search.Enable = 'on';

function rd_search(varargin)
global ne_ui;
hFig = ne_ui.renamedicom.hFig;
hTag = ne_ui.renamedicom.hTag;

% extension/depth limitation
if hTag.CB_renamedicom_dcmonly.Value > 0
    ext = {'*.dcm', '*.DCM'};
else
    ext = '*';
end
if hTag.CB_renamedicom_subdirs.Value > 0
    depth = {};
else
    depth = {'depth=1'};
end

% try search
try
    hTag.TX_renamedicom_files.String = 'Files found:';
    hFig.Pointer = 'watch';
    drawnow;
    files = findfiles(hTag.ED_renamedicom_folder.String, ext, depth{:});

    % temporary bug fix!!
    rfiles = findfiles(hTag.ED_renamedicom_folder.String, ext, depth{:});
    if ~isequal(files, rfiles)
        files = union(files(:), rfiles(:));
    end
catch ne_eo;
    warning( ...
        'neuroelf:CallbackError', ...
        'Error finding files in ''%s'': %s.', ...
        hTag.ED_renamedicom_folder.String, ne_eo.message ...
    );
    hTag.LB_renamedicom_files.String = {};
    hTag.LB_renamedicom_files.Value = [];
    hTag.LB_renamedicom_files.ListBoxTop = 1;
    hTag.BT_renamedicom_getinfo.Enable = 'Off';
    hFig.Pointer = 'arrow';
    return;
end

% empty
if isempty(files)
    warning( ...
        'neuroelf:FileNotFound', ...
        'No DICOM files found in ''%s''.', ...
        hTag.ED_renamedicom_folder.String ...
    );
    newidx = [];
    hTag.BT_renamedicom_getinfo.Enable = 'Off';
else
    newidx = 1;
    hTag.BT_renamedicom_getinfo.Enable = 'On';
end

% set to listbox (unique and sorted)
files = unique(files(:));
hTag.TX_renamedicom_files.String = sprintf('%d files found:', numel(files));
hTag.LB_renamedicom_files.String = files;
hTag.LB_renamedicom_files.Value = newidx;
hTag.LB_renamedicom_files.ListBoxTop = 1;
hFig.Pointer = 'arrow';

function rd_getinfo(varargin)
global ne_ui;
hFig = ne_ui.renamedicom.hFig;
hTag = ne_ui.renamedicom.hTag;

% single entry ?
sidx = hTag.LB_renamedicom_files.Value;
if numel(sidx) ~= 1
    return;
end
files = hTag.LB_renamedicom_files.String;
if ~iscell(files)
    files = cellstr(files);
end
file = files{sidx};

% try to load file
try
    hFig.Pointer = 'watch';
    drawnow;
    dcmfile = cell(1, 1);
    dcmfile{1} = xff(file);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects(dcmfile);
    warning( ...
        'neuroelf:CallbackError', ...
        'Not a DICOM file: ''%s''; removing from list...', ...
        file ...
    );
    files(sidx) = [];
    hTag.LB_renamedicom_files.String = files;
    hTag.LB_renamedicom_files.Value = [];
    hFig.Pointer = 'arrow';
    return;
end

% put information into further listboxes
dcmfile = dcmfile{1};

% get dictionary and data
dict = dicom_dic(4, dcmfile.DataDictionary);
k2tg = dict.KeyToTag;
keys = sort(dcmfile.DataKeys);
data = dcmfile.Data;
keyl = dcmfile.DataKeyLookup;
ne_ui.renamedicom.info = struct( ...
    'dict', dict, ...
    'file', file, ...
    'k2tg', k2tg, ...
    'keys', {keys}, ...
    'data', data, ...
    'keyl', keyl);

% clear object
dcmfile.ClearObject;

% put keys and values into list
fields = cell(size(keys));
for fc = 1:numel(keys)
    value = data(keyl.(keys{fc})).Value;
    if isnumeric(value)
        if numel(value) == 1
            fields{fc} = sprintf('%g', value);
        elseif isempty(value)
            fields{fc} = '[]';
        elseif numel(value) < 7
            fields{fc} = ...
                sprintf('[%g%s]', value(1), sprintf(', %g', value(2:end)));
        else
            sz = size(value);
            fields{fc} = ...
                sprintf('%d%s value', sz(1), sprintf('-by-%d', sz(2:end)));
        end
    elseif ischar(value)
        if numel(value) <= 64
            fields{fc} = value;
        else
            fields{fc} = [value(1:28) ' ... ' value(end-27:end)];
        end
    else
        sz = size(value);
        fields{fc} = sprintf( ...
            '%d%s %s', sz(1), sprintf('-by-%d', sz(2:end)), class(value));
    end
    if isfield(k2tg, keys{fc})
        keys{fc} = sprintf('%s (k%s)', k2tg.(keys{fc}), keys{fc}([3:6, 8:11]));
    else
        keys{fc} = keys{fc}([1, 3:6, 8:11]);
    end
end
hTag.LB_renamedicom_fields.String = keys;
hTag.LB_renamedicom_fields.Value = 1;
hTag.LB_renamedicom_fields.ListBoxTop = 1;
hTag.LB_renamedicom_fieldc.String = fields;
hTag.LB_renamedicom_fieldc.Value = 1;
hTag.LB_renamedicom_fieldc.ListBoxTop = 1;

% enable test button
hTag.BT_renamedicom_testpat.Enable = 'On';

% test pattern
rd_testpat

% reset pointer
hFig.Pointer = 'arrow';

function rd_fields2c(varargin)
global ne_ui;
hTag = ne_ui.renamedicom.hTag;
hTag.LB_renamedicom_fieldc.Value = hTag.LB_renamedicom_fields.Value;
hTag.LB_renamedicom_fieldc.ListBoxTop = hTag.LB_renamedicom_fields.ListBoxTop;

function rd_fieldc2s(varargin)
global ne_ui;
hTag = ne_ui.renamedicom.hTag;
hTag.LB_renamedicom_fields.Value = hTag.LB_renamedicom_fieldc.Value;
hTag.LB_renamedicom_fields.ListBoxTop = hTag.LB_renamedicom_fieldc.ListBoxTop;

function rd_delfiles(varargin)
global ne_ui;
hTag = ne_ui.renamedicom.hTag;

% remove any selected values
sidx = hTag.LB_renamedicom_files.Value;
if ~isempty(sidx)
    files = hTag.LB_renamedicom_files.String;
    if ~iscell(files)
        files = cellstr(files);
    end
    files(sidx) = [];
    hTag.LB_renamedicom_files.String = files;
    hTag.LB_renamedicom_files.Value = [];
end

function rd_testpat(varargin)
global ne_ui;
hTag = ne_ui.renamedicom.hTag;

% only valid if info is valid
info = ne_ui.renamedicom.info;
if ~isstruct(info)
    return;
end

% get required information
fpat = hTag.ED_renamedicom_pattern.String;
keyl = info.keyl;
data = info.data;
file = info.file;

% set output to parsed string
[fpat, warn] = parsepattern(fpat, keyl, data, file);
hTag.LB_renamedicom_pattest.String = fpat;
if ~isempty(warn)
    hTag.BT_renamedicom_copy.Enable = 'Off';
    hTag.BT_renamedicom_pack.Enable = 'Off';
    hTag.BT_renamedicom_rename.Enable = 'Off';
    uiwait(warndlg(warn, 'Rename DICOM files - warning', 'modal'));
else
    hTag.BT_renamedicom_copy.Enable = 'On';
    hTag.BT_renamedicom_pack.Enable = 'On';
    hTag.BT_renamedicom_rename.Enable = 'On';
end

function rd_rename(varargin)
global ne_ui;
hTag = ne_ui.renamedicom.hTag;

% copy flag
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmpi(varargin{3}(:)', {'copy', 'pack', 'rename'}))
    actn = 'rename';
else
    actn = lower(varargin{3}(:)');
end

% call testpat and check anyway!
rd_testpat;

% if button rename is still enabled, it worked
if ~strcmpi(hTag.BT_renamedicom_rename.Enable, 'on')
    return;
end

% set files and pattern
files = hTag.LB_renamedicom_files.String;
if ~iscell(files)
    files = cellstr(files);
end
ne_ui.renamedicom.actn = actn;
ne_ui.renamedicom.files = files;
ne_ui.renamedicom.patt = hTag.ED_renamedicom_pattern.String;

% close dialog
ne_ui.renamedicom.hFig.Delete;

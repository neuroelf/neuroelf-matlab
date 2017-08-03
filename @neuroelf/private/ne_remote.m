% FUNCTION ne_remote: handle remote configuration
function varargout = ne_remote(varargin)

% Version:  v0.9d
% Build:    14072516
% Date:     Jul-25 2014, 4:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, 2014, Jochen Weber
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

% global variable
global ne_gcfg;

% initialize output
varargout = cell(1, nargout);

% requires a valid command
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmpi(varargin{3}(:)', {'listen', 'unlisten'}))
    return;
end

% unlisten
if strcmpi(varargin{3}(:)', 'unlisten')

    % simply disable listener (most of the work done in timer callback!)
    ne_gcfg.c.remote = false;

    % then delete files that are no longer needed
    ifmt = lower(ne_gcfg.c.ini.Remote.ImageFormat);
    ipath = [neuroelf_path('remote') '/images'];
    for dimgs = {'corslice', 'render', 'sagslice', 'surface', 'traslice', 'zoomslice'}
        dimg = sprintf('%s/%s.%s', ipath, dimgs{1}, ifmt);
        if exist(dimg, 'file') > 0
            try
                delete(dimg);
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end
    end
    for dimgs = {'tcplot.csv'}
        dimg = sprintf('%s/%s', ipath, dimgs{1});
        if exist(dimg, 'file') > 0
            try
                delete(dimg);
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end
    end

    % then return
    return;
end

% CODE FROM HERE: initialize listener (and setup everything!)

% no folder set up for listening
if nargin > 3 && ...
    ischar(varargin{4}) && ...
   ~isempty(varargin{4}) && ...
    exist(varargin{4}(:)', 'dir') > 0
    [ia, ne_gcfg.c.remotecfg.scanpath] = isabsolute(varargin{4}(:)');
end
if isempty(ne_gcfg.c.remotecfg.scanpath) || ...
    exist(ne_gcfg.c.remotecfg.scanpath, 'dir') ~= 7
    rp = uigetdir(pwd, 'Please select the folder containing the remote commands...');
    if ~ischar(rp) || ...
        isempty(rp) || ...
        exist(rp, 'dir') ~= 7
        ne_gcfg.c.remote = false;
        return;
    end
    [ia, ne_gcfg.c.remotecfg.scanpath] = isabsolute(rp);
end
scanpath = ne_gcfg.c.remotecfg.scanpath;
ne_gcfg.c.ini.Remote.ScanFolder = scanpath;

% remove old content
fs = filesep;
rpc = dir([scanpath fs '*.nro']);
try
    for fc = 1:numel(rpc)
        delete([scanpath fs rpc(fc).name]);
    end
    rpc = dir([scanpath fs '*.nro']);
    if ~isempty(rpc)
        error( ...
            'neuroelf:FolderNotWritable', ...
            'NeuroElf listener folder not writable.' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end
rpc = dir([scanpath fs '*.nrt']);
try
    for fc = 1:numel(rpc)
        delete([scanpath fs rpc(fc).name]);
    end
    rpc = dir([scanpath fs '*.nrt']);
    if ~isempty(rpc)
        error( ...
            'neuroelf:FolderNotWritable', ...
            'NeuroElf listener folder not writable.' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% get configuration
rini = ne_gcfg.c.ini.Remote;

% initialize basic content
ne_gcfg.c.remotecfg.cmdcount = 0;
ne_gcfg.c.remotecfg.cmdfiles = cell(0, 6);
ne_gcfg.c.remotecfg.commandid = zeros(0, 2);
ne_gcfg.c.remotecfg.commands = cell(0, 6);
ne_gcfg.c.remotecfg.gcwcount = rini.GCWaitCount;
ne_gcfg.c.remotecfg.lastcmd = -1;
ne_gcfg.c.remotecfg.lastscan = -1;
ne_gcfg.c.remotecfg.logfile = -1;
ne_gcfg.c.remotecfg.scanning = false;
ne_gcfg.c.remotecfg.scantimer = [];
ne_gcfg.c.remotecfg.stopping = false;
ne_gcfg.c.remotecfg.stoptimer = false;

% try catch around rest
try

    % open log file
    ne_gcfg.c.remotecfg.logfile = fopen([scanpath fs rini.LogFile], 'a');

    % log starting the remote
    ne_remote_log('NE_REMOTE%internal%Starting remote...');
    ne_gcfg.c.remotecfg.scantimer = timer;

    % set up scanner
    set(ne_gcfg.c.remotecfg.scantimer,  ...
        'ExecutionMode', 'fixedSpacing', ...
        'Period',        0.001 * max(10, ceil(1000 / rini.ScanFrequency)), ...
        'StartDelay',    0.2, ...
        'StopFcn',       @ne_remote_cleanup, ...
        'TimerFcn',      @ne_remote_listen);

    % update controls
    ts = ne_gcfg.h.MainFig.TagStruct;
    ts.UIM_NeuroElf_ListenerStart.Enable = 'off';
    ts.UIM_NeuroElf_ListenerStop.Enable = 'on';
    ts.BT_NeuroElf_RListener.Value = 1;
    ts.BT_NeuroElf_RListener.Callback = {@ne_remote, 'unlisten'};

    % start timer
    ne_gcfg.c.remote = true;
    start(ne_gcfg.c.remotecfg.scantimer);

% handle errors
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;

    % clean up
    try
        ne_remote_cleanup;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    rethrow(ne_eo);
end

% update once
drawnow;



% SUB-functions



% clean up function
function ne_remote_cleanup(varargin)
global ne_gcfg;

% don't run twice
if ne_gcfg.c.remotecfg.stopping
    return;
end
ne_gcfg.c.remotecfg.stopping = true;

% while processing a command or scanning
if ~isempty(ne_gcfg.c.remotecfg.commands) || ...
    ne_gcfg.c.remotecfg.scanning
    warning( ...
        'neuroelf:InternalTimeout', ...
        'Cleaning up, but jobs still running(?).' ...
    );
    for cc = size(ne_gcfg.c.remotecfg.commands, 1):-1:1
        tdata = ne_gcfg.c.remotecfg.commands(cc, :);
        try
            stop(tdata{2});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        try
            delete(tdata{2});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ne_gcfg.c.remotecfg.commands(cc, :) = [];
    end
end

% clean up timers
if numel(ne_gcfg.c.remotecfg.scantimer) == 1 && ...
    isa(ne_gcfg.c.remotecfg.scantimer, 'timer')
    try
        stop(ne_gcfg.c.remotecfg.scantimer);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    try
        delete(ne_gcfg.c.remotecfg.scantimer);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
ne_gcfg.c.remotecfg.scantimer = [];

% close and delete all files
for fc = 1:size(ne_gcfg.c.remotecfg.cmdfiles, 1)
    if ne_gcfg.c.remotecfg.cmdfiles{fc, 2} > 1
         try
             fclose(ne_gcfg.c.remotecfg.cmdfiles{fc, 2});
         catch ne_eo;
             ne_gcfg.c.lasterr = ne_eo;
         end
    end
    if exist(ne_gcfg.c.remotecfg.cmdfiles{fc, 1}, 'file') > 0
        try
            delete(ne_gcfg.c.remotecfg.cmdfiles{fc, 1});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
end
ne_gcfg.c.remotecfg.cmdfiles = cell(0, 6);
pause(0.1);
fs = filesep;
scanpath = ne_gcfg.c.remotecfg.scanpath;
rpc = dir([scanpath fs '*.nro']);
try
    for fc = 1:numel(rpc)
        delete([scanpath fs rpc(fc).name]);
    end
    rpc = dir([scanpath fs '*.nro']);
    if ~isempty(rpc)
        error( ...
            'neuroelf:FolderNotWritable', ...
            'NeuroElf listener folder not writable (couldn''t delete outputs).' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
if ne_gcfg.c.remotecfg.logfile > 1
    try
        fclose(ne_gcfg.c.remotecfg.logfile);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
ne_gcfg.c.remotecfg.logfile = -1;

% set controls
ts = ne_gcfg.h.MainFig.TagStruct;
ts.UIM_NeuroElf_ListenerStart.Enable = 'on';
ts.UIM_NeuroElf_ListenerStop.Enable = 'off';
ts.BT_NeuroElf_RListener.Value = 0;
ts.BT_NeuroElf_RListener.Callback = {@ne_remote, 'listen'};

% store config
ne_gcfg.c.remotecfg.stopping = false;

% disable listener
ne_gcfg.c.remote = false;


% find files (with restrictions)
function [f, folder] = ne_remote_findfiles(uid, folder)
global ne_gcfg;
rini = ne_gcfg.c.ini.Remote;
ubf = rini.UserPrivBase{uid};
upf = rini.UserPrivFolders{uid};

% folder must be absolute doesn't exist
folder = strrep(folder, '\', '/');
while ~isempty(strfind(folder, '//'))
    folder = strrep(folder, '//', '/');
end
if isempty(folder) || ...
    folder(1) ~= '/' || ...
    exist([ubf folder(2:end)], 'dir') ~= 7 || ...
    isempty(upf)
    f = {};
    folder = '/';
    return;
end

% see if folder is in list of upf
folderparts = splittocellc(folder, '/');
if numel(folderparts) > 1
    fc = 3;
    while fc <= numel(folderparts)
        if strcmp(folderparts{fc}, '..')
            if fc > 2
                folderparts(fc-1:fc) = [];
                fc = fc - 2;
            else
                folderparts(fc:fc) = [];
                fc = fc - 1;
            end
        end
        fc = fc + 1;
    end
    if numel(folderparts) > 1
        folder = gluetostringc(folderparts, '/');
    else
        folder = '/';
    end
end
if isempty(folderparts) || ...
    strcmp(folderparts{1}, '..')
    [f, folder] = re_remote_findfiles(uid, '/');
    return;
end
nfu = zeros(numel(upf), 1);
for fc = 1:numel(upf)
    upf{fc} = splittocellc(upf{fc}, '/');
    nfu(fc) = numel(upf{fc});
end
nfp = numel(folderparts);
fpfound = any(strcmp(upf, '/'));
fpmatch = fpfound;
if ~fpfound
    for fc = 1:numel(upf)
        nft = min(nfu(fc), nfp);
        if all(strcmpi(folderparts(1:nft), upf{fc}(1:nft)))
            fpfound = true;
            if nfu(fc) <= nfp
                fpmatch = true;
            end
        end
    end
end
if ~fpfound
    f = {};
    return;
end

% find entries in folder
d = findfiles([ubf folder(2:end)], '*', 'depth=1', 'dirs=1');
if fpmatch
    f = findfiles([ubf folder(2:end)], '*.*', 'depth=1', 'dirs=0');
else
    f = {};
end

% remove bad entries
if ~isempty(d)
    d(~cellfun('isempty', regexpi(d, '.*\/\.[^\/]*$'))) = [];
    d = regexprep(d, ['^' strrep(ubf, '/', '\/')], '/');
end
if ~isempty(f)
    f(~cellfun('isempty', regexpi(f, '.*\/\.[^\/]*$'))) = [];
    f = regexprep(f, ['^' strrep(ubf, '/', '\/')], '/');
    ff = false(size(f));
    for fc = 1:numel(rini.FindFiles)
        ff = ff | (~cellfun('isempty', regexpi(f, rini.FindFiles{fc})));
    end
    f = f(ff);
end

% only present sub-dirs if matching
nfp = nfp + 1;
nfu = min(nfu, nfp);
for dc = numel(d):-1:1
    [fpp, fpn, fpe] = fileparts(d{dc});
    folderparts{nfp} = [fpn fpe];
    fpfound = false;
    for fc = 1:numel(upf)
        if all(strcmpi(folderparts(1:nfu(fc)), upf{fc}(1:nfu(fc))))
            fpfound = true;
            break;
        end
    end
    if ~fpfound
        d(dc) = [];
    else
        d{dc} = [d{dc} '/'];
    end
end

% return value
if strcmp(folder, '/')
    f = [d(:); f(:)];
else
    f = [{'../'}; d(:); f(:)];
end


% log entry
function ne_remote_log(entry)
global ne_gcfg;

% file valid
if ne_gcfg.c.remotecfg.logfile > 1

    % get info
    logtime = datestr(now, 'dd-mmm-yyyy HH:MM:SS.FFF');
    logfmem = java.lang.Runtime.getRuntime.freeMemory;
    logtmem = java.lang.Runtime.getRuntime.totalMemory;
    fprintf(ne_gcfg.c.remotecfg.logfile, '%s %09d/%09d %s\n', ...
        logtime, logfmem, logtmem, entry);
end


% authentication
function ne_remote_auth(fid, auth)
global ne_gcfg;

% split auth request
auths = splittocellc(auth, '%');
if numel(auths) < 4 || ...
   ~strcmpi(auths{2}, 'auth')
    return;
end

% check against database
outfile = sprintf('%s/%s-%s.nrot', ne_gcfg.c.remotecfg.scanpath, fid, auths{1});
rini = ne_gcfg.c.ini.Remote;
uid = find(strcmp(rini.User, auths{3}));
if numel(uid) ~= 1 || ...
   ~strcmp(rini.UserPasswd{uid}, auths{4})
    pass = 'failed';
else
    % session
    if numel(auths) > 4 && ...
        strcmp(['S' auths{5}], makelabel(['S' auths{5}])) && ...
        isfield(ne_gcfg.c.remotecfg.session, ['S' auths{5}])
        sess = auths{5};
    else
        sess = sprintf('%06x', floor(1 + 16777215 * rand(1, 1)));
        while isfield(ne_gcfg.c.remotecfg.session, ['S' sess])
            sess = sprintf('%06x', floor(1 + 16777215 * rand(1, 1)));
        end
        try
            ne_remote_initsess(sess);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            binwrite(outfile, uint8(ne_eo.message));
            movefile(outfile, outfile(1:end-1));
            return;
        end
    end
    pass = sprintf('ok%s', sess);
end

% log entry
ne_remote_log([auth '%' pass]);

% write out file
binwrite(outfile, uint8(pass));
movefile(outfile, outfile(1:end-1));


% session cleanup
function ne_remote_cleanupsess(sess)
global ne_gcfg;

% get session
sess = ne_gcfg.c.remotecfg.session.(sess);


% session initialization
function ne_remote_initsess(sess)
global ne_gcfg;
ch = ne_gcfg.h;
ts = ch.MainFig.TagStruct;

% populate tag struct
tf = fieldnames(ts);
for fc = 1:numel(tf)
    ts.(tf{fc}) = get(ts.(tf{fc}).MLHandle);
end

% initialize some special controls' fields
ts.DD_NeuroElf_varlist.String = ch.SliceVar.String(1);
if ~isempty(ch.SliceVar.UserData)
    ts.DD_NeuroElf_varlist.UserData = ch.SliceVar.UserData(1, :);
else
    ts.DD_NeuroElf_varlist.UserData = cell(0, 4);
end
ts.DD_NeuroElf_varlist.Value = 1;

% copy over to ch
tf = fieldnames(ch);
for fc = 1:numel(tf)
    if isxfigure(ch.(tf{fc}))
        try
            ch.(tf{fc}) = ts.(ch.(tf{fc}).Tag);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            ch.(tf{fc}) = struct;
        end
    elseif isstruct(ch.(tf{fc})) && ...
        numel(ch.(tf{fc})) == 1
        stf = fieldnames(ch.(tf{fc}));
        for sfc = 1:numel(stf)
            if isxfigure(ch.(tf{fc}).(stf{sfc}))
                try
                    ch.(tf{fc}).(stf{fc}) = ts.(ch.(tf{fc}).(stf{fc}).Tag);
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end
            end
        end
    end
end
disp(ch);

% create session struct
ne_gcfg.c.remotecfg.session.(['S' sess]) = struct( ...
    'fcfg', ne_gcfg.fcfg, ...
    'h',    ch, ...
    'lacc', now, ...
    'uics', ts, ...
    'wspc', struct);


% TIMER functions



% folder listener
function varargout = ne_remote_listen(varargin)
global ne_gcfg;

% initialize output
varargout = cell(1, nargout);

% no longer listening
if ~ne_gcfg.c.remote

    % clean up and exit
    try
        stop(varargin{1});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    return;
end

% we ARE listening NOW
if ne_gcfg.c.remotecfg.scanning
    return;
end
ne_gcfg.c.remotecfg.scanning = true;

% get config
ne_gcfg.c.remotecfg.lastscan = now;
stfi = ne_gcfg.c.ini.Remote.StaleFileTime / 86400;

% read folder contents
fs = filesep;
rpc = dir([ne_gcfg.c.remotecfg.scanpath fs '*.nrt']);
rpc = {rpc.name};

% add files to listener
for fc = 1:numel(rpc)

    % open file
    try
        fnid = rpc{fc}(1:end-4);
        rpc{fc} = [ne_gcfg.c.remotecfg.scanpath fs rpc{fc}];
        fid = fopen(rpc{fc}, 'r');
        if fid < 1
            continue;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        continue;
    end

    % rename file
    try
        movefile(rpc{fc}, [rpc{fc}(1:end-1) 'l']);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        fclose(fid);
        continue;
    end

    % store information in config
    ne_gcfg.c.remotecfg.cmdfiles(end+1, :) = ...
        {[rpc{fc}(1:end-1) 'l'], fid, 0, now + stfi, fnid, 0};
end

% handle listened-to files
for fc = size(ne_gcfg.c.remotecfg.cmdfiles, 1):-1:1

    % check whether file is larger than before
    fid = ne_gcfg.c.remotecfg.cmdfiles{fc, 2};
    fnid = ne_gcfg.c.remotecfg.cmdfiles{fc, 5};
    try
        fseek(fid, 0, 1);
        fsize = ftell(fid);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        ne_gcfg.c.remotecfg.cmdfiles(fc, :) = [];
        continue;
    end
    if fsize > ne_gcfg.c.remotecfg.cmdfiles{fc, 3}

        % read additional commands
        try
            fseek(fid, ne_gcfg.c.remotecfg.cmdfiles{fc, 3}, -1);
            cmdadd = fread(fid, [1, Inf], 'uint8=>char');
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            continue;
        end
        ne_gcfg.c.remotecfg.cmdfiles{fc, 3} = ftell(fid);

        % split commands
        cmds = splittocellc(cmdadd, char([10, 13]), true, true);
        cmds(cellfun('isempty', cmds)) = [];

        % create new commands
        for cc = 1:numel(cmds)

            % invalid command
            if ~any(cmds{cc} == '%') || ...
                isempty(regexpi(cmds{cc}, '%[a-z]+'))
                continue;
            end
            ppos = findfirst(cmds{cc} == '%');
            cmdid = cmds{cc}(1:ppos-1);

            % authentication request
            if ~isempty(regexpi(cmds{cc}(ppos:end), '^%auth%[^%]+%[^%]+$'))
                try
                    ne_remote_auth(fnid, cmds{cc});
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end
                continue;
            end

            % create new command ID
            cid = floor(16777216 * rand(1, 1));
            while any(ne_gcfg.c.remotecfg.commandid(:, 1)) == cid
                cid = floor(16777216 * rand(1, 1));
            end

            % create new timer object
            try
                ct = timer;
                ne_gcfg.c.remotecfg.commandid(end+1, :) = [cid, 0];
                ne_gcfg.c.remotecfg.commands(end+1, :) = ...
                    {cid, ct, fnid, cc, cmdid, cmds{cc}};
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end

            % set timer properties
            set(ct, ...
                'ExecutionMode', 'fixedSpacing', ...
                'Period',        0.5, ...
                'StartDelay',    0.005 + 0.001 * floor(6 * rand(1, 1)), ...
                'TimerFcn',      {@ne_remote_process, cid});
            start(ct);
        end

    % file unchanged, stale?
    elseif now > ne_gcfg.c.remotecfg.cmdfiles{fc, 4}

        % close file
        try
            fclose(fid);
            delete(ne_gcfg.c.remotecfg.cmdfiles{fc, 1});
            ne_gcfg.c.remotecfg.cmdfiles(fc, :) = [];
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        continue;
    end
end

% reduce garbage-collection wait counter
ne_gcfg.c.remotecfg.gcwcount = ne_gcfg.c.remotecfg.gcwcount - 1;

% perform garbage collection?
if ne_gcfg.c.remotecfg.gcwcount < 1

    % re-set counter
    ne_gcfg.c.remotecfg.gcwcount = ne_gcfg.c.ini.Remote.GCWaitCount;

    % drop stale sessions
    sf = fieldnames(ne_gcfg.c.remotecfg.session);
    sf(strcmp(sf, 'S000000')) = [];
    coff = now - 1 / 24;
    for fc = 1:numel(sf)
        if ne_gcfg.c.remotecfg.session.(sf{fc}).lacc < coff
            try
                ne_remote_cleanupsess(sf{fc});
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
            ne_gcfg.c.remotecfg.session = ...
                rmfield(ne_gcfg.c.remotecfg.session, sf{fc});
        end
    end

    % run Java GC
    java.lang.Runtime.getRuntime.gc;
    ne_remote_log('NE_REMOTE%internal%Java heap memory garbage-collector run...');
end

% done scanning
ne_gcfg.c.remotecfg.scanning = false;


% process command
function varargout = ne_remote_process(varargin)
global ne_gcfg;
ch = ne_gcfg.h;
rini = ne_gcfg.c.ini.Remote;

% setup output
varargout = cell(1, nargout);

% only valid if timer found
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1
    return;
end

% get timer data
cid = varargin{3};
cpos = find(ne_gcfg.c.remotecfg.commandid(:, 1) == cid);

% already being processed
if ne_gcfg.c.remotecfg.commandid(cpos, 2) > 0
    return;
end

% find other commands from same file
tdata = ne_gcfg.c.remotecfg.commands(cpos, :);
ocmds = (strcmp(ne_gcfg.c.remotecfg.commands(:, 3), tdata{3}) & ...
    strcmp(ne_gcfg.c.remotecfg.commands(:, 5), tdata{5}));

% if another command
if sum(ocmds) > 1

    % with lower command number?
    ocmdi = cat(1, ne_gcfg.c.remotecfg.commands{ocmds, 4});
    if any(ocmdi < tdata{4})

        % then don't do anything (and try again on next iteration)
        return;
    end
end

% set command to being processed
stop(tdata{2});
ne_gcfg.c.remotecfg.commandid(cpos, 2) = 1;
ne_gcfg.c.remotecfg.lastcmd = now;

% split command into particles
cmd = splittocellc(deblank(tdata{6}), '%');

% reject invalid commands outright
if numel(cmd) < 5
    try
        delete(tdata{2});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    cpos = find(ne_gcfg.c.remotecfg.commandid(:, 1) == cid);
    ne_gcfg.c.remotecfg.commandid(cpos, :) = [];
    ne_gcfg.c.remotecfg.commands(cpos, :) = [];
    return;
end

% create output filename
outfile = sprintf('%s/%s-%s.nrot', ne_gcfg.c.remotecfg.scanpath, tdata{3}, cmd{1});

% invalid authentication
uname = cmd{2};
uid = find(strcmp(rini.User, uname));
if numel(uid) ~= 1 || ...
   ~strcmp(rini.UserPasswd{uid}, cmd{3})

    % log bad attempt
    ne_remote_log(['UNAUTHORIZED:' tdata{6}]);
    binwrite(outfile, 'bad_auth');
    movefile(outfile, outfile(1:end-1));

    % remove timer and command
    try
        delete(tdata{2});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    cpos = find(ne_gcfg.c.remotecfg.commandid(:, 1) == cid);
    ne_gcfg.c.remotecfg.commandid(cpos, :) = [];
    ne_gcfg.c.remotecfg.commands(cpos, :) = [];

    % return early
    return;
end
if ~strcmp(['S' cmd{4}], makelabel(['S' cmd{4}])) || ...
   ~isfield(ne_gcfg.c.remotecfg.session, ['S' cmd{4}])
    ne_remote_log(['BADSESSION:' tdata{6}]);
    binwrite(outfile, 'session_timeout');
    movefile(outfile, outfile(1:end-1));
    try
        delete(tdata{2});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    cpos = find(ne_gcfg.c.remotecfg.commandid(:, 1) == cid);
    ne_gcfg.c.remotecfg.commandid(cpos, :) = [];
    ne_gcfg.c.remotecfg.commands(cpos, :) = [];
    return;
else
    sess = ['S' cmd{4}];
    ne_gcfg.c.remotecfg.session.(sess).lacc = now;
end
cmd(2:4) = [];

% log command
ne_gcfg.c.remotecfg.cmdcount = ne_gcfg.c.remotecfg.cmdcount + 1;
ne_remote_log(tdata{6});

% prepare output settings
rini = ne_gcfg.c.ini.Remote;
imf = rini.ImageFormat;
if any(strcmpi(imf, {'jpeg', 'jpg'}))
    imq = {'Quality', rini.ImageJPGQuality};
else
    imq = {};
end
tio = ne_gcfg.tio;
uics = ne_gcfg.h.MainFig.TagStruct;
uistate = false;

% depending on command
switch lower(cmd{2})

    % update the position
    case {'cpos'}

        % requires 3 more arguments
        if numel(cmd) > 4
            try
                xp = str2double(cmd{3});
                yp = str2double(cmd{4});
                zp = str2double(cmd{5});
                if numel(xp) == 1 && ...
                    numel(yp) == 1 && ...
                    numel(zp) == 1 && ...
                   ~isinf(xp) && ...
                   ~isnan(xp) && ...
                    xp >= -128 && ...
                    xp <= 128 && ...
                   ~isinf(yp) && ...
                   ~isnan(yp) && ...
                    yp >= -128 && ...
                    yp <= 128 && ...
                   ~isinf(zp) && ...
                   ~isnan(zp) && ...
                    zp >= -128 && ...
                    zp <= 128
                    ne_gcfg.fcfg.cpos = [xp, yp, zp];
                    ne_setslicepos;
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end
        svals = sprintf('%g%%', ne_gcfg.c.svals(:)');
        if ~isempty(svals)
            svals(end) = [];
        end
        out = sprintf('%d%%%d%%%d%%%s\n', ne_gcfg.fcfg.cpos, svals);
        binwrite(outfile, uint8(out));

    % coronal slice
    case {'corslice'}
        im = tio.imCor.Rendered;
        cp = round(min(256, max(1, 128 - ne_gcfg.fcfg.cpos)));
        im(cp(3), :, :) = min(255, double(im(cp(3), :, :)) + 64);
        im(:, 257 - cp(1), :) = min(255, double(im(:, 257 - cp(1), :)) + 64);
        imwrite(im, outfile, imf, imq{:});

    % saggital slice
    case {'sagslice'}
        im = tio.imSag.Rendered;
        cp = round(min(256, max(1, 128 - ne_gcfg.fcfg.cpos)));
        im(cp(3), :, :) = min(255, double(im(cp(3), :, :)) + 64);
        im(:, cp(2), :) = min(255, double(im(:, cp(2), :)) + 64);
        imwrite(im, outfile, imf, imq{:});

    % transversal slice
    case {'traslice'}
        im = tio.imTra.Rendered;
        cp = round(min(256, max(1, 128 - ne_gcfg.fcfg.cpos)));
        im(cp(2), :, :) = min(255, double(im(cp(2), :, :)) + 64);
        im(:, 257 - cp(1), :) = min(255, double(im(:, 257 - cp(1), :)) + 64);
        imwrite(im, outfile, imf, imq{:});

    % zoomed slice
    case {'zoomslice'}
        imwrite(tio.imSlZ.Rendered, outfile, imf, imq{:});

    % time course plot
    case {'tcplot'}

        % get data
        if isempty(ch.TCPlotChildren)
            plotdata = get(ch.TCPlotChild);
            xdata = plotdata.XData(:);
            ydata = plotdata.YData(:);
            cdata = round(255.001 .* plotdata.Color);
            ldata = plotdata.LineWidth;
        else
            plotdata = cell(1, numel(ch.TCPlotChildren));
            plottype = plotdata;
            for pc = 1:numel(plotdata)
                plotdata{pc} = get(ch.TCPlotChildren(pc));
                plottype{pc} = plotdata{pc}.Type;
            end
            plotdata(~strcmp(plottype, 'line')) = [];
            xdata = plotdata{1}.XData;
            ydata = NaN .* zeros(numel(xdata), numel(plotdata));
            cdata = zeros(numel(plotdata), 3);
            ldata = zeros(numel(plotdata), 1);
            for pc = 1:numel(plotdata)
                try
                    ydata(:, pc) = plotdata{pc}.YData(:);
                    cdata(pc, :) = round(255.001 .* plotdata{pc}.Color);
                    ldata(pc) = plotdata{pc}.LineWidth;
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end
            end
        end

        % create output
        validx = (0:(numel(xdata)-1))';
        xout = sprintf('x,%d,%d\n', lsqueeze(([validx, xdata])'));
        yout = cell(1, size(ydata, 2));
        for pc = 1:numel(yout)
            ysym = sprintf('y%d', pc - 1);
            yout{pc} = sprintf([ysym ',%d,%d\n'], lsqueeze(([validx, ydata(:, pc)])'));
        end
        yout = cat(2, yout{:});
        cout = sprintf('c%d,0,%x%x%x\n', lsqueeze(([(0:(size(cdata, 1)-1))', cdata])'));
        lout = sprintf('l%d,0,%d\n', lsqueeze(([(0:(numel(ldata)-1))', ldata])'));
        tout = sprintf('symbol,index,value\n%s%s%s%s', xout, yout, cout, lout);
        binwrite(outfile, uint8(tout));
        
    % click on button
    case {'click'}

        % requires button
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3}) && ...
           ~isempty(regexpi(cmd{3}, '^bt_neuroelf')) && ...
            isfield(uics, cmd{3})

            % perform callback
            uic = uics.(cmd{3});
            cb = uic.Callback;
            jb = [];

            % background mode
            cbbground = any(strcmpi(cmd(3:end), 'background'));

            % override arguments
            if cbbground
                cmd(strcmpi(cmd, 'background')) = [];
            end
            if isa(cb, 'function_handle')
                cb = {cb};
            end
            if iscell(cb)
                for cpc = 4:numel(cmd)
                    if ~isempty(regexpi(cmd{cpc}, '^[\-\+]?\d+(\.\d+)?(e[\-\+]?\d+)?$'))
                        cmd{cpc} = str2double(cmd{cpc});
                    elseif strcmpi(cmd{cpc}, 'true')
                        cmd{cpc} = true;
                    elseif strcmpi(cmd{cpc}, 'false')
                        cmd{cpc} = false;
                    end
                end
                if numel(cmd) > 3
                    cb(2:(numel(cmd)-2)) = cmd(4:end);
                end
            end

            % background mode
            if cbbground
                try
                    jb = javacomponent(javax.swing.JButton('REMOTE'), ...
                        [1, 1, 100, 20], ne_gcfg.h.MainFigMLH);
                    set(jb, 'ActionPerformedCallback', cb);
                    jb.doClick();
                    delete(jb);
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                    if ~isempty(jb)
                        try
                            delete(jb);
                        catch ne_eo;
                            ne_gcfg.c.lasterr = ne_eo;
                        end
                    end
                end

            % foreground mode
            else
                if iscell(cb)
                    feval(cb{1}, uic.MLHandle, 0, cb{2:end});
                elseif ~isempty(cb)
                    if ischar(cb)
                        try
                            evalin('base', cb);
                        catch ne_eo;
                            ne_gcfg.c.lasterr = ne_eo;
                        end
                    else
                        feval(cb, uic.MLHandle, 0);
                    end
                end
            end
        end

        % read UI state
        uistate = true;

    % list directory
    case {'dir'}

        % requires folder
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3})
            folder = cmd{3};
        else
            folder = '/';
        end

        % get list
        try
            [flist, folder] = ne_remote_findfiles(uid, folder);
            flist = regexprep(flist, ['^' strrep(folder, '/', '\/') '\/'], '');
            flist = regexprep(flist, ['^' strrep(folder, '/', '\/')], '');
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            flist = {'../'};
        end

        % write out file
        binwrite(outfile, uint8([folder '/\n' gluetostring(flist, '\n')]));

    % open file
    case {'openfile'}

        % requires file
        s = 'BADFILE';
        ubf = rini.UserPrivBase{uid};
        if numel(cmd) > 2 && ...
            exist([ubf(1:end-1) cmd{3}], 'file') > 0
            try
                obj = fopen([ubf(1:end-1) cmd{3}], 'r');
                if obj < 1
                    error( ...
                        'neuroelf:remote:FileNotReadable', ...
                        'File %s not readable.', ...
                        cmd{3} ...
                    );
                end
                fclose(obj);
                obj = ne_openfile(0, 0, [ubf(1:end-1) cmd{3}]);
                drawnow;
                objtype = lower(obj.Filetype);
                switch (objtype)
                    case {'hdr'}
                        if isequal(obj, ne_gcfg.fcfg.SliceVar)
                            tvar = 'SliceVar';
                        else
                            tvar = 'StatsVar';
                        end
                    case {'dmr', 'fmr', 'head', 'msk', 'vmr', 'vtc'}
                        tvar = 'SliceVar';
                    case {'cmp', 'glm', 'vmp'}
                        tvar = 'StatsVar';
                    case {'srf'}
                        tvar = 'SurfVar';
                    case {'smp'}
                        tvar = 'SurfStatsVar';
                end
                tvud = ne_gcfg.h.(tvar).UserData;
                tdum = ne_gcfg.c.ini.Remote.(sprintf('Max%ss', tvar));
                if size(tvud, 1) == (tdum + 1)
                    try
                        tvud{1, 4}.ClearObject;
                    catch ne_eo;
                        ne_gcfg.c.lasterr = ne_eo;
                    end
                    try
                        obj.Browse;
                    catch ne_eo;
                        ne_gcfg.c.lasterr = ne_eo;
                    end
                end
                s = 'OK';
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                s = ne_eo.message;
            end
        end

        % write status
        binwrite(outfile, uint8(s));

    % set a text field
    case {'setstring'}

        % requires field (and string)
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3}) && ...
           ~isempty(regexpi(cmd{3}, '^(ed|txt?)_neuroelf')) && ...
            isfield(uics, cmd{3})

            % string to set
            if numel(cmd) > 3
                newstring = cmd{4};
            else
                newstring = '';
            end

            % set string
            uics.(cmd{3}).String = newstring;

            % perform callback
            cb = uics.(cmd{3}).Callback;
            try
                if iscell(cb)
                    feval(cb{1}, uics.(cmd{3}).MLHandle, 0, cb{2:end});
                elseif ~isempty(cb)
                    if ischar(cb)
                        try
                            evalin('base', cb);
                        catch ne_eo;
                            ne_gcfg.c.lasterr = ne_eo;
                        end
                    else
                        feval(cb, uics.(cmd{3}).MLHandle, 0);
                    end
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end

        % read UI state ?
        if numel(cmd) > 4 && ...
            strcmp(cmd{5}, '1')
            uistate = true;
        else
            binwrite(outfile, uint8('OK'));
        end

    % select from dropdown/listbox
    case {'setvalue'}

        % requires valid controls
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3}) && ...
            isfield(uics, cmd{3})
            uic = uics.(cmd{3});

            % requires new value
            if numel(cmd) > 3
                try

                    % parse new value
                    newidx = u8str2double(cmd{4});
                    newidx = newidx(:);

                    % type of control
                    if any(strcmpi(uic.Style, {'listbox', 'popupmenu'}))
                        ustr = uic.String;
                        if ~iscell(ustr)
                            ustr = cellstr(ustr);
                        end
                        minidx = 1;
                        maxidx = numel(ustr);
                    else
                        minidx = 0;
                        maxidx = 1;
                    end
                    newidx(isinf(newidx) | isnan(newidx) | newidx < minidx | ...
                        newidx > maxidx | newidx ~= fix(newidx)) = [];

                    % update
                    if strcmpi(uic.Style, 'listbox') || ...
                        numel(newidx) == 1
                        uistate = true;
                    end

                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                    newidx = uic.Value(:);
                end
            else
                newidx = uic.Value(:);
            end

            % update value
            uic.Value = newidx;
        end

        % update UI state ?
        if uistate

            % callback
            cb = uic.Callback;
            try
                if iscell(cb)
                    feval(cb{1}, uic.MLHandle, 0, cb{2:end});
                elseif ~isempty(cb)
                    if ischar(cb)
                        try
                            evalin('base', cb);
                        catch ne_eo;
                            ne_gcfg.c.lasterr = ne_eo;
                        end
                    else
                        feval(cb, uic.MLHandle, 0);
                    end
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        else
            binwrite(outfile, uint8('NOTHING'));
        end

    % set volume number
    case {'setvolume'}

        % requires valid volume
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3}) && ...
            all(cmd{3} >= '0' & cmd{3} <= '9') && ...
            strcmpi(uics.SL_NeuroElf_TempPos.Enable, 'on')
            newvol = str2double(cmd{3});
            if newvol >= 1 && ...
                newvol <= uics.SL_NeuroElf_TempPos.Max
                uics.SL_NeuroElf_TempPos.Value = newvol;
                feval(uics.SL_NeuroElf_TempPos.Callback, ...
                    uics.SL_NeuroElf_TempPos.MLHandle, 0);
            end
        end

        % update state
        uistate = true;

    % get list of slicevars
    case {'dd_neuroelf_varlist'}

        % update slicevar?
        svar = ne_gcfg.h.SliceVar;
        svaru = svar.UserData;
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3}) && ...
            all(cmd{3} >= '0' && cmd{3} <= '9')
            cmd{3} = str2double(cmd{3});
            if cmd{3} > 0 && ...
                cmd{3} <= size(svaru, 1)
                svar.Value = cmd{3};
                ne_setcvar;
            end
        end

        % get list of all loaded objects in workspace
        svari = svar.Value;
        svars = svaru(:, [1, 3]);

        % create %-separated string
        out = uint8(sprintf('%d%%%d%%%s', size(svars, 1), ...
            svari, gluetostringc(lsqueeze(svars'), '%')));
        binwrite(outfile, out);

    % get list of statsvars
    case {'dd_neuroelf_statlist'}

        % update statsvar?
        stvar = ne_gcfg.h.StatsVar;
        stvaru = stvar.UserData;
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3}) && ...
            all(cmd{3} >= '0' & cmd{3} <= '9')
            cmd{3} = str2double(cmd{3});
            if cmd{3} > 0 && ...
                cmd{3} <= size(stvaru, 1)
                stvar.Value = cmd{3};
                ne_setcstats;
            end
        end

        % get list of all loaded objects in workspace
        stvari = stvar.Value;
        stvars = stvaru(:, [1, 3]);

        % create %-separated string
        out = uint8(sprintf('%d%%%d%%%s', size(stvars, 1), ...
            stvari, gluetostringc(lsqueeze(stvars'), '%')));
        binwrite(outfile, out);

    % get list of statsvar maps
    case {'lb_neuroelf_statmaps'}

        % update statsvar?
        stmap = ne_gcfg.h.StatsVarMaps;
        stmaps = stmap.String;
        if ~iscell(stmaps)
            stmaps = cellstr(stmaps);
        end
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3})
            try
                newmapi = u8str2double(cmd{3});
                if ~isempty(newmapi) && ...
                   ~any(isinf(newmapi(:)) | isnan(newmapi(:)) | ...
                        newmapi(:) < 1 | newmapi(:) > numel(stmaps))
                    stmap.Value = newmapi(:);
                    ne_setcstatmap;
                elseif isequal(newmapi, 0)
                    stmap.Value = [];
                    ne_setcstatmap;
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        elseif numel(cmd) > 2
            try
                stmap.Value = [];
                ne_setcstatmap;
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end

        % get list of all loaded objects in workspace
        stmapi = stmap.Value;

        % create %-separated string
        out = uint8(sprintf('%d%%%s%%%s', numel(stmaps), ...
            sprintf('%d,', stmapi(:)'), gluetostringc(stmaps, '%')));
        binwrite(outfile, out);

    % set p-value
    case {'dd_neuroelf_statsetp'}

        % only if control is enabled
        if numel(cmd) > 2 && ...
           ~isempty(cmd{3}) && ...
            all(cmd{3} >= '0' & cmd{3} <= '9') && ...
            strcmpi(ne_gcfg.h.Stats.PThresh.Enable, 'on')
            pth = ne_gcfg.h.Stats.PThresh;
            pti = str2double(cmd{3});
            pts = pth.String;
            if ~iscell(pts)
                pts = cellstring(pts);
            end
            if pti > 0 && ...
                pti <= numel(pts)
                try
                    pth.Value = pti;
                    feval(pth.Callback);
                    out = 'ok';
                catch ne_eo;
                    out = ne_eo.message;
                end
            end
        end
        binwrite(outfile, uint8(out));

    % get UI screenshot
    case {'uishot'}

        % set new position
        if ne_gcfg.fcfg.page < 3 && ...
            numel(cmd) > 4
            try
                xcrd = str2double(cmd{3});
                ycrd = str2double(cmd{4});
                zcrd = str2double(cmd{5});
                if numel(xcrd) == 1 && ...
                   ~isinf(xcrd) && ...
                   ~isnan(xcrd) && ...
                    xcrd >= -128 && ...
                    xcrd <= 128 && ...
                    numel(ycrd) == 1 && ...
                   ~isinf(ycrd) && ...
                   ~isnan(ycrd) && ...
                    ycrd >= -128 && ...
                    ycrd <= 128 && ...
                    numel(zcrd) == 1 && ...
                   ~isinf(zcrd) && ...
                   ~isnan(zcrd) && ...
                    zcrd >= -128 && ...
                    zcrd <= 128
                    ne_gcfg.fcfg.cpos = [xcrd, ycrd, zcrd];
                    ne_setslicepos;
                    drawnow;
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end

        % write out screenshot (only works correctly if on top)
        mfh = ne_gcfg.h.MainFigMLH;
        ss = getframe(mfh);
        imwrite(ss.cdata, outfile, imf, imq{:});

    % UI state
    case {'uistate'}

        % handle later
        uistate = true;

    % no valid command
    otherwise

        % inform user
        if numel(cmd) > 2
            opts = sprintf('%s,', cmd{3:end});
            opts(end) = [];
        else
            opts = '';
        end
        binwrite(outfile, uint8(sprintf('Unknown command: %s(%s)', cmd{2}, opts)));
end

% handle UI state
if uistate

    % get UI components
    uicf = fieldnames(uics);
    uicf(cellfun('isempty', regexpi(uicf, '^(BT_|CB_|DD_|ED_|LB_|TX)'))) = [];
    state = cell(4, numel(uicf));
    state(1, :) = uicf';
    state(2:4, :) = {''};
    for uc = 1:numel(uicf)
        uic = uics.(uicf{uc});
        switch (uicf{uc}(1:2))
            case {'CB'}
                uicval = uic.Value;
                if ~isempty(uicval)
                    uicval = sprintf('%d,', uicval(:)');
                    state{3, uc} = uicval(1:end-1);
                end
            case {'DD', 'LB'}
                uicstr = uic.String;
                if ~iscell(uicstr)
                    uicstr = cellstr(uicstr);
                end
                uicstr = gluetostringc(uicstr, '\n');
                state{2, uc} = uicstr;
                uicval = uic.Value;
                if ~isempty(uicval)
                    uicval = sprintf('%d,', uicval(:)');
                    state{3, uc} = uicval(1:end-1);
                end
            case {'ED', 'TX'}
                uicstr = char(uint8(char(uic.String)));
                if size(uicstr, 1) > 1
                    uicstr = gluetostringc(cellstr(uicstr), '\n');
                end
                uicstr(uicstr < 32 | uicstr > 127) = ' ';
                state{2, uc} = uicstr;
        end
        state{4, uc} = uic.Enable;
    end
    state = state(:)';
    state = sprintf('%s%%', state{:});
    binwrite(outfile, uint8(state));
end

% make filename as requested
movefile(outfile, outfile(1:end-1));

% remove timer and command
try
    delete(tdata{2});
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
cpos = find(ne_gcfg.c.remotecfg.commandid(:, 1) == cid);
ne_gcfg.c.remotecfg.commandid(cpos, :) = [];
ne_gcfg.c.remotecfg.commands(cpos, :) = [];

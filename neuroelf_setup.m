function neuroelf_setup
% neuroelf_setup  - NeuroElf post installation setup
%
% FORMAT:       neuroelf_setup
%
% No input/output fields.
%
% Using: ddeblank, gluetostring, importvmrfromanalyze, license,
%        neuroelf_makelibs, neuroelf_splash, splittocell, tdlocal2.

% Version:  v1.1
% Build:    17080709
% Date:     Aug-07 2017, 9:33 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2017, Jochen Weber
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

% global variables
global ne_methods ne_ui;

% close all files
fclose all;

% load neuroelf library (for other functions)
nelf = neuroelf;

% potentially remove older BVQXtools version from path first
if exist('BVQXfile', 'file') > 0
    bf = which('BVQXfile');
    if ~isempty(strfind(bf, '@BVQXfile'))
        disp('Removing BVQXtools from MATLAB path!');
        while ~isempty(strfind(bf, '@BVQXfile'))
            bf = fileparts(bf);
        end
        try
            rmpath(bf);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
        rehash path;
    end
end

% remove existing xff cache file
xfc = [neuroelf_path('cache') filesep 'xffcache.mat'];
if exist(xfc, 'file') == 2
    try
        delete(xfc);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% get current pwd and switch..
opwd = pwd;
npwd = neuroelf_path;
cd(npwd);
if any(npwd(end) == '0123456789')
    npbuild = regexprep(npwd, '^.*_(\d+)$', '$1');
else
    npbuild = 'release';
end

% write out "done" flag file
ne_methods.asciiwrite([neuroelf_path('cache') filesep 'setup_' npbuild '.done'], 'done');

try
    % get screen size
    rs = get(0, 'ScreenSize');

    % if screen is large enough (and available)
    if all(rs(3:4) >= [640, 480])

        % open splash figure
        hFig = xfigure([neuroelf_path('tfg') '/splash.tfg']);

        % get handles to progress text, bar, and children
        hPrg = hFig.TagStruct.LB_Progress;
        hBar = hFig.TagStruct.PB_Progress;
        imgh = hFig.TagStruct.IM_Splash.Children;
        imgh = imgh(strcmpi('image', get(imgh, 'Type')));

        % initialize splash animation
        imgt = ne_methods.neuroelf_splash(imgh);

        % show figure
        hFig.HandleVisibility = 'callback';
        hFig.Visible = 'on';
        drawnow;

        % and start animation
        start(imgt);

    % with too small (or non existing) screen
    else
        hPrg = [];
        hBar = [];
    end

% deal with an error here
catch ne_eo;
    warning('neuroelf:xfigureError:errorOpeningSplashScreen', ...
        'Error with xfigure class: %s.', ne_eo.message);
    hPrg = [];
    hBar = [];
end

% print banner
splittocell = ne_methods.splittocell;
nelftv = neuroelf_version;
nelftb = sprintf('%d', neuroelf_build);
banner = [ ...
    char(10), ...
    '==========================================================', char(10), ...
    char(10), ...
    'NeuroElf v' nelftv ', build ' nelftb ' (unified edition)', char(10), ...
    char(10), ...
    'Please direct any bug reports and/or feature requests directly to:', char(10), ...
    char(10), ...
    'Jochen Weber <jw2661@columbia.edu>', char(10), ...
    char(10), ...
    'Thanks must go to Chih-Jen Lin for allowing the re-use of the', char(10), ...
    'libSVM code repository :)', char(10), ...
    char(10), ...
    'Another thank-you is going to Aapo Hyvarinen for letting me integrate the', char(10), ...
    'FastICA algorithm into the toolbox!', char(10), ...
    char(10), ...
    'All original source code (which contains this license statement) is', char(10), ...
    char(10), ...
    '==========================================================', char(10), ...
    ne_methods.license(), ...
    '==========================================================', char(10), ...
    char(10)];
instprog(hPrg, hBar, banner, 0.05, '', splittocell);

% build libs (MEX)
mxt = lower(mexext);
mxs = strrep(mxt, 'mex', '');
try
    instprog(hPrg, hBar, ...
        ['Checking/compiling MEX files for current platform (' upper(mxs) ')...'], ...
         0.05, 'Checking MEX files...', splittocell);
    ne_methods.neuroelf_makelibs(hPrg, hBar, [0.05, 0.5]);

% if an error occurred
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    disp(['There was some error compiling MEX files and NeuroElf does', char(10), ...
          'not come with pre-compiled files for your platform: ', upper(mxs)]);
end

% make sure path and toolbox cache are updated
try
    rehash path;
    rehash toolboxcache;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% check "helper classes"
try
    instprog(hPrg, hBar, 'Checking xini class...', 0.5, 'Checking classes...', splittocell);
    v = xini;
    if ~isxini(v)
        error('CLASS_FAILURE');
    end
    try

        % copy "clean" configuration file if file doesn't exist
        cfgpath = [neuroelf_path('config') '/'];
        if exist([cfgpath 'neuroelf.ini'], 'file') ~= 2
            copyfile([cfgpath 'neuroelf_clean.ini'], [cfgpath 'neuroelf.ini']);
        end

        % load configuration file and original file
        v = xini([cfgpath 'neuroelf.ini'], 'convert');
        vc = xini([cfgpath 'neuroelf_clean.ini'], 'convert');

        % copy values over (ensure that erroneous values don't break things)
        vs = vc.GetSections;
        for sc = 1:numel(vs)
            ss = vc.(vs{sc});
            sn = fieldnames(ss);
            for nc = 1:numel(sn)
                v.(vs{sc}).(sn{nc}) = vc.(vs{sc}).(sn{nc});
            end
        end

        % get information on tools
        vtools = v.Tools;
        
        % save (update) ini file
        v.Save;
        
        % release both objects
        v.Release;
        vc.Release;
    catch ne_eo;
        vtools = struct('dcm2nii', struct( ...
            'a64',    'dcm2nii.a64', ...
            'glx',    'dcm2nii.glx', ...
            'maci',   'dcm2nii.osx', ...
            'maci64', 'dcm2nii.osx', ...
            'w32',    'dcm2nii.exe', ...
            'w64',    'dcm2nii.exe'));
        rethrow(ne_eo);
    end

    % check xff.ini
    try
        v = xini([neuroelf_path('config') '/xff.ini'], 'convert');

        % correctly request new temporary folder name (needs fix for 2017b+)
        tdir = tempdir;
        if numel(tdir) > 8 && strcmpi(tdir(1:8), '/private') && exist('/tmp')
            tdir = '/tmp';
        end
        if ~isempty(hPrg)
            vui = inputdlg({'Use the following folder as temporary disk space:'}, ...
                'NeuroElf - temp folder config', 1, {['  ' tdir]});
        else
            vui = {input(sprintf( ...
                'Please enter the folder used for temporary disk space (%s): ', ...
                tdir), 's')};
        end

        % with new temp dir entry
        if iscell(vui) && numel(vui) == 1 && ~isempty(ne_methods.ddeblank(vui{1}))

            % remove leading/trailing spaces
            vui = ne_methods.ddeblank(vui{1});

            % if directory exists
            if ~isempty(vui) && ~strcmpi(vui, tempdir) && exist(vui, 'dir') == 7

                % make sure it ends in filesep
                if vui(end) ~= '/' && vui(end) ~= '\'
                    vui(end+1) = filesep;
                end

                % update
                v.GZip.TempDir = vui;
                v.Save;
            end
        end

        % release object
        v.Release;
    catch ne_eo;
        rethrow(ne_eo);
    end
catch ne_eo;
    warning('neuroelf:xiniError:objectMethodFailed', ...
        'Error with xini class: %s.', ne_eo.message);
end

% check xff class
try
    instprog(hPrg, hBar, 'Checking xff class...', 0.52, '', splittocell);

    % instantiate factory
    v = xff();

    % get available extensions (file types)
    x = sort(fieldnames(v.Extensions));

    % and print out a Rx8 matrix with extensions
    nx = numel(x);
    if mod(nx, 8) > 0
        x(end+1:8 * ceil(nx / 8)) = {''};
    end
    x8 = cell(1, numel(x) / 8);
    for xc = 1:numel(x8)
        x8{xc} = sprintf('%5s ', x{(xc - 1) * 8 + 1:xc*8});
    end

    % add to output
    info = [char(10), ...
        sprintf('Current release supports %d xff filetypes:', nx), char(10), ...
    	char(10), ...
        upper(ne_methods.gluetostring(x8, char(10))), char(10)];

    % and keep track of which types fail (new:XXX)
    xok = false(size(x));
    for xc = 1:numel(xok)
        try
            xobj = {[]};
            if ~isempty(x{xc})
                xobj{1} = xff(sprintf('new:%s', x{xc}));
                xobj{1}.ClearObject;
            end
            xok(xc) = true;
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            clearxffobjects(xobj);
        end
    end

    % set ROOT object to true (will always fail!)
    xok(strcmpi('root', x(:)')) = true;

    % print out failed types (code needs fixing!)
    if sum(xok) < numel(xok)
        info = [info, char(10), sprintf('xff(''new:*'') failed on:%s.', ...
                sprintf(' %s', x{find(~xok)})), char([10, 10])];

    % or report full success
    else
        info = [info, char(10), ...
            'xff(''new:*'') working for all object types.', char([10, 10])];
    end
    instprog(hPrg, hBar, info, 0.7, '', splittocell);
catch ne_eo;
    warning('neuroelf:xff:generalError', 'Error checking xff class: %s.', ne_eo.message);
end

% create colin brain vmr from colin.img and mask
% see http://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach for info
try

    % path
    cpath = [neuroelf_path('colin') filesep];

    % check if VMR already exists
    if exist([cpath 'colin.vmr'], 'file') ~= 2
        instprog(hPrg, hBar, 'Creating colin.vmr from colin.hdr/img...', 0.70, ...
            'Patching/creating required binary files', splittocell);

        % create if necessary
        colin = ne_methods.importvmrfromanalyze([cpath 'colin.img'], 'lanczos3', [0.01, 0.99]);
        colin.SaveAs([cpath 'colin.vmr']);
    else

        % or load
        colin = xff([cpath 'colin.vmr']);
    end

    % clear colin
    colin.ClearObject;
    colin = [];
catch ne_eo;
    warning('neuroelf:xffError:errorCreatingFile', ...
        'Error creating colin_brain.vmr: %s.', ne_eo.message);
    colin = [];
end

% test tdlocal2 w/ files
try
    instprog(hPrg, hBar, 'Testing tdlocal2 (local TD database)...', ...
        0.9, 'Testing tdlocal2', splittocell);
    ne_methods.tdlocal2(2, 0, 0, 0);
catch ne_eo;
    warning('neuroelf:TDLocalError', 'Error using tdlocal2: %s.', ne_eo.message);
end

% make sure dcm2nii binary is marked executable
try
    if ~ispc
        instprog(hPrg, hBar, 'Setting executable permissions...', ...
            0.9, 'Setting executable permissions.', splittocell);
        system(sprintf('chmod a+x ''%s''', [neuroelf_path('contrib') filesep ...
            vtools.dcm2nii.(strrep(lower(mexext), 'mex', ''))]));
    end
catch ne_eo;
    warning('neuroelf:chmodError', ...
        'Error changing (executable) permission of file(s): %s', ne_eo.message);
end

% unload colin
if numel(colin) == 1 && isxff(colin, true)
    try
        colin.ClearObject;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% try to add path
try
    addpath(neuroelf_path);
    savepath;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% stop splash after another second of pause
if isfield(ne_ui, 'splash') && isstruct(ne_ui.splash) && isfield(ne_ui.splash, 'timer')
    try
        for xc = 1:25
            pause(0.04);
            drawnow;
        end
        stop(ne_ui.splash.timer);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% clean up...
instprog(hPrg, hBar, 'Done.', 1, 'Cleaning up...', splittocell);
pause(0.2);
if ~isempty(hPrg)
    disp(ne_methods.gluetostring(hPrg.String, char(10)));
end
try
    xfigure(xfigure, 'DeleteAllFigures');
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% call neuroelf_makefiles?
try
    if ~isempty(hPrg)
        vui = questdlg('Create all files with neuroelf_makefiles?', ...
            'NeuroElf - makefiles ?', 'Yes');
    else
        vui = input('Create all files with neuroelf_makefiles (Y)es/(N)o: ', 's');
    end
    if ~isempty(vui) && lower(vui(1)) == 'y'
        ne_methods.neuroelf_makefiles('all');
        xfigure(xfigure, 'DeleteAllFigures');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% switch back
try
    cd(opwd);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% rehash everything
evalin('base', 'rehash toolboxcache;', '');
evalin('base', 'rehash;', '');
evalin('base', 'clear xff;', '');


% sub function for progress
function instprog(hedt, hbar, txt, p, ptxt, splittocell)
if ~isempty(hedt)
    if ~isempty(txt)
        str = hedt.String;
        if ischar(str) && size(str, 1) > 1
            str = cellstr(str);
        elseif ischar(str)
            str = splittocell(str, char(10));
        end
        stxt = splittocell(txt, char(10));
        str = [str(:); stxt(:)];
        hedt.String = str;
        hedt.Value = numel(str);
        hedt.ListboxTop = max(1, numel(str) - 15);
    end
    if nargin > 4 && ~isempty(ptxt)
        hbar.Progress(p, ptxt);
    elseif nargin > 3
        hbar.Progress(p);
    end
    drawnow;
else
    disp(txt);
    pause(0.001);
end

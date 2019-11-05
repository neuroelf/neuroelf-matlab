function ne_initwindow(MainFig)
% FUNCTION ne_initwindow: initialize main window UI

% Version:  v1.1
% Build:    16060916
% Date:     Jun-09 2016, 4:45 PM EST
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

% global variable
global ne_gcfg;

% reject if not a valid figure
if ~isxfigure(MainFig)
    return;
end

% don't re-init
if numel(ne_gcfg) == 1 && isstruct(ne_gcfg) && isfield(ne_gcfg, 'c') && isstruct(ne_gcfg.c) && ...
    numel(ne_gcfg.c) == 1 && isfield(ne_gcfg.c, 'init') && islogical(ne_gcfg.c.init) && ne_gcfg.c.init
    return;
end

% set root property (to work properly) and get all tags (by name)
set(0, 'ShowHiddenHandles', 'off');
tags = MainFig.TagStruct;

% create and populate NeuroElf's global configuration
fprintf(' core config...');
pause(0.001);
c = struct;

% atlas labels (also initializes local TAL Daemon Database)
c.atlas.tal = struct('labels', {tdlocal2(8, 'labels')});
c.atlas.tal.shorts = regexprep(lower(c.atlas.tal.labels), '\s+', '_');

% list of blocking callbacks used to reject double calls and keep the
% window open as long as needed
c.blockcb = {};
c.breakcb = {};

% keep track of mouse button down
c.btdown = [];
c.btdoup = false;
c.btdwnf = '';
c.btupact = {};

% create list of public callbacks (begining in "ne_")
c.callbacks = struct;
flist = findfiles([neuroelf_path filesep '@neuroelf' filesep 'private'], 'ne_*.m', 'depth=1');
for cbc = 1:numel(flist)
    [fpath, fname] = fileparts(flist{cbc});
    c.callbacks.(fname(4:end)) = eval(['@' fname]);
end

% manually add those that don't fit the pattern
c.callbacks.slicevar = @cv_slicevar;
c.callbacks.statsvar = @cv_statsvar;
c.callbacks.varlist = @cv_varlist;

% echo calls to prompt (is set from ini.MainFig.EchoCalls later)
c.echo = false;

% extended map names (is set from ini.Statistics.ExtendedMapNames later)
c.extmapnames = false;

% used for simple test if currently in a modal callback
c.incb = false;

% load configuration file
try
    c.ini = xini([neuroelf_path('config') '/neuroelf.ini'], 'convert');
    if ~isxini(c.ini)
        error('neuroelf:GUI:configError', 'Error loading neuroelf.ini');
    end

    % get complete content (for faster access during this function)
    cini = c.ini.GetComplete;

    % quick content-check for necessary sections
    if ~isfield(cini, 'MainFig') || ~isfield(cini, 'Satellites') || ...
       ~isfield(cini, 'BetaPlot') || ~isfield(cini, 'Children') || ...
       ~isfield(cini, 'Drawing') || ~isfield(cini, 'Mediation') || ...
       ~isfield(cini, 'MKDA') || ~isfield(cini, 'PLPPlot') || ...
       ~isfield(cini, 'RecentFiles') || ~isfield(cini, 'Remote') || ...
       ~isfield(cini, 'Render') || ~isfield(cini, 'Statistics') || ...
       ~isfield(cini, 'Surface') || ~isfield(cini, 'SurfMontage') || ...
       ~isfield(cini, 'SurfMontageConfigs') || ~isfield(cini, 'SurfMontageElems') || ...
       ~isfield(cini, 'TCMovie') || ~isfield(cini, 'Tools') || ~isfield(cini, 'VMP')
        error('neuroelf:GUI:configError', ...
            'Configuration file (_core/config/neuroelf.ini) corrupted.');
    end
catch ne_eo;
    abt.Delete;
    rethrow(ne_eo);
end

% initialization status (will be set to true once run through)
c.init = false;

% last error caught (for diagnostics)
c.lasterr = [];

% last update (useful to sync across remote-linked windows)
c.lastupd = -1;

% linked browsing
c.linked = false;

% matlab version
c.mlversion = str2double(splittocell(regexprep(version, '\s+.*$', ''), '.'));
c.mlversion = sum([100, 1] .* c.mlversion(1:2));

% get OS type and if it's Mac
c.ostype = ostype;
c.ostypemac = strcmpi(c.ostype.machine, 'mac');

% physio file-browsing filter configuration (order of file types)
c.physio = struct('filter', {{ ...
     '*.acq', 'Acknowledge files (*.acq)'; ...
     '*.mat', 'MAT files (*.mat)'; ...
     '*.ntt;*.log;*.txt;*.csv', 'Numeric text files (*.ntt, *.log, *.txt, *.csv)'; ...
     '*.*',   'All files (*.*)'}});

% progress information (visible flag, amount [0 .. 1], text)
c.progress = {false, 0, 'NeuroElf - progress'};
 
% remote (folder, status, and configuration)
if isempty(cini.Remote.ScanFolder) || exist(cini.Remote.ScanFolder, 'dir') ~= 7
    c.ini.Remote.ScanFolder = neuroelf_path('cache');
end
c.remote = false;
c.remotecfg = struct( ...
    'cmdcount',  0, ...
    'cmdfiles',  {cell(0, 6)}, ...
    'commandid', zeros(0, 2), ...
    'commands',  {cell(0, 6)}, ...
    'gcwcount',  cini.Remote.GCWaitCount, ...
    'lastcmd',   -1, ...
    'lastscan',  -1, ...
    'logfile',   -1, ...
    'scanning',  false, ...
    'scanpath',  c.ini.Remote.ScanFolder, ...
    'scantimer', [], ...
    'session',   struct('S000000', struct), ...
    'stopping',  false, ...
    'stoptimer', false);

% resize and render preview calls in progress
c.resize = false;
c.rpreview = false;
c.satresize = false;

% sampled values (e.g. used for sending data via remote)
c.svals = [];

% title information (anatomical, stats, and maps information, same for SRF)
c.title = {'anatomical.vmr', '', ''; 'surface.srf', '', ''};

% get reference to XFF factory object (to avoid further calls)
c.xff = xff();

% initialize global storage structure
% .c       - core configuration (already initialized)
% .cc      - children config (one field per UI, such as BetaPlots BP######)
% .fcfg    - figure configuration (reflecting controls, etc.)
% .h       - handles (used for increased performance and readability)
% .lut     - main LUT object (for non-VMP StatsVar display)
% .poi     - POI object reflecting the cluster list
% .tio     - transimg objects for main window slice display
% .voi     - VOI object reflecting the cluster list
% .w       - workspace (variables in NeuroElf control)
% .wc      - workspace control (variables to be cleared upon exit)
fprintf(' UI config...');
pause(0.001);
ne_gcfg = struct('c', c, 'cc', struct, 'fcfg', struct, 'h', struct, ...
    'lut', [], 'poi', [], 'tio', struct, 'voi', [], 'w', struct, 'wc', struct);

% begin with figure default configuration
fcfg = ne_gcfg.fcfg;

% configs for ANCOVA, AP, conman, console, MDM, MKDA, mediation, render, surf+vis montage
fcfg.AC = [];
fcfg.AP = [];
fcfg.CM = [];
fcfg.Console = [];
fcfg.MDM = [];
fcfg.MKDA = [];
fcfg.RM = [];
fcfg.Render = [];
fcfg.SurfMontage = [];
fcfg.VisMontage = [];

% alphasim thresholds
fcfg.asimthr = cini.Tools.alphasim.Thresholds;

% crosshair visible and color
fcfg.chair = cini.MainFig.Crosshairs;
fcfg.chcol = cini.MainFig.CrosshairColor;

% initialize position of children
fcfg.chpos = {{}, {}, zeros(0, 4)};

% cluster connectivity setting, limit (radius), and sorting
fcfg.clconn = cini.Statistics.ClusterConnectivity;
fcfg.clim = cini.Statistics.ExtractRadius;
fcfg.clsort = cini.Statistics.Sorting;

% current position (in TAL coordinates and order)
fcfg.cpos = [0, 0, 0];

% cursor stepsize (depending on dataset)
fcfg.cstep = 1;

% current drawing "direction" (slice) and direction order (cycling)
fcfg.ddir = [1, 2];
fcfg.dirorder = {'sag', 'cor', 'tra'};

% current DTI object
fcfg.dti = [];

% full figure size
fcfg.fullsize = MainFig.Position(3:4);
fcfg.fullsized = true;
fcfg.fullsizes = tags.TX_NeuroElf_SValues.Position(1:2) + [-10, 4];

% get all handles required for resize events
set(0, 'ShowHiddenHandles', 'on');
fcfg.fullsizexh = lsqueeze(get(MainFig.MLHandle, 'Children'));
set(0, 'ShowHiddenHandles', 'off');

% make sure to remove uimenu objects
fctype = get(fcfg.fullsizexh, 'Type');
fcfg.fullsizexh(strcmpi(fctype, 'uimenu') | strcmpi(fctype, 'annotationpane')) = [];

% prepare sizes array with settings
fcfg.fullsizex = zeros(numel(fcfg.fullsizexh), 10);

% record current position of all controls
fcpos = get(fcfg.fullsizexh(:, 1), 'Position');
fcpos = cat(1, fcpos{:});

% store in columns 3 through 6
fcfg.fullsizex(:, 3:6) = fcpos;

% and compute small (swapped) size
fcfg.fullsizex(:, 7:10) = fcpos - repmat([fcfg.fullsizes, 0, 0], size(fcpos, 1), 1);

% central controls, shift a quarter up + a quarter right
fcfg.fullsizex( ...
    fcpos(:, 1) > (tags.FR_NeuroElf_vertdivide.Position(1) + 4) & ...
    fcpos(:, 1) < (tags.CB_NeuroElf_Interpolate.Position(1) + 2) & ...
    fcpos(:, 2) > (tags.RB_NeuroElf_LUTColor.Position(2) - 4) & ...
    fcpos(:, 2) < (tags.ED_NeuroElf_BVSX.Position(2) + 16), 2) = 10;

% surface-space stats controls
fcfg.fullsizex( ...
    fcpos(:, 1) > (tags.FR_NeuroElf_vertdivide.Position(1) + 4) & ...
    fcpos(:, 2) < (tags.CB_NeuroElf_SrfNegStat.Position(2) + 16), 2) = 13;

% top left and left-side buttons (shift up)
fcfg.fullsizex( ...
    fcpos(:, 1) < (tags.FR_NeuroElf_vertdivide.Position(1)) & ...
    fcpos(:, 2) > (tags.LB_NeuroElf_clusters.Position(2) - 2), 2) = 1;
fcfg.fullsizex( ...
    fcpos(:, 1) < (tags.BT_NeuroElf_slvartrf.Position(1) - 2) & ...
    fcpos(:, 2) > (tags.LB_NeuroElf_clusters.Position(2) - 2), 2) = 11;
fcfg.fullsizex( ...
    fcpos(:, 1) > (tags.FR_NeuroElf_vertdivide.Position(1)) & ...
    fcpos(:, 1) < (tags.BT_NeuroElf_draw0.Position(1) + 2) & ...
    fcpos(:, 2) > (tags.BT_NeuroElf_showv16.Position(2) - 1), 2) = 1;

% right-side buttons (shift right, up)
fcfg.fullsizex( ...
    fcpos(:, 1) > (tags.BT_NeuroElf_undock.Position(1) - 2) & ...
    fcpos(:, 2) > (tags.BT_NeuroElf_render.Position(2) - 2), 2) = 3;

% right-bottom button (shift right)
fcfg.fullsizex( ...
    fcpos(:, 1) > (tags.BT_NeuroElf_undock.Position(1) - 2) & ...
    fcpos(:, 2) < (tags.BT_NeuroElf_render.Position(2) - 2), 2) = 5;

% cluster table (text) and dividing frame (size up) 1
fcfg.fullsizex(fcfg.fullsizexh(:, 1) == tags.ED_NeuroElf_clusters.MLHandle, 2) = 12;
fcfg.fullsizex(fcfg.fullsizexh(:, 1) == tags.FR_NeuroElf_vertdivide.MLHandle, 2) = 2;

% axes (full slices, surface, render) resize +width +height
fcfg.fullsizex( ...
    fcfg.fullsizexh(:, 1) == tags.IM_NeuroElf_Slice_Zoom.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.AX_NeuroElf_Slice_Zoom.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.IM_NeuroElf_Slice_Rend.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.AX_NeuroElf_Slice_Rend.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.AX_NeuroElf_Surface.MLHandle, 2) = 4;

% some controls (resize +width)
fcfg.fullsizex( ...
    fcfg.fullsizexh(:, 1) == tags.LB_NeuroElf_Scenery.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.ED_NeuroElf_SrfViewPnt.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.TX_NeuroElf_SValues.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.AX_NeuroElf_TC_Plot.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.AX_NeuroElf_STC_Plot.MLHandle, 2) = 6;

% left slice, shift half up, resize half (both)
fcfg.fullsizex( ...
    fcfg.fullsizexh(:, 1) == tags.IM_NeuroElf_Slice_SAG.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.AX_NeuroElf_Slice_SAG.MLHandle, 2) = 7;

% right-top slice, shift hal
fcfg.fullsizex( ...
    fcfg.fullsizexh(:, 1) == tags.IM_NeuroElf_Slice_COR.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.AX_NeuroElf_Slice_COR.MLHandle, 2) = 8;

% right-top slice, shift hal
fcfg.fullsizex( ...
    fcfg.fullsizexh(:, 1) == tags.IM_NeuroElf_Slice_TRA.MLHandle | ...
    fcfg.fullsizexh(:, 1) == tags.AX_NeuroElf_Slice_TRA.MLHandle, 2) = 9;

% current GLM
fcfg.glm = [];

% gray-scale LUT
fcfg.graylut = [];

% histogram position
fcfg.histpos = tags.AX_NeuroElf_Slice_Hist.Position;
fcfg.histpos(3:4) = fcfg.histpos(1:2) + fcfg.histpos(3:4);
fcfg.histset = 0;
fcfg.histval = [0, 1];

% interpolation method (for statistical vars, SliceVar always linear)
fcfg.imethod = cini.Statistics.InterpMethod;

% instantaneous seed-correlations
fcfg.instscorr = cini.Statistics.InstantSeedCorr;
fcfg.instscrad = cini.Statistics.InstantSeedRadius;
fcfg.instscvar = [];

% join stats mode (vs. overlaying/overriding of later maps)
fcfg.join = cini.Statistics.JoinMaps;
fcfg.joinmd2 = cini.Statistics.JoinMapsMaxDist;
fcfg.joinulay = 5;

% split into local maxima
fcfg.localmax = cini.Statistics.LocalMax;
fcfg.localmaxsrf = cini.Statistics.LocalMaxSrfNeigh;
fcfg.localmaxsz = cini.Statistics.LocalMaxSizes;

% keyboard modifiers pressed at present
fcfg.mods = {};

% mouse position
fcfg.mpos = struct('cur', [0, 0], 'ddat', {{}}, 'down', [-1, -1], 'last', [0, 0], 'mods', {{}});

% neurosynth terms
fcfg.nsynth = struct( ...
    'termm', {findfiles([neuroelf_path('nsynth') '/terms'], '*.nii.gz', 'depth=1', 'relative=')}, ...
    'terms', {deblank(splittocellc(asciiread([neuroelf_path('nsynth') '/terms.txt']), ','))});

% no update flag
fcfg.noupdate = true;

% voxel-space orientation
if isempty(cini.MainFig.Orientation) || lower(cini.MainFig.Orientation(1)) ~= 'n'
    fcfg.orient = 'r';
else
    fcfg.orient = 'n';
end

% currently displayed page (xfigure property of figure object)
fcfg.page = 1;

% paint color code, mode (and settings; VMRs only)
fcfg.paint = struct('bbox', [-128, -128, -128; 128, 128, 128], 'code', 240, ...
    'mode', 1, 'over', [0, 32767], 'rad', 0, 'shap2', [0, 0], 'shap2w', 1, ...
    'shap3', [0, 0, 0], 'shap3w', 1, 'shape', 's', 'smooth', 0, 'smootk', ones(1001, 1));

% PLP handle (current PLP object, reference only)
fcfg.plp = [];

% p-values range factor (so, for p<0.05 what's the upper threshold?)
fcfg.prange = 0.0002;

% progress counter for tasks
fcfg.progress = struct;

% renderer
fcfg.renderer = cini.MainFig.Renderer;

% sampling frames (normal/zoom)
fcfg.sframe = [128, 128, 128; -127.9999, -127.9999, -127.9999];
fcfg.sframez = [96, 80, 104; -95.9999, -111.9999, -87.9999];

% show V16 content
fcfg.showv16 = false;

% position of 3-slice images
fcfg.slicepos = [tags.IM_NeuroElf_Slice_SAG.Position; ...
    tags.IM_NeuroElf_Slice_COR.Position; tags.IM_NeuroElf_Slice_TRA.Position];
fcfg.slicepos(:, 3:4) = fcfg.slicepos(:, 1:2) + fcfg.slicepos(:, 3:4);

% surface "position"
fcfg.spos = {[], 0, [0, 0, 0]};

% surface configuration
fcfg.srfcfg = struct('anglex', 180, 'angley', 0, 'time', 0, 'trans', [0, 0, 0], 'zoom', 1);
fcfg.srfcfgload = struct('anglex', 0, 'angley', 0, 'time', 0, 'trans', [0, 0, 0], 'zoom', 1);

% sampling stepsize (default: 1mm)
fcfg.sstep = 1;
fcfg.sstepz = 0.75;

% statistics-on-anatomical alpha-reduction factor
fcfg.stalphared = 2;

% surface-time-course plot visible
fcfg.stcplot = false;
fcfg.stcplotdata = [];

% position of surface-time-course display and variable handle
fcfg.stcpos = tags.AX_NeuroElf_STC_Plot.Position;
fcfg.stcpos(3:4) = fcfg.stcpos(1:2) + fcfg.stcpos(3:4);
fcfg.stcvar = [];

% data extracts from stats
fcfg.stext = struct('cons', {cell(0, 2)}, 'data', [], 'vois', {{}});

% transformation matrix (can be used to rotate display)
fcfg.strans = eye(4);
fcfg.strrot = [0, 0, 0];
fcfg.strscl = [1, 1, 1];
fcfg.strtra = [0, 0, 0];
fcfg.strzoom = false;

% position of surface axes
fcfg.surfpos = tags.AX_NeuroElf_Surface.Position;
fcfg.surfpos(:, 3:4) = fcfg.surfpos(:, 1:2) + fcfg.surfpos(:, 3:4);

% sampling zoom (MNI brain only)
fcfg.szoom = false;

% time-course plot visible
fcfg.tcplot = false;
fcfg.tcplotdata = [];
fcfg.tcplotylim = [-Inf, Inf];

% position of time-course display
fcfg.tcpos = tags.AX_NeuroElf_TC_Plot.Position;
fcfg.tcpos(3:4) = fcfg.tcpos(1:2) + fcfg.tcpos(3:4);

% default transimg object size
fcfg.tioosz = 256;

% text position objects
fcfg.txtpos = [tags.ED_NeuroElf_TALX.MLHandle, tags.ED_NeuroElf_TALY.MLHandle, tags.ED_NeuroElf_TALZ.MLHandle, ...
    tags.ED_NeuroElf_BVSX.MLHandle, tags.ED_NeuroElf_BVSY.MLHandle, tags.ED_NeuroElf_BVSZ.MLHandle];

% update GLM beta plots?
fcfg.updglmbp = true;

% which zoom view (0: 3 slices, 1...3: sag, cor, tra)
fcfg.zoom = 0;

% position of zoomed slices (times 3, for three "objects")
fcfg.zslicepos = tags.IM_NeuroElf_Slice_Zoom.Position;
fcfg.zslicepos(3:4) = fcfg.zslicepos(1:2) + fcfg.zslicepos(3:4);
fcfg.zslicepos = fcfg.zslicepos([1, 1, 1], :);

% initialiaze SliceVar and StatsVar to empty
fcfg.SliceUnder = struct('Filetype', 'NONE', 'RunTimeVars', struct('Trf', eye(4)));
fcfg.SliceVar = struct('Filetype', 'NONE', 'RunTimeVars', struct('Trf', eye(4)));
fcfg.StatsVar = struct('Filetype', 'NONE', 'RunTimeVars', struct('Trf', eye(4)));
fcfg.StatsVarIdx = [];

% display threshold, parameters and alpha
fcfg.StatsVarPar = {'t', 1, 1};

% stats var references object
fcfg.StatsVarRefObj = [];

% and repeat for surface files
fcfg.SurfBackColor = [0, 0, 0];
fcfg.SurfBarSize = [256, 64];
fcfg.SurfVar = struct('Filetype', 'NONE');
fcfg.SurfStatsVar = struct('Filetype', 'NONE');
fcfg.SurfStatsVarIdx = [];
fcfg.SurfStatsVarPar = {'t', 1, 1};
fcfg.SurfStatsVarRefObj = [];

% try to load standard LUT
try
    ne_gcfg.lut = xff([neuroelf_path('lut') '/Standard_extended.olt']);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ne_gcfg.lut = xff('new:olt');
end

% create layered image objects
tioosz = 256;
try
    if cini.MainFig.HiResSlicing
        tiofig = [];
        tiofig = figure;
        set(tiofig, 'Position', [20, 20, 100, 100]);
        set(tiofig, 'Visible', 'off');
        tiofrm = getframe(tiofig);
        delete(tiofig);
        if all(round([size(tiofrm.cdata, 1), size(tiofrm.cdata, 2)] ./ 100) == 2)
            tioosz = 2 * tioosz;
            fcfg.tioosz = tioosz;
        end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    delete(tiofig);
end
ne_gcfg.tio = struct( ...
    'imSag', transimg(tioosz, tioosz), ...
    'imCor', transimg(tioosz, tioosz), ...
    'imTra', transimg(tioosz, tioosz), ...
    'imSlZ', transimg(2 * tioosz, 2 * tioosz), ...
    'imRnd', transimg(512, 512), ...
    'satSag', [], 'satCor', [], 'satTra', [], 'satSlZ', [], 'satRnd', []);

% and set handles of images into transimg object (to allow use of display)
sethandle(ne_gcfg.tio.imSag, get(tags.IM_NeuroElf_Slice_SAG.MLHandle, 'Children'));
sethandle(ne_gcfg.tio.imCor, get(tags.IM_NeuroElf_Slice_COR.MLHandle, 'Children'));
sethandle(ne_gcfg.tio.imTra, get(tags.IM_NeuroElf_Slice_TRA.MLHandle, 'Children'));
sethandle(ne_gcfg.tio.imSlZ, get(tags.IM_NeuroElf_Slice_Zoom.MLHandle, 'Children'));
sethandle(ne_gcfg.tio.imRnd, get(tags.IM_NeuroElf_Slice_Rend.MLHandle, 'Children'));

% create default (new) POI/VOI objects
ne_gcfg.poi = xff('new:poi');
ne_gcfg.voi = xff('new:voi');

% put important ones into internal structure with short names
ch = ne_gcfg.h;

% main figure (xfigure and MLHandle)
ch.MainFig = MainFig;
ch.MainFigMLH = MainFig.MLHandle;
ch.MainFigTags = tags;

% direct children (contrast manager and vismontage UIs, etc.)
ch.AC = [];
ch.AP = [];
ch.CM = [];
ch.Console = [];
ch.MDM = [];
ch.MKDA = [];
ch.RM = [];
ch.Render = [];
ch.SurfMontage = [];
ch.VisMontage = [];

% initialize Children
ch.Children = struct;

% specific DTI menu entries
ch.DTIPlotFibers = tags.UIM_NeuroElf_DTIPlotFib;
ch.DTITrackFibers = tags.UIM_NeuroElf_HDRDTIFibTrack;

% echo calls menu item
ch.EchoCalls = tags.UIM_NeuroElf_EchoCalls;

% linked browsing menu item
ch.LinkedBrowse = tags.UIM_NeuroElf_LinkedBrowse;
ch.LinkedBrowseBT = tags.BT_NeuroElf_togglelink;

% recent files menu(s)
rfnum = cini.RecentFiles.Number;
ch.RecentFiles = struct( ...
    'slc', {cell(rfnum, 1)}, 'stat',  {cell(rfnum, 1)}, ...
    'srf', {cell(rfnum, 1)}, 'srfst', {cell(rfnum, 1)});
rffs = {'slc', 'stat', 'srf', 'srfst'};
rfms = {tags.UIM_NeuroElf_recentslice.MLHandle, tags.UIM_NeuroElf_recentstats.MLHandle, ...
    tags.UIM_NeuroElf_recentsrf.MLHandle, tags.UIM_NeuroElf_recentsrfst.MLHandle};
for rfc = 1:rfnum
    for rff = 1:4
        ch.RecentFiles.(rffs{rff}){rfc} = uimenu( ...
            'Enable',   'on', 'Callback', {@ne_openfile, ''}, ...
            'Label',    sprintf('%s_%04d', rffs{rff}, rfc), ...
            'Parent',   rfms{rff}, 'Visible',  'off');
    end
end

% listener toggle
ch.Listener = tags.BT_NeuroElf_RListener;

% neurosynth menu
ch.NeuroSynth = tags.UIM_NeuroElf_neurosynth;

% SVC menu entries
ch.SVCEntries = [tags.UIM_NeuroElf_VMPSVCVOI.MLHandle, tags.UIM_NeuroElf_VMPSVCMask.MLHandle, ...
    tags.UIM_NeuroElf_VMPSVCVMR.MLHandle, tags.UIM_NeuroElf_VMPSVCColin.MLHandle];

% slicing and statistics variable selection (dropdown)
ch.SliceVar = tags.DD_NeuroElf_varlist;
ch.SliceVar.UserData = cell(0, 4);
ch.StatsVar = tags.DD_NeuroElf_statlist;
ch.StatsVar.UserData = cell(0, 4);
ch.SurfVar = tags.DD_NeuroElf_varlistsrf;
ch.SurfVar.UserData = cell(0, 4);
ch.SurfStatsVar = tags.DD_NeuroElf_statlistsrf;
ch.SurfStatsVar.UserData = cell(0, 4);
ch.Scenery = tags.LB_NeuroElf_Scenery;
ch.Scenery.UserData = cell(0, 4);
ch.Scenery.Value = [];
ch.SceneryProps = tags.BT_NeuroElf_SceneProps;
ch.SceneryViewPoint = tags.ED_NeuroElf_SrfViewPnt;

% available maps from StatsVar
ch.StatsVarMaps = tags.LB_NeuroElf_statmaps;
ch.StatsVarMaps.Value = [];
ch.StatsVarRefs = tags.DD_NeuroElf_statsref;
ch.StatsVarRefs.UserData = {[]};
ch.StatsVarRefRegs = tags.BT_NeuroElf_statsrefreg;
ch.StatsVarRefNuis = tags.BT_NeuroElf_statsrefixx;
ch.SurfStatsVarMaps = tags.LB_NeuroElf_statmapssrf;
ch.SurfStatsVarMaps.Value = [];
ch.SurfStatsVarRefs = tags.DD_NeuroElf_statsrefsrf;
ch.SurfStatsVarRefs.UserData = {[]};

% projection button
ch.StatsVarProject = tags.BT_NeuroElf_statmproj;

% TC Movie menu
ch.TCMovieMenu = tags.UIM_NeuroElf_VTCCondAvgMovie;

% ToolTip
ch.ToolTip = tags.UIM_NeuroElf_TTIP;

% show V16 button
ch.VMRShowV16 = tags.BT_NeuroElf_showv16;

% list of clusters
ch.Clusters = tags.LB_NeuroElf_clusters;
ch.Clusters.Value = [];
ch.Clusters.String = {};
ch.ClustersSrf = tags.LB_NeuroElf_Srfclust;
ch.ClustersSrf.Value = [];
ch.ClustersSrf.String = {};

% cluster output table
ch.ClusterTable = tags.ED_NeuroElf_clusters;

% cluster zoom toggle button
ch.ClusterZoom = tags.BT_NeuroElf_clustzoom;

% progress bar
ch.Progress = tags.PB_NeuroElf_mainpbar;
ch.Progress.Visible = 'off';
ch.Progress.HandleVisibility = 'off';
ch.Progress.UserData = {@neuroelf_gui, 'progress', []};

% surface axes
ch.Surface = tags.AX_NeuroElf_Surface.MLHandle;
srf = ch.Surface;
srfcfg = cini.Surface;
srfbcl = (1 / 255) .* srfcfg.BackgroundColor(:)';
fcfg.SurfBackColor = srfbcl;
set(srf, 'Color', srfbcl);
set(srf, 'View', [90, 0]);
slim = [-128, 128];
set(srf, 'XLim', 4 * slim, 'YLim', slim, 'ZLim', slim);
set(srf, 'XTick', [], 'YTick', [], 'ZTick', []);
for lc = 1:numel(srfcfg.Lights)
    light('Parent', srf, 'Position', srfcfg.Lights{lc}, 'Color', (1 / 255) .* srfcfg.LightColors{lc});
end
set(srf, 'XColor', [0, 0, 0], 'YColor', [0, 0, 0], 'ZColor', [0, 0, 0]);
ssbp = cini.Statistics.ThreshBarPos;
meshsize = floor(512 .* (ssbp(1, [4, 3]) - ssbp(1, [2, 1])));
fcfg.SurfBarSize = meshsize;
[ssbv, ssbf] = mesh3d(meshsize, struct('orient', 4, ...
    'xline', [0.5, round(256 * (ssbp(2) - 0.5))], ...
    'yline', [0.5, round(256 * (ssbp(1) - 0.5))], 'zvalue', -256));
ch.SurfaceStatsBar = patch(ssbv(:, 1), ssbv(:, 2), ssbv(:, 3), zeros(size(ssbv, 1), 1), ...
    'FaceColor', 'none', 'EdgeColor', 'none', 'Parent', ch.Surface, 'Visible', 'off');
set(ch.SurfaceStatsBar, 'Faces', ssbf, 'FaceVertexCData', repmat(srfbcl, size(ssbf, 1), 1), ...
    'FaceColor', 'flat', 'Visible', 'off');
stextxpos = 0.015 * slim(1) + 0.985 * slim(2);
stextypos1 = 0.78 * slim(1) + 0.22 * slim(2);
stextypos2 = 0.22 * slim(1) + 0.78 * slim(2);
ch.SurfaceStatsText = [text(0, stextxpos, stextypos1, ' ', 'Parent', ch.Surface, ...
    'Color', [1, 1, 1], 'FontSize', 12, 'HorizontalAlignment', 'right'), ...
    text(0, stextxpos, stextypos2, ' ', 'Parent', ch.Surface, 'Color', [1, 1, 1], ...
    'FontSize', 12, 'HorizontalAlignment', 'right')];

% create SurfaceHGTransform to hold transformable patches
ch.SurfaceTransform = hgtransform('Parent', ch.Surface);
set(ch.SurfaceTransform, 'Matrix', btc_meshtrf(fcfg.srfcfg), 'Visible', 'on');

% add crosshair lines to images axes objects -> SAG
chax = tags.AX_NeuroElf_Slice_SAG.MLHandle;
set(chax, 'Units', 'pixels');
ch.SagLineX = line([0; 0.999], [0.5; 0.5], 'Color', fcfg.chcol, 'Parent', chax);
ch.SagLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', fcfg.chcol, 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');

% -> COR
chax = tags.AX_NeuroElf_Slice_COR.MLHandle;
set(chax, 'Units', 'pixels');
ch.CorLineX = line([0; 0.999], [0.5; 0.5], 'Color', fcfg.chcol, 'Parent', chax);
ch.CorLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', fcfg.chcol, 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');
ch.CorStatsText = [text(0.99, 0.22, ' ', 'Parent', chax, 'Color', [1, 1, 1], ...
    'FontSize', 10, 'HorizontalAlignment', 'right'), ...
    text(0.99, 0.77, ' ', 'Parent', chax, 'Color', [1, 1, 1], ...
    'FontSize', 10, 'HorizontalAlignment', 'right')];

% -> TRA
chax = tags.AX_NeuroElf_Slice_TRA.MLHandle;
set(chax, 'Units', 'pixels');
ch.TraLineX = line([0; 0.999], [0.5; 0.5], 'Color', fcfg.chcol, 'Parent', chax);
ch.TraLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', fcfg.chcol, 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');

% -> Zoom
chax = tags.AX_NeuroElf_Slice_Zoom.MLHandle;
set(chax, 'Units', 'pixels');
ch.ZoomLineX = line([0; 0.999], [0.5; 0.5], 'Color', fcfg.chcol, 'Parent', chax);
ch.ZoomLineY = line([0.5; 0.5], [0.001; 0.999], 'Color', fcfg.chcol, 'Parent', chax);
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');
ch.ZoomStatsText = [text(0.99, 0.22, ' ', 'Parent', chax, 'Color', [1, 1, 1], ...
    'FontSize', 12, 'HorizontalAlignment', 'right'), ...
    text(0.99, 0.78, ' ', 'Parent', chax, 'Color', [1, 1, 1], ...
    'FontSize', 12, 'HorizontalAlignment', 'right')];

% -> Render
chax = tags.AX_NeuroElf_Slice_Rend.MLHandle;
set(chax, 'Units', 'pixels');
set(chax, 'Units', 'pixels', 'XTick', [], 'YTick', [], 'Visible', 'off');

% -> Histogram
chax = tags.AX_NeuroElf_Slice_Hist.MLHandle;
set(chax, 'Units', 'pixels');
ch.HistImage = tags.IM_NeuroElf_Slice_Hist.Children;
ch.HistLine1 = line([0; 0.5], [0.002; 0.002], 'Color', [0.5, 0.5, 0.5], ...
    'LineWidth', 3, 'Parent', chax);
ch.HistLine2 = line([0.5; 1], [0.998; 0.998], 'Color', [0.5, 0.5, 0.5], ...
    'LineWidth', 3, 'Parent', chax);
ch.HistPlot = line(0.1 * ones(256, 1), (1/512:1/256:511/512)', ...
    'Color', [0.25, 1, 0.25], 'LineWidth', 2, 'Parent', chax);
set(ch.HistImage, 'CData', repmat(uint8(0:255)', [1, 16, 3]));
set(chax, 'Units', 'pixels', 'YDir', 'normal', 'XTick', [], 'YTick', [], 'Visible', 'off');
tags.IM_NeuroElf_Slice_Hist.YDir = 'normal';

% current position
ch.Coord.TEdX = mlhandle(tags.ED_NeuroElf_TALX);
ch.Coord.TEdY = mlhandle(tags.ED_NeuroElf_TALY);
ch.Coord.TEdZ = mlhandle(tags.ED_NeuroElf_TALZ);
ch.Coord.VEdX = mlhandle(tags.ED_NeuroElf_BVSX);
ch.Coord.VEdY = mlhandle(tags.ED_NeuroElf_BVSY);
ch.Coord.VEdZ = mlhandle(tags.ED_NeuroElf_BVSZ);

% time-dim controls
ch.Coord.Temp = tags.ED_NeuroElf_TempPos;
ch.Coord.TempSlider = tags.SL_NeuroElf_TempPos;
ch.Coord.TempSlider.Max = 120;
ch.Coord.TempSlider.Value = 1;
ch.Coord.TempSlider.Min = 1;

% undo paint toggle button
ch.DrawUndo = tags.BT_NeuroElf_drawu;

% menu item fast access
ch.Menu.CloseFile = tags.UIM_NeuroElf_closefile;
ch.Menu.CloseFile.Visible = 'off';
ch.Menu.LimitVMR = tags.UIM_NeuroElf_VMRLimitVMR;
ch.Menu.SelectFile = tags.UIM_NeuroElf_selectfile;
ch.Menu.SelectFile.Visible = 'off';
ch.Menu.Stats = tags.UIM_NeuroElf_STATS;
ch.Menu.Stats.Visible = 'off';
ch.Menu.VOI = tags.UIM_NeuroElf_VOI;
ch.Menu.VOI.Visible = 'off';

% add some text properties
ch.Stats.LThresh = tags.ED_NeuroElf_LowerThresh;
ch.Stats.UThresh = tags.ED_NeuroElf_UpperThresh;
ch.Stats.PosTail = tags.CB_NeuroElf_PositivStat;
ch.Stats.NegTail = tags.CB_NeuroElf_NegativStat;
ch.Stats.PThresh = tags.DD_NeuroElf_StatSetP;
ch.Stats.PThresh.String = lsqueeze(splittocellc(sprintf('%g ', ...
    cini.Statistics.Thresholds), ' ', true, true));
ch.Stats.kThresh = tags.ED_NeuroElf_kExtThresh.MLHandle;
ch.Stats.UsekThr = tags.CB_NeuroElf_kExtThresh;
ch.Stats.ICBM2TAL = tags.CB_NeuroElf_ICBM2TAL;
ch.Stats.TDClient = tags.CB_NeuroElf_TDClient;
ch.Stats.UseLUT = tags.RB_NeuroElf_LUTColor;
ch.Stats.UseRGB = tags.RB_NeuroElf_RGBColor;
ch.Stats.RGBLPos = tags.BT_NeuroElf_RGBLowPos;
ch.Stats.RGBUPos = tags.BT_NeuroElf_RGBUppPos;
ch.Stats.RGBLNeg = tags.BT_NeuroElf_RGBLowNeg;
ch.Stats.RGBUNeg = tags.BT_NeuroElf_RGBUppNeg;
ch.SrfStats.LThresh = tags.ED_NeuroElf_SrfLowerThr;
ch.SrfStats.UThresh = tags.ED_NeuroElf_SrfUpperThr;
ch.SrfStats.PosTail = tags.CB_NeuroElf_SrfPosStat;
ch.SrfStats.NegTail = tags.CB_NeuroElf_SrfNegStat;
ch.SrfStats.PThresh = tags.DD_NeuroElf_SrfStatSetP;
ch.SrfStats.PThresh.String = {'0.05'; '0.02'; '0.01'; '0.005'; '0.002'; ...
    '0.001'; '0.0005'; '0.0001'; '1e-5'; '1e-6'};
ch.SrfStats.PThreshTxt = tags.TX_NeuroElf_SrfSetP;
ch.SrfStats.ClusterTable = tags.BT_NeuroElf_SrfClusterT;
ch.SrfStats.kThresh = tags.ED_NeuroElf_SrfkExtThr.MLHandle;
ch.SrfStats.UsekThr = tags.CB_NeuroElf_SrfkExtThr;
ch.SrfStats.UseLUT = tags.RB_NeuroElf_SrfLUTColor;
ch.SrfStats.UseRGB = tags.RB_NeuroElf_SrfRGBColor;
ch.SrfStats.RGBLPos = tags.BT_NeuroElf_SrfRGBLPos;
ch.SrfStats.RGBUPos = tags.BT_NeuroElf_SrfRGBUPos;
ch.SrfStats.RGBLNeg = tags.BT_NeuroElf_SrfRGBLNeg;
ch.SrfStats.RGBUNeg = tags.BT_NeuroElf_SrfRGBUNeg;

% split local max
ch.SplitLocalMax = tags.CB_NeuroElf_SplitMaxima;
ch.SplitLocalMax.Value = double(fcfg.localmax);

% interpolation checkbox
ch.Interpolate = tags.CB_NeuroElf_Interpolate;

% label with currently sampled values
ch.SampledValues = tags.TX_NeuroElf_SValues;

% time-course splot (plot a flat line as default)
ch.STCPlot = tags.AX_NeuroElf_STC_Plot;
ch.STCPlotChild = plot(ch.STCPlot.MLHandle, (1:120)', zeros(120, 1), ...
    'LineWidth', 2, 'Color', [0.8, 0.8, 0.8]);
set(ch.STCPlot.MLHandle, 'HandleVisibility', 'off');
ch.STCPlotChildren = [];
ch.TCPlot = tags.AX_NeuroElf_TC_Plot;
ch.TCPlotChild = plot(ch.TCPlot.MLHandle, (1:120)', zeros(120, 1));
set(ch.TCPlot.MLHandle, 'HandleVisibility', 'off');
ch.TCPlotChildren = [];
hold(ch.TCPlot.MLHandle, 'on');
ch.TCPlotDiscards = image([1, 120], [0, 0], ...
    uint8(repmat(reshape([255, 0, 0], [1, 1, 3]), 1, 120)), 'Parent', ch.TCPlot.MLHandle);
set(ch.TCPlotDiscards, 'AlphaData', zeros(1, 120));
ch.TCPlotUndock = tags.BT_NeuroElf_tcundock;

% minimize/maximize buttons
ch.Maximize = tags.BT_NeuroElf_maximize;
ch.Minimize = tags.BT_NeuroElf_minimize;

% for children ob objects (handles will be deleted upon exit)
ch.fChild = [];

% add a child to SelectFile list
uimenu(mlhandle(tags.UIM_NeuroElf_selectfile), 'Label', 'none', 'Enable', 'off');

% populate Colin-27 list
try
    mmh = mlhandle(tags.UIM_NeuroElf_opencolin);
    sep = 'off';
    c27p = neuroelf_path('colin');
    shnp = neuroelf_path('shen');
    c27l = splittocell(asciiread( ...
        [neuroelf_path('config') '/colin27fileorder.txt']), char([10, 13]), 1, 1);
    c27sl = grep(c27l, '-x', '\.vmr$');
    if ~isempty(c27sl)
        mmhc = uimenu(mmh, 'Label', 'VMRs');
        for lc = 1:numel(c27sl)
            if strcmpi(c27sl{lc}, 'separator.vmr')
                sep = 'on';
                continue;
            elseif exist([c27p '/' c27sl{lc}], 'file') == 2
                uimenu(mmhc, 'Label', ['    ' strrep(c27sl{lc}, '.vmr', '')], ...
                    'Callback', {@ne_openfile, [c27p '/' c27sl{lc}]}, ...
                    'Separator', sep);
            elseif isempty(regexpi(c27sl{lc}, '^colin'))
                uimenu(mmhc, 'Label', strrep(c27sl{lc}, '.vmr', ''), ...
                    'Enable', 'off');
            end
            sep = 'off';
        end
        sep = 'off';
    end
    c27sl = grep(c27l, '-x', '\.srf$');
    if ~isempty(c27sl)
        mmhc = uimenu(mmh, 'Label', 'Surfaces', 'Separator', 'on');
        for lc = 1:numel(c27sl)
            if strcmpi(c27sl{lc}, 'separator.srf')
                sep = 'on';
                continue;
            elseif exist([c27p '/' c27sl{lc}], 'file') == 2
                uimenu(mmhc, 'Label', ['    ' strrep(c27sl{lc}, '.srf', '')], ...
                    'Callback', {@ne_openfile, [c27p '/' c27sl{lc}]}, 'Separator', sep);
            elseif isempty(regexpi(c27sl{lc}, '^colin'))
                uimenu(mmhc, 'Label', strrep(c27sl{lc}, '.srf', ''), 'Enable', 'off');
            end
            sep = 'off';
        end
    end
    c27sl = grep(c27l, '-x', '\.special$');
    if ~isempty(c27sl)
        mmhc = uimenu(mmh, 'Label', 'Specialized surfaces');
        for lc = 1:numel(c27sl)
            if exist([c27p '/' strrep(c27sl{lc}, '.special', '.srf')], 'file') == 2
                uimenu(mmhc, 'Label', ['    ' strrep(c27sl{lc}, '.special', '')], ...
                    'Callback', {@ne_openfile, [c27p '/' strrep(c27sl{lc}, '.special', '.srf')]}, ...
                    'Separator', 'off');
            elseif exist([shnp '/' strrep(c27sl{lc}, '.special', '.srf')], 'file') == 2
                uimenu(mmhc, 'Label', ['    ' strrep(c27sl{lc}, '.special', '')], ...
                    'Callback', {@ne_openfile, [shnp '/' strrep(c27sl{lc}, '.special', '.srf')]}, ...
                    'Separator', 'off');
            elseif isempty(regexpi(c27sl{lc}, '^colin'))
                uimenu(mmhc, 'Label', strrep(c27sl{lc}, '.special', ''), 'Enable', 'off');
            end
        end
    end
    c27sl = grep(c27l, '-x', '\.subcort$');
    if ~isempty(c27sl)
        mmhc = uimenu(mmh, 'Label', 'Sub-cortical surfaces');
        for lc = 1:numel(c27sl)
            if exist([c27p '/' strrep(c27sl{lc}, '.subcort', '.srf')], 'file') == 2
                uimenu(mmhc, 'Label', ['    ' strrep(c27sl{lc}, '.subcort', '')], ...
                    'Callback', {@ne_openfile, [c27p '/' strrep(c27sl{lc}, '.subcort', '.srf')]}, ...
                    'Separator', 'off');
            end
        end
        uimenu(mmhc, 'Label', 'Load all subcortical surfaces', 'Callback', ...
            {@ne_srf_tools, 'loadsubcort', 0}, 'Separator', 'on');
        if ~isempty(findfiles(c27p, 'colin_subcort_*ICBMnorm.srf'))
        uimenu(mmhc, 'Label', 'Load all subcortical ICBMnorm surfaces', 'Callback', ...
            {@ne_srf_tools, 'loadsubcort', 'ICBMnorm'}, 'Separator', 'off');
        end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    tags.UIM_NeuroElf_opencolin.Visible = 'off';
end

% populate LUT list
try
    mmh = mlhandle(tags.UIM_NeuroElf_LUTList);
    luts = findfiles(neuroelf_path('lut'), '*.olt', 'depth=1');
    for lc = 1:numel(luts)
        [lutp, lutf] = fileparts(luts{lc});
        lutf = strrep(lutf, '_', ' ');
        uimenu(mmh, 'Label', [upper(lutf(1)) lutf(2:end)], ...
            'Callback', {@ne_openfile, luts{lc}});
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    tags.UIM_NeuroElf_LUTList.Visible = 'off';
end
if exist([neuroelf_path('files') '/shenparcel/shen_parcels.voi'], 'file') ~= 2
    tags.UIM_NeuroElf_LoadShenVOIs.Visible = 'off';
end

% coordinate selection (support cursor keys in text boxes)
tags.ED_NeuroElf_TALX.KeyPressFcn = {@ne_keypress_uic, 0, @ne_setwindowpos};
tags.ED_NeuroElf_TALY.KeyPressFcn = {@ne_keypress_uic, 0, @ne_setwindowpos};
tags.ED_NeuroElf_TALZ.KeyPressFcn = {@ne_keypress_uic, 0, @ne_setwindowpos};
tags.ED_NeuroElf_BVSX.KeyPressFcn = {@ne_keypress_uic, 0, @ne_setwindowpos};
tags.ED_NeuroElf_BVSY.KeyPressFcn = {@ne_keypress_uic, 0, @ne_setwindowpos};
tags.ED_NeuroElf_BVSZ.KeyPressFcn = {@ne_keypress_uic, 0, @ne_setwindowpos};

% hide tooltip (menu)
tags.UIM_NeuroElf_TTIP.Visible = 'off';

% main UI Fcn callbacks
MainFig.CloseRequestFcn = @ne_closewindow;
MainFig.KeyPressFcn = @ne_keypress;
MainFig.KeyReleaseFcn = @ne_keyrelease;
MainFig.WindowButtonDownFcn = @ne_btdown;
MainFig.WindowButtonMotionFcn = @ne_btmove;
MainFig.WindowButtonUpFcn = @ne_btup;

% disable OpenGL hardware support on Windows?
if strcmpi(fcfg.renderer, 'opengl') && ispc && cini.Surface.OpenGLHWAccelOnWindows ~= 1
    try
        opengl('software', true);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        warning('neuroelf:GUI:openGLError', ...
            'Error setting OpenGL to SoftwareAcceleration for PCs.');
        fcfg.renderer = 'zbuffer';
    end
end

% set group states
MainFig.SetGroupEnabled('SLoaded', 'off');
MainFig.SetGroupEnabled('SLdNVMP', 'off');
MainFig.SetGroupVisible('FMRMenu', 'off');
MainFig.SetGroupVisible('GLMMenu', 'off');
MainFig.SetGroupVisible('VMPMenu', 'off');
MainFig.SetGroupVisible('VMRMenu', 'off');
MainFig.SetGroupVisible('VTCMenu', 'off');
MainFig.SetGroupVisible('SRFMenu', 'off');
MainFig.SetGroupVisible('SMPMenu', 'off');
MainFig.SetGroupVisible('HDRMenu', 'off');
MainFig.SetGroupEnabled('SLdNSMP', 'off');

% set back to global structure
ne_gcfg.fcfg = fcfg;
ne_gcfg.h = ch;

% pre-set button colors
lutc = (1 / 255) * ne_gcfg.lut.Colors;
lutn = size(lutc, 1);
ch.Stats.RGBLPos.BackgroundColor = lutc(1, :);
ch.Stats.RGBUPos.BackgroundColor = lutc(0.5 * lutn, :);
ch.Stats.RGBLNeg.BackgroundColor = lutc(0.5 * lutn + 1, :);
ch.Stats.RGBUNeg.BackgroundColor = lutc(lutn, :);
ch.SrfStats.RGBLPos.BackgroundColor = lutc(1, :);
ch.SrfStats.RGBUPos.BackgroundColor = lutc(0.5 * lutn, :);
ch.SrfStats.RGBLNeg.BackgroundColor = lutc(0.5 * lutn + 1, :);
ch.SrfStats.RGBUNeg.BackgroundColor = lutc(lutn, :);

% set options
ne_setoption(0, 0, 'ctableadd', cini.Statistics.ClusterTableAdd);
ne_setoption(0, 0, 'ctablelupcrd', cini.Statistics.LookupCoord);
ne_setoption(0, 0, 'ctablescsizes', fcfg.localmaxsz);
ne_setoption(0, 0, 'ctablesort', fcfg.clsort);
ne_setoption(0, 0, 'drawon', 'Mouse', cini.Drawing.OnMouse);
ne_setoption(0, 0, 'drawon', 'Cursor', cini.Drawing.OnCursor);
ne_setoption(0, 0, 'drawon', 'Position', cini.Drawing.OnPosition);
ne_setoption(0, 0, 'drawon', 'Linked', cini.Drawing.OnLinked);
ne_setoption(0, 0, 'echocalls', cini.MainFig.EchoCalls);
ne_setoption(0, 0, 'extmapnames', cini.Statistics.ExtendedMapNames);
ne_setoption(0, 0, 'extonselect', cini.Statistics.ExtractOnSelect);
ne_setoption(0, 0, 'extsepchars', cini.Statistics.ExtractSepChars);
ne_setoption(0, 0, 'exttransio', cini.Statistics.ExtractTransIO);
ne_setoption(0, 0, 'extwithsids', cini.Statistics.ExtractWithSubIDs);
ne_setoption(0, 0, 'instseedcorr', cini.Statistics.InstantSeedCorr);
ne_setoption(0, 0, 'joinmd2', cini.Statistics.JoinMapsMaxDist);
ne_setoption(0, 0, 'mkda', 'lookuponcursor', cini.MKDA.LookupOnCursor);
ne_setoption(0, 0, 'mkda', 'lookuponcluster', cini.MKDA.LookupOnCluster);
ne_setoption(0, 0, 'newvmrres', cini.MainFig.NewVMRRes);
ne_setoption(0, 0, 'orientation', fcfg.orient, [], false);
ne_setoption(0, 0, 'showthreshbars', cini.Statistics.ShowThreshBars, false);
ne_setoption(0, 0, 'showthreshtext', cini.Statistics.ShowThreshText, false);
ne_setoption(0, 0, 'surfrecoonesurf', cini.Surface.RecoOneSurfaceOnly);
ne_setoption(0, 0, 'surfreco4tps', cini.Surface.RecoTriPerVoxelFace == 4);
ne_setoption(0, 0, 'surfreusexsm', cini.Surface.RecoOneSurfaceOnly);
ne_setoption(0, 0, 'vmpusefdr', cini.Statistics.VMPUseFDR);

% load objects from xff class (includes setcvar, setcstats, etc.)
fprintf(' loading objects...');
pause(0.001);
ne_loadobjects(1);

% fill neurosynth menu
ne_neurosynth(0, 0, 'menu');

% show the first page (3-slices view, includes call to setslicepos)
ne_gcfg.fcfg.noupdate = false;
ne_showpage(0, 0, 1);

% size adaptation
ss0 = get(0, 'ScreenSize');
lastsize = cini.MainFig.Size;
if ~cini.MainFig.FullSize || any(ss0(3:4) < [1024, 720])
    ne_swapfullsize(0, 0, 'swap');
elseif all(lastsize >= fcfg.fullsize)
    ne_swapfullsize(0, 0, lastsize);
end

% set to last known position
try
    lastpos = cini.MainFig.Position;
    if any(lastpos ~= -1)
        ch.MainFig.Position(1:2) = lastpos;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% force buttons to correct colors
pause(0.01);
drawnow;
try
    getframe(ch.MainFigMLH);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% ensure buttons look OK for older versions
if c.mlversion < 900
    bts = fieldnames(tags);
    bts(cellfun('isempty', regexp(bts, '^BT_'))) = [];
    for cc = 1:numel(bts)
        if strcmpi(tags.(bts{cc}).Type, 'axes')
            bpos = tags.(bts{cc}).Position;
            set(tags.(bts{cc}).MLHandle, 'XLim', [0, bpos(3) + 1], 'YLim', [0, bpos(4) + 1]);
        end
    end
end

% get positions
[ne_gcfg.fcfg.chpos{1:3}] = ch.MainFig.ChildPositions;
ne_gcfg.fcfg.chpos{3}(:, 3:4) = ne_gcfg.fcfg.chpos{3}(:, 1:2) + ne_gcfg.fcfg.chpos{3}(:, 3:4);

% set correct renderer and make figure visible
ch.Progress.HandleVisibility = 'off';
ch.MainFig.HandleVisibility = 'callback';
set(ch.MainFigMLH, 'Resize', 'on', 'ResizeFcn', @ne_swapfullsize);
ch.MainFig.Visible = 'on';
ne_setoption(0, 0, 'renderer', fcfg.renderer);
figure(ch.MainFigMLH);

% force color buttons to display correctly
ch.Stats.RGBLPos.BackgroundColor = 1 - ch.Stats.RGBLPos.BackgroundColor;
ch.Stats.RGBUPos.BackgroundColor = 1 - ch.Stats.RGBUPos.BackgroundColor;
ch.Stats.RGBLNeg.BackgroundColor = 1 - ch.Stats.RGBLNeg.BackgroundColor;
ch.Stats.RGBUNeg.BackgroundColor = 1 - ch.Stats.RGBUNeg.BackgroundColor;
ch.SrfStats.RGBLPos.BackgroundColor = 1 - ch.SrfStats.RGBLPos.BackgroundColor;
ch.SrfStats.RGBUPos.BackgroundColor = 1 - ch.SrfStats.RGBUPos.BackgroundColor;
ch.SrfStats.RGBLNeg.BackgroundColor = 1 - ch.SrfStats.RGBLNeg.BackgroundColor;
ch.SrfStats.RGBUNeg.BackgroundColor = 1 - ch.SrfStats.RGBUNeg.BackgroundColor;
drawnow;
ch.Stats.RGBLPos.BackgroundColor = 1 - ch.Stats.RGBLPos.BackgroundColor;
ch.Stats.RGBUPos.BackgroundColor = 1 - ch.Stats.RGBUPos.BackgroundColor;
ch.Stats.RGBLNeg.BackgroundColor = 1 - ch.Stats.RGBLNeg.BackgroundColor;
ch.Stats.RGBUNeg.BackgroundColor = 1 - ch.Stats.RGBUNeg.BackgroundColor;
ch.SrfStats.RGBLPos.BackgroundColor = 1 - ch.SrfStats.RGBLPos.BackgroundColor;
ch.SrfStats.RGBUPos.BackgroundColor = 1 - ch.SrfStats.RGBUPos.BackgroundColor;
ch.SrfStats.RGBLNeg.BackgroundColor = 1 - ch.SrfStats.RGBLNeg.BackgroundColor;
ch.SrfStats.RGBUNeg.BackgroundColor = 1 - ch.SrfStats.RGBUNeg.BackgroundColor;

% start up remote
if cini.Remote.OnStartup
    fprintf(' remote...');
    ne_remote(0, 0, 'listen');
end

% done
ne_gcfg.c.init = true;
fprintf(' done.\n');

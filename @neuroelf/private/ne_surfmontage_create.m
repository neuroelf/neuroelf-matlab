function varargout = ne_surfmontage_create(varargin)
% ne_surfmontage_create  - create surface-based montage image
%
% FORMAT:       ne_surfmontage_create(SRC, EVT [, config])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       config      optional 1x1 struct with two fields
%        .elems     Elements-by-15 cell specification of stitched images
%                   - name of SRF (is looked up in the Scenery UI widget)
%                   - number of srf.Smooth operations (after sampling)
%                   - Smooth operation factor
%                   - number of srf.Inflate operations (after Smoothing)
%                   - Inflate operation factor
%                   - X-translation for screenshot
%                   - Y-translation for screenshot
%                   - azimuth angle for screenshot
%                   - zenith angle for screenshot
%                   - zoom factor for screenshot
%                   - time index for screenshot
%                   - resize X before stitching
%                   - resize Y before stitching
%                   - stitched X offset (relative 0...1)
%                   - stiched Y offset (relative 0...1)
%        .output    1x10 cell with general output options
%                   - image size (1x2 double, width and height)
%                   - RGB background color
%                   - copy stats bars flag
%                   - write output (set to false/0 or true/1)
%                   - output filename
%                   - sample VMP flag (set to false/0 or true/1)
%                   - restrict VMP flag (set to false/0 or true/1)
%                   - multiply threshold flag (set to false/0 or true/1)
%                   - threshold multiplication factor (set to [LowF, UppF])
%                   - save smoothed/inflated SRF (set to false/0 or true/1)
%
% No output fields.
%
% Example:
%
%     ne_surfmontage_create(0, 0, struct('elems', {{ ...
%         leftsrf.FilenameOnDisk('r'),  0, 0.5, 0, 0.5, -25, -22, 180, 0, 1.375, 0, 0.5, 0.5, 0, 0;
%         leftsrf.FilenameOnDisk('r'),  0, 0.5, 0, 0.5,  25, -22,   0, 0, 1.375, 0, 0.5, 0.5, 0, 0.5;
%         rightsrf.FilenameOnDisk('r'), 0, 0.5, 0, 0.5,  25, -22,   0, 0, 1.375, 0, 0.5, 0.5, 0.5, 0;
%         rightsrf.FilenameOnDisk('r'), 0, 0.5, 0, 0.5, -25, -22, 180, 0, 1.375, 0, 0.5, 0.5, 0.5, 0.5}}, ...
%         'output', {{[3074, 2304], [255, 255, 255], 0, 1, './surfmontage_out.png', 0, 0, 0, [1, 1], 0}}));

% Version:  v1.1
% Build:    18080413
% Date:     Aug-04 2016, 1:32 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015 - 2018, Jochen Weber
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
ch = ne_gcfg.h;
sm = ne_gcfg.h.SurfMontage;
cc = ne_gcfg.fcfg.SurfMontage;
cini = ne_gcfg.c.ini;

% output
if nargout > 0
    varargout = cell(1, nargout);
end

% only allow one concurrent call
if any(strcmp('surfmontagecreate', ne_gcfg.c.blockcb))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'surfmontagecreate';

% either options must be given all
if nargin > 2 && isstruct(varargin{3}) && numel(varargin{3}) == 1 && ...
    isfield(varargin{3}, 'elems') && iscell(varargin{3}.elems) && ...
    ndims(varargin{3}.elems) == 2 && size(varargin{3}.elems, 2) == 15 && ...
    isfield(varargin{3}, 'output') && iscell(varargin{3}.output) && ...
    numel(varargin{3}.output) == 10
    o = varargin{3};

% or taken from UI
else

    % initialize
    o = struct('elems', {{}}, 'output', {cell(1, 10)});

    % then fill in selected elements
    if ~isstruct(cc) || ~isfield(cc, 'cc') || ~isfield(cc, 'ec') || ~isfield(cc, 'mc') || ...
       ~isstruct(sm) || ~isfield(sm, 'h') || ~isstruct(sm.h)
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'surfmontagecreate')) = [];
        return;
    end
    mc = cc.mc;
    ec = cc.ec;
    cc = cc.cc;
    tags = sm.h;
    fcc = fieldnames(cc);
    if ~isempty(fcc)
        cci = tags.DD_surfmontage_configs.Value;
        ccc = cc.(fcc{cci});
    else
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'surfmontagecreate')) = [];
        return;
    end
    ecf = ccc{5};
    o.elems = cell(numel(ecf), 15);
    for ecc = 1:numel(ecf)
        o.elems(ecc, :) = ec.(ecf{ecc});
    end
    o.output = [ccc(2:4), {mc.WriteToFile, mc.WriteFilename, mc.SampleVMP, ...
        mc.RestrictVMPsToVOIs, mc.MultiplyThresholds, mc.ThresholdFactors, ...
        mc.SaveAndReuseSmoothed}];
end
outsize = o.output{1};
outcolor = o.output{2};
cpstbars = o.output{3};
writeout = o.output{4};
writefile = o.output{5};
samplevmp = o.output{6};
restrictvmp = o.output{7};
multithresh = o.output{8};
threshfactor = o.output{9};
savesmoothed = o.output{10};

% sampling VMP requires selected maps
if samplevmp
    vmp = ne_gcfg.fcfg.StatsVar;
    vmpi = ne_gcfg.h.StatsVarMaps.Value;
    if numel(vmp) ~= 1 || ~isxff(vmp, {'hdr', 'head', 'vmp'}) || isempty(vmp.Map) || ...
        isempty(vmpi)
        uiwait(warndlg('No VMP (or compatible) maps loaded/selected.', ...
            'NeuroElf - error', 'modal'));
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'surfmontagecreate')) = [];
        return;
    end
    vmpt = lower(vmp.FileType);
end

% prepare process
mfp = ch.MainFig.Pointer;
objs = cell(1, 3);
opwd = pwd;
[tdir, tname] = fileparts(tempname);
tpwd = [tdir, filesep, tname];
try
    mkadir(tpwd);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
if exist(tpwd, 'dir') == 0
    uiwait(warndlg('Error creating temporary folder for files.', ...
        'NeuroElf - warning', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'surfmontagecreate')) = [];
    return;
end
cd(tpwd);
outputok = false;
elems = o.elems;
colinpath = neuroelf_path('colin');
scenery = ne_gcfg.h.Scenery;
upsampling = cini.Surface.SatelliteUpsampling;

% store stats bar setting
showbars = cini.Statistics.ShowThreshBars;
if showbars
    ne_setoption(0, 0, 'showthreshbars', false);
end

% generate output image (in memory)
oi = uint8(0);
oi(outsize(2), outsize(1), 3) = 0;

% max offsets
maxoffx = max(cat(1, elems{:, 14}));
maxoffy = max(cat(1, elems{:, 15}));

% make the call
cprog = ne_progress(0, 0, {true, 0, 'Creating surface montage'});
try

    % set up progress bar
    ch.MainFig.Pointer = 'watch';
    drawnow;

    % set background color
    for pc = 1:3
        oi(:, :, pc) = uint8(outcolor(pc));
    end

    % if not sampling from VMP, all SRFs must be loaded (to be configured!)
    if ~samplevmp

        % get Scenery data
        sc = scenery.UserData(:, 4);
        rsc = sc;
        for oc = 1:numel(sc)
            sc{oc} = sc{oc}.FilenameOnDisk;
            [epath, ename, eext] = fileparts(sc{oc});
            rsc{oc} = [ename, eext];
        end
        esidx = zeros(size(elems, 1), 1);
        for oc = 1:size(elems, 1)
            [epath, ename, eext] = fileparts(elems{oc, 1});
            if isempty(epath)
                tidx = findfirst(strcmpi(rsc, [ename, eext]));
            else
                tidx = findfirst(strcmpi(sc, [epath, filesep, ename, eext]));
            end
            if ~isempty(tidx)
                esidx(oc) = tidx;
            end
        end
        if any(esidx == 0)
            error('neuroelf:guierror:missingSurface', ...
                'Some surface(s) not loaded and configured for montage.');
        end

    % otherwise they must be loadable
    else

        % check files
        for oc = 1:size(elems, 1)
            [epath, ename, eext] = fileparts(elems{oc, 1});
            if isempty(epath)
                if exist([opwd, filesep, ename, eext], 'file') > 0
                    elems{oc, 1} = [opwd, filesep, ename, eext];
                elseif exist([colinpath, filesep, ename, eext], 'file') > 0
                    elems{oc, 1} = [colinpath, filesep, ename, eext];
                else
                    error('neuroelf:guierror:surfaceMissing', ...
                        'Some element surface(s) cannot be found/loaded.');
                end
            else
                if exist(elems{oc, 1}, 'file') == 0
                    error('neuroelf:guierror:surfaceMissing', ...
                        'Some element surface(s) cannot be found/loaded.');
                end
            end
        end

        % copy VMP, keep only relevant maps
        objs{3} = vmp.CopyObject;
        switch (vmpt)
            case 'hdr'
                objs{3}.VoxelData = objs{3}.VoxelData(:, :, :, vmpi);
                objs{3}.ImgDim.Dim(5) = numel(vmpi);
                objs{3}.RunTimeVars.Map = objs{3}.RunTimeVars.Map(vmpi);
            case 'head'
                objs{3}.Brick = objs{3}.Brick(vmpi);
                objs{3}.NrOfVolumes = numel(vmpi);
                objs{3}.RunTimeVars.Map = objs{3}.RunTimeVars.Map(vmpi);
            case 'vmp'
                objs{3}.Map = objs{3}.Map(vmpi);
                objs{3}.NrOfMaps = numel(vmpi);
        end

        % add to main UI
        objs{3}.Browse(1:numel(vmpi));

        % restrict each map to VOI clusters
        if restrictvmp && strcmp(vmpt, 'vmp')
            for mapc = 1:numel(objs{3}.Map)
                if objs{3}.Map(mapc).EnableClusterCheck > 0
                    objs{3}.ClusterTable(mapc);
                    objs{3}.Map(mapc).VMPData = objs{3}.Map(mapc).VMPData .* ...
                        single(objs{3}.Map(mapc).VMPDataCT);
                end
            end
            if ~isempty(ne_gcfg.h.Clusters.Value)
                ne_maskstatswithcls;
                nvmpi = numel(vmpi) + 1;
                switch (vmpt)
                    case 'hdr'
                        objs{3}.VoxelData = objs{3}.VoxelData(:, :, :, nvmpi:end);
                        objs{3}.ImgDim.Dim(5) = numel(vmpi);
                        objs{3}.RunTimeVars.Map = objs{3}.RunTimeVars.Map(1:numel(vmpi));
                    case 'head'
                        objs{3}.Brick = objs{3}.Brick(nvmpi:end);
                        objs{3}.NrOfVolumes = numel(vmpi);
                        objs{3}.RunTimeVars.Map = objs{3}.RunTimeVars.Map(1:numel(vmpi));
                    case 'vmp'
                        objs{3}.Map = objs{3}.Map(nvmpi:end);
                        objs{3}.NrOfMaps = numel(vmpi);
                end
                objs{3}.Browse(1:numel(vmpi));
            end
        end
    end

    % iterate over elements
    pec = 1 / size(elems, 1);
    for oc = 1:size(elems, 1)

        % progress
        ch.Progress.Progress((oc - 1) * pec, ...
            sprintf('Surface montage element %d (prep)...', oc));

        % element offsets
        soffx = elems{oc, 14};
        soffy = elems{oc, 15};
        
        % select element
        if ~samplevmp
            
            % set in Scenery index
            scenery.Value = esidx(oc);
            doubleclick(scenery);

        % load SRF, sample VMP, configure SMP
        else

            % smoothed/inflated
            if elems{oc, 2} > 0 || elems{oc, 4} > 0

                % reuse/save SRFs
                if savesmoothed

                    % filename
                    smoothedname = sprintf('%s_sm%d_%s_inf%d_%s%s', ...
                        elems{oc, 1}(1:end-4), elems{oc, 2}, ...
                        strrep(sprintf('%g', elems{oc, 3}), '.', ''), elems{oc, 4}, ...
                        strrep(sprintf('%g', elems{oc, 5}), '.', ''), elems{oc, 1}(end-3:end));

                    % already exists
                    if exist(smoothedname, 'file') > 0

                        % load
                        objs{1} = xff(smoothedname);

                    % create
                    else
                        objs{1} = xff(elems{oc, 1});
                        ch.Progress.Progress((oc - 1) * pec, ...
                            sprintf('Surface montage element %d (SRF prep)...', oc));
                        if elems{oc, 2} > 0
                            objs{1}.Smooth(elems{oc, 2}, elems{oc, 3}, ...
                                struct('show', false, 'distwsq', true));
                        end
                        if elems{oc, 4} > 0
                            objs{1}.Inflate(elems{oc, 4}, elems{oc, 5}, ...
                                struct('show', false));
                        end

                        % and save
                        objs{1}.SaveAs(smoothedname);
                    end

                % don't re-use
                else

                    % load and smooth/inflate
                    objs{1} = xff(elems{oc, 1});
                    if elems{oc, 2} > 0
                        objs{1}.Smooth(elems{oc, 2}, elems{oc, 3}, ...
                            struct('show', false, 'distwsq', true));
                    end
                    if elems{oc, 4} > 0
                        objs{1}.Inflate(elems{oc, 4}, elems{oc, 5}, ...
                            struct('show', false));
                    end
                end

            % no smoothing/inflation required
            else

                % just load
                objs{1} = xff(elems{oc, 1});
            end

            % add to UI (updates scenery)
            objs{1}.Browse;

            % get SMP (already in GUI then!)
            objs{2} = ne_vmp_createsmp;

            % alter SMP?
            if multithresh
                for mapc = 1:numel(objs{2}.Map)
                    objs{2}.Map(mapc).LowerThreshold = ...
                        threshfactor(1) * objs{2}.Map(mapc).LowerThreshold;
                    objs{2}.Map(mapc).UpperThreshold = ...
                        threshfactor(2) * objs{2}.Map(mapc).UpperThreshold;
                end
            end

            % update with all maps
            objs{2}.Browse(1:numel(objs{2}.Map));
        end

        % progress
        ch.Progress.Progress((oc - 0.75) * pec, ...
            sprintf('Surface montage element %d (create image)...', oc));

        % configure bars
        if cpstbars && soffx == maxoffx && soffy == maxoffy
            ne_setoption(0, 0, 'showthreshbars', true);
        end

        % pop out
        [sat, stags, satid] = ne_undock;
        ne_satsetcolor(0, 0, satid, outcolor);
        ne_setcsrfstatbars(0, 0, satid);
        ne_setsurfpos(0, 0, satid, ...
            {elems{oc, 8:9}, [0, elems{oc, 6}, elems{oc, 7}], elems{oc, 10:11}});

        % resize
        szw = elems{oc, 12};
        szh = elems{oc, 13};
        if szw < 1
            szw = outsize(1) * szw;
        end
        if szh < 1
            szh = outsize(2) * szh;
        end
        srsize = [szw, szh];
        sssize = ones(1, 2) * max(ceil(srsize ./ upsampling));
        ne_satresize(0, 0, satid, sssize);
        sssize = upsampling .* sssize;

        % generate screenshot
        ssname = [tempname '.png'];
        ne_screenshot(0, 0, sat.MLHandle, ssname, 'high-q');

        % close undocked window
        ne_closesatwindow(0, 0, satid);

        % progress
        ch.Progress.Progress((oc - 0.1) * pec, ...
            sprintf('Surface montage element %d (stitch image)...', oc));

        % turn bars off again
        if cpstbars && soffx == maxoffx && soffy == maxoffy
            ne_setoption(0, 0, 'showthreshbars', false);
        end

        % read screenshot
        sscont = imread(ssname);
        delete(ssname);

        % cut out from image
        if ~isequal(srsize, sssize)
            sroff = 1 + floor(0.5 .* (sssize - srsize));
            sscont = sscont(sroff(2):sroff(2)+srsize(2)-1, sroff(1):sroff(1)+srsize(1)-1, :);
        end

        % store in output
        if soffx < 1
            soffx = round(soffx * outsize(1));
        end
        if soffy < 1
            soffy = round(soffy * outsize(2));
        end
        soff = [soffx, soffy];
        oi(soff(2)+1:soff(2)+srsize(2), soff(1)+1:soff(1)+srsize(1), :) = sscont;

        % clear SRF/SMP if needed
        if samplevmp
            objs{2}.ClearObject;
            objs{1}.ClearObject;
            objs(1:2) = {[], []};
        end
    end
    
    % all went well
    outputok = true;

% deal with errors
catch ne_eo;
    uiwait(warndlg(ne_eo.message, 'NeuroElf - warning', 'modal'));
end

% re-set bars
if showbars
    ne_setoption(0, 0, 'showthreshbars', true);
end

% remove objects from memory
clearxffobjects(objs);

% browse VMP again?
if samplevmp && numel(vmp) == 1 && isxff(vmp)
    vmp.Browse(vmpi);
end

% delete all files in this folder
tfiles = findfiles(pwd, '*.*');
if ~isempty(tfiles)
    mdelete(tfiles);
end

% change up
cd('..');

% delete temporary folder
rmdir(tname);

% change back to original folder
cd(opwd);

% output ok
ax = [];
if outputok

    % write output
    if writeout
        imwrite(oi, writefile);

    % show in figure
    else

        % create figure
        nf = figure;
        figure(nf);
        rs = get(0, 'ScreenSize');

        % create axes
        ax = axes('Parent', nf);

        % compute figure size from image
        mims = size(oi);
        mims = mims([2, 1]);

        % compare to screen size
        rc = floor(0.5 * rs(3:4));
        rs = 2 * floor(0.45 * rs(3:4));

        % if image (either width/height) larger than screen
        if any(mims > rs)

            % reduce size to match available space
            di = max(mims ./ rs);
            np = [rc - ceil(0.5 * (mims / di)), ceil(mims / di)];

        % or simply use image size
        else
            np = [rc - ceil(0.5 * mims), mims];
        end

        % figure settings
        set(nf, 'Units', 'pixels');
        set(nf, 'Position', np);
        set(nf, 'Units', 'normalized', 'NumberTitle', 'off', 'Name', ...
            'Surface montage');
        drawnow;
        pause(0.1);
        drawnow;

        % create image
        set(0, 'CurrentFigure', nf);
        set(nf, 'CurrentAxes', ax);
        image(oi);

        % axes settings
        set(ax, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
        set(ax, 'Visible', 'off');

        % set up screen-shot keypress
        set(nf, 'WindowKeyPressFcn', @smc_shot);

        % make sure to update screen
        drawnow;
    end
end

% allow further montage creations
ch.MainFig.Pointer = mfp;
ne_progress(0, 0, cprog);
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'surfmontagecreate')) = [];
if nargout > 0
    varargout{1} = oi;
    if nargout > 1
        varargout{2} = ax;
    end
end



% keypress screenshot
function smc_shot(src, ke, varargin)

% get Key and Modifier from keyboard event (see Matlab docu!)
kk = ke.Key;
mn = ke.Modifier;

% screenshot key combination
if numel(mn) == 1 && strcmpi(mn{1}, 'shift') && strcmpi(kk, 's')

    % create screenshot
    if strcmpi(get(src, 'Type'), 'figure')
        ne_screenshot(0, 0, src, '', 'high-q');
    else
        ne_screenshot(0, 0, src);
    end
end

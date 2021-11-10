function varargout = ne_openfile(varargin)
% ne_openfile  - add file to workspace and, possibly, controls
%
% FORMAT:       [obj = ] ne_openfile(SRC, EVT, filespec [, toolloaded])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       filespec    either filename or object reference
%       toolloaded  1x1 boolean flag to set whether loaded by the GUI
%
% Output fields:
%
%       obj         object reference (e.g. for filenames)
%
% Examples:
%
%    vmr = ne_openfile(0, 0, 'subject123_IHC.vmr');
%    ne_openfile(0, 0, srf_object, false);

% Version:  v1.1
% Build:    21111013
% Date:     Nov-10 2021, 1:15 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2021, Jochen Weber
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% not while in a callback
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% current status of pointer
mfp = ne_gcfg.h.MainFig.Pointer;

% get handles
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
ci = ne_gcfg.c.ini;
try

    % a file (name or handle) was given
    if nargin > 2

        % if not is a handle
        if numel(varargin{3}) ~= 1 || ~isxff(varargin{3}, true)

            % set pointer to watch
            ne_gcfg.h.MainFig.Pointer = 'watch';
            drawnow;

            % pass on to xff
            loadedds = ne_gcfg.c.xff.Documents();
            f = ne_gcfg.c.xff.Document(varargin{3}, true);
            toolloaded = ~isequal(loadedds, ne_gcfg.c.xff.Documents());
            if nargout > 0
                varargout{1} = f;
            end
            try
                ne_gcfg.c.incb = false;
                ne_openfile(0, 0, f, toolloaded, true);
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
            ne_gcfg.h.MainFig.Pointer = mfp;
            return;

        % otherwise
        else

            % get handle into f and object's struct
            f = varargin{3};
            if nargout > 0
                varargout{1} = f;
            end

            % update?
            if nargin < 5 || ~islogical(varargin{5}) || numel(varargin{5}) ~= 1
                updatecvar = true;
            else
                updatecvar = varargin{5};
            end

            % get the list of loaded objects (SliceVar)
            svars = ch.SliceVar.UserData;

            % compare to list, and if found, simply set index
            for loc = 1:size(svars, 1)
                if f == svars{loc, 4}
                    ne_gcfg.c.incb = false;
                    if updatecvar
                        ch.SliceVar.Value = loc;
                        ne_setcvar;
                    end
                    ne_gcfg.h.MainFig.Pointer = mfp;
                    return;
                end
            end

            % get the list of loaded objects (StatsVar)
            svars = ch.StatsVar.UserData;

            % compare to list, and if found, simply set index
            for loc = 1:size(svars, 1)
                if f == svars{loc, 4}
                    ne_gcfg.c.incb = false;
                    if updatecvar
                        ch.StatsVar.Value = loc;
                        ne_setcstats;
                    end
                    ne_gcfg.h.MainFig.Pointer = mfp;
                    return;
                end
            end

            % get the list of loaded objects (SurfVar)
            svars = ch.SurfVar.UserData;

            % compare to list, and if found, simply set index
            for loc = 1:size(svars, 1)
                if f == svars{loc, 4}
                    ne_gcfg.c.incb = false;
                    if updatecvar
                        ch.SurfVar.Value = loc;
                        if numel(ch.Scenery.Value) < 2
                            svar = ch.Scenery.UserData;
                            for sloc = 1:size(svar, 1)
                                if f == svars{sloc, 4}
                                    ch.Scenery.Value = sloc;
                                    ne_sceneryselect;
                                end
                            end
                        end
                        ne_setcsrf;
                    end
                    ne_gcfg.h.MainFig.Pointer = mfp;
                    return;
                end
            end

            % get the list of loaded objects (SurfStatsVar)
            if ~isempty(svars)
                svars = ch.SurfStatsVar.UserData;
                for loc = 1:size(svars, 1)
                    if f == svars{loc, 4}
                        ne_gcfg.c.incb = false;
                        if updatecvar
                            ch.SurfStatsVar.Value = loc;

                            % also get time
                            rtv = svars{loc, 4}.RunTimeVars;
                            if isfield(rtv, 'AvgMTC') && islogical(rtv.AvgMTC) && ...
                                numel(rtv.AvgMTC) == 1 && rtv.AvgMTC
                                ne_gcfg.fcfg.srfcfg.time = 0.001 * ...
                                    ((rtv.SubMapVol - 1) * rtv.AvgWindowStep + rtv.AvgWindowFrom);
                                ne_setsurfpos;
                            end

                            % then update
                            ne_setcsrfstats;
                        end
                        ne_gcfg.h.MainFig.Pointer = mfp;
                        return;
                    end
                end
            end

            % otherwise, get workspace variables
            lfld = fieldnames(ne_gcfg.w);

            % see if already loaded
            for loc = 1:numel(lfld)
                if f == ne_gcfg.w.(lfld{loc});
                    ne_gcfg.c.incb = false;
                    ne_gcfg.h.MainFig.Pointer = mfp;
                    return;
                end
            end

            % default: we did NOT load this (was a handle!)
            if nargin < 4 || ~islogical(varargin{4}) || numel(varargin{4}) ~= 1
                toolloaded = false;
            else
                toolloaded = varargin{4};
            end
        end

    % without argument
    else

        % set pointer to watch
        ne_gcfg.h.MainFig.Pointer = 'watch';
        drawnow;

        % and let xff do the rest
        f = xff('*.*', 0);
        if isempty(f)
            ne_gcfg.c.incb = false;
            ne_gcfg.h.MainFig.Pointer = mfp;
            drawnow;
            return;
        end
        toolloaded = true;
        if nargout > 0
            varargout{1} = f;
        end

        % try to switch to loaded file's dir
        try
            if ~isempty(f)
                cd(fileparts(f.FilenameOnDisk(true)));
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        
        % pass on
        try
            ne_gcfg.c.incb = false;
            ne_openfile(0, 0, f, toolloaded, true);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ne_gcfg.h.MainFig.Pointer = mfp;
        return;
    end

    % check validity of handle
    if numel(f) ~= 1 || ~isxff(f, true)
        ne_gcfg.c.incb = false;
        ne_gcfg.h.MainFig.Pointer = mfp;
        return;
    end

% error handling
catch ne_eo;
    ne_gcfg.c.incb = false;
    ne_gcfg.h.MainFig.Pointer = mfp;
    uiwait(warndlg(ne_eo.message, 'NeuroElf GUI - error message', 'modal'));
    return;
end

% echo
fname = f.FilenameOnDisk(2);
if isempty(fname)
    try
        fname = f.RunTimeVars.xffID;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
if ne_gcfg.c.echo && ...
   ~isempty(fname)
    if nargin > 3 && ...
       ~isempty(varargin{4})
        vargs = any2ascii(varargin(4:end));
        ne_echo({'neuroelf_gui(''openfile'', xff(''%s''), %s)', fname, vargs(2:end-1)});
    else
        ne_echo({'neuroelf_gui(''openfile'', xff(''%s''))', fname});
    end
end

% get filename
[fn{1:3}] = fileparts(fname);
ftype = lower(f.Filetype);
frtv = f.RunTimeVars;

% and create variable name
vn = makelabel([upper(ftype) '_' fn{2}]);
vnl = numel(vn);
if vnl > 26
    vn = vn(1:26);
    vnl = 26;
end

% in case this name is already in the, recreate until uniquely defined
while isfield(ne_gcfg.w, vn)
    vn = sprintf('%s_%04x', vn(1:vnl), floor(65535.999 * rand(1)));
end

% then add to workspace and workspace control
ne_gcfg.w.(vn) = f;
ne_gcfg.wc.(vn) = toolloaded;

% for hdr (Anayze/NIftI) or head (AFNI/BRIK/HEAD) files, check for stats
if (strcmp(ftype, 'hdr') && ...
    (~isempty(regexpi(f.DataHist.Description, ...
          'spm\{[ft]_\[\d+(\.\d+)?(,\s*[1-9][0-9\.]*)?\]\}')) || ...
     ~isempty(regexpi(f.DataHist.Description, 'spm_spm\:beta\s*\(\d+\)')) || ...
     (isfield(frtv, 'StatsObject') && islogical(frtv.StatsObject) && ...
      numel(frtv.StatsObject) == 1 && frtv.StatsObject))) || ...
   (strcmp(ftype, 'head') && ~isempty(strfind(f.TypeOfVolumes, '_FUNC')))

    % parse data?
    if ~isfield(frtv, 'StatsObject') || ~islogical(frtv.StatsObject) || ...
        numel(frtv.StatsObject) ~= 1 || ~frtv.StatsObject

        % get Map field (transparent from RunTimeVars!)
        map = f.Map;

        % and ensure that certain fields are set
        switch ftype

            % for HDR
            case 'hdr'

                % test whether map is good to go
                if numel(map) ~= size(f.VoxelData, 4)
                    map = map(ones(1, size(f.VoxelData, 4)));
                end

                % get description string
                stdesc = regexpi(f.DataHist.Description, ...
                    'spm\{([ft])_(\[[^\}]+\])\}', 'tokens');
                if ~isempty(stdesc)
                    sttype = stdesc{1}{1};
                    stdf = eval(stdesc{1}{2}, '[10, 100]');
                else
                    sttype = 'b';
                    stdf = [240, 1];
                end

                % for each volume
                for vc = 1:numel(map)

                    % check type
                    if numel(map(vc).Type) ~= 1 || map(vc).Type < 1
                        map(vc).Type = find(strcmpi(sttype, {'t', 'r', 'c', 'f'}));
                        if isempty(map(vc).Type)
                            map(vc).Type = 15;
                        end
                    end

                    % check d.f.
                    if numel(map(vc).DF1) ~= 1 || map(vc).DF1 < 1
                        map(vc).DF1 = stdf(1);
                    end
                    if map(vc).Type == 4 && (numel(map(vc).DF2) ~= 1 || map(vc).DF2 < 1)
                        map(vc).DF2 = stdf(2);
                    end
                end

            % for HEAD
            case 'head'

                % test whether RunTimeVars.Map is good to go
                if numel(map) ~= numel(f.Brick)
                    map = map(ones(1, numel(f.Brick)));
                end

                % some definitions
                brik_ftypes = [15, 15, 2, 1, 4, 12, 15, 15, 15, 15, 15, 15];

                % check bricks
                for vc = 1:numel(map)
                    % check type
                    if numel(map(vc).Type) ~= 1 || map(vc).Type < 1
                        map(vc).Type = brik_ftypes(f.Brick(vc).FuncType + 1);
                    end

                    % check d.f.
                    stdf = [f.Brick(vc).FuncParams, 0];
                    if numel(stdf) == 1
                        stdf = [10, 1, 1];
                    end
                    if numel(map(vc).DF1) ~= 1 || map(vc).DF1 < 1
                        if map(vc).Type ~= 2
                            map(vc).DF1 = stdf(1);
                        else
                            map(vc).DF1 = stdf(1) - sum(stdf(2:end));
                        end
                    end
                    if map(vc).Type == 4 && (numel(map(vc).DF2) ~= 1 || map(vc).DF2 < 1)
                        map(vc).DF2 = stdf(2);
                    end
                end
        end

        % reassign updated map
        f.Map = map;
        f.RunTimeVars.StatsObject = true;
        try
            f.SaveRunTimeVars;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end

    % and then replace it for later
    ftype = 'hfunc';

% for component files, check documenttype
elseif strcmp(ftype, 'cmp')

    % what type...
    switch (f.DocumentType)

        % FMR
        case 0
            ftype = 'fcmp';

        % VTC
        case 1
            ftype = 'vcmp';

        % MTC
        case 2
            ftype = 'mcmp';
    end

% for average VTCs
elseif strcmp(ftype, 'vtc') && isfield(frtv, 'AvgVTC') && frtv.AvgVTC
    ftype = 'atc';
end

% recent files number
cir = ci.RecentFiles;
rfn = cir.Number - 1;

% depending on file type
switch (ftype)

    % for supported stats maps
    case {'atc', 'ava', 'glm', 'hfunc', 'map', 'vcmp', 'vmp'}

        % ensure the object has a map selection field
        if numel(frtv.MapSelection) < 2 || isempty(frtv.MapSelection{2})
            mnames = f.MapNames;
            if ~isempty(mnames)
                f.RunTimeVars.MapSelection = {mnames(1), 1};
            else
                f.RunTimeVars.MapSelection = {{}, []};
            end
        end

        % for VMPs, make sure the colors are set
        if any(strcmp(ftype, {'atc', 'ava', 'glm', 'hfunc', 'vcmp', 'vmp'}))
            f.SetColors([], 'xauto');
        end

        % then add to StatsVar list
        vtype = 'vox';
        if ~any(strcmp(ftype, {'ava', 'glm'})) || ...
           (strcmp(ftype, 'glm') && f.ProjectType < 2) || ...
           (strcmp(ftype, 'ava') && f.ProjectType == 2)
            cv_addtolist(ch.StatsVar, f, vn);
        elseif (strcmp(ftype, 'glm') && f.ProjectType == 2) || ...
           (strcmp(ftype, 'ava') && f.ProjectType == 3)
            smaps = f.RunTimeVars.Map;
            for vc = 1:numel(smaps)
                if isempty(smaps(vc).LowerThreshold)
                    smaps(vc).LowerThreshold = 0.0001;
                end
                if isempty(smaps(vc).UpperThreshold)
                    smaps(vc).UpperThreshold = max(smaps(vc).LowerThreshold + 0.5, 1);
                end
            end
            f.RunTimeVars.Map = smaps;
            cv_addtolist(ch.SurfStatsVar, f, vn);
            vtype = 'surf';
        end
        ne_gcfg.c.incb = false;

        % for toolloaded, add to recent files
        if toolloaded && ~isempty(fname)
            cir.StatsVar(strcmp(cir.StatsVar, fname)) = [];
            ci.RecentFiles.StatsVar = [{fname}; ...
                cir.StatsVar(1:min(rfn, numel(cir.StatsVar)))];
        end

        % and set current object
        if updatecvar
            if vtype(1) == 'v'
                ne_setcstats;
            else
                ne_setcsrfstats;
                if isempty(ch.Scenery.Value) && ~isempty(ch.Scenery.UserData)
                    nvert = f.NrOfVertices;
                    scu = ch.Scenery.UserData;
                    for scuc = 1:size(scu, 1)
                        if isxff(scu{scuc, 4}, true) && scu{scuc}.NrOfVertices == nvert
                            ch.Scenery.Value = scuc;
                            ne_sceneryselect;
                            break;
                        end
                    end
                end
            end
        end

        % for GLMs, test if the contrast manager is open
        if strcmp(ftype, 'glm') && f.ProjectTypeRFX > 0 && isfield(ch, 'CM') && ...
            isstruct(ch.CM) && isfield(ne_gcfg.fcfg, 'CM') && isstruct(ne_gcfg.fcfg.CM)

            % check if GLM is already in the list
            ff = [];
            for glmc = 1:numel(ne_gcfg.fcfg.CM.GLMs)
                if ne_gcfg.fcfg.CM.GLMs{glmc} == f
                    ff = f;
                    break;
                end
            end

            % not found
            if ~isxff(ff)

                % then add to lists
                ne_gcfg.fcfg.CM.GLMs{end+1} = f;
                glms = ch.CM.h.GLMs.String;
                if ~iscell(glms)
                    glms = cellstr(glms);
                end
                if isempty(fname)
                    fnname = sprintf('<#%d: %d subs, %d preds>', ...
                        f.Filenumber, f.NrOfSubjects, f.NrOfSubjectPredictors);
                else
                    [fnpath, fnname] = fileparts(fname);
                end
                glms{end+1} = fnname;
                ch.CM.h.GLMs.String = glms;
                ch.CM.h.GLMs.Value = numel(glms);
                if updatecvar
                    ne_cm_setglm;
                end
            end
        end

    % for OLT/LUT objects
    case 'olt'

        % replace currently loaded object
        try
            ne_gcfg.lut.ClearObject;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ne_gcfg.lut = f;

        % and check the loaded VMP/MAP, if replacement is in order
        if cc.page < 3
            stvar = cc.StatsVar;
            stvix = cc.StatsVarIdx;
            typt = {'hdr', 'vmp'};
        else
            stvar = cc.SurfStatsVar;
            stvix = cc.SurfStatsVarIdx;
            typt = {'fsmf', 'smp'};
        end
        if numel(stvar) == 1 && isxff(stvar, typt) && numel(stvix) == 1

            % replace old LUTName
            try
                stvar.Map(stvix).LUTName = f.FilenameOnDisk;

                % and set OverlayColors, if requested
                if stvar.Map(stvix).UseRGBColor <= 0
                    stvar.Map(stvix).OverlayColors = f.Colors;
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end
        ne_gcfg.c.incb = false;

        % update output
        if updatecvar
            if any(strcmp(typt, 'vmp'))
                ne_setslicepos;
            else
                ne_setcsrfstatmap;
            end
        end

    % for anything NeuroElf GUI can display as "anatomical" information
    case {'dmr', 'fmr', 'hdr', 'head', 'mgh', 'msk', 'vmr', 'vtc'}

        % for VTC
        if strcmp(ftype, 'vtc') && isempty(f.VTCData)
            if toolloaded
                f.ClearObject;
            end
            if nargin < 3
                uiwait(warndlg('Cannot show empty VTC.', 'NeuroElf - info', 'modal'));
            end
            ne_gcfg.h.MainFig.Pointer = mfp;
            ne_gcfg.c.incb = false;
            return;
        elseif strcmp(ftype, 'fmr') && istransio(f.Slice(1).STCData)
            f.LoadSTC;
        elseif strcmp(ftype, 'hdr') && istransio(f.VoxelData)
            f.LoadVoxelData;
        elseif (strcmp(ftype, 'head') && istransio(f.Brick(1).Data)) || ...
            strcmp(ftype, 'mgh') && istransio(f.MGHData)
            f.LoadTransIOData;
        end

        % add to SliceVar
        cv_addtolist(ch.SliceVar, f, vn);
        ne_gcfg.c.incb = false;

        % for toolloaded, add to recent files
        if toolloaded
            fname = f.FilenameOnDisk(2);
            if ~isempty(fname)
                cir.SliceVar(strcmp(cir.SliceVar, fname)) = [];
                ci.RecentFiles.SliceVar = [{fname}; ...
                    cir.SliceVar(1:min(rfn, numel(cir.SliceVar)))];
            end
        end

        % make sure we define the gradient type
        if ~isfield(frtv, 'DisplayType') || ~ischar(frtv.DisplayType)
            f.RunTimeVars.DisplayType = 'gray';
        end

        % make sure a good scaling window is available
        if ~isfield(frtv, 'ScalingWindow') || ...
           ~isfield(frtv, 'ScalingWindowLim') || ...
           ~isfield(frtv, 'ScalingHist')
            f.SetScalingWindow;
        end

        % add slicing range bounding box
        if ~isfield(frtv, 'RenderBBox')
            if any(strcmp(ftype, {'dmr', 'fmr'}))
                f.RunTimeVars.RenderBBox = ones(3, 2);
                f.RunTimeVars.RenderBBox(:, 2) = ...
                    [f.ResolutionX; f.ResolutionY; f.NrOfSlices];
            else
                f.RunTimeVars.RenderBBox = [ones(3, 1), f.BoundingBox.DimXYZ(1:3)'];
            end
            f.RunTimeVars.RenderBBoxFull = f.RunTimeVars.RenderBBox;
        end

        % and update output
        if updatecvar
            ne_setcvar;
        end

    % for surface files
    case {'fsbf', 'srf', 'tom'}

        % empty?
        if isempty(f.VertexCoordinate)
            if toolloaded
                f.ClearObject;
            end
            ne_gcfg.w = rmfield(ne_gcfg.w, vn);
            ne_gcfg.wc = rmfield(ne_gcfg.wc, vn);
            uiwait(warndlg('Cannot display empty surface.', 'NeuroElf - information', 'modal'));
            ne_gcfg.c.incb = false;
            ne_gcfg.h.MainFig.Pointer = mfp;
            return;
        end

        hnd = handles(f);
        if ~isfield(hnd, 'Surface') || ~ishandle(hnd.Surface)

            % create transimg for coloring
            scfg = ci.Surface;
            numv = size(f.VertexCoordinate, 1);
            ptio = transimg(1, numv, scfg.BackgroundColor);
            render(ptio);

            % set initial handles
            if ~isfield(hnd, 'MTC') || numel(hnd.MTC) ~= 1 || ~isxff(hnd.MTC, 'mtc')
                f.SetHandle('MTC', []);
            end
            if ~isfield(hnd, 'Stats') || ~iscell(hnd.Stats) || numel(hnd.Stats) ~= 2 || ...
                numel(hnd.Stats{1}) ~= 1 || ~isxff(hnd.Stats{1}, {'fsmf', 'smp'})
                f.SetHandle('Stats', {[], []});
            else
                f.SetHandle('Stats', hnd.Stats(:)');
            end
            f.SetHandle('SurfProps', {[0, 0, 0], [0, 0, 0], [1, 1, 1], 1, 'f', [], 'none'});
            f.SetHandle('SurfTIO', ptio);
            f.SetHandle('VertexMorphIndex', 0);
            
            % prep: Coordinate, Normal, MeshCenter, TIO.Rendered, SurfProps
            if ~isfield(hnd, 'VertexMorphMeshes') || ...
                isempty(hnd.VertexMorphMeshes) || ...
               ~iscell(hnd.VertexMorphMeshes) || ...
                size(hnd.VertexMorphMeshes, 2) ~= 5 || ...
                size(hnd.VertexMorphMeshes{1, 1}, 1) ~= numv
                f.SetHandle('VertexMorphMeshes', cell(0, 5));
                hnd.VertexMorphMeshes = cell(0, 5);
            end

            % get coordinates and normals
            [p, pn] = btc_meshcn(f, ne_gcfg.fcfg.srfcfg, true);
            np = size(p, 1);

            % prepare axes for adding the patch
            hold(ne_gcfg.h.Surface, 'on');
            
            % create sub-group
            hgtrf = hgtransform('Parent', ne_gcfg.h.SurfaceTransform);

            % regular (with faces) patch
            if ~isempty(f.TriangleVertex)

                % create surface patch
                hsrf = patch('Faces', f.TriangleVertex(:, [1, 3, 2]), ...
                    'Vertices', p, 'VertexNormals', pn, ...
                    'FaceVertexCData', [0, 0, 0], 'FaceColor', 'none', ...
                    'EdgeColor', 'none', 'Parent', hgtrf);

                % and make some more initial settings
                set(hsrf, 'AmbientStrength', scfg.LightAmbient, ...
                    'BackFaceLighting', 'unlit', ...
                    'ButtonDownFcn', @ne_btdown, ...
                    'DiffuseStrength', scfg.LightDiffuse, ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', scfg.Alpha, 'FaceColor', 'interp', ...
                    'EdgeLighting', 'gouraud', 'FaceLighting', 'gouraud', ...
                    'LineStyle', 'none', ...
                    'SpecularStrength', scfg.LightSpecular, ...
                    'SpecularExponent', scfg.LightSpecularExponent, ...
                    'SpecularColorReflectance', scfg.LightSpecularReflectance, ...
                    'UserData', struct('SRF', f), 'Visible', 'off');

            % fiber tracking (no faces)
            else

                % create surface patch
                hsrf = patch(p(:, 1), p(:, 2), p(:, 3), 0, 'Parent', hgtrf);
                set(hsrf, 'CData', (1 / 255) .* double(reshape(f.VertexColor(:, 2:4), [np, 1, 3])), ...
                    'EdgeColor', 'interp', 'EdgeAlpha', 1, 'FaceAlpha', 0, 'FaceColor', 'none', ...
                    'AmbientStrength', scfg.LightAmbient, ...
                    'DiffuseStrength', scfg.LightDiffuse, ...
                    'SpecularStrength', scfg.LightSpecular, ...
                    'SpecularExponent', scfg.LightSpecularExponent, ...
                    'SpecularColorReflectance', scfg.LightSpecularReflectance);
                set(hsrf, 'Visible', 'off');
                f.SetHandle('SurfProps', {[0, 0, 0], [0, 0, 0], [1, 1, 1], 1, 'w', [], 'none'});
            end

            % set handles
            f.SetHandle('Surface', hsrf);
            f.SetHandle('SurfaceTransform', hgtrf);
            f.SetHandle('SUpdate', {@ne_srfupdatecoords, f, hsrf});

            % and colorize
            btc_meshcolor(f);
        end
        if ~isfield(hnd, 'SurfStatsBars') || ~iscell(hnd.SurfStatsBars)
            f.SetHandle('SurfStatsBars', {zeros(cc.SurfBarSize(1), 0, 3)});
        end

        % add to SurfVar *and* Scenary
        cv_addtolist(ch.SurfVar, f, vn);
        if isempty(hnd.VertexMorphMeshes)
            cv_addtolist(ch.Scenery, f, vn, {'NrOfVertices'}, 'Scene compilation (empty)');
        else
            cv_addtolist(ch.Scenery, f, vn, {'NrOfVertices', {'NrOfMorphStates: %d', ...
                size(hnd.VertexMorphMeshes, 1) + 1}}, 'Scene compilation (empty)');
        end

        % for toolloaded, add to recent files
        if toolloaded
            fname = f.FilenameOnDisk(2);
            if ~isempty(fname)
                cir.SurfVar(strcmp(cir.SurfVar, fname)) = [];
                ci.RecentFiles.SurfVar = [{fname}; ...
                    cir.SurfVar(1:min(rfn, numel(cir.SurfVar)))];
            end
        end

        % and update output
        if updatecvar
            ch.Scenery.Value = size(ch.Scenery.UserData, 1);
            ne_sceneryselect;
        end

    % for surface statistics
    case {'fsmf', 'mtc', 'smp', 'scmp'}

        % ensure the object has a map selection field
        if ~isfield(f.RunTimeVars, 'MapSelection') || ...
           ~iscell(f.RunTimeVars.MapSelection ) || ...
            numel(f.RunTimeVars.MapSelection) ~= 2
            mnames = f.MapNames;
            if ~isempty(mnames)
                f.RunTimeVars.MapSelection = {mnames(1), 1};
            else
                f.RunTimeVars.MapSelection = {{}, []};
            end
        end

        % ensure colors are set
        f.SetColors([], 'xauto');

        % add to list
        cv_addtolist(ch.SurfStatsVar, f, vn);

        % for toolloaded, add to recent files
        if toolloaded
            fname = f.FilenameOnDisk(2);
            if ~isempty(fname)
                cir.SurfStatsVar(strcmp(cir.SurfStatsVar, fname)) = [];
                ci.RecentFiles.SurfStatsVar = [{fname}; ...
                    cir.SurfStatsVar(1:min(rfn, numel(cir.SurfStatsVar)))];
            end
        end

        % and set current object
        if updatecvar
            ne_setcsrfstats;
            if isempty(ch.Scenery.Value) && ~isempty(ch.Scenery.UserData)
                nvert = f.NrOfVertices;
                scu = ch.Scenery.UserData;
                for scuc = 1:size(scu, 1)
                    if isxff(scu{scuc, 4}, true) && scu{scuc, 4}.NrOfVertices == nvert
                        ch.Scenery.Value = scuc;
                        ne_sceneryselect;
                        break;
                    end
                end
            end
        end
end

% and don't forget the pointer !
ne_gcfg.c.incb = false;
ne_gcfg.h.MainFig.Pointer = mfp;

% and put into output
if nargout > 0
    varargout{1} = f;
end

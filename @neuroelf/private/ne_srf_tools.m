function varargout = ne_srf_tools(varargin)
% ne_srf_tools  - tools for currently selected SurfVar
%
% FORMAT:       [result = ] ne_srf_tools(SRC, EVT, action [, arguments, ...])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       action      one of the supported actions (see below)
%       arguments   additional (usually optional) arguments
%
% Output fields:
%
%       result      resulting object (if any)
%
% Notes: calls that do not specify a surface object will always use the
%        currently selected one: ne_gcfg.fcfg.SurfVar --if no surface is
%        open in the GUI, those calls are then invalid (no action taken)
%
% Actions:
%
%   'backproject'   - project SRF object back into a ISOmetric VMR
%       syntax:     ne_srf_tools(0, 0, 'backproject' [, VMR [, CONFIG]])
%       VMR         optional VMR, use instead of ne_gcfg.fcfg.SliceVar
%       CONFIG      1x3 cell array, {COLOR_CODE, TOS_FACTOR, SAMPLING}
%                   each value MUST be a string (parsed), whereas 
%                   COLOR_CODE is the value set at each nearest coordinate
%                   TOS_FACTOR is a triangle-oversampling factor
%                   SAMPLING is the from:step:to along the normal vectors
%       example:    ne_srf_tools(0, 0, 'backproject', [], {'240','3','-1:1:1'})
%       notes:      this uses the SRF::BackToVMR method and sets the
%                   VMR.RunTimeVars.UndoBuffer with the current content if
%                   it is not already set (like drawing does)
%
%   'cancelmorph'   - cancel the current morphing (GUI queued only)
%       syntax:     ne_srf_tools(0, 0, 'cancelmorph')
%
%   'clone'         - clone existing SRF (AFT::CopyObject) and add to GUI
%       syntax:     ne_srf_tools(0, 0, 'clone' [, SRF])
%       SRF         optional SRF, use instead of ne_gcfg.fcfg.SurfVar;
%
%   'closefiles'    - close selected files in Scenery list (incl. SMPs)
%       syntax:     ne_srf_tools(0, 0, 'closefiles')
%
%   'clustersrf'    - generate surface from VOI clusters
%       syntax:     ne_srf_tools(0, 0, 'clustersrf' [, VOI [, VIDX]])
%       VOI         optional VOI object, use instead of ne_gcfg.voi
%       VIDX        optional indices of VOIs within VOI object
%
%   'createsphere'  - create a sphere (icosahedron) mesh
%       syntax:     ne_srf_tools(0, 0, 'createsphere' [, CONFIG])
%       CONFIG      1x5 cell array (must all be given) with string (!) fields
%                   {radius_mm, center, triangulation_factor, rgb1, rgb2}
%       example:    ne_srf_tools(0, 0, 'createsphere', {'120', '0 0 0', '6','180 180 180', '210 210 210'})
%
%   'curvsmp'       - create and open a curvature SMP for current SRF
%       syntax:     ne_srf_tools(0, 0, 'curvsmp')
%
%   'denssmp'       - create and open a density SMP for current SRF
%       syntax:     ne_srf_tools(0, 0, 'denssmp')
%
%   'findintensity' - morph (spherical) SRF towards intensity in slicing var
%       syntax:     ne_srf_tools(0, 0, 'findintensity' [, CONFIG [, NSC]])
%       CONFIG      1x4 cell array (must all be given) with string (!) fields
%                   {target_intensity, pull_force, iterations, in_or_out}
%                   whereas in_or_out is either 'in' or 'out'
%       NSC         if given and set to 'nospherecheck' skip checking SRF
%       example:    ne_srf_tools(0, 0, 'findintensity', {'80', '0.5', '100', 'in'})
%
%   'inflate'       - inflate currently selected SRF
%       syntax:     ne_srf_tools(0, 0, 'inflate' [, CONFIG])
%       CONFIG      1x2 cell array (must all be given) with string (!) fields
%                   {iterations, smoothing_force}
%       example:    ne_srf_tools(0, 0, 'inflate', {'5000', '0.5'})
%
%   'loadscenery'   - load a previously created scenery file
%       syntax:     ne_srf_tools(0, 0, 'loadscenery' [, SCNFILE])
%       SCNFILE     MAT-filename with scenery content
%
%   'loadsubcort'   - load all sub-cortical meshes (if available)
%       syntax:     ne_srf_tools(0, 0, 'loadsubcort' [, ICBMNORM])
%       ICBMNORM    if set to 'icbmnorm', load ICBMnorm versions
%
%   'maptoico'      - map spherized surface to icosahedron mesh
%       syntax:     ne_srf_tools(0, 0, 'maptoico' [, CONFIG])
%       CONFIG      1x3 cell array (must all be given) with string (!) fields
%                   {triangulation_factor, ssm_or_tsm, inverse_also}
%                   whereas ssm_or_tsm should be either 'ssm' or 'tsm'
%       example:    ne_srf_tools(0, 0, 'maptoico', {'6', 'tsm', 'y'})
%
%   'morph'         - general interface to SRF::Morph
%       syntax:     ne_srf_tools(0, 0, 'morph' [, CONFIG])
%       CONFIG      1x9 cell array (must all be given) with string (!) fields
%                   {niter, force, areac, areaw, distc, distw, norm, nr, sph},
%                   whereas the fields are:
%                   niter - number of iterations (e.g. '1000')
%                   force - general smoothing force (e.g. '0.1')
%                   areac - keep the area constant ('yes' or 'no')
%                   areaw - area-based weighting of smoothing ('yes' or 'no')
%                   distc - distortion correction force (e.g. '0.5')
%                   distw - distance weighting (either 'no', 'log', or 'sq')
%                   norm  - force along normal vector (e.g. '0.01')
%                   nr    - ramp up along-normal force ('yes' or 'no')
%                   sph   - to-sphere force (e.g. '0.001')
%
%   'savescenery'   - save all meshes with configuration in a MAT file
%       syntax:     ne_srf_tools(0, 0, 'savescenery' [, SCNFILE [, REFS]])
%       SCNFILE     MAT-filename to use for writing scenery MAT-file
%                   if the filename has a ".ply" extension, write PLY file
%       REFS        if set to true (and using a MAT-file), store references
%                   to existing SRF/SMP file objects instead of data
%
%   'setcolors'     - set two main colors of currently selected surface
%       syntax:     ne_srf_tools(0, 0, 'setcolors' [, COLORS])
%       COLORS      2x3 RGB colors (0..255), will be scaled to 0..1
%
%   'setmorphtarget' - set morphing target in current SRF
%       syntax:     ne_srf_tools(0, 0, 'setmorphtarget' [, TSRF])
%       TSRF        target SRF
%
%   'smooth'        - smooth current surface
%       syntax:     ne_srf_tools(0, 0, 'smooth' [, CONFIG])
%       CONFIG      1x3 cell array (must all be given) with string (!) fields
%                   {iterations, smoothing_force, weighting}, whereas the
%                   weighting can be one of 'n' (none), 'd' (distance-based
%                   weighting), or 'sq' (square-of-distance weighting)
%
%   'smoothsmp'     - smooth currently selected SMP map
%       syntax:     ne_srf_tools(0, 0, 'smoothsmp' [, SMP [, MAPS [, CONFIG]]])
%       SMP         indicate alternative SMP object if not current
%       MAPS        indicate other than selected maps
%       CONFIG      1x2 cell array (must all be given) with string (!) fields
%                   {smooth_kernel_mm, iterations}
%
%   'smpformula'    - SMP::ComputeFormula interface
%       syntax:     ne_srf_tools(0, 0, 'smpformula')
%
%   'srfinfo'       - print information about surface to the console
%       syntax:     ne_srf_tools(0, 0, 'srfinfo')
%
%   'tosphere'      - morph currently selected SRF to spherical surface
%       syntax:     ne_srf_tools(0, 0, 'tosphere' [, CONFIG])
%       CONFIG      1x3 cell array (must all be given) with string (!) fields
%                   {iterations_array, smoothing_forces_array,
%                   to_sphere_forces_array} whereas array sizes must match
%
%   'undomorph'     - revert current surface to state prior to last morph
%       syntax:     ne_srf_tools(0, 0, 'undomorph')

% Version:  v1.1
% Build:    16060813
% Date:     Jun-08 2016, 1:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
ci = ne_gcfg.c.ini.Surface;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% active SRF loaded
srf = cc.SurfVar;
if (numel(srf) ~= 1 || ...
   ~isxff(srf, {'fsbf', 'srf'})) && ...
   (nargin < 3 || ...
    ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
    ~any(strcmpi(varargin{3}(:)', ...
     {'clustersrf', 'createsphere', 'findintensity', 'loadscenery', ...
      'loadsubcort', 'recosmsph'})))
    return;
end
if isxff(srf, {'fsbf', 'srf'})
    srfh = handles(srf);
else
    srfh = struct( ...
        'xff',          -1, ...
        'ShownInGUI',   true, ...
        'MTC',          [], ...
        'Stats',        {{struct('Filetype', 'NONE'), []}}, ...
        'SurfProps',    {{[0, 0, 0], [0, 0, 0], [1, 1, 1], 1, 'f', [], 'none'}}, ...
        'Surface',      -1, ...
        'CancelMorph',  false, ...
        'Morphing',     false);
end

% blocked?
if any(strcmpi(ne_gcfg.c.blockcb, 'ne_srf_tools')) && ...
   ~strcmpi(varargin{3}(:)', 'cancelmorph')
    return;
end

% what to do
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3})
    return;
else
    action = lower(varargin{3}(:)');
end

% decide
switch (action)

    % apply TSM (with dialog if not provided)
    case 'applytsm'

        % provided
        if nargin < 5 || ~ischar(varargin{4}) || isempty(regexpi(varargin{4}(:)', '\.tsm')) || ...
           ~ischar(varargin{5}) || isempty(regexpi(varargin{5}(:)', '\.srf'))
            [tsmfile, tsmpath] = uigetfile({'*.tsm', ...
                'Triangle-based Sphere-Mapping files (*.tsm)'}, ...
                'Please choose the TSM file to apply...');
            if isequal(tsmfile, 0) || isequal(tsmpath, 0) || ...
                isempty(regexpi(tsmfile, '\.tsm$'))
                return;
            end
            if isempty(tsmpath)
                tsmpath = pwd;
            end
            tsm = [tsmpath filesep tsmfile];
            try
                tsmc = [];
                tsmc = xff(tsm);
                if ~isxff(tsmc, 'tsm')
                    error('NOT_A_TSM');
                end
                tsmv = tsmc.NrOfTargetVertices;
                tsmc.ClearObject;
                if any(tsmv == (2 + 10 .* (4 .^ (0:7))))
                    sph = '';
                else
                    [sphfile, sphpath] = uigetfile({'.srf', ...
                        'Spherical Surface file (*.srf)'}, ...
                        'Please select the matching surface file...');
                    if isequal(sphfile, 0) || isequal(sphpath, 0) || ...
                        isempty(regexpi(sphfile, '\.srf$'))
                        return;
                    end
                    if isempty(sphpath)
                        sphpath = pwd;
                    end
                    sph = [sphpath filesep sphfile];
                end
            catch ne_eo;
                clearxffobjects({tsmc});
                ne_gcfg.c.lasterr = ne_eo;
                return;
            end
        else
            tsm = varargin{4}(:)';
            sph = varargin{5}(:)';
        end
        ne_srf_tools(0, 0, 'maptoico', tsm, sph);

    % back-project to VMR
    case {'backproject'}

        % only valid for VMRs
        if nargin < 4 || ...
            numel(varargin{4}) ~= 1 || ...
           ~isxff(varargin{4}, 'vmr')
            slvar = cc.SliceVar;
        else
            slvar = varargin{4};
        end
        if ~isxff(slvar, 'vmr') || ...
            slvar.VoxResX ~= slvar.VoxResY || ...
            slvar.VoxResX ~= slvar.VoxResZ || ...
           ~any([0.5, 1] == slvar.VoxResX)
            uiwait(warndlg('Back-projection only valid for 1mm/0.5mm ISO-VMRs.', ...
                'NeuroElf - info', 'modal'));
            return;
        end

        % request a few things
        if nargin < 5 || ...
           ~iscell(varargin{5}) || ...
            numel(varargin{5}) ~= 3
            bpconf = inputdlg({'Target color code', ...
                'Triangle oversampling factor', 'Sampling along normal vectors'}, ...
                'NeuroElf - SRF backprojection configuration', 1, ...
                {'  235', '  5', '  -0.25:0.25:0.25'});
        else
            bpconf = varargin{5};
        end

        % assess config
        try
            tcode = str2double(bpconf{1});
            if isinf(tcode) || ...
                isnan(tcode) || ...
                tcode < 0 || ...
                tcode > 255
                return;
            else
                tcode = round(tcode);
            end
            tosfactor = str2double(bpconf{2});
            if isinf(tosfactor) || ...
                isnan(tosfactor) || ...
                tosfactor < 0 || ...
                tosfactor > 6
                return;
            else
                tosfactor = round(tosfactor);
            end
            nsrange = eval(['[' bpconf{3} ']']);
            if ~isa(nsrange, 'double') || ...
                isempty(nsrange) || ...
                any(isinf(nsrange(:)) | isnan(nsrange(:)))
                return;
            else
                nsfrom = nsrange(1);
                if numel(nsrange) > 1
                    nsto = nsrange(end);
                    nsstep = nsrange(2) - nsrange(1);
                else
                    nsto = nsfrom;
                    nsstep = 1;
                end
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % run backprojection
        bopt = struct('nfrom', nsfrom, 'nstep', nsstep, 'nto', nsto, ...
            'res', slvar.VoxResX, 'tcode', tcode, 'triovsmp', tosfactor, 'vmr', slvar);
        if tcode < 226
            bopt.fillmode = 'linear';
        end
        mfp = ch.MainFig.Pointer;
        ch.MainFig.Pointer = 'watch';
        drawnow;
        try
            srf.BackToVMR(bopt);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ch.MainFig.Pointer = mfp;

        % switch to first page and update
        ne_showpage(0, 0, 1);
        drawnow;

    % back-project to new VMR
    case {'backtovmr'}

        % settings
        bpconf = inputdlg({'VMR resolution (mm, 1 or 0.5)', 'Target color code', ...
            'Triangle oversampling factor', 'Sampling along normal vectors'}, ...
            'NeuroElf - SRF backprojection configuration', 1, ...
            {'  1', '  235', '  5', '  -0.25:0.25:0.25'});
        if ~iscell(bpconf) || numel(bpconf) ~= 4 || isempty(bpconf{1}) || ...
           ~ischar(bpconf{1}) || isnan(str2double(ddeblank(bpconf{1})))
            return;
        end
        vmrres = str2double(ddeblank(bpconf{1}));
        if vmrres ~= 0.5
            vmrres = 1;
        end

        % generate VMR
        try
            vmr = [];
            oldres = ne_gcfg.c.ini.MainFig.NewVMRRes;
            ne_gcfg.c.ini.MainFig.NewVMRRes = vmrres;
            vmr = ne_newfile(0, 0, 'vmr');
            ne_gcfg.c.ini.MainFig.NewVMRRes = oldres;
            ne_srf_tools(0, 0, 'backproject', vmr, bpconf(2:4));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            ne_gcfg.c.ini.MainFig.NewVMRRes = oldres;
            if isxff(vmr, true)
                vmr.ClearObject;
            end
            return;
        end

    % cancel morph
    case {'cancelmorph'}

        % is morphing?
        if isfield(srfh, 'Morphing') && islogical(srfh.Morphing) && ...
            numel(srfh.Morphing) == 1 && srfh.Morphing

            % cancel
            srf.SetHandle('CancelMorph', true);

            % update
            ne_setsurfpos(0, 0, 'upshape');
        end

    % cloning
    case {'clone'}
        
        % object
        if nargin > 3 && ...
            numel(varargin{4}) == 1 && ...
            isxff(varargin{4}, 'srf')
            srf = varargin{4};
        end
        
        % copy
        srfcpy = srf.CopyObject;

        % open copy
        ne_openfile(0, 0, srfcpy, true);

        % get corresponding variable name
        srfname = cv_varlist(0, 0, srfcpy);

        % set as tool-loaded
        if ~isempty(srfname)
            ne_gcfg.wc.(srfname) = true;
        end

        % improve name a little
        [fp, fn] = fileparts(srf.FilenameOnDisk(true));
        if isempty(fn)
            fn = 'untitled';
        end
        fe = '.srf';
        fn = sprintf('<copy of %s%s>', fn, fe);

        % return
        srf = srfcpy;

        % find in controls
        ud = ch.SurfVar.UserData;
        udf = false;
        for udc = 1:size(ud, 1)
            if ud{udc, 4} == srfcpy
                udf = true;
                break;
            end
        end
        if udf
            ch.SurfVar.String{udc} = fn;
        end
        ud = ch.Scenery.UserData;
        udf = false;
        for udc = 1:size(ud, 1)
            if ud{udc, 4} == srfcpy
                udf = true;
                break;
            end
        end
        if udf
            srfh = handles(srf);
            if isempty(srfh.VertexMorphMeshes)
                ch.Scenery.String{udc} = ...
                    sprintf('%s (NrOfVertices: %d)', fn, srfcpy.NrOfVertices);
            else
                ch.Scenery.String{udc} = ...
                    sprintf('%s (NrOfVertices: %d, NrOfMorphStates: %d)', ...
                    fn, srfcpy.NrOfVertices, size(srfh.VertexMorphMeshes, 1) + 1);
            end
        end

    % close selected files
    case {'closefiles'}

        % which files
        fi = ch.Scenery.Value;
        fi = sort(fi(:));
        fu = ch.Scenery.UserData;
        if isempty(fi) || isempty(fu)
            return;
        end

        % which to show after
        ti = unique(max(1, min(fi - 1)));

        % delete
        for fc = numel(fi):-1:1
            try

                % first get surface object
                srf = fu{fi(fc), 4};

                % then check in handles
                srfh = handles(srf);

                % whether stats are set and valid
                if isfield(srfh, 'Stats') && ...
                    iscell(srfh.Stats) && ...
                    numel(srfh.Stats) == 2 && ...
                    isxff(srfh.Stats{1}, {'mtc', 'smp'})

                    % and close that file as well
                    ne_closefile(0, 0, srfh.Stats{1});
                end

                % then remove the corresponding file
                ne_closefile(0, 0, srf);
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end

        % set new index
        if ~isempty(ch.Scenery.UserData)
            ch.Scenery.Value = ti;

            % and update
            ne_sceneryselect;
            ch.Scenery.Visible = 'off';
            drawnow;
            ch.Scenery.ListboxTop = ...
                max(1, min(size(ch.Scenery.UserData, 1), ch.Scenery.ListboxTop));
            ch.Scenery.Visible = 'on';
            drawnow;
        else
            ch.Scenery.Value = [];
        end


    % cluster-table based SRF
    case {'clustersrf'}

        % check VOI is good
        if nargin < 4 || ...
            numel(varargin{4}) ~= 1 || ...
           ~isxff(varargin{4}, 'voi') || ...
            isempty(varargin{4}.VOI)
            voi = ne_gcfg.voi;
            if ~isxff(voi, 'voi') || ...
                isempty(voi.VOI)
                return;
            end
        else
            voi = varargin{4};
        end
        if nargin < 5 || ...
           ~isa(varargin{5}, 'double') || ...
            isempty(varargin{5}) || ...
            any(isinf(varargin{5}(:)) | isnan(varargin{5}(:)) | varargin{5}(:) < 1)
            voii = 1:numel(voi.VOI);
        else
            voii = unique(round(min(numel(voi.VOI), varargin{5}(:)')));
        end
        vois = voi.VOI(voii);

        % block further access
        ne_gcfg.c.blockcb{end+1} = 'ne_srf_tools';

        % prepare progress bar
        tcnum = 0;
        for vc = 1:numel(vois)
            if ~isempty(regexpi(vois(vc).Name, 'SC\d+_\d+')) || ...
                (isfield(vois, 'IsLocalMax') && ...
                 islogical(vois(vc).IsLocalMax) && ...
                 numel(vois(vc).IsLocalMax) == 1 && ...
                 vois(vc).IsLocalMax)
                continue;
            end
            tcnum = tcnum + 1;
        end
        cprog = ne_progress(0, 0, ...
            {true, 0, sprintf('Creating SRF from %d clusters...', tcnum)});

        % check resolution
        if isfield(voi.RunTimeVars, 'ClusterTableResXYZ')
            res = ceil(max(voi.RunTimeVars.ClusterTableResXYZ));
        else
            res = 1;
        end

        % create VMR for reco
        try
            slvar = [];
            srf = [];
            csrf = [];
            nsrf = [];
            slvar = xff('new:vmr');
            cnum = 1;
            for vc = 1:numel(vois)
                if ~isempty(regexpi(vois(vc).Name, 'SC\d+_\d+')) || ...
                    (isfield(vois, 'IsLocalMax') && ...
                     islogical(vois(vc).IsLocalMax) && ...
                     numel(vois(vc).IsLocalMax) == 1 && ...
                     vois(vc).IsLocalMax)
                    continue;
                end
                if res == 1
                    slvar.VMRData(bvcoordconv(vois(vc).Voxels, 'tal2bvx', slvar.BoundingBox)) = 240;
                else
                    voxi = bvcoordconv(vois(vc).Voxels, 'tal2bvc', slvar.BoundingBox);
                    from = ceil(-0.5 * res);
                    to = ceil(0.5 * (res - 1));
                    szvmr = size(slvar.VMRData);
                    for xc = from:to
                        for yc = from:to
                            for zc = from:to
                                voxt = max(1, [min(szvmr(1), voxi(:, 1) + xc), ...
                                    min(szvmr(2), voxi(:, 2) + yc), ...
                                    min(szvmr(3), voxi(:, 3) + zc)]);
                                slvar.VMRData(sub2ind(szvmr, voxt(:, 1), voxt(:, 2), voxt(:, 3))) = 240;
                            end
                        end
                    end
                end

                % reconstruct (without handle warnings)
                csrf = slvar.DBReco(struct('tps', 4, 'warn', false));
                if csrf.NrOfVertices < 3
                    error( ...
                        'neuroelf:BadCluster', ...
                        'Error creating surface for cluster %d.', ...
                        vc ...
                    );
                end

                % set color
                csrf.VertexColor = ones(size(csrf.VertexCoordinate, 1), 1) * [NaN, vois(vc).Color];

                % either copy for first cluster
                if cnum == 1
                    srf = csrf;

                % or combine and delete particles
                else
                    nsrf = srf.Combine(csrf);
                    srf.ClearObject;
                    csrf.ClearObject;
                    srf = nsrf;
                end
                ch.Progress.Progress(cnum / tcnum);
                cnum = cnum + 1;

                % prepare for next cluster
                if vc < numel(vois)
                    slvar.VMRData(:) = 0;
                end
            end

            % clear VMR
            slvar.ClearObject;

            % then open surface
            ne_openfile(0, 0, srf, true);

        % deal with issues
        catch ne_eo;
            ne_progress(0, 0, cprog);
            if isxff(slvar, true)
                slvar.ClearObject;
            end
            if isxff(srf, true)
                srf.ClearObject;
            end
            if isxff(csrf, true)
                csrf.ClearObject;
            end
            if isxff(ncsrf, true)
                nsrf.ClearObject;
            end
            uiwait(warndlg(ne_eo.message, 'NeuroElf - warning', 'modal'));
            ne_gcfg.c.blockcb(strcmpi(ne_gcfg.c.blockcb, 'ne_srf_tools')) = [];
            return;
        end

        % patch name
        fname = sprintf('<clusters from %d VOIs.srf>', tcnum);
        ud = ch.SurfVar.UserData;
        udf = false;
        for udc = 1:size(ud, 1)
            if ud{udc, 4} == srf
                udf = true;
                break;
            end
        end
        if udf
            ch.SurfVar.String{udc} = fname;
        end
        ud = ch.Scenery.UserData;
        udf = false;
        for udc = 1:size(ud, 1)
            if ud{udc, 4} == srf
                udf = true;
                break;
            end
        end
        if udf
            ch.Scenery.String{udc} = ...
                sprintf('%s (NrOfVertices: %d)', fname, srf.NrOfVertices);
        end
        ne_progress(0, 0, cprog);
        ne_gcfg.c.blockcb(strcmpi(ne_gcfg.c.blockcb, 'ne_srf_tools')) = [];

    % create sphere
    case {'createsphere'}

        % inquire about radius, center, and triangle degree, and colors?
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 5
            srfcfg = inputdlg({'Sphere radius (mm)', 'Center (TAL coordinate)', ...
                'Triangulation degree (5 = 20480, 6 = 81920, 7 = 327680)', ...
                'Color 1', 'Color 2'}, 'NeuroElf - create sphere config', 1, ...
                {'  120', '  [0, 0, 0]', '  6', '  [85, 173, 250]', '  [26, 61, 85]'});
            if ~iscell(srfcfg) || ...
                numel(srfcfg) ~= 5
                return;
            end
        else
            srfcfg = varargin{4};
        end

        % check config
        try
            radius = str2double(srfcfg{1});
            if ~isa(radius, 'double') || ...
                numel(radius) ~= 1 || ...
                isinf(radius) || ...
                isnan(radius) || ...
                radius < 1 || ...
                radius > 256
                return;
            end
            center = 128 - lsqueeze(eval(srfcfg{2}))';
            center = center(1, [2, 3, 1]);
            if ~isa(center, 'double') || ...
                any(isinf(center) | isnan(center) | center < 0 | center > 256)
                return;
            end
            tfactor = str2double(srfcfg{3});
            if ~isa(tfactor, 'double') || ...
                numel(tfactor) ~= 1 || ...
                isinf(tfactor) || ...
                isnan(tfactor) || ...
                tfactor < 2 || ...
                tfactor > 7
                return;
            else
                tfactor = round(tfactor);
            end
            color1 = lsqueeze(eval(srfcfg{4}))';
            color1 = color1(1, [1, 2, 3]);
            if ~isa(color1, 'double') || ...
                any(isinf(color1) | isnan(color1) | color1 < 0 | color1 > 255)
                return;
            else
                color1 = (1 / 255) .* [color1, 255];
            end
            color2 = lsqueeze(eval(srfcfg{5}))';
            color2 = color2(1, [1, 2, 3]);
            if ~isa(color2, 'double') || ...
                any(isinf(color2) | isnan(color2) | color2 < 0 | color2 > 255)
                return;
            else
                color2 = (1 / 255) .* [color2, 255];
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % create sphere
        srf = spheresrf(radius, tfactor);
        center = center - srf.MeshCenter;
        srf.VertexCoordinate = srf.VertexCoordinate + ones(srf.NrOfVertices, 1) * center;
        srf.ConvexRGBA = color1;
        srf.ConcaveRGBA = color2;

        % open
        ne_openfile(0, 0, srf, true);

        % patch name
        fname = sprintf('<%d-triangle-sphere.srf>', size(srf.TriangleVertex, 1));
        ud = ch.SurfVar.UserData;
        udf = false;
        for udc = 1:size(ud, 1)
            if ud{udc, 4} == srf
                udf = true;
                break;
            end
        end
        if udf
            ch.SurfVar.String{udc} = fname;
        end
        ud = ch.Scenery.UserData;
        udf = false;
        for udc = 1:size(ud, 1)
            if ud{udc, 4} == srf
                udf = true;
                break;
            end
        end
        if udf
            ch.Scenery.String{udc} = ...
                sprintf('%s (NrOfVertices: %d)', fname, srf.NrOfVertices);
        end

    % curvature map
    case {'curvsmp'}

        % call and open
        try
            smp = [];
            smp = srf.CurvatureMap(struct('medrem', true));
            ne_openfile(0, 0, smp, true);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            if isxff(smp, true)
                smp.ClearObject;
            end
        end

        % reassign for output
        srf = smp;

    % curvature map
    case {'denssmp'}

        % call and open
        try
            smp = [];
            smp = srf.DensityMap(struct('dist', true, 'mindist', true, 'stddist', true));
            ne_openfile(0, 0, smp, true);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            if isxff(smp, true)
                smp.ClearObject;
            end
        end

        % reassign for output
        srf = smp;

    % find intensity
    case {'findintensity'}

        % only valid with SliceVar
        slvar = ne_gcfg.fcfg.SliceVar;
        if ~isxff(slvar, {'hdr', 'head', 'vmr'})
            return;
        end

        % request for intensity
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 4
            v = slvar.GetVolume(1);
            vi = round(0.6666 * mean(v(v > 0)));
            ivalues = inputdlg({'Target intensity', 'Pulling force', 'Iterations', ...
                'Default direction, "in" or "out"'}, ...
                'NeuroElf - mesh-to-intensity morph', 1, {sprintf('  %d', vi), ...
                '  0.75', '  75', '  in'});
        else
            ivalues = varargin{4};
        end
        if ~iscell(ivalues) || ...
            numel(ivalues) ~= 4
            return;
        end

        % resolve inputs
        try
            ivalue = str2double(ivalues{1});
            force = str2double(ivalues{2});
            niter = str2double(ivalues{3});
            defdir = ddeblank(lower(ivalues{4}(:)'));
            if isinf(ivalue) || ...
                isnan(ivalue) || ...
                ivalue <= 0 || ...
                isinf(force) || ...
                isnan(force) || ...
                force <= 0 || ...
                force > 1 || ...
                isinf(niter) || ...
                isnan(niter) || ...
                niter < 1 || ...
               ~any(strcmp(defdir, {'in', 'out'}))
                return;
            end
            niter = min(1000, round(niter));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % mesh available already
        if isxff(srf, 'srf') && ...
           (nargin < 5 || ...
            ~ischar(varargin{5}) || ...
            isempty(varargin{5}) || ...
            ~strcmpi(varargin{5}(:)', 'nospherecheck'))

            % test whether it's a sphere
            r = srf.VertexCoordinate;
            rm = mean(r, 1);
            r = r - repmat(rm, size(r, 1), 1);
            r = sqrt(sum(r .* r, 2));

            % too much variability?
            if std(r) > (0.01 * mean(r))

                % don't use it
                srf = [];
            end
        end

        % create sphere
        if ~isxff(srf, 'srf')

            % sample slice var to estimate center
            slvals = slvar.SampleTalBox(struct( ...
                'BBox',   [-127.5, -127.5, -127.5; 127.5, 127.5, 127.5], ...
                'ResXYZ', [2, 2, 2]));

            % compute histogram
            slh = histcount(slvals, 0, 2 * ivalue, 1);

            % don't count 0 values
            slh(1) = 0;

            % scale to 1
            slhs = cumsum(slh);
            slhs = slhs ./ slhs(end);

            % middle value is what we are looking for
            ivalp = slhs(1 + round(ivalue));

            % add some portion to the histogram
            ival1 = findfirst(slhs >= ivalp - 0.1);
            ival2 = findfirst(slhs <= ivalp + 0.1, -1);

            % then mask
            slvals = (slvals >= ival1 & slvals <= ival2);

            % convert coordinates
            [slx, sly, slz] = ind2sub(size(slvals), find(slvals(:)));

            % estimate required center
            rcenter = mean([slx(:), sly(:), slz(:)]);
            slx = slx(:) - rcenter(1);
            sly = sly(:) - rcenter(2);
            slz = slz(:) - rcenter(3);
            sll = sqrt(slx .* slx + sly .* sly + slz .* slz);
            rcenter = 2 .* rcenter - 129;
            rcenter = sprintf('[%.1f, %.1f, %.1f]', rcenter);
            rradius = sprintf('%d', ceil(2 * mean(sll) + 3 * std(sll)));

            tfac = sprintf('%d', ci.HeadMeshTriFactor);
            srf = ne_srf_tools(0, 0, 'createsphere', ...
                {rradius, rcenter, tfac, '[248, 240, 216]', '[252, 248, 236]'});
        end

        % make sure it's selected
        ud = ch.Scenery.UserData;
        udf = false;
        for udc = 1:size(ud, 1)
            if ud{udc, 4} == srf
                udf = true;
                break;
            end
        end
        if udf
            ch.Scenery.Value = udc;
            ne_sceneryselect;
        end

        % then block this
        ne_gcfg.c.blockcb{end+1} = 'ne_srf_tools';

        % then morph the sphere
        cprog = ne_progress(0, 0, {true, 0, 'Finding intensity...'});
        try
            srf.MorphTo(slvar, ivalue, niter, struct( ...
                'default', defdir, 'force', force, 'pbar', ch.Progress, ...
                'snmat', true));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        sv = ch.Scenery.Value;
        if numel(sv) == 1 && sv == udf
            ne_setsurfpos(0, 0, 'upshape');
        end

        % unblock
        ne_progress(0, 0, cprog);
        ne_gcfg.c.blockcb(strcmpi(ne_gcfg.c.blockcb, 'ne_srf_tools')) = [];

    % inflate
    case {'inflate'}

        % ask for number of iterations and force
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 2 || ...
           ~ischar(varargin{4}{1}) || ...
           ~ischar(varargin{4}{2}) || ...
            isempty(varargin{4}{1}) || ...
            isempty(varargin{4}{2})
            uio = inputdlg( ...
                {'Number of inflation iterations:', 'Smoothing force ([0 ... 1]):'}, ...
                'NeuroElf GUI - user input', 1, {'  5000', '  0.5'});
            if numel(uio) ~= 2 || ...
               ~iscell(uio) || ...
               ~ischar(uio{1}) || ...
               ~ischar(uio{2})
                return;
            end
        else
            uio = varargin{4};
        end
        uio{1}(uio{1} == ' ') = [];
        uio{2}(uio{2} == ' ') = [];
        try
            uio{1} = str2double(uio{1});
            uio{2} = str2double(uio{2});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        if isinf(uio{1}) || ...
            isnan(uio{1}) || ...
            uio{1} <= 0 || ...
            isinf(uio{2}) || ...
            isnan(uio{2}) || ...
            uio{2} <= 0 || ...
            uio{2} > 1
            return;
        end

        % echo
        if ne_gcfg.c.echo
            ne_echo('srf', 'Inflate', uio{:});
        end

        % make call
        cprog = ne_progress(0, 0, {true, 0, 'Inflating surface...'});
        ne_gcfg.c.blockcb{end+1} = 'ne_srf_tools';
        try
            srf.Inflate(uio{:}, struct('pbar', ch.Progress));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ne_gcfg.c.blockcb(strcmpi(ne_gcfg.c.blockcb, 'ne_srf_tools')) = [];
        ne_progress(0, 0, cprog);

        % update
        ne_setsurfpos(0, 0, 'upshape');

    % load scenery
    case {'loadscenery'}

        % request filename
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{4}) || ...
            exist(varargin{4}(:)', 'file') ~= 2
            [scfile, scpath] = uigetfile( ...
                {'*.mat', 'NeuroElf Scenery MAT files (*.mat)'}, ...
                'Reload NeuroElf Scenery file...');
            if isequal(scfile, 0) || ...
                isequal(scpath, 0)
                return;
            end
            if isempty(scpath)
                scpath = pwd;
            end
            scfile = [scpath, filesep, scfile];
        else
            scfile = varargin{4}(:)';
        end

        % get current scenery content
        cscu = ch.Scenery.UserData;

        % load file
        try
            scdata = load(scfile);
            if ~isstruct(scdata) || ...
                numel(scdata) ~= 1 || ...
               ~isfield(scdata, 'scdata') || ...
                numel(fieldnames(scdata)) ~= 1 || ...
               ~iscell(scdata.scdata) || ...
                numel(scdata.scdata) ~= 4 || ...
               ~isa(scdata.scdata{1}, 'double') || ...
               ~iscell(scdata.scdata{2}) || ...
               ~iscell(scdata.scdata{3}) || ...
               ~isequal(size(scdata.scdata{2}), size(scdata.scdata{3})) || ...
                numel(scdata.scdata{4}) ~= 1 || ...
               ~isstruct(scdata.scdata{4}) || ...
               ~isfield(scdata.scdata{4}, 'anglex') || ...
               ~isfield(scdata.scdata{4}, 'angley') || ...
               ~isfield(scdata.scdata{4}, 'time') || ...
               ~isfield(scdata.scdata{4}, 'trans') || ...
               ~isfield(scdata.scdata{4}, 'zoom')
                error( ...
                    'neuroelf:BadFileContent', ...
                    'Invalid file content for scenery file.' ...
                );
            end
            scidx = false(numel(scdata.scdata{2}) + size(cscu, 1));
            scsrf = scdata.scdata{2};
            scsrfh = scdata.scdata{3};
        catch ne_eo;
            uiwait(warndlg(['Invalid Scenery file: ' ne_eo.message], ...
                'NeuroElf - warning', 'modal'));
        end

        % iterate over surfaces
        remsurf = false(numel(scsrf), 1);
        for oc = 1:numel(scsrf)

            % start with re-loading stats
            if isfield(scsrfh{oc}, 'Stats') && ...
                iscell(scsrfh{oc}.Stats) && ...
               ~isempty(scsrfh{oc}.Stats) && ...
                isstruct(scsrfh{oc}.Stats{1}) && ...
                isfield(scsrfh{oc}.Stats{1}, 'FILE') && ...
                ischar(scsrfh{oc}.Stats{1}.FILE) && ...
                isfield(scsrfh{oc}.Stats{1}, 'CONT') && ...
                isstruct(scsrfh{oc}.Stats{1}.CONT) && ...
                numel(scsrfh{oc}.Stats{1}.CONT) == 1 && ...
               ((isfield(scsrfh{oc}.Stats{1}.CONT, 'Map') && ...
                 ~isempty(scsrfh{oc}.Stats{1}.CONT.Map)) || ...
                (isfield(scsrfh{oc}.Stats{1}.CONT, 'MTCData') && ...
                 ~isempty(scsrfh{oc}.Stats{1}.CONT.MTCData)))

                % re-generate
                stfile = scsrfh{oc}.Stats{1}.FILE;
                stcont = scsrfh{oc}.Stats{1}.CONT;

                % cannot re-create?
                fromfile = ~isempty(stfile);
                if fromfile
                    if exist(stfile, 'file') ~= 2
                        scsrfh{oc}.Stats{1} = struct('Filetype', 'NONE');
                        scsrfh{oc}.Stats{2} = [];
                        continue;
                    end

                    % try loading
                    try
                        oldsmp = xff(stfile);
                        oldcont = getcont(oldsmp);
                    catch ne_eo;
                        ne_gcfg.c.lasterr = ne_eo;
                        scsrfh{oc}.Stats{1} = struct('Filetype', 'NONE');
                        scsrfh{oc}.Stats{2} = [];
                        continue;
                    end
                end

                % create or reload SMP
                if fromfile
                    newsmp = oldsmp;
                else
                    if isfield(stcont, 'MTCData')
                        newsmp = xff('new:mtc');
                    else
                        newsmp = xff('new:smp');
                    end
                end

                % copy content
                smpflds = fieldnames(stcont);
                for fc = 1:numel(smpflds)
                    if ~strcmp(smpflds{fc}, 'RunTimeVars')
                        if fromfile && ...
                            ischar(stcont.(smpflds{fc})) && ...
                            isequal(size(stcont.(smpflds{fc})), [1, 5]) && ...
                            strcmp(stcont.(smpflds{fc}), 'CCFRF') && ...
                            isfield(oldcont, smpflds{fc});
                            newsmp.(smpflds{fc}) = oldcont.(smpflds{fc});
                        else
                            newsmp.(smpflds{fc}) = stcont.(smpflds{fc});
                        end
                    else
                        if fromfile
                            rtvflds = fieldnames(oldcont.RunTimeVars);
                            for rfc = 1:numel(rtvflds)
                                if ~strcmp(rtvflds{rfc}, 'xffID')
                                    newsmp.RunTimeVars.(rtvflds{rfc}) = ...
                                        oldcont.RunTimeVars.(rtvflds{rfc});
                                end
                            end
                        end
                        if isstruct(stcont.RunTimeVars)
                            rtvflds = fieldnames(stcont.RunTimeVars);
                            for rfc = 1:numel(rtvflds)
                                if ~strcmp(rtvflds{rfc}, 'xffID')
                                    newsmp.RunTimeVars.(rtvflds{rfc}) = ...
                                        stcont.RunTimeVars.(rtvflds{rfc});
                                end
                            end
                        end
                    end
                end

                % open file
                ne_openfile(0, 0, newsmp, true);

                % then replace in Stats
                scsrfh{oc}.Stats{1} = newsmp;

            % no stats
            else
                scsrfh{oc}.Stats = {struct('Filetype', 'NONE'), []};
            end

            % SRF ilename given
            srffile = scsrf{oc}.FILE;
            srfcont = scsrf{oc}.CONT;
            fromfile = ~isempty(srffile);

            % try loading
            try
                if fromfile
                    newsrf = xff(srffile);
                    if ~isxff(newsrf, 'srf')
                        if isxff(newsrf, 'true')
                            newsrf.ClearObject;
                        end
                        error( ...
                            'neuroelf:BadRefObject', ...
                            'Invalid surface object file.' ...
                        );
                    end
                else
                    newsrf = emptysrf;
                end
                oldcont = getcont(newsrf);

                % update content
                srfflds = fieldnames(srfcont);
                for fc = 1:numel(srfflds)
                    if ~strcmp(srfflds{fc}, 'RunTimeVars')
                        if ~fromfile || ...
                           ~ischar(srfcont.(srfflds{fc})) || ...
                           ~isequal(size(srfcont.(srfflds{fc})), [1, 5]) || ...
                           ~strcmp(srfcont.(srfflds{fc}), 'CCFRF') || ...
                           ~isfield(oldcont, srfflds{fc});
                            newsrf.(srfflds{fc}) = srfcont.(srfflds{fc});
                        end
                    else
                        if isstruct(srfcont.RunTimeVars)
                            rtvflds = fieldnames(srfcont.RunTimeVars);
                            for rfc = 1:numel(rtvflds)
                                if ~strcmp(rtvflds{rfc}, 'xffID')
                                    newsrf.RunTimeVars.(rtvflds{rfc}) = ...
                                        srfcont.RunTimeVars.(rtvflds{rfc});
                                end
                            end
                        end
                    end
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                if isxff(scsrfh{oc}.Stats{1}, true)
                    scsrfh{oc}.Stats{1}.ClearObject;
                end
                remsurf(oc) = true;
                continue;
            end

            % then set handles
            hflds = fieldnames(scsrfh{oc});
            for hc = 1:numel(hflds)
                newsrf.SetHandle(hflds{hc}, scsrfh{oc}.(hflds{hc}));
            end

            % show surface
            ne_openfile(0, 0, newsrf, true);

            % we got here, activate in index
            if any(scdata.scdata{1} == oc)
                scidx(oc + size(cscu, 1)) = true;
            end
        end

        % update index
        scidx(find(remsurf) + size(cscu, 1)) = [];
        scval = find(scidx);
        if isempty(scval)
            ch.Scenery.ListboxTop = 1;
        else
            ch.Scenery.ListboxTop = scval(1);
        end
        ch.Scenery.Value = scval;

        % then position
        ne_gcfg.fcfg.srfcfg = scdata.scdata{4};
        ne_sceneryselect;
        drawnow;

    % load subcortical areas
    case {'loadsubcort'}

        % filename pattern
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
           ~strcmpi(varargin{4}(:)', 'icbmnorm')
            icbmnorm = false;
            spat = 'colin_subcort_*.srf';
        else
            icbmnorm = true;
            spat = 'colin_subcort_*_ICBMnorm.srf';
        end

        % find files
        srfs = findfiles(neuroelf_path('colin'), spat, 'depth=1');
        if ~icbmnorm
            srfs(~cellfun('isempty', regexpi(srfs, 'icbmnorm'))) = [];
        end
        srfs(~cellfun('isempty', regexpi(srfs, '_combined'))) = [];

        % reorder?
        if mod(numel(srfs), 2) == 1
            nsrfs = srfs(:);
            nsrfs(2:2:end) = srfs(2:ceil(numel(srfs)/2));
            nsrfs(3:2:end) = srfs(ceil(numel(srfs)/2)+1:end);
            srfs = nsrfs;
        end

        % load SRFs
        for fc = 1:numel(srfs)
            ne_openfile(0, 0, srfs{fc});
        end

    % map to icosahedron
    case {'maptoico'}

        % make sure SRF is saved
        srfc = getcont(srf);
        srff = srf.FilenameOnDisk;
        if isempty(srff)
            vans = questdlg('Surface not yet saved. Would you like to save it now?', ...
                'NeuroElf - user request');
            if isempty(vans) || ...
                lower(vans(1)) ~= 'y'
                return;
            end
            srf.SaveAs;
            srff = srf.FilenameOnDisk;
            if isempty(srff)
                return;
            end
        end
        
        % TSM given
        if nargin > 3 && ischar(varargin{4}) && ~isempty(regexpi(varargin{4}(:)', '\.tsm$'))
            try
                tsm = [];
                sph = [];
                tsph = [];
                tsm = xff(varargin{4}(:)');
                if ~isxff(tsm, 'tsm')
                    error('NO_TSM_FILE');
                end
                f = (2 + 10 .* (4 .^ (0:7))) ./ tsm.NrOfTargetVertices;
                if ~any(f == 1)
                    if nargin > 4 && ischar(varargin{5}) && ...
                       ~isempty(regexpi(varargin{5}(:)', '\.srf$'))
                        sph = xff(varargin{5}(:)');
                    else
                        sph = xff('*.srf', 'Please select matching (sphere) surface...');
                    end
                else
                    sph = spheresrf(100, find(f == 1) - 1);
                end
                if ~isxff(sph, 'srf')
                    clearxffobjects([tsm, sph, tsph]);
                    uiwait(warndlg('No matching surface available.', 'NeuroElf - info', 'modal'));
                    return;
                end
                tsph = srf.ApplyTSM(tsm, sph);
                tsm.ClearObject;
                sph.ClearObject;
                ne_openfile(0, 0, tsph, true);
                return;
            catch ne_eo;
                clearxffobjects([tsm, sph, tsph]);
                uiwait(warndlg(['Error applying TSM: ' ne_eo.message], ...
                    'NeuroElf - info', 'modal'));
                return;
            end
        end
        
        % already mapped?
        [srfp, srff] = fileparts(srff);
        if ci.ReuseSphereMapping
            tsms = findfiles(srfp, [srff '_to*kSPH.tsm'], 'depth=1');
        else
            tsms = {};
        end

        % TSM
        sph = {[]};
        if ~isempty(tsms)
            try
                tsm = {[]};
                tsm{1} = xff(tsms{maxpos(cellfun('prodofsize', tsms))});
                sph{1} = spheresrf(100, round(log2((tsm{1}.NrOfTargetVertices - 2) / 10) / 2));
                srf = srf.ApplyTSM(tsm{1}, sph{1});
                clearxffobjects([tsm, sph]);
                ne_openfile(0, 0, srf, true);
                return;
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                clearxffobjects([tsm, sph]);
            end
        elseif ci.ReuseSphereMapping
            ssms = findfiles(srfp, [srff '_to*kSPH.ssm'], 'depth=1');
            if ~isempty(ssms)
                try
                    ssm = {[]};
                    ssm{1} = xff(ssms{maxpos(cellfun('prodofsize', ssms))});
                    sph{1} = spheresrf(100, round(log2((ssm{1}.NrOfTargetVertices - 2) / 10) / 2));
                    srf = srf.ApplySSM(ssm{1}, sph{1});
                    clearxffobjects([ssm, sph]);
                    ne_openfile(0, 0, srf, true);
                    return;
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                    clearxffobjects([ssm, sph]);
                end
            end
        end

        % only valid if a reference SRF is set
        checksphere = true;
        srfals = srfc.AutoLinkedSRF;
        if isempty(srfals)
            if isempty(srff)
                return;
            end
            vans = questdlg('No reference surface defined. Use this one instead?', ...
                'NeuroElf - user request');
            if ~isequal(lower(vans), 'yes')
                return;
            end
            srfals = srf.FilenameOnDisk;
            checksphere = false;
        end

        % load original mesh
        try
            osrf = xff(srfals);
            if ~isxff(osrf, 'srf') || ...
               ~isequal(srfc.TriangleVertex, osrf.TriangleVertex)
                if isxff(osrf, true)
                    osrf.ClearObject;
                end
                error( ...
                    'neuroelf:GUI:NotValidLinkFile', ...
                    'Not a valid SRF linked.' ...
                );
            end
        catch ne_eo;
            uiwait(warndlg(ne_eo.message, 'NeuroElf - info', 'modal'));
            return;
        end
        
        % only valid with spherical SRF
        ocoords = srfc.VertexCoordinate;
        center = mean(ocoords);
        coords = ocoords - ones(size(ocoords, 1), 1) * center;
        lengths = sqrt(sum(coords .* coords, 2));
        mlength = mean(lengths);
        lengths = abs((1 / mlength) .* (lengths - mlength));
        if any(abs(lengths) > 0.01) && ...
            checksphere
            osrf.ClearObject;
            uiwait(warndlg('Surface not spherical (enough).', 'NeuroElf - info', 'modal'));
            return;
        end
        ncoords = numel(lengths);
        
        % config
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 3
            uio = inputdlg({ ...
                'Triangulation factor (5 = 20480, 6 = 81920, 7 = 327680)', ...
                'SSM or TSM', 'Also create inverse mapping (yes/no)'}, ...
                'NeuroElf - Icosahedron mapping configuration', 1, ...
                {'  6', '  tsm', '  yes'});
            if ~iscell(uio) || ...
                numel(uio) ~= 3
                osrf.ClearObject;
                return;
            end
        else
            uio = varargin{4};
        end
        try
            tfac = round(str2double(ddeblank(uio{1}(:)')));
            if numel(tfac) ~= 1 || ...
                isinf(tfac) || ...
                isnan(tfac) || ...
                tfac < 3 || ...
                tfac > 8
                osrf.ClearObject;
                return;
            end
            type = lower(ddeblank(uio{2}));
            if type(1) ~= 's' && ...
                type(1) ~= 't'
                osrf.ClearObject;
                return;
            end
            type = type(1);
            ialso = (lower(uio{3}(1)) == 'y');
        catch ne_eo;
            osrf.ClearObject;
            uiwait(warndlg(ne_eo.message, 'NeuroElf - error', 'modal'));
            return;
        end
        
        % create sphere first
        sph = spheresrf(ceil(mlength), tfac, [128, 128, 128]);
        
        % smooth sphere so the vertices are better distributed
        sph.Inflate(1000, 0.5, struct('show', false));
        
        % then create mappings
        ofln = osrf.FilenameOnDisk;
        if type == 's'
            ssm = sph.SphereToSphereMapCoord(srf);
            ssm.SaveAs(sprintf('%s_to%dkSPH.ssm', ...
                ofln(1:end-4), round(0.001 * ncoords)));
            if ialso
                issm = srf.SphereToSphereMapCoord(sph);
                issm.SaveAs(sprintf('%s_to%dkSPH_inverse.ssm', ...
                    ofln(1:end-4), round(0.001 * ncoords)));
                issm.ClearObject;
            end
            srf = osrf.ApplySSM(ssm, sph);
            ssm.ClearObject;
        else
            tsm = sph.SphereToSphereMapTri(srf, struct('fwithin', true));
            tsm.SaveAs(sprintf('%s_to%dkSPH.tsm', ...
                ofln(1:end-4), round(0.001 * ncoords)));
            if ialso
                itsm = srf.SphereToSphereMapTri(sph, struct('fwithin', true));
                itsm.SaveAs(sprintf('%s_to%dkSPH_inverse.ssm', ...
                    ofln(1:end-4), round(0.001 * ncoords)));
                itsm.ClearObject;
            end
            srf = osrf.ApplyTSM(tsm, sph);
            tsm.ClearObject;
        end
        
        % remove non-needed object
        osrf.ClearObject;
        sph.ClearObject;
        
        % then show output
        ne_openfile(0, 0, srf, true);

    % general morphing
    case {'morph'}

        % config not given
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 9 || ...
           ~all(cellfun(@ischar, varargin{4}(:)))
            uio = inputdlg({ ...
                'Number of morphing operations:', ...
                'General smoothing force (0 .. 1):', ...
                'Keep surface area constant (yes/no):', ...
                'Weight smoothing force by area associated with vertex (yes/no):', ...
                'Apply distortion-correction force (0 .. 6):', ...
                'Weight smoothing force by distance (''no'', ''log'', or ''sq''):', ...
                'Apply along-normal force (-1 .. 1):', ...
                'Ramp up along-normal force (yes/no)', ...
                'Apply to-sphere force (0 .. 1):'}, ...
                'NeuroElf - surface morphing config', 1, ...
                {'  5000', '  0.25', '  yes', '  no', '  0', '  no', '  0', '  yes', '  0'});
            if ~iscell(uio) || ...
                numel(uio) ~= 9 || ...
               ~all(cellfun(@ischar, uio)) || ...
                any(cellfun('isempty', ddeblank(uio)))
                return;
            end
        else
            uio = varargin{4};
        end

        % parse options
        try
            uio = ddeblank(uio(:));
            niter = str2double(uio{1});
            if isinf(niter) || ...
                isnan(niter) || ...
                niter < 1
                return;
            else
                niter = min(50000, round(niter));
            end
            force = str2double(uio{2});
            if isinf(force) || ...
                isnan(force) || ...
                force <= 0 || ...
                force > 1
                return;
            end
            areac = (lower(uio{3}(1)) == 'y');
            areaw = (lower(uio{4}(1)) == 'y');
            distc = str2double(uio{5});
            if isinf(distc) || ...
                isnan(distc) || ...
                distc < 0
                return;
            end
            distc = min(6, distc);
            distw = (lower(uio{6}(1)) ~= 'n');
            distwl = (lower(uio{6}(1)) == 'l');
            distwsq = (lower(uio{6}(1)) == 's');
            normf = str2double(uio{7});
            if isinf(normf) || ...
                isnan(normf) || ...
                normf < -1 || ...
                normf > 1
                return;
            end
            normramp = (lower(uio{8}(1)) == 'y');
            tosphere = str2double(uio{9});
            if isinf(tosphere) || ...
                isnan(tosphere) || ...
                tosphere < 0 || ...
                tosphere > 1
                return;
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % apply morphing
        ne_gcfg.c.blockcb{end+1} = 'ne_srf_tools';
        cprog = ne_progress(0, 0, {true, 0, 'Morphing surface...'});
        try
            srf.Morph(niter, force, 'smooth', struct( ...
                'areac',     areac, ...
                'areaw',     areaw, ...
                'distc',     distc, ...
                'distw',     distw, ...
                'distwl',    distwl, ...
                'distwsq',   distwsq, ...
                'norm',      normf, ...
                'normramp',  normramp, ...
                'sphere',    tosphere, ...
                'pbar',      ch.Progress));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ne_gcfg.c.blockcb(strcmpi(ne_gcfg.c.blockcb, 'ne_srf_tools')) = [];
        ne_progress(0, 0, cprog);

        % update
        ne_setsurfpos(0, 0, 'upshape');
        
    % save scenery
    case {'savescenery'}

        % request output filename
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{4})
            [scfile, scpath] = uiputfile( ...
                {'*.mat', 'NeuroElf Scenery MAT files (*.mat)'; ...
                 '*.ply', 'PoLYgon file format (*.ply)'}, ...
                'Save NeuroElf Scenery as...');
            if isequal(scfile, 0) || ...
                isequal(scpath, 0)
                return;
            end
            if isempty(scpath)
                scpath = pwd;
            end
            scfile = [scpath, filesep, scfile];
        else
            scfile = varargin{4}(:)';
        end

        % grap Scenery UserData
        scu = ch.Scenery.UserData;

        % prepare output
        scdata = cell(1, 4);
        scdata{1} = ch.Scenery.Value;
        scdata{2} = cell(size(scu, 1), 1);
        scdata{3} = scdata{2};
        scdata{4} = ne_gcfg.fcfg.srfcfg;

        % store file references instead of full content
        ismat = true;
        if (nargin < 5 || ...
            ~islogical(varargin{5}) || ...
            numel(varargin{5}) ~= 1) && ...
            strcmpi(scfile(end-3:end), '.mat')
            stref = questdlg('Store file references instead of content?', ...
                'NeuroElf - user input request', 'Yes', 'No', 'Cancel', 'Yes');
            if ~ischar(stref) || ...
               ~any(strcmpi(stref, {'no', 'yes'}))
                return;
            end
            stref = (lower(stref(1)) == 'y');
        elseif ~strcmpi(scfile(end-3:end), '.mat')
            stref = false;
            ismat = false;
        else
            stref = varargin{5};
        end

        % update
        ne_setsurfpos;
        mfp = ch.MainFig.Pointer;
        ch.MainFig.Pointer = 'watch';
        drawnow;

        % for each object in the scenery
        for oc = 1:size(scu, 1)

            % get objects contents
            scucont = getcont(scu{oc, 4});
            if stref

                % get fieldnames and filename
                scucfld = fieldnames(scucont);
                scufile = scu{oc, 4}.FilenameOnDisk;

                % if file exists
                if ~isempty(scufile) && ...
                    exist(scufile, 'file') > 0

                    % load file and get content
                    try
                        orgobj = xff(scufile);
                        orgcont = getcont(orgobj);
                        orgobj.ClearObject;

                    % with an error after all
                    catch ne_eo;
                        ne_gcfg.c.lasterr = ne_eo;
                        scufile = '';

                        % replace all content with error fields
                        orgcont = scucont;
                        for fc = 1:numel(scucfld)
                            orgcont.(scucfld{fc}) = 'ERROR';
                        end
                    end

                    % for each field
                    for fc = 1:numel(scucfld)

                        % check if content is the same
                        if isequaln(scucont.(scucfld{fc}), orgcont.(scucfld{fc}))

                            % and mark as to "copy content from ref file"
                            scucont.(scucfld{fc}) = 'CCFRF';
                        end
                    end

                % invalid file
                else
                    scufile = '';
                end

                % store
                scdata{2}{oc} = struct('FILE', scufile, 'CONT', scucont);

            % not interested in file reference anyway
            else
                scdata{2}{oc} = struct('FILE', '', 'CONT', scucont);
            end

            % get handles as well
            h = handles(scu{oc, 4});
            
            % add special handle for PLY format
            if ~ismat
                h.SurfTIORendered = h.SurfTIO.Rendered;
            end

            % then work on handles
            for fc = {'xff', 'CleanUp', 'ShownInGUI', 'SurfTIO', 'SUpdate', 'Surface', ...
                    'SourceObject', 'VertexCoordinateTal', 'VertexNormalTal'}
                if isfield(h, fc{1})
                    h = rmfield(h, fc{1});
                end
            end
            if isfield(h, 'Stats') && ...
                iscell(h.Stats) && ...
               ~isempty(h.Stats) && ...
                isxff(h.Stats{1}, true)

                % same as above for surfaces
                smpcont = getcont(h.Stats{1});
                if stref

                    % get fieldnames and filename
                    smpcfld = fieldnames(smpcont);
                    smpfile = h.Stats{1}.FilenameOnDisk;

                    % if file exists
                    if ~isempty(smpfile) && ...
                        exist(smpfile, 'file') > 0

                        % load file and get content
                        try
                            orgobj = xff(smpfile);
                            orgcont = getcont(orgobj);
                            orgobj.ClearObject;

                        % with an error after all
                        catch ne_eo;
                            ne_gcfg.c.lasterr = ne_eo;
                            smpfile = '';

                            % replace all content with error fields
                            orgcont = smpcont;
                            for fc = 1:numel(smpcfld)
                                orgcont.(smpcfld{fc}) = 'ERROR';
                            end
                        end

                        % for each field
                        for fc = 1:numel(smpcfld)

                            % check if content is the same
                            if isequaln(smpcont.(smpcfld{fc}), orgcont.(smpcfld{fc}))

                                % and mark as to "copy content from ref file"
                                smpcont.(smpcfld{fc}) = 'CCFRF';
                            end
                        end
                    end

                % no reference file
                else
                    smpfile = '';
                end

                % store
                h.Stats{1} = struct('FILE', smpfile, 'CONT', smpcont);
            end

            % and store handles
            scdata{3}{oc} = h;
        end

        % try saving the file
        try
            
            % MAT file
            fid = [];
            if ismat
                
                % save
                save(scfile, 'scdata');
                
                % then return early
                ch.MainFig.Pointer = mfp;
                return;
            end
            
            % open PLY file for writing
            fid = fopen(scfile, 'w');
            if fid < 1
                error( ...
                    'neuroelf:FileOpenError', ...
                    'Error opening Scenery (PLY) file for writing.' ...
                );
            end
            
            % count total number of vertices and triangles
            vcs = zeros(numel(scdata{2}), 1);
            tcs = vcs;
            for oc = 1:numel(scdata{2})
                vcs(oc) = size(scdata{2}{oc}.CONT.VertexCoordinate, 1);
                tcs(oc) = size(scdata{2}{oc}.CONT.TriangleVertex, 1);
            end
            
            % then compile data
            vdata = zeros(sum(vcs), 7);
            tdata = zeros(sum(tcs), 4);
            tdata(:, 1) = 3;
            vc = 1;
            tc = 1;
            for oc = 1:numel(scdata{2})
                vdata(vc:vc+vcs(oc)-1, 1:3) = ...
                    128 - scdata{2}{oc}.CONT.VertexCoordinate(:, [3, 1, 2]);
                vdata(vc:vc+vcs(oc)-1, 4:6) = ...
                    double(squeeze(scdata{3}{oc}.SurfTIORendered));
                vdata(vc:vc+vcs(oc)-1, 7) = ...
                    limitrangec(round(255 * scdata{3}{oc}.SurfProps{4}), 0, 255, 0);
                tdata(tc:tc+tcs(oc)-1, 2:4) = ...
                    scdata{2}{oc}.CONT.TriangleVertex + (vc - 2);
                vc = vc + vcs(oc);
                tc = tc + tcs(oc);
            end
            
            % then write header and data
            lf = char(10);
            fwrite(fid, ['ply' lf]);
            fwrite(fid, ['format ascii 1.0' lf]);
            fprintf(fid, 'comment NeuroElf GUI scenery with %d objects\n', numel(scdata{2}));
            fprintf(fid, 'element vertex %d\n', size(vdata, 1));
            fprintf(fid, ['property float x\nproperty float y\nproperty float z\n' ...
                'property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\n']);
            fprintf(fid, 'element face %d\n', size(tdata, 1));
            fwrite(fid, ['property list uchar int vertex_index' lf]);
            fwrite(fid, ['end_header' lf]);

            % write data
            fprintf(fid, '%g %g %g %d %d %d %d\n', vdata');
            fprintf(fid, '%d %d %d %d\n', tdata');

            % close file
            fclose(fid);
            
            % print header
        catch ne_eo;
            uiwait(warndlg(ne_eo.message, 'NeuroElf - error', 'modal'));
            if ~isempty(fid)
                try
                    fclose(fid);
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end
            end
        end
        ch.MainFig.Pointer = mfp;

    % set colors
    case {'setcolors'}

        % get current colors
        cols = round(255 .* [srf.ConvexRGBA(1:3); srf.ConcaveRGBA(1:3)]);

        % get new colors
        if nargin < 4 || ...
           ~isa(varargin{4}, 'double') || ...
           ~isequal(size(varargin{4}), [2, 3]) || ...
            any(isinf(varargin{4}(:)) | isnan(varargin{4}(:)) | varargin{4}(:) < 0 | varargin{4}(:) > 255)
            newcols = colorpicker(cols, {'Color 1'; 'Color 2'});
        else
            newcols = round(varargin{4});
        end

        % compare colors
        if ~isequal(cols, newcols)

            % set colors
            srf.ConvexRGBA(1:3) = (1 / 255) .* newcols(1, :);
            srf.ConcaveRGBA(1:3) = (1 / 255) .* newcols(2, :);

            % udpate coloring
            btc_meshcolor(srf, true);
        end

    % set morphing target
    case {'setmorphtarget'}

        % get list of all SRFs that are available for a match
        srftri = srf.TriangleVertex;
        srfs = ne_gcfg.h.SurfVar.UserData;
        srfi = false(size(srfs, 1), 1);
        for vc = 1:numel(srfi)
            if srfs{vc, 4} ~= srf && isequal(srfs{vc, 4}.TriangleVertex, srftri)
                srfi(vc) = true;
            end
        end

        % no match
        if ~any(srfi)
            return;
        end
        
        % input given
        tsrf = [];
        if nargin > 3 && ...
            numel(varargin{4}) == 1 && ...
            isxff(varargin{4}, 'srf') && ...
            isequal(varargin{4}.TriangleVertex, srftri)
            tsrfh = handles(tsrf);
            if isfield(tsrfh, 'SurfProps') && ...
                iscell(tsrfh.SurfProps) && ...
                numel(tsrfh.SurfProps) >=5 && ...
                isfield(tsrfh, 'SurfTIO') && ...
                numel(tsrfh.SurfTIO) == 1 && ...
                isa(tsrfh.SurfTIO, 'transimg') && ...
                isfield(tsrfh, 'ShownInGUI') && ...
                islogical(tsrfh.ShownInGUI) && ...
                numel(tsrfh.ShownInGUI) == 1 && ...
                tsrfh.ShownInGUI
                tsrf = varargin{4};
            end
        end
        
        % display matches
        if isempty(tsrf)
            srfn = srfs(srfi, 1);
            srfo = srfs(srfi, 4);
            for vc = 1:numel(srfn)
                srffn = srfo{vc}.FilenameOnDisk;
                if isempty(srffn)
                    srfn{vc} = '<unsaved.srf>';
                else
                    srfn{vc} = srffn;
                end
            end

            % show a selector
            [csel, csok] = listdlg( ...
                'PromptString',  'Please select SRF to set as morphing target...', ...
                'ListString',    srfn, ...
                'InitialValue',  1, ...
                'SelectionMode', 'single', ...
                'ListSize',      [min(600, max(300, 10 .* size(char(srfn), 2))), 200]);
            if ~isequal(csok, 1) || ...
                isempty(csel)
                srf.SetHandle('VertexMorphMeshes', cell(0, 5));
                ne_setsurfpos;
                return;
            end
            tsrf = srfo{csel};
        end

        % set as target
        mtsc = getscont(tsrf);
        morphmesh = {mtsc.C.VertexCoordinate, mtsc.C.VertexNormal, ...
             mtsc.C.MeshCenter, double(mtsc.H.SurfTIO.Rendered), mtsc.H.SurfProps};
        if isfield(mtsc.C.RunTimeVars, 'VertexAlpha') && ...
            isa(mtsc.C.RunTimeVars.VertexAlpha, 'double') && ...
            numel(mtsc.C.RunTimeVars.VertexAlpha) == size(morphmesh{1}, 1)
            morphmesh{5}{4} = limitrangec(mtsc.C.RunTimeVars.VertexAlpha, 0, 1, 0);
        end
        if ~isempty(mtsc.H.VertexMorphMeshes)
            morphmesh = cat(1, morphmesh, mtsc.H.VertexMorphMeshes);
        end
        srf.SetHandle('VertexMorphMeshes', morphmesh);
        
        % update filename in scenery
        scu = ch.Scenery.UserData;
        for sci = 1:size(scu, 1)
            if srf == scu{sci, 4}
                scu{sci, 1} = regexprep(scu{sci, 1}, 'NrOfVertices\:.*', ...
                    sprintf('NrOfVertices: %d, NrOfMorphStates: %d)', ...
                    size(morphmesh{1}, 1), size(morphmesh, 1) + 1));
                scs = ch.Scenery.String;
                if ~iscell(scs)
                    scs = cellstr(scs);
                end
                scs{sci} = scu{sci, 1};
                ch.Scenery.String = scs;
                ch.Scenery.UserData = scu;
            end
        end

        % update
        ne_setsurfpos(0, 0, 'upshape');

    % smoothing
    case {'smooth'}

        % ask for number of iterations and force
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 3 || ...
           ~ischar(varargin{4}{1}) || ...
           ~ischar(varargin{4}{2}) || ...
           ~ischar(varargin{4}{3})
            uio = inputdlg( ...
                {'Number of smoothing iterations:', 'Smoothing force ([0 ... 1]):', ...
                 'Weighting (n)one, (d)istance, (sq)uared distance:'}, ...
                 'NeuroElf GUI - user input', 1, {'  100', '  0.07', '  n'});
            if numel(uio) ~= 3 || ...
               ~iscell(uio) || ...
               ~ischar(uio{1}) || ...
               ~ischar(uio{2}) || ...
               ~ischar(uio{3})
                return;
            end
        else
            uio = varargin{4};
        end
        uio{1}(uio{1} == ' ') = [];
        uio{2}(uio{2} == ' ') = [];
        uio{3} = lower(ddeblank(uio{3}));
        try
            uio{1} = str2double(uio{1});
            uio{2} = str2double(uio{2});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        if isinf(uio{1}) || ...
            isnan(uio{1}) || ...
            uio{1} <= 0 || ...
            isinf(uio{2}) || ...
            isnan(uio{2}) || ...
            uio{2} <= 0 || ...
            uio{2} > 1 || ...
            isempty(uio{3}) || ...
           ~any('dns' == uio{3}(1))
            return;
        end
        if uio{3}(1) == 'n'
            xopts = {};
        elseif uio{3}(1) == 'd'
            xopts = {'distw', 1};
        elseif uio{3}(1) == 's'
            xopts = {'distw', 1, 'distwsq', 1};
        end

        % echo
        if ne_gcfg.c.echo
            if isempty(xopts)
                ne_echo('srf', 'Smooth', uio{1:2});
            else
                ne_echo('srf', 'Smooth', uio{1:2}, struct(xopts{:}));
            end
        end

        % make call
        cprog = ne_progress(0, 0, {true, 0, 'Smoothing surface...'});
        srf.Smooth(uio{1:2}, struct('pbar', ch.Progress, xopts{:}));
        ne_progress(0, 0, cprog);

        % update
        ne_setsurfpos(0, 0, 'upshape');

    % smooth SMP
    case {'smoothsmp'}

        % valid SMP object
        if nargin < 4 || ...
            numel(varargin{4}) ~= 1 || ...
           ~isxff(varargin{4}, 'smp')
            smp = cc.SurfStatsVar;
        else
            smp = varargin{4};
        end
        if ~isxff(smp, 'smp') || ...
            isempty(smp.Map) || ...
            numel(smp.Map(1).SMPData) ~= size(srf.VertexCoordinate, 1)
            return;
        end

        % map numbers
        if nargin < 5 || ...
           ~isa(varargin{5}, 'double') || ...
            any(isinf(varargin{5}(:)) | isnan(varargin{5}(:)) | varargin{5}(:) < 1)
            mapidx = ne_gcfg.h.SurfStatsVarMaps.Value;
            if isempty(mapidx)
                return;
            end
        else
            mapidx = unique(round(varargin{5}(:)))';
        end

        % config
        if nargin < 6 || ...
           ~iscell(varargin{6}) || ...
            numel(varargin{6}) ~= 2
            smcfg = inputdlg({'Smoothing kernel (mm)', 'Number of iterations'}, ...
                'NeuroElf - SMP smoothing configuration', 1, {'  6', '  1'});
        else
            smcfg = varargin{6};
        end

        % assess config
        try
            smk = str2double(smcfg{1});
            if isinf(smk) || ...
                isnan(smk) || ...
                smk <= 0
                return;
            end
            niter = str2double(smcfg{2});
            if isinf(niter) || ...
                isnan(niter) || ...
                niter < 1
                return;
            else
                niter = floor(niter);
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % run smoothing, then re-load SMP
        smp.Smooth(srf, niter, mapidx, smk);
        if nargin < 4 || ...
            numel(varargin{4}) ~= 1 || ...
           ~isxff(varargin{4}, 'smp')
            ne_setcsrfstats;
        end

        % return smp
        srf = smp;

    % clone current SMP
    case {'smpclone'}

        % get current selection
        stvar = srfh.Stats{1};
        if ~isxff(stvar, 'smp')
            return;
        end
        
        % copy object
        stvar = stvar.CopyObject;
        
        % add to interface (with same map selection)
        ne_openfile(0, 0, stvar, true);
        ne_setcsrfstatmap(0, 0, srfh.Stats{2});
        
        % return
        srf = stvar;
        
    % SMP-based formula
    case {'smpformula'}

        % get current selection
        stvar = srfh.Stats{1};
        if ~isxff(stvar, 'smp')
            return;
        end
        stvix = srfh.Stats{2}(:);
        formdesc = '#i := i-th map';
        if ~isempty(stvix)
            formdesc = [formdesc ', $i := i-th selected map'];
        end

        % get formula and other values
        formset = {'  #1', sprintf('  %d', numel(stvar.Map) + 1), ...
            '  ', '  ', '  ', '  ', '  n'};
        try
            formset = ddeblank(inputdlg({ ...
                ['Formula: ' formdesc], ...
                'Target map number:', ...
                'Target map name (leave blank to use formula):', ...
                'Map type (either of t, r, CC, or F; leave blank to copy from 1st):', ...
                'DF1 (leave blank to copy from first):', ...
                'DF2 (leave blank to copy from first):', ...
                'Use p-values (conversion required) instead of stats, (y)es/(n)o:'}, ...
                'NeuroElf GUI - Map settings', 1, formset));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % check formula fields (as much as possible)
        if ~iscell(formset) || ...
            numel(formset) ~= 7
            return;
        end
        maptypes = {'t', 'r', 'CC', 'F', '', '', '', '', '', '', ...
            '%', 'z_ica', 'TH', '', 'beta', 'prob', '', '', '', ...
            'MD', 'FA'};
        if isempty(formset{1}) || ...
           ~any(formset{1} == '$' | formset{1} == '#') || ...
            isempty(formset{2}) || ...
            any(formset{2} < '0' | formset{2} > '9') || ...
            numel(formset{4}) > 1 || ...
           (~isempty(formset{4}) && ...
            ~any(strcmpi(formset{4}, maptypes))) || ...
            isempty(formset{7}) || ...
            ~any(lower(formset{7}(1)) == 'ny')
            return;
        end
        formset{2} = str2double(formset{2});
        if isempty(formset{3})
            formset{3} = formset{1};
        end
        if formset{2} < 1 || ...
            formset{2} > (numel(stvar.Map) + 1)
            return;
        end

        % perform command
        opts = struct( ...
            'mapsel',  stvix(:)', ...
            'name',    formset{3}, ...
            'pvalues', (lower(formset{7}(1)) == 'y'), ...
            'target',  formset{2});
        stvar.ComputeFormula(formset{1}, opts);

        % echo
        if ne_gcfg.c.echo
            ne_echo('smp', 'ComputeFormula', opts);
        end

        % adapt type and DF if requested
        if ~isempty(formset{4})
            stvar.Map(formset{2}).Type = findfirst(strcmpi(maptypes, formset{4}));
        end
        if ~isempty(formset{5}) && ...
           ~isempty(regexpi(formset{5}, '^\d+$')) && ...
            formset{5}(1) > '0'
            stvar.Map(formset{2}).DF1 = str2double(formset{5});
        end
        if ~isempty(formset{6}) && ...
           ~isempty(regexpi(formset{6}, '^\d+$')) && ...
            formset{6}(1) > '0'
            stvar.Map(formset{2}).DF2 = str2double(formset{6});
        end

        % bring up new map
        stvar.SetColors(formset{2}, 'xauto');
        ne_openfile(0, 0, stvar);
        ne_gcfg.h.SurfStatsVarMaps.Value = formset{2};
        ne_setcsrfstatmap;

    % SRF info
    case {'srfinfo'}

        % get (super-) content
        srfsc = getscont(srf);

        % computations
        c = srfsc.C.VertexCoordinate;
        t = srfsc.C.TriangleVertex;
        [vn, tn] = mesh_normals(c, t);
        vn = -vn;
        vnd = sqrt(sum((vn - srfsc.C.VertexNormal) .^ 2, 2));
        gn = 0.5 * (2 - mean(vnd));
        mc = srfsc.C.MeshCenter;
        [a, ta] = srf.Area;
        mmta = minmaxmean(ta);
        rc = mean(c);
        minc = min(c) - 0.04999;
        maxc = max(c) + 0.04999;
        c = c - ones(size(c, 1), 1) * rc;
        ct = cat(3, c(t(:, 1), :), c(t(:, 2), :), c(t(:, 3), :));
        mct = mean(ct, 3);
        tnal = sqrt(sum(mct .* mct, 2));
        tnal = tnal .* min(1, sum(mct .* tn, 2) ./ (tnal .* sqrt(sum(tn .* tn, 2))));
        tnv = (1 / 3) .* (ta(:) .* tnal(:));
        ttnv = sum(tnv);

        % print out stuff
        disp(' ');
        disp('NeuroElf GUI - Surface information');
        disp('------------------------------------------------------------------------');
        if isempty(srfsc.F)
            srff = '<unsaved.srf>';
        else
            srff = srfsc.F;
            if numel(srff) > 62
                srff = ['''' srff(1:20) ' ... ' srff(end-34:end) ''''];
            end
        end
        disp(['Filename: ' srff]);
        fprintf('NrOfVertices:  %d\n', size(c, 1));
        fprintf('NrOfTriangles: %d\n', size(t, 1));
        fprintf('MeshCenter:    [%f, %f, %f]\n', 128 - mc(1, [3, 1, 2]));
        disp('------------------------------------------------------------------------');
        fprintf('Real center coordinate:      [%.1f, %.1f, %.1f]\n', 128 - rc(1, [3, 1, 2]));
        fprintf('Data fits into box:          [%.1f, %.1f, %.1f; %.1f, %.1f, %.1f]\n', ...
            128 - [maxc(1, [3, 1, 2]), minc(1, [3, 1, 2])]);
        disp('------------------------------------------------------------------------');
        fprintf('Total area of triangles:     %.2fmm^2\n', a);
        fprintf('Smallest triangle covers:    %.4fmm^2\n', mmta(1));
        fprintf('Largest triangle covers:     %.2fmm^2\n', mmta(2));
        disp('------------------------------------------------------------------------');
        fprintf('Total enclosed volume:       %.1fmm^3\n', ttnv);
        fprintf('Corresponding sphere radius: %.2fmm\n', (3 * ttnv / (4 * pi)) ^ (1 / 3));
        disp('------------------------------------------------------------------------');
        fprintf('Goodness of normals:         %.1f\n\n', 100 * gn);
        drawnow;

    % tosphere
    case {'tosphere'}

        % ask for number of iterations and force
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 3 || ...
           ~ischar(varargin{4}{1}) || ...
           ~ischar(varargin{4}{2}) || ...
           ~ischar(varargin{4}{3})

            % estimate required steps
            c = srf.VertexCoordinate;
            c = c - ones(size(c, 1), 1) * mean(c);
            c = sqrt(sum(c .* c, 2));
            c = (1 / mean(c)) .* c;
            c = min(3, 2 * std(c));
            rsteps = 5000 * min(25, max(0.4, round(c * 0.05 * sqrt(srf.NrOfVertices))));
            rsteps = rsteps .* [0.4, 0.24, 0.16, 0.1, 0.1];
            rforce = max([1e-5, 2e-5, 1e-4, 2e-4, 2e-3], [0.4, 0.6, 0.8, 2.5, 5] ./ rsteps);
            rsteps = sprintf('  [%d, %d, %d, %d, %d]', rsteps);
            rforce = sprintf('  [%.5f, %.5f, %.4f, %.4f, %.3f]', rforce);
            uio = inputdlg( ...
                {'Number of to-sphere-morphing iterations (array, called separatly):', ...
                'Smoothing forces ([0 ... 1]):',  'To-sphere forces ([0 ... 0.1])'}, ...
                'NeuroElf GUI - user input', 1, ...
                {rsteps, '  [0.500, 0.250, 0.150, 0.125, 0.100]', rforce});
            if numel(uio) ~= 3 || ...
               ~iscell(uio) || ...
               ~ischar(uio{1}) || ...
               ~ischar(uio{2}) || ...
               ~ischar(uio{3})
                return;
            end
        else
            uio = varargin{4};
        end
        try
            uio{1} = eval(uio{1});
            uio{2} = eval(uio{2});
            uio{3} = eval(uio{3});
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        if ~isa(uio{1}, 'double') || ...
           ~isa(uio{2}, 'double') || ...
           ~isa(uio{3}, 'double') || ...
            numel(uio{1}) ~= size(uio{1}, 2) || ...
           ~isequal(size(uio{1}), size(uio{2})) || ...
           ~isequal(size(uio{1}), size(uio{3})) || ...
            any(isinf(uio{1}) | isnan(uio{1}) | uio{1} < 1  | uio{1} > 50000 | ...
                isinf(uio{2}) | isnan(uio{2}) | uio{2} <= 0 | uio{2} > 1 | ...
                isinf(uio{3}) | isnan(uio{3}) | uio{3} < 0  | uio{3} > 0.1)
            return;
        end

        % echo
        if ne_gcfg.c.echo
            ne_echo('srf', 'ToSphere', uio{:});
        end

        % make call
        ne_gcfg.c.blockcb{end+1} = 'ne_srf_tools';
        cprog = ne_progress(0, 0, {true, 0, 'To-sphere-morphing...'});
        try
            srf.ToSphere(uio{:}, struct('pbar', ch.Progress));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        ne_progress(0, 0, cprog);
        ne_gcfg.c.blockcb(strcmpi(ne_gcfg.c.blockcb, 'ne_srf_tools')) = [];

        % update
        ne_setsurfpos(0, 0, 'upshape');

    % undo last morph
    case {'undomorph'}

        % check
        if isfield(srf.RunTimeVars, 'UndoBuffer') && ...
            isequal(size(srf.VertexCoordinate), size(srf.RunTimeVars.UndoBuffer))

            % put back
            srf.VertexCoordinate = srf.RunTimeVars.UndoBuffer;
            srf.SetHandle('VertexCoordinateTal', []);

            % recalc normals
            srf.RecalcNormals;
        end

        % update
        ne_setsurfpos(0, 0, 'upshape');

    % update colors/shape
    case {'update'}

        % is surface displayed
        scu = ch.Scenery.UserData;
        scv = ch.Scenery.Value;
        srfd = false;
        for srfc = 1:numel(scv)
            if srf == scu{scv(srfc), 4}
                srfd = true;
                break;
            end
        end

        % update
        if srfd
            if isfield(srf.Handles, 'VertexCoordinateTal')
                srf.DeleteHandle('VertexCoordinateTal');
            end
            ne_setsurfpos(0, 0, 1, 'upshape');
        end

% nothing after all
    otherwise
        return;
end

% and put into output
if nargout > 0
    varargout{1} = srf;
end

% and update one last time
ne_setsurfpos;

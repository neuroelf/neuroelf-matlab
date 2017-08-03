function varargout = ne_hdr_dti(varargin)
% ne_hdr_dti  - DTI tools for currently selected SliceVar (HDR)
%
% FORMAT:       [result = ] ne_hdr_dti(SRC, EVT, action [, arguments, ...])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       action      one of the supported actions (see below)
%       arguments   additional (usually optional) arguments
%
% Output fields:
%
%       result      resulting object or variable (if any)
%
% Notes: calls that do not specify a surface object will always use the
%        currently selected one: ne_gcfg.fcfg.SurfVar --if no surface is
%        open in the GUI, those calls are then invalid (no action taken)
%
% Actions:
%
%   'dtcalc'        - compute diffusion tensor and eigenvectors/-values
%       syntax:     [DTI, RGB = ] ne_hdr_dti(0, 0, 'dtcalc')
%       outputs:
%       DTI         NIftI object with 26 volumes:
%           ( 1- 3) 4D, 3-axes FA-weighted first eigenvector map
%           ( 4- 6) 3 eigenvalue maps
%           ( 7)    FA (fractional anisotropy) map
%           ( 8)    B0 estimated from full dataset
%           ( 9-11) 4D, 3-axes first eigenvector map (original units)
%           (12-14) 4D, 3-axes second eigenvector map
%           (15-17) 4D, 3-axes third eigenvector map
%           (18-26) 4D, 9-element (full) DTI (tensor) map
%       RGB         same as first three volumes in DTI dataset, but as
%                   uint8 and 5-dim display (RGB-dataset hack in NeuroElf)
%
%   Note: both objects will be opened by the GUI, so it is not strictly
%         necessary to receive outputs to gain access
%
%   'flipbvecs'     - flip BVecs in RunTimeVars in some direction
%       syntax:     ne_hdr_dti(0, 0, 'flipbvecs', DIRS)
%       DIRS        string containing 'x', 'y', and/or 'z'
%       example:    ne_hdr_dti(0, 0, 'flipbvecs', 'x')
%
%   'motcorr'       - apply (3-parameter, 2-pass) motion correction
%       syntax:     ne_hdr_dti(0, 0, 'motcorr')
%
%   'plotfibers'    - plot detected fibers (from whole brain or VOIs)
%       syntax      
%
%   'trackfibers'   - perform fiber-tracking algorithm

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:26 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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

% part of this code was derived from DTISearch / Streamline found at
% http://www.mathworks.com/matlabcentral/fileexchange/34008
% and courtesy of Chang Chia-Hao (also under BSD license)

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% active SRF loaded
hdr = cc.SliceVar;
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   (~strcmpi(varargin{3}(:)', 'plotfibers') && ...
    (numel(hdr) ~= 1 || ...
     ~isxff(hdr, 'hdr')))
    return;
end
action = lower(varargin{3}(:)');
if strcmp(action, 'plotfibers')
    srf = cc.SurfVar;
    if numel(srf) ~= 1 || ...
       ~isxff(srf, 'srf') || ...
       ~isfield(srf.RunTimeVars, 'SourceDTIxffID')
        xffroot = xff;
        srfs = xffroot.Documents('srf');
        if isempty(srfs)
            uiwait(warndlg('Plotting fibers requires a DTI-Tracking-SRF object to be selected.', ...
                'NeuroElf - info', 'modal'));
            return;
        end
        srfo = srfs;
        for sc = numel(srfs):-1:1
            try
                srfo{sc} = xffroot.Document(srfs{sc});
                if numel(srfo{sc}) ~= 1 || ...
                   ~isxff(srfo{sc}, 'srf') || ...
                   ~isfield(srfo{sc}.RunTimeVars, 'SourceDTIxffID') || ...
                   ~isfield(srfo{sc}.RunTimeVars, 'FiberTrackTitle')
                    srfs(sc) = [];
                    srfo(sc) = [];
                else
                    srfs{sc} = sprintf('%s (%d vertices)', ...
                        srfo{sc}.RunTimeVars.FiberTrackTitle, srfo{sc}.NrOfVertices);
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                srfs(sc) = [];
                srfo(sc) = [];
            end
        end
        if isempty(srfo)
            uiwait(warndlg('Plotting fibers requires a DTI-Tracking-SRF object to be selected.', ...
                'NeuroElf - info', 'modal'));
            return;
        end
        if numel(srfo) > 1
            [srfsel, selok] = listdlg('ListString', srfs, 'SelectionMode', 'single', ...
                'ListSize', [640, 200], 'InitialValue', 1, ...
                'Name', 'Please select the DTI-trackts object to plot from...');
        else
            srfsel = 1;
            selok = 1;
        end
        if isequal(selok, 0) || ...
            numel(srfsel) ~= 1 || ...
           ~any(srfsel == (1:numel(srfs)))
            return;
        end
        srf = srfo{srfsel};
    end
    rtv = srf.RunTimeVars;
else
    rtv = hdr.RunTimeVars;
end

% certain actions don't require BVals/BVecs to match in size!
if any(strcmp(action, {'plotfibers', 'trackfibers'}))
    if numel(cc.dti) == 1 && ...
        isxff(cc.dti, 'hdr')
        hdr = cc.dti;
    end
    bvmatch = false;
else
    bvmatch = true;
end

% no double calls
if ~any(strcmp(action, {'canceldti'})) && ...
    any(strcmp(ne_gcfg.c.blockcb, 'hdr_dti'))
    return;
end
mfp = ch.MainFig.Pointer;

% no BVals/BVecs
vxs = size(hdr.VoxelData);
if ~isfield(rtv, 'BVals') || ...
   ~isa(rtv.BVals, 'double') || ...
   (bvmatch && ~isequal(size(rtv.BVals), [vxs(4), 1])) || ...
    any(isinf(rtv.BVals) | isnan(rtv.BVals) | rtv.BVals < 0) || ...
   ~isfield(rtv, 'BVecs') || ...
   ~isa(rtv.BVecs, 'double') || ...
   (bvmatch && ~isequal(size(rtv.BVecs), [vxs(4), 3])) || ...
    any(isinf(rtv.BVecs(:)) | isnan(rtv.BVecs(:)))
    return;
end

% B0 images
bvals = rtv.BVals;
bvecs = rtv.BVecs;
B0 = all([bvals, bvecs] == 0, 2);
if ~any(B0) || ...
    all(B0)
    return;
end
B1 = find(~B0);
B0 = find(B0);

% blocked?
if any(strcmpi(ne_gcfg.c.blockcb, 'ne_hdr_dti'))
    return;
end

% resolution
hdrcfr = hdr.CoordinateFrame(1);
hdrtrf = hdrcfr.Trf;
hdrres = hdrcfr.Resolution;

% decide
switch (action)

    % diffusion tensor calculation
    case {'dtcalc'}

        % block further access
        ne_gcfg.c.blockcb{end+1} = 'hdr_dti';
        ch.MainFig.Pointer = 'watch';
        drawnow;

        % apply some ad-hoc masking to reduce computational effort
        dwi = double(hdr.VoxelData(:, :, :, :));
        if all([0, 1] ~= hdr.ImgDim.ScalingSlope)
            if hdr.ImgDim.ScalingIntercept ~= 1
                dwi = hdr.ImgDim.ScalingSlope .* dwi + hdr.ImgDim.ScalingIntercept;
            else
                dwi = hdr.ImgDim.ScalingSlope .* dwi;
            end
        elseif hdr.ImgDim.ScalingIntercept ~= 0
            dwi = dwi + hdr.ImgDim.ScalingIntercept;
        end
        
        % run function
        cprog = ne_progress(0, 0, {true, 0, 'dwitensorcalc'});
        [dtidata, dwib0, fa, evals, evecs, faev1, b2] = dwitensorcalc(dwi, ...
            bvals, bvecs, struct('pbar', ch.Progress, 'res', hdrres, 'trf', hdrtrf));
        ne_progress(0, 0, cprog);

        % create 4D dataset containing all of these
        dtivol = repmat(single(NaN), [vxs(1:3), 26]);
        dtivol(:, :, :, 1:3) = faev1;
        dtivol(:, :, :, 4:6) = evals;
        dtivol(:, :, :, 7) = fa;
        dtivol(:, :, :, 8) = dwib0;
        dtivol(:, :, :, 9:17) = reshape(evecs, [vxs(1:3), 9]);
        dtivol(:, :, :, 18:26) = dtidata;

        % create copy of object and patch some information
        dti = hdr.CopyObject;
        dti.ImgDim.Dim(5) = 25;
        dti.ImgDim.DataType = 16;
        dti.ImgDim.BitsPerPixel = 32;
        dti.ImgDim.ScalingSlope = 1;
        dti.ImgDim.ScalingIntercept = 0;

        % store in volume
        dti.VoxelData = dtivol;

        % then adapt RunTimeVars
        dti.RunTimeVars.AutoSave = true;
        dti.RunTimeVars.BVals2 = b2;
        dti.RunTimeVars.Map = dti.RunTimeVars.Map(ones(1, 25));
        mapnames = { ...
            'FAwEV1(1)', 'FAwEV1(2)', 'FAwEV1(3)', ...
            'EVal(1)', 'EVal(2)', 'EVal(3)', ...
            'FA', ...
            'DWIb0Est', ...
            'EVec1(1)', 'EVec1(2)', 'EVec1(3)', ...
            'EVec2(1)', 'EVec2(2)', 'EVec2(3)', ...
            'EVec3(1)', 'EVec3(2)', 'EVec3(3)', ...
            'DTI(1,1)', 'DTI(2,1)', 'DTI(3,1)', ...
            'DTI(1,2)', 'DTI(2,2)', 'DTI(3,2)', ...
            'DTI(1,3)', 'DTI(2,3)', 'DTI(3,3)'};
        for mc = 1:numel(mapnames)
            dti.RunTimeVars.Map(mc).Name = mapnames{mc};
        end

        % re-scale to -1 .. 1
        dti.RunTimeVars.ScalingWindow = [-1, 1];
        dti.RunTimeVars.ScalingWindowLim = [-1, 1];
        dti.SetScalingWindow([-1, 1], true);
        dti.RunTimeVars.ScalingWindow = [0, 1];

        % open in GUI
        ne_openfile(0, 0, dti, true);

        % and return this object
        hdr = dti;

        % also create colored version with FAwEV1 only
        hc = dti.CopyObject;
        hc.VoxelData = hc.VoxelData(:, :, :, 1:3);
        hc.RunTimeVars.Map(4:end) = [];
        hc.VoxelData = reshape(hc.VoxelData, [vxs(1:3), 1, 3]);
        hc.ImgDim.Dim(1) = 5;
        hc.ImgDim.Dim(5:6) = [1,3];
        hc.VoxelData(isnan(hc.VoxelData(:))) = 0;
        hc.VoxelData = limitrangec(hc.VoxelData, 0, 1, 0);
        hc.VoxelData = uint8(round(255 .* hc.VoxelData));
        hc.ImgDim.DataType = 2;
        hc.ImgDim.BitsPerPixel = 8;
        hc.RunTimeVars.ScalingWindow = [0, 255];
        hc.RunTimeVars.ScalingWindowLim = [0, 255];
        hc.SetScalingWindow([0, 255], true);
        hc.RunTimeVars.ScalingWindow = [0, 128];
        hc.VoxelDataRGBA = hc.VoxelData;
        hc.VoxelData = [];
        hc.ImgDim.DataType = 128;
        hc.ImgDim.BitsPerPixel = 24;
        ne_openfile(0, 0, hc, true);
        varargout{2} = hc;

        % unblock access
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];
        ch.MainFig.Pointer = mfp;
        drawnow;

    % flip BVecs
    case {'flipbvecs'}

        % required argument
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{4}) || ...
           ~any(lower(varargin{4}(:)) == 'x' | lower(varargin{4}(:)) == 'y' | lower(varargin{4}(:)) == 'z')
            return;
        end
        flipdirs = lower(varargin{4}(:));
        if any(flipdirs == 'x')
            hdr.RunTimeVars.BVecs(:, 1) = -rtv.BVecs(:, 1);
        end
        if any(flipdirs == 'y')
            hdr.RunTimeVars.BVecs(:, 2) = -rtv.BVecs(:, 2);
        end
        if any(flipdirs == 'z')
            hdr.RunTimeVars.BVecs(:, 3) = -rtv.BVecs(:, 3);
        end

    % apply motion correction
    case {'motcorr'}

        % block further access
        ne_gcfg.c.blockcb{end+1} = 'hdr_dti';

        % progress bar
        cprog = ne_progress(0, 0, {true, 0, 'Setting up DTI motion correction...'});
        drawnow;

        % alignment of B0 images needed
        if numel(B0) > 1

            % but not yet implemented
            uiwait(warndlg('Alignment of B0 images not yet implemented.', ...
                'NeuroElf - warning', 'modal'));
        end

        % compute mean of non-B0 images
        meanimage = single(smoothdata3(mean(hdr.VoxelData(:, :, :, B1), 4), ...
            [3, 3, 3] ./ hdrres));

        % align images and keep track of parameters
        [t1, t2, vxdata] = rbalign(meanimage, single(hdr.VoxelData(:, :, :, B1)), ...
            struct('interpe', 'cubic', 'params', [true(1, 3), false(1, 3)], ...
            'pbar', ch.Progress, 'prange', [0, 0.4]));

        % recompute mean
        meanimage = single(smoothdata3(mean(vxdata, 4), [3, 3, 3] ./ hdrres));

        % second alingment step
        t2 = rbalign(meanimage, vxdata, struct('interpe', 'cubic', ...
            'params', [true(1, 3), false(1, 3)], 'pbar', ch.Progress, ...
            'prange', [0.4, 0.8]));

        % combine alignments, then subtract first volume
        t1(1:3, 4, :) = -(t1(1:3, 4, :) + t2(1:3, 4, :));
        t1(1:3, 4, 2:end) = t1(1:3, 4, 2:end) - repmat(t1(1:3, 4, 1), [1, 1, numel(B1) - 1]);

        % apply
        for vc = 1:numel(B1);
            ch.Progress.Progress(0.8 + 0.2 * vc / numel(B1), ...
                sprintf('Final resampling volume %d.', vc));
            hdr.VoxelData(:, :, :, B1(vc)) = flexinterpn_method( ...
                hdr.VoxelData(:, :, :, B1(vc)), [Inf, Inf, Inf; ones(2, 3); vxs(1:3)], ...
                0, t1(:, :, vc), 'lanczos3');
        end

        % progress bar again
        ne_progress(0, 0, cprog);
        drawnow;

        % unblock access
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];

    % plot colored fibers
    case {'plotfibers'}

        % only valid if fields set
        if ~isfield(rtv, 'FiberLookupMask') || ...
           ~isfield(rtv, 'FiberLookupCells') || ...
           ~isinteger(rtv.FiberLookupMask) || ...
           ~iscell(rtv.FiberLookupCells) || ...
            numel(rtv.FiberLookupCells) ~= sum(rtv.FiberLookupMask(:) > 0)
            uiwait(warndlg('Invalid fiber tracks.', 'NeuroElf - info', 'modal'));
            return;
        end
        fibers = srf.VertexCoordinate;
        fevec2 = srf.VertexNormal;
        fcolor = srf.VertexColor;
        flookup = cell(size(rtv.FiberLookupMask));
        vxs = size(flookup);
        vxx = vxs(1);
        vxy = vxs(2);
        vxz = vxs(3);
        flookup(rtv.FiberLookupMask > 0) = rtv.FiberLookupCells;
        itrf = rtv.TransTalToVoxel';

        % configuration
        if nargin < 5 || ...
           ~iscell(varargin{5}) || ...
            numel(varargin{5}) ~= 4
            uio = inputdlg({ ...
                'Alpha transparency of fibers:', ...
                'Number of steps that fiber endpoints must be within:', ...
                'Draw how many segments:', ...
                'Plot target: (o)bject, (s)atellite, or (n)ew'}, ...
                'NeuroElf - user input', 1, ...
                {'  1.0', '  Inf', '  Inf', '  object'});
            if ~iscell(uio) || ...
                numel(uio) ~= 4
                return;
            end
        else
            uio = varargin{5};
        end
        try
            uio = ddeblank(uio(:));
            falpha = str2double(uio{1});
            if isinf(falpha) || ...
                isnan(falpha) || ...
                falpha <= 0
                return;
            else
                falpha = min(1, falpha);
            end
            endsteps = str2double(uio{2});
            if isnan(endsteps) || ...
                endsteps < 1
                return;
            end
            drawsteps = str2double(uio{3});
            if isnan(drawsteps) || ...
                drawsteps <= 0
                return;
            end
            drawstepsp = drawsteps + 1;
            drawstepst = ceil(drawstepsp);
            drawstepw1 = drawstepst - drawstepsp;
            drawstepw2 = 1 - drawstepw1;
            drawstepsf = (drawstepsp ~= round(drawstepsp));
            if isempty(uio{4}) || ...
               ~any(lower(uio{4}(1)) == 'nos')
                return;
            end
            ptarget = lower(uio{4}(1));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end

        % spatial specification: fibers must cross through ...
        if nargin > 5 && ...
            numel(varargin{6}) == 1 && ...
            isxff(varargin{6}, 'voi')
            voi = varargin{6};
            if nargin > 6 && ...
                isa(varargin{7}, 'double') && ...
               ~isempty(varargin{7}) && ...
               ~any(isinf(varargin{7}(:)) | isnan(varargin{7}(:)) | ...
                    varargin{7}(:) < 1 | varargin{7}(:) > numel(voi.VOI))
                voiidx = unique(round(varargin{7}(:)))';
            else
                voiidx = 1:numel(voi.VOI);
            end
            voi = voi.VOI(voiidx);
            tvox = cat(1, voi.Voxels);
            cvox = round(mean(tvox));
            tvox(:, 4) = 1;
            vox = tvox * itrf;
            vox(:, 4) = [];
            vox = unique(round(vox), 'rows');
            vox = sub2ind(vxs, vox(:, 1), vox(:, 2), vox(:, 3));
        elseif numel(ne_gcfg.voi) == 1 && ...
            isxff(ne_gcfg.voi, 'voi') && ...
           ~isempty(ne_gcfg.voi.VOI) && ...
           ~isempty(ne_gcfg.h.Clusters.Value)
            try
                voi = ne_gcfg.voi.VOI(ne_gcfg.h.Clusters.Value);
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                return;
            end
            tvox = cat(1, voi.Voxels);
            cvox = round(mean(tvox));
            tvox(:, 4) = 1;
            vox = tvox * itrf;
            vox(:, 4) = [];
            vox = unique(round(vox), 'rows');
            vox = sub2ind(vxs, vox(:, 1), vox(:, 2), vox(:, 3));
        else
            vox = find(cellfun('prodofsize', flookup(:)) > 0);
            [tvox, tvoy, tvoz] = ind2sub(vxs, vox);
            cvox = round(mean([tvox(:), tvoy(:), tvoz(:)]));
        end

        % get all tracts we want to plot
        numvox = numel(vox);
        tidx = lsqueezecells(flookup(vox));
        tidx = unique(cat(1, tidx{:}))';

        % empty
        if isempty(tidx)
            uiwait(warndlg('No fibers selected to plot.', 'NeuroElf - info', 'modal'));
            return;
        end

        % then reconstruct indices
        fromidx = rtv.FiberStarts(tidx);
        toidx = rtv.FiberStarts(tidx+1) - 1;
        flength = double(toidx - fromidx);

        % sub-select
        if ~isinf(endsteps) && ...
            endsteps < max(flength)

            % match fiber ends to mask
            tfibers = [128 - fibers(:, [3, 1, 2]), ones(size(fibers, 1), 1)] * itrf;
            endmask = false(vxs);
            endmask(vox) = true;
            fmatch = true(numel(tidx), 1);
            for fc = 1:numel(fmatch)
                f = tfibers(fromidx(fc):toidx(fc)-1, 1:3);
                if size(f, 1) <= endsteps
                    continue;
                end
                tcrd = round(double(f(1:endsteps, :)));
                tcrd(any(tcrd < 1, 2) | tcrd(:, 1) > vxx | tcrd(:, 2) > vxy | tcrd(:, 3) > vxz, :) = [];
                if any(endmask(sub2ind(vxs, tcrd(:, 1), tcrd(:, 2), tcrd(:, 3))))
                    continue;
                end
                tcrd = round(double(f(end+1-endsteps:end, :)));
                tcrd(any(tcrd < 1, 2) | tcrd(:, 1) > vxx | tcrd(:, 2) > vxy | tcrd(:, 3) > vxz, :) = [];
                if any(endmask(sub2ind(vxs, tcrd(:, 1), tcrd(:, 2), tcrd(:, 3))))
                    fibers(fromidx(fc):toidx(fc)-1, :) = ...
                        fibers(toidx(fc)-1:-1:fromidx(fc), :);
                    fevec2(fromidx(fc):toidx(fc)-1, :) = ...
                        fevec2(toidx(fc)-1:-1:fromidx(fc), :);
                    fcolor(fromidx(fc):toidx(fc)-1, :) = ...
                        fcolor(toidx(fc)-1:-1:fromidx(fc), :);
                else
                    fmatch(fc) = false;
                end
            end

            % mask
            fromidx = fromidx(fmatch);
            toidx = toidx(fmatch);
            flength = flength(fmatch);
            if isempty(fromidx)
                uiwait(warndlg('No fibers remaining after endpoint selection.', ...
                    'NeuroElf - info', 'modal'));
                return;
            end
        end

        % additional mask (endpoints)
        if ~isinf(endsteps) && ...
            nargin > 6 && ...
            numel(varargin{end-1}) == 1 && ...
            isxff(varargin{end-1}, 'voi') && ...
           ~isempty(varargin{end-1}.VOI) && ...
            isa(varargin{end}, 'double') && ...
           ~any(isinf(varargin{end}(:)) | isnan(varargin{end}(:)) | ...
                varargin{end}(:) < 1 | varargin{end}(:) > numel(varargin{end-1}.VOI))
            voi = varargin{end-1};
            if ~isempty(varargin{end})
                voiidx = unique(round(varargin{end}(:)))';
            else
                voiidx = 1:numel(voi.VOI);
            end
            voi = voi.VOI(voiidx);
            vox = 128 - cat(1, voi.Voxels);
            vox = unique(round(vox), 'rows');
            vox = unique(sub2ind([256, 256, 256], vox(:, 2), vox(:, 3), vox(:, 1)));
            tfibers = [128 - fibers(:, [3, 1, 2]), ones(size(fibers, 1), 1)] * itrf;
            endmask = false(vxs);
            endmask(vox) = true;
            fmatch = true(numel(fromidx), 1);
            for fc = 1:numel(fmatch)
                f = tfibers(fromidx(fc):toidx(fc)-1, :);
                tcrd = round(double(f(max(1, end+1-endsteps):end, :)));
                tcrd(any(tcrd < 1, 2) | tcrd(:, 1) > vxx | tcrd(:, 2) > vxy | tcrd(:, 3) > vxz, :) = [];
                if ~any(endmask(sub2ind(vxs, tcrd(:, 1), tcrd(:, 2), tcrd(:, 3))))
                    fmatch(fc) = false;
                end
            end
            fromidx = fromidx(fmatch);
            toidx = toidx(fmatch);
            flength = flength(fmatch);
            if isempty(fromidx)
                uiwait(warndlg('No fibers remaining after endpoint selection.', ...
                    'NeuroElf - info', 'modal'));
                return;
            end
        end

        % actual plot
        if ptarget ~= 'o'

            % create plot
            hFig = [];
            if ptarget == 's'
                sats = fieldnames(ne_gcfg.cc);
                for satc = 1:numel(sats)
                    if isstruct(ne_gcfg.cc.(sats{satc})) && ...
                        isfield(ne_gcfg.cc.(sats{satc}), 'Config') && ...
                        isstruct(ne_gcfg.cc.(sats{satc}).Config) && ...
                        isfield(ne_gcfg.cc.(sats{satc}).Config, 'sattag') && ...
                        ischar(ne_gcfg.cc.(sats{satc}).Config.sattag) && ...
                        isfield(ne_gcfg.cc.(sats{satc}).Config, 'sattype') && ...
                        ischar(ne_gcfg.cc.(sats{satc}).Config.sattype) && ...
                        strcmp(ne_gcfg.cc.(sats{satc}).Config.sattype(:)', 'surf') && ...
                        isfield(ne_gcfg.cc.(sats{satc}), 'Surface') && ...
                        ishandle(ne_gcfg.cc.(sats{satc}).Surface)
                        htag = ne_gcfg.cc.(sats{satc});
                        hFig = htag.Satellite;
                        iSat = htag.Config.sattag;
                        break;
                    end
                end
            end
            if isempty(hFig)
                cpage = ne_gcfg.fcfg.page;
                if cpage ~= 3
                    ne_showpage(0, 0, 3);
                end
                [hFig, htag, iSat] = ne_undock;
                if cpage ~= 3
                    ne_showpage(0, 0, cpage);
                end
            end
            x = ne_gcfg.cc.(iSat).Surface;
            ps = zeros(numel(fromidx), 1);
            hold(x, 'on');
            hdrtrft = hdrtrf';
            for fc = 1:numel(ps)
                t = double(fibers{fc});
                if ~isinf(drawsteps) && ...
                    size(t, 1) > drawstepsp
                    t = t(1:drawstepst, :);
                    if drawstepsf
                        t(end, :) = ...
                            drawstepw1 .* t(end-1, :) + drawstepw2 .* t(end, :);
                    end
                end
                crd = [t(:, 1:3); NaN, NaN, NaN];
                crd(:, 4) = 1;
                crd = crd * hdrtrft;
                st = size(crd, 1);
                ps(fc) = patch(crd(:, 1), crd(:, 2), crd(:, 3), 0, 'Parent', x);
                cdata = [abs(t(:, 5:7)); NaN, NaN, NaN];
                cdata = reshape(cdata, [st, 1, 3]);
                set(ps(fc), 'CData', cdata, 'EdgeColor', 'interp', ...
                    'EdgeAlpha', falpha, 'FaceAlpha', 0, ...
                    'SpecularColorReflectance', 1, 'UserData', struct('pcrd', crd));
            end

            % update position (coordinates of patches!)
            ne_setsurfpos(0, 0, iSat);

            % return patch handles, as well as the tag portion and figure
            hdr = ps;
            if nargout > 1
                varargout{2} = iSat;
                if nargout > 2
                    varargout{3} = hFig;
                    if nargout > 3
                        varargout{4} = htag;
                    end
                end
            end

        % SRF object
        else

            % compile new SRF data
            if ~isinf(drawsteps)
                flength = min(flength, drawsteps);
            end
            tlength = sum(flength) + numel(flength);
            flength = flength - 1;
            toidx = fromidx + uint32(flength);
            newfibers = nan(tlength, 3);
            newfevecs = nan(tlength, 3);
            newcolors = nan(tlength, 4);
            tfrom = 1;
            flookup = cell(vxs);
            fcount = uint32(0);
            fcount(vxx, vxy, vxz) = 0;
            for fc = 1:numel(flength)
                fl = flength(fc);
                tto = tfrom + fl;
                fidx = fromidx(fc);
                eidx = toidx(fc);
                newfibers(tfrom:tto, :) = fibers(fidx:eidx, :);
                newfevecs(tfrom:tto, :) = fevec2(fidx:eidx, :);
                newcolors(tfrom:tto, :) = fcolor(fidx:eidx, :);
                tfibers = round([128 - fibers(fidx:eidx, [3, 1, 2]), ones(eidx + 1 - fidx, 1)] * itrf);
                tfibers(any(tfibers < 1, 2) | tfibers(:, 1) > vxx | tfibers(:, 2) > vxy | tfibers(:, 3) > vxz, :) = [];
                tfibers = unique(sub2ind(vxs, tfibers(:, 1), tfibers(:, 2), tfibers(:, 3)));
                tcount = fcount(tfibers) + 1;
                for lc = 1:numel(tfibers)
                    flookup{tfibers(lc)}(tcount(lc), 1) = fc;
                end
                fcount(tfibers) = tcount;
                tfrom = tto + 2;
            end

            % create new SRF (copy)
            srf = srf.CopyObject;
            srf.NrOfVertices = tlength;
            srf.VertexCoordinate = newfibers;
            srf.VertexNormal = newfevecs;
            srf.VertexColor = newcolors;
            srf.Neighbors = repmat({0, []}, tlength, 1);
            flmask = cellfun('prodofsize', flookup);
            mflmask = max(flmask(:));
            if mflmask < 256
                srf.RunTimeVars.FiberLookupMask = uint8(flmask);
            elseif mflmask < 65536
                srf.RunTimeVars.FiberLookupMask = uint16(flmask);
            else
                srf.RunTimeVars.FiberLookupMask = uint32(flmask);
            end
            srf.RunTimeVars.FiberLookupCells = flookup(flmask(:) > 0);
            srf.RunTimeVars.FiberStarts = uint32(1 + [0; cumsum(2 + flength(:))]);
            srf.RunTimeVars.FiberTrackTitle = ...
                sprintf('%s (masked with %d voxels around [%d, %d, %d])', ...
                rtv.FiberTrackTitle, numvox, cvox(1), cvox(2), cvox(3));

            % open in GUI
            ne_openfile(0, 0, srf, true);
            hdr = srf;
        end

    % track fibers
    case {'trackfibers'}

        % only if 25-volume matrix
        mapnames = {rtv.Map.Name};
        if ~any(strcmp(mapnames, 'FA')) || ...
           ~any(strcmp(mapnames, 'EVec1(1)')) || ...
           ~any(strcmp(mapnames, 'EVec1(2)')) || ...
           ~any(strcmp(mapnames, 'EVec1(3)')) || ...
           ~any(strcmp(mapnames, 'EVec2(1)')) || ...
           ~any(strcmp(mapnames, 'EVec2(2)')) || ...
           ~any(strcmp(mapnames, 'EVec2(3)'))
            return;
        end

        % some settings
        if nargin < 4 || ...
           ~iscell(varargin{4}) || ...
            numel(varargin{4}) ~= 9 || ...
           ~all(cellfun(@ischar, varargin{4}(:))) || ...
            any(cellfun('isempty', varargin{4}(:)))
            uio = inputdlg({ ...
                'Name for this fiber tracking:', ...
                'FA seed-threshold:', ...
                'Minimum fiber-length threshold (mm):', ...
                'Threshold of fiber deviation angle (degrees):', ...
                'FA stopping-criterion:', ...
                'Sampling step size (mm):', ...
                'Maximum number of steps:', ...
                'Spatial over-sampling factor: (1 - 3)', ...
                'Track from selected VOI clusters only:'}, ...
                'NeuroElf - user input', 1, ...
                {'  whole-brain', '  0.5', '  25', '  20', '  0.2', '  1.0', '  350', '  1', '  no'});
            if ~iscell(uio) || ...
                numel(uio) ~= 9
                ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];
                return;
            end
        else
            uio = varargin{4}(:);
        end
        try
            uio = ddeblank(uio(:));
            fname = uio{1};
            ftag = makelabel(fname);
            if ftag(1) == 'V'
                return;
            end
            seedfa = str2double(uio{2});
            if isinf(seedfa) || ...
                isnan(seedfa) || ...
                seedfa <= 0 || ...
                seedfa > 1
                return;
            end
            flmin = str2double(uio{3});
            if isinf(flmin) || ...
                isnan(flmin) || ...
                flmin < 1
                ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];
                return;
            end
            tangle = str2double(uio{4});
            if isinf(tangle) || ...
                isnan(tangle) || ...
                tangle <= 0 || ...
                tangle >= 90
                ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];
                return;
            end
            tfa = str2double(uio{5});
            if isinf(tfa) || ...
                isnan(tfa) || ...
                tfa <= 0 || ...
                tfa >= 1
                ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];
                return;
            end
            stepsize = str2double(uio{6});
            if isinf(stepsize) || ...
                isnan(stepsize) || ...
                stepsize <= 0 || ...
                stepsize > (2 * harmmean(hdrres))
                ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];
                return;
            end
            maxtiter = round(max(2, min(2000, str2double(uio{7}))));
            if isinf(maxtiter) || ...
                isnan(maxtiter)
                return;
            end
            oversmp = str2double(uio{8});
            if isinf(oversmp) || ...
                isnan(oversmp) || ...
               ~any((1:11) == oversmp)
                return;
            end
            if ~isempty(uio{9}) && ...
                lower(uio{9}(1)) == 'y'
                clonly = true;
            else
                clonly = false;
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];
            return;
        end

        % block further access
        ne_gcfg.c.blockcb{end+1} = 'hdr_dti';

        % seed VOI
        if clonly && ...
            numel(ne_gcfg.voi) == 1 && ...
            isxff(ne_gcfg.voi, 'voi') && ...
           ~isempty(ne_gcfg.voi.VOI) && ...
           ~isempty(ch.Clusters.Value)
            seedvois = ne_gcfg.voi.VOI(ch.Clusters.Value);
        elseif nargin > 4 && ...
            numel(varargin{5}) == 1 && ...
            isxff(varargin{5}, 'voi') && ...
           ~isempty(varargin{5}.VOI)
            seedvoi = varargin{5};
            if nargin > 5 && ...
                isa(varargin{6}, 'double') && ...
               ~isempty(varargin{6}) && ...
               ~any(isinf(varargin{6}(:)) | isnan(varargin{6}(:)) | varargin{6}(:) < 1)
                seedvois = unique(min(numel(seedvoi.VOI), round(varargin{6}(:)')));
            else
                seedvois = 1:numel(seedvoi.VOI);
            end
            seedvois = seedvoi.VOI(seedvois);
        else
            seedvois = [];
        end
        if ~isempty(seedvois)
            vox = cat(1, seedvois.Voxels);
            vox(:, 4) = 1;
            vox = vox * inv(hdrtrf)';
            vox(:, 4) = [];
            vox = unique(round(vox), 'rows');
            vox(any(vox < 1, 2) | vox(:, 1) > vxs(1) | vox(:, 2) > vxs(2) | vox(:, 3) > vxs(3), :) = [];
            vox = sub2ind(vxs(1:3), vox(:, 1), vox(:, 2), vox(:, 3));
            seedvois = false(vxs(1:3));
            seedvois(vox) = true;
        end

        % pointer
        ch.MainFig.Pointer = 'watch';
        drawnow;

        % track
        cprog = ne_progress(0, 0, {true, 0, 'DTI tracking...'});
        hdr = hdr.DTITrackEVs(struct( ...
            'angthresh', tangle, 'fathresh', tfa, 'flmin', flmin, 'faweight', true, ...
            'maxtiter', maxtiter, 'oversmp', oversmp, 'pbar', ch.Progress, ...
            'seedfa', seedfa, 'seedmask', seedvois, 'stepsize', stepsize));
        ne_progress(0, 0, cprog);
        hdr.RunTimeVars.FiberTrackTitle = fname;

        % unblock access
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'hdr_dti')) = [];
        ne_openfile(0, 0, hdr, true);

        % alpha level
        if oversmp > 1 && ...
            isxff(hdr, true)
            hdr.SetHandle('SurfProps', ...
                {[0, 0, 0], [0, 0, 0], [1, 1, 1], 1 / oversmp, 'w', [], 'none'});
            btc_meshcolor(hdr, true);
            ne_setsurfpos;
        end
        ch.MainFig.Pointer = mfp;
        drawnow;

    % nothing after all
    otherwise
        return;
end

% and put into output
if nargout > 0
    varargout{1} = hdr;
end

% and update one last time
if cc.page ~= 3
    ne_setslicepos;
end;

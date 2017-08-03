function [xo, densm] = srf_Morph(xo, niter, force, type, opts)
% SRF::Morph  - apply vertex morphing
%
% FORMAT:       [srf, densm] = srf.Morph([niter, force [, type [, opts]])
%
% Input fields:
%
%       niter       number of iterations (default 1)
%       force       smoothing force applied to morphing (default: 0.07)
%       type        either of 'even', {'smooth'}
%       opts        optional settings
%        .areac     if given and true, keep area constant
%        .areaw     area-weighted smoothing (higher precedence, false)
%        .distc     additional distortion correction force (>0 ... 3)
%        .distw     distance-weighted smoothing (lower precedence, false)
%        .distwsq   square-of-distance weighting (also sets distw, false)
%        .norm      force along normal vector (default: 0)
%        .normramp  ramp-up normal force from 0 to .norm (default: false)
%        .pbar      optional progress bar object (xfigure/xprogress, [])
%        .show      show during morphing (default: true)
%        .sphere    additionally applied to-sphere force
%        .title     display title for progress bar
%
% Output fields:
%
%       srf         altered object
%       densm       density map SMP
%
% Using: mesh_morph.

% Version:  v1.1
% Build:    16060713
% Date:     Jun-07 2016, 1:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% neuroelf library and global configuration
global ne_methods xffsngl;
mesh_morph = ne_methods.mesh_morph;
stc = xffsngl.CONF.settings;
if ~isfield(stc, 'Morphing')
    stc = struct('Morphing', struct);
end
stc = stc.Morphing;
if isfield(stc, 'ShowProgress')
    spbar = stc.ShowProgress;
else
    spbar = true;
end
if isfield(stc, 'ShowSRF')
    spsrf = stc.ShowSRF;
else
    spsrf = false;
end

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 5 || numel(opts) ~= 1 || ~isstruct(opts)
    opts = struct;
end
if nargin < 2 || numel(niter) ~= 1 || ~isa(niter, 'double') || isnan(niter) || ...
    niter < 1 || niter > 1e5
    niter = 1;
else
    niter = floor(niter);
end
if nargin < 3 || numel(force) ~= 1 || ~isa(force, 'double') || isnan(force) || ...
    force < -16 || force > 16 || force == 0
    force = 0.07;
end
opts.force = force;
if nargin < 4 || isempty(type) || ~ischar(type) || ...
   ~any(strcmpi(type(:)', {'even', 'smooth'}))
    opts.type = 'smooth';
else
    opts.type = lower(type(:)');
end
if isfield(opts, 'areac')
    if numel(opts.areac) == 1
        if (islogical(opts.areac) && opts.areac) || ...
           (isnumeric(opts.areac) && opts.areac ~= 0)
            opts.areac = double(0 + (opts.areac && true));
        end
    else
        opts.areac = 1;
    end
end
if isfield(opts, 'areaw') && islogical(opts.areaw) && numel(opts.areaw) == 1
    opts.areaw = double(opts.areaw);
else
    opts.areaw = 0;
end
if isfield(opts, 'distw') && islogical(opts.distw) && numel(opts.distw) == 1
    opts.distw = double(opts.distw);
else
    opts.distw = 0;
end
if isfield(opts, 'distwsq') && islogical(opts.distwsq) && numel(opts.distwsq) == 1
    opts.distwsq = double(opts.distwsq);
else
    opts.distwsq = 0;
end
if ~isfield(opts, 'norm') || numel(opts.norm) ~= 1 || ~isa(opts.norm, 'double') || ...
    isinf(opts.norm) || isnan(opts.norm) || opts.norm < -1 || opts.norm > 1
    opts.norm = 0;
end
if isfield(opts, 'normramp') && (islogical(opts.normramp) || ...
    (isa(opts.normramp, 'double') && isequal(opts.normramp, 1))) && numel(opts.normramp) == 1
    opts.normramp = double(opts.normramp);
else
    opts.normramp = 0;
end
if ~isfield(opts, 'pbar') || numel(opts.pbar) ~= 1 || ...
   (~isxfigure(opts.pbar, true) && ~isa(opts.pbar, 'xprogress'))
    opts.pbar = [];
else
    spbar = true;
end
if ~isfield(opts, 'show') || ~islogical(opts.show) || numel(opts.show) ~= 1
    opts.show = true;
end
if ~opts.show
    spsrf = false;
end
if ~isfield(opts, 'stepsize') || ~isa(opts.stepsize, 'double') || numel(opts.stepsize) ~= 1 || ...
    isinf(opts.stepsize) || isnan(opts.stepsize) || opts.stepsize < 1
    opts.stepsize = min(250, max(50, ceil(niter ^ (2 / 3))));
else
    opts.stepsize = min(niter, floor(opts.stepsize));
end
stepsize = opts.stepsize;
if ~isfield(opts, 'title') || ~ischar(opts.title) || isempty(opts.title)
    opts.title = [upper(opts.type(1)) opts.type(2:end) 'ing'];
end

% get content
xo.H.CancelMorph = false;
xo.H.Morphing = true;
xo.H.VertexCoordinateTal = [];
xo.H.VertexNormalTal = [];
if ~isfield(xo.H, 'VertexCoordinateOrig')
    xo.H.VertexCoordinateOrig = xo.C.VertexCoordinate;
    xo.H.VertexNormalOrig = xo.C.VertexNormal;
end

% set UndoBuffer
xo.C.RunTimeVars.UndoBuffer = xo.C.VertexCoordinate;

% also set auto-linked file (if not already set)
if any(strcmpi(xo.C.AutoLinkedSRF, {'', 'none', '<none>'})) && ~isempty(xo.F)
    xo.C.AutoLinkedSRF = xo.F;
end

% update
canceled = false;
sh = xo.H;
if isfield(sh, 'SUpdate') && ...
   (~isa(sh.SUpdate, 'function_handle') || numel(sh.SUpdate) ~= 1) && ...
   (~iscell(sh.SUpdate) || isempty(sh.SUpdate) || ~isa(sh.SUpdate{1}, 'function_handle'))
    xo.H = rmfield(xo.H, 'SUpdate');
    sh = xo.H;
end
bc = xo.C;
c = bc.VertexCoordinate;
n = bc.Neighbors;
[fn{1:3}] = fileparts(xo.F);
fn{2} = [fn{2} fn{3}];
if numel(fn{2}) > 28
    fn{2} = [fn{2}(1:13) '...' fn{2}(end-12:end)];
end

% determine progress bar capabilities
pbar = opts.pbar;
if spbar && isempty(pbar)
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 120, 640, 36]);
        xprogress(pbar, 'settitle', sprintf('%s %s (%d steps)...', opts.title, fn{2}, niter));
        xprogress(pbar, 0, 'Morphing...', 'visible', 0, 1);
        drawnow;
    catch xfferror
        neuroelf_lasterr(xfferror);
        spbar = false;
    end
elseif spbar
    pvis = pbar.Visible;
    pbar.Progress(0, sprintf('%s %s (%d steps)...', opts.title, fn{2}, niter));
    pbar.Visible = 'on';
    drawnow;
end

% also show SRF
if spsrf
    if ~isfield(sh, 'SUpdate') && xo.S.Extensions{1}(1) == 's'
        srf_Show(xo);
    end
    bc = xo.C;
end

% ramp normal force
if opts.norm ~= 0 && opts.normramp ~= 0
    rnorm = true;
    rnormto = opts.norm;
    rniter = niter;
else
    rnorm = false;
end

% if pbar, loop while niter >= stepsize
done = 0;
titer = niter + 40;
if ~isempty(pbar) || spsrf
    while niter >= stepsize
        opts.niter = stepsize;
        if rnorm
            opts.rampfrom = rnormto * done / rniter;
            opts.norm = rnormto * (done + stepsize) / rniter;
        end
        try
            if niter == stepsize && nargout > 1
                [c, densmv] = mesh_morph(c, n, bc.TriangleVertex, opts);
            else
                c = mesh_morph(c, n, bc.TriangleVertex, opts);
            end
        catch xfferror
            if ~isempty(pbar)
                closebar(pbar);
            end
            error('neuroelf:xff:MEXError', ...
                'Error morphing mesh (after %d pre-iterations): %s.', done, xfferror.message);
        end
        niter = niter - stepsize;
        done = done + stepsize;
        if spbar
            pbar.Progress(done / titer);
        end
        xo.H.VertexCoordinateTal = [];
        xo.H.VertexNormalTal = [];
        xo.C.VertexCoordinate = c;
        if spsrf
            srf_RecalcNormals(xo);
            if ~isfield(sh, 'SUpdate') && xo.S.Extensions{1}(1) == 's'
                srf_Show(xo);
            elseif iscell(sh.SUpdate)
                feval(sh.SUpdate{1}, 0, 0, sh.SUpdate{2:end});
            else
                feval(sh.SUpdate);
            end
            drawnow;
        end

        % cancel?
        sh = xo.H;
        if isfield(sh, 'CancelMorph') && islogical(sh.CancelMorph) && ...
            numel(sh.CancelMorph) == 1 && sh.CancelMorph

            % no more processing
            niter = 0;
            canceled = true;

            % and cancel only once!
            xo.H.CancelMorph = false;
        end
    end
end
opts.niter = niter;

% last call to morphmesh
try
    if opts.niter > 0
        if rnorm
            opts.rampfrom = rnormto * done / rniter;
            opts.norm = rnormto * (done + stepsize) / rniter;
        end
        if nargout > 1
            [c, densmv] = mesh_morph(c, n, bc.TriangleVertex, opts);
        else
            c = mesh_morph(c, n, bc.TriangleVertex, opts);
        end
    end
catch xfferror
    if ~isempty(pbar) && isempty(opts.pbar)
        closebar(pbar);
    end
    error('neuroelf:xff:MEXError', ...
        'Error morphing mesh (after %d pre-iterations): %s.', done, xfferror.message);
end

% put back
if ~canceled
    xo.C.VertexCoordinate = c;
end

% progress bar
if spbar
    pbar.Progress((titer - 40) / titer, 'Post-morph re-calculcating normals');
end

% and recalc normals after morphing
srf_RecalcNormals(xo);

% and disable morphing flag
xo.H.CancelMorph = false;
xo.H.Morphing = false;
xo.H.VertexCoordinateTal = [];
xo.H.VertexNormalTal = [];

% was canceled
if canceled && spsrf && isfield(xo.H, 'ShownInGUI') && xo.H.ShownInGUI

    % update UI
    neuroelf_gui('setsurfpos');
end

% delete figure with surface
if spsrf
    if isfield(sh, 'Surface') && numel(sh.Surface) == 1 && ishandle(sh.Surface) && isvalid(sh.Surface)
        h = sh.Surface;
        while ~strcmpi(get(h, 'Type'), 'figure')
            h = get(h, 'Parent');
            if strcmpi(get(h, 'Type'), 'root')
                warning('neuroelf:xff:internalError', ...
                    'Invalid UI handle structure. Leaving surface.');
                break;
            end
        end
        if strcmpi(get(h, 'Type'), 'figure') && strcmpi(get(h, 'Tag'), 'surface')
            delete(h);
            aft_DeleteHandle(xo, 'Surface');
        elseif isfield(sh, 'SUpdate')
            if iscell(sh.SUpdate)
                feval(sh.SUpdate{1}, 0, 0, sh.SUpdate{2:end});
            else
                feval(sh.SUpdate);
            end
            drawnow;
        end
    end
end

% produce true output map
if nargout > 1

    % progress bar
    if spbar
        pbar.Progress((titer - 25) / titer, 'Post-morph computing maps...');
    end

    % create map
    densmv = densmv.DensitySMPData;
    md = mean(densmv);
    sd = std(densmv);
    densm = xff('new:smp');
    dbc = densm.C;
    dbc.NrOfVertices = numel(densmv);
    dbc.NameOfOriginalSRF = xo.F;
    dbc.Map.LowerThreshold = max(0.01, md - 1.5 * sd);
    dbc.Map.UpperThreshold = md + 1.5 * sd;
    dbc.Map.ShowPositiveNegativeFlag = 1;
    dbc.Map.UseValuesAboveThresh = 0;
    dbc.Map.DF1 = 6;
    dbc.Map.BonferroniValue = numel(densmv);
    dbc.Map.Name = 'Mesh density map (average area per vertex)';
    dbc.Map.SMPData = densmv;
    dv = (densmv - md) / sd;
    dbc.Map(2) = dbc.Map(1);
    dbc.Map(2).LowerThreshold = 1;
    dbc.Map(2).UpperThreshold = 2.5;
    dbc.Map(2).UseValuesAboveThresh = 1;
    dbc.Map(2).DF1 = 6;
    dbc.Map(2).BonferroniValue = numel(densmv);
    dbc.Map(2).Name = 'Mesh density map (normalized average area per vertex)';
    dbc.Map(2).SMPData = dv;
    densm.C = dbc;
end

% clean up progress bar
if ~isempty(pbar) && isempty(opts.pbar)
    closebar(pbar);
elseif spbar
    pbar.Visible = pvis;
end

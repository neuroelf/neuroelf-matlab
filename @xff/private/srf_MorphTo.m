function xo = srf_MorphTo(xo, target, ivalue, niter, opts)
% SRF::MorphTo  - apply vertex morphing towards target
%
% FORMAT:       [srf] = srf.Morph(target, [ivalue, [niter, [, opts]]])
%
% Input fields:
%
%       target      xff object of type HDR, HEAD or VMR
%       ivalue      intensity value to find, default: 100
%       niter       number of iterations (default 1)
%       opts        optional settings
%        .default   either of {'in'} or 'out'
%        .force     force, value > 0 and <= 1 (default: 0.5)
%        .pbar      optional progress bar object (xfigure/xprogress, [])
%        .show      show during morphing (default: true)
%        .snmat     apply SPM-based spatial normalization
%
% Output fields:
%
%       srf         altered object
%
% Using: limitrangec, mesh_morph, minmaxmean, skew.

% Version:  v1.1
% Build:    16051916
% Date:     May-19 2016, 4:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library and global configuration
global ne_methods xffsngl;
limitrangec = ne_methods.limitrangec;
mesh_morph  = ne_methods.mesh_morph;
minmaxmean  = ne_methods.minmaxmean;
skew        = ne_methods.skew;
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
    numel(target) ~= 1 || ~xffisobject(target, true, {'hdr', 'head', 'vmr'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if isempty(target.F)
    tfn = 'anatomical dataset';
else
    [tfp, tfn] = fileparts(target.F);
end
if nargin < 3 || ~isa(ivalue, 'double') || numel(ivalue) ~= 1 || isinf(ivalue) || isnan(ivalue)
    ivalue = 100;
end
if nargin < 4 || numel(niter) ~= 1 || ~isa(niter, 'double') || isnan(niter) || ...
    niter < 1 || niter > 1e5
    niter = 1;
else
    niter = floor(niter);
end
if nargin < 5 || numel(opts) ~= 1 || ~isstruct(opts)
    opts = struct;
end
opts.areac = 0;
opts.type = 'smooth';
if ~isfield(opts, 'default') || ~ischar(opts.default) || isempty(opts.default) || ...
   ~any(strcmpi(opts.default(:)', {'in', 'out'}))
    opts.default = 'in';
else
    opts.default = lower(opts.default(:)');
end
if ~isfield(opts, 'force') || ~isa(opts.force, 'double') || numel(opts.force) ~= 1 || ...
    isinf(opts.force) || isnan(opts.force)
    opts.force = 0.5;
else
    opts.force = max(sqrt(eps), min(1, opts.force));
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
if isfield(opts, 'snmat') && islogical(opts.snmat) && numel(opts.snmat) == 1 && ...
    opts.snmat && isfield(target.C.RunTimeVars, 'SPMsn') && ...
    numel(target.C.RunTimeVars.SPMsn) == 1 && isstruct(target.C.RunTimeVars.SPMsn)
    opts.snmat = target.C.RunTimeVars.SPMsn;
end
if ~isfield(opts, 'snmat') || numel(opts.snmat) ~= 1 || ~isstruct(opts.snmat) || ...
   ~isfield(opts.snmat, 'VF') || ~isfield(opts.snmat, 'VG') || ...
   ~isfield(opts.snmat, 'Tr') || ~isfield(opts.snmat, 'Affine') || ...
   ~isstruct(opts.snmat.VF) || numel(opts.snmat.VF) ~= 1 || ...
   ~isfield(opts.snmat.VF, 'mat') || ~isa(opts.snmat.VF.mat, 'double') || ...
   ~isequal(size(opts.snmat.VF.mat), [4, 4]) || ~isstruct(opts.snmat.VG) || ...
    isempty(opts.snmat.VG) || ~isfield(opts.snmat.VG, 'dim') || ~isfield(opts.snmat.VG, 'mat') || ...
   ~isa(opts.snmat.VG(1).dim, 'double') || ~isequal(size(opts.snmat.VG(1).dim), [1, 3]) || ...
   ~isa(opts.snmat.VG(1).mat, 'double') || ~isequal(size(opts.snmat.VG(1).mat), [4, 4]) || ...
   ~isa(opts.snmat.Tr, 'double') || ndims(opts.snmat.Tr) ~= 4 || ...
   ~isa(opts.snmat.Affine, 'double') || ~isequal(size(opts.snmat.Affine), [4, 4])
    opts.snmat = [];
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
c = bc.VertexCoordinate(:, [3, 1, 2]);
cc = bc.MeshCenter(1, [3, 1, 2]);
if all(cc == cc(1))
    c = cc(1) - c;
else
    c = [cc(1) - c(:, 1), cc(2) - c(:, 2), cc(3) - c(:, 3)];
end
n = bc.VertexNormal(:, [3, 1, 2]);
[nn, ne] = srf_NeighborsNDegree(xo, 1);
nn = full(sum(ne > 0, 2));
nn3 = repmat(nn, 1, 3);
tri = bc.TriangleVertex;
oc1 = ones(size(c, 1), 1);

% determine progress bar capabilities
pbar = opts.pbar;
if spbar && isempty(pbar)
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 120, 640, 36]);
        xprogress(pbar, 'settitle', sprintf('finding intensity %d...', ivalue));
        xprogress(pbar, 0, ...
            sprintf('finding intensity %d (%d steps)...', ivalue, niter), 'visible', 0, 1);
    catch xfferror
        neuroelf_lasterr(xfferror);
        spbar = false;
    end
elseif spbar
    pvis = pbar.Visible;
    pbar.Progress(0, sprintf('Finding intensity %d in %s (%d steps)...', ivalue, tfn, niter));
    pbar.Visible = 'on';
end

% also show SRF
if spsrf
    if ~isfield(sh, 'SUpdate')
        srf_Show(xo);
    end
    bc = xo.C;
end

% if pbar, loop while niter >= stepsize
done = 0;
titer = niter;
iforce = min(10, 1 / opts.force);
maxsmooth = 10;
rcmaxsmooth = ceil([(4 / 3) .^ (2:20), titer]);
cvval = aft_GetVolume(target, 1);
cvval = minmaxmean(cvval(cvval > 0), 5);
cvval = max(1, min(sqrt(cvval(2) - cvval(1)), sqrt(cvval(6))));
while niter > 0

    % restrict normals to length 1 (in the desired direction)
    n = n ./ repmat(sqrt(sum(n .* n, 2)), 1, 3);
    if opts.default(1) ~= 'i'
        n = -n;
    end

    % smooth normal vectors
    n = full(ne * n) ./ nn3;

    % set to length 1 again (for sampling)
    nl = repmat(sqrt(sum(n .* n, 2)), 1, 3);
    nh = 0.5 .* (n ./ nl);

    % but slow down motion if not in the same direction
    n = n .* nl;

    % sample target at coordinate and +/- normal, ne1, ne2
    cval = aft_SampleData3D(target, c, struct('method', 'lanczos3', 'snmat', opts.snmat)) - ivalue;
    macval = abs(cval);
    maxval = median(macval) + std(macval);
    macval = mean(macval);
    maxval = max(maxval, macval);
    cval = (1 ./ max(1, macval)) .* cval;
    cvaln = aft_SampleData3D(target, c + nh, struct('method', 'lanczos3', 'snmat', opts.snmat));
    cvalnm = aft_SampleData3D(target, c - nh, struct('method', 'lanczos3', 'snmat', opts.snmat));
    cvaln = max(-iforce, min(iforce, cval)) .* limitrangec(abs(cvaln - cvalnm), opts.force, 1, 0);

    % smooth the forces like the normals
    cvaln = full(ne * cvaln) ./ nn;

    % determine which way to go and how far
    v = n .* cvaln(:, [1, 1, 1]);
    v(any(isinf(v) | isnan(v), 2), :) = 0;

    % then add to c (as well as inversely weighted n)
    c = c + opts.force .* v;

    % for out direction
    if opts.default(1) ~= 'i'

        % estimate volume as overall std
        vol = std(c, [], 1);
    end

    % slight smoothing between iterations
    sforce = max(0.05, 0.25 * opts.force * 0.1 * maxsmooth);
    sniter = ceil(min(opts.force .* niter, ...
        max([10 / (opts.force ^ 2), maxsmooth, 2 * maxval / cvval])));
    c = mesh_morph(c, bc.Neighbors, tri, struct('areac', 1, ...
        'distwsq', 1, 'distc', 1 + opts.force, 'force', sforce, 'niter', sniter));

    % match again for center and area?
    if opts.default(1) ~= 'i'
        nvol = std(c, [], 1);

        % then patch coordinates
        c = c .* (oc1 * (vol ./ nvol));
    end

    % set back in original mesh
    if all(cc == cc(1))
        cn = cc(1) - c(:, [2, 3, 1]);
    else
        cn = [cc(2) - c(:, 2), cc(3) - c(:, 3), cc(1) - c(:, 1)];
    end
    xo.C.VertexCoordinate = cn;
    xo.H.VertexCoordinateTal = [];
    xo.H.VertexNormalTal = [];

    % more to do
    srf_RecalcNormals(xo);
    n = xo.C.VertexNormal(:, [3, 1, 2]);

    % also recompute smoothing factor
    if any(rcmaxsmooth == niter)
        smp = srf_CurvatureMap(xo);
        smpc = smp.C;
        delete(smp);
        smpc = double(smpc.Map.SMPData);
        smp = srf_DensityMap(xo, struct('dist', false, 'mindist', false));
        smpd = smp.C;
        delete(smp);
        smpd = double(smpd.Map(1).SMPData);
        maxsmooth = min(15, ceil(mean(abs(smpc)) + 1.5 .* std(smpc) + 2 .* abs(skew(smpc))));
        maxsmooth = maxsmooth + ...
            min(15, ceil(mean(abs(smpd)) + 1.5 .* std(smpd) + 2 .* abs(skew(smpd))));
    end

    % decrease inter and increase done
    niter = niter - 1;
    done = done + 1;

    % report progress
    if spbar
        pbar.Progress(done / titer, ...
            sprintf('Finding intensity %d in %s (%d/%d steps, %.2f, %d)...', ...
            ivalue, tfn, done, titer, sforce, sniter), false);
    end

    % show update
    if spsrf
        if ~isfield(sh, 'SUpdate')
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

% clean up progress bar
if ~isempty(pbar) && isempty(opts.pbar)
    closebar(pbar);
elseif spbar
    pbar.Visible = pvis;
end

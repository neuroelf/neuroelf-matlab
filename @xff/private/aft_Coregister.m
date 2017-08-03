function m44 = aft_Coregister(xo, xo2, opts)
% AFT::Coregister  - coregister two objects
%
% FORMAT:       [m44 = ] obj.Coregister(obj2 [, opts])
%
% Input fields:
%
%       obj2        second xff object
%       opts        optional settings
%        .bbox      bounding box (TAL notation, default: -127 .. 128)
%        .grad      use gradient information instead (default: false)
%        .norm      normalize volumes to mid range (default: true)
%        .normrng   normalization range (default: 0.95)
%        .grid      resolution gridsize (mm, default: [4, 2])
%        .smoothsrc pre-smoothing kernel of source data in mm (default: 8)
%        .smoothtrg pre-smoothing kernel of target data in mm (default: 8)
%        .srcvol    source volume (default: 1)
%        .trgvol    target volume (default: 1)
%
% Output fields:
%
%       m44         4x4 transformation matrix
%
% TYPES: HDR, HEAD, MGH, VMR, VTC
%
% Using: findfirst, flexinterpn, histcount, minmaxmean, rbalign, smoothdata3.

% Version:  v1.1
% Build:    16031615
% Date:     Mar-16 2016, 3:32 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011 - 2016, Jochen Weber
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

% neuroelf library
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'hdr', 'head', 'vmr', 'vtc'}) || ...
    numel(xo2) ~= 1 || ~xffisobject(xo2, true, {'hdr', 'head', 'vmr', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bbox') || ~isa(opts.bbox, 'double') || ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:)) | opts.bbox(:) < -256 | opts.bbox(:) > 256)
    opts.bbox = [-127, -127, -127; 128, 128, 128];
end
if ~isfield(opts, 'grad') || ~islogical(opts.grad) || numel(opts.grad) ~= 1
    opts.grad = false;
end
if ~isfield(opts, 'grid') || ~isa(opts.grid, 'double') || isempty(opts.grid) || numel(opts.grid) > 4 || ...
    any(isinf(opts.grid(:)) | isnan(opts.grid(:)) | opts.grid(:) <= 0 | opts.grid(:) > 64)
    opts.grid = [4, 2];
end
if ~isfield(opts, 'norm') || ~islogical(opts.norm) || numel(opts.norm) ~= 1
    opts.norm = true;
end
if ~isfield(opts, 'normrng') || ~isa(opts.normrng, 'double') || numel(opts.normrng) ~= 1 || ...
    isinf(opts.normrng) || isnan(opts.normrng) || opts.normrng < 0.25 || opts.normrng > 1
    opts.normrng = 0.95;
end
if ~isfield(opts, 'smoothsrc') || ~isa(opts.smoothsrc, 'double') || ~any(numel(opts.smoothsrc) == [1, numel(opts.grid)]) || ...
    isinf(opts.smoothsrc) || isnan(opts.smoothsrc) || opts.smoothsrc < 0
    opts.smoothsrc = 8;
else
    opts.smoothsrc = min(opts.smoothsrc, 32);
end
if numel(opts.smoothsrc) ~= numel(opts.grid)
    opts.smoothsrc = opts.smoothsrc .* ones(1, numel(opts.grid));
end
if ~isfield(opts, 'smoothtrg') || ~isa(opts.smoothtrg, 'double') || ...
   ~any(numel(opts.smoothtrg) == [1, numel(opts.grid)]) || ...
    isinf(opts.smoothtrg) || isnan(opts.smoothtrg) || opts.smoothtrg < 0
    opts.smoothtrg = 8;
else
    opts.smoothtrg = min(opts.smoothtrg, 32);
end
if numel(opts.smoothtrg) ~= numel(opts.grid)
    opts.smoothtrg = opts.smoothtrg .* ones(1, numel(opts.grid));
end
if ~isfield(opts, 'srcvol') || ~isa(opts.srcvol, 'double') || numel(opts.srcvol) ~= 1 || ...
    isinf(opts.srcvol) || isnan(opts.srcvol) || opts.srcvol < 1 || opts.srcvol ~= fix(opts.srcvol)
    opts.srcvol = 1;
end
if ~isfield(opts, 'trgvol') || ~isa(opts.trgvol, 'double') || numel(opts.trgvol) ~= 1 || ...
    isinf(opts.trgvol) || isnan(opts.trgvol) || opts.trgvol < 1 || opts.trgvol ~= fix(opts.trgvol)
    opts.trgvol = 1;
end

% get data handle
bc1 = xo.C;

% ensure main object has RunTimeVars.TrfPlus field
if ~isfield(bc1.RunTimeVars, 'TrfPlus') || ~isa(bc1.RunTimeVars.TrfPlus, 'double') || ...
   ~isequal(size(bc1.RunTimeVars.TrfPlus), [4, 4]) || ...
    any(isinf(bc1.RunTimeVars.TrfPlus(:)) | isnan(bc1.RunTimeVars.TrfPlus(:))) || ...
    any(bc1.RunTimeVars.TrfPlus(4, :) ~= [0, 0, 0, 1])
    bc1.RunTimeVars.TrfPlus = eye(4);
    xo.C = bc;
end

% initialize output
m44 = eye(4);
seps = sqrt(eps);

% iterate over each grid size
for cc = 1:numel(opts.grid)

    % bounding box
    bbox = struct('BBox', opts.bbox, 'ResXYZ', opts.grid(cc) .* ones(1, 3));

    % get volumes
    vs = aft_SampleTalBox(xo, bbox, opts.srcvol);
    vt = aft_SampleTalBox(xo2, bbox, opts.trgvol);

    % smoothing
    if opts.smoothsrc > 0
        vs = ne_methods.smoothdata3(vs, (opts.smoothsrc(cc) / opts.grid(cc)) .* ones(1, 3));
    end
    if opts.smoothtrg > 0
        vt = ne_methods.smoothdata3(vt, (opts.smoothtrg(cc) / opts.grid(cc)) .* ones(1, 3));
    end

    % gradient
    if opts.grad
        flexinterpn = ne_methods.flexinterpn;
        vsg = vs;
        vtg = vt;
        fib = [inf, inf, inf; 1, 1, 1; 1, 1, 1; size(vsg)];
        fx = zeros(4, 3);
        fy = zeros(4, 3);
        fz = zeros(4, 3);
        fx([2, 4], 1) = 0.125;
        fy([2, 4], 2) = 0.125;
        fz([2, 4], 3) = 0.125;
        for zc = 1:size(vsg, 3)
            fib([2, 4], 3) = zc;
            vsg(:, :, zc) = (flexinterpn(vs, fib + fx) - flexinterpn(vs, fib - fx)) .^ 2 + ...
                (flexinterpn(vs, fib + fy) - flexinterpn(vs, fib - fy)) .^ 2 + ...
                (flexinterpn(vs, fib + fz) - flexinterpn(vs, fib - fz)) .^ 2;
            vtg(:, :, zc) = (flexinterpn(vt, fib + fx) - flexinterpn(vt, fib - fx)) .^ 2 + ...
                (flexinterpn(vt, fib + fy) - flexinterpn(vt, fib - fy)) .^ 2 + ...
                (flexinterpn(vt, fib + fz) - flexinterpn(vt, fib - fz)) .^ 2;
        end

        % copy back
        vs = vsg;
        vt = vtg;
    end

    % normalize volumes
    if opts.norm

        % get values that are good
        gs = ~(isinf(vs) | isnan(vs) | vs == 0);
        gt = ~(isinf(vt) | isnan(vt) | vt == 0);
        ngs = sum(gs(:));
        ngt = sum(gt(:));

        % compute histograms
        mms = ne_methods.minmaxmean(vs(gs));
        mms(2) = mms(2) + seps;
        hs = ne_methods.histcount(vs(gs), mms(1) - seps, mms(2), (mms(2) - mms(1)) * 0.0001);
        hsc = cumsum(hs);
        mmt = ne_methods.minmaxmean(vt(gt));
        mmt(2) = mmt(2) + seps;
        ht = ne_methods.histcount(vt(gt), mmt(1) - seps, mmt(2), (mmt(2) - mmt(1)) * 0.0001);
        htc = cumsum(ht);

        % get mid range percentile values
        lr = 0.5 * opts.normrng;
        ur = 1 - lr;
        rs1 = mms(1) + (ne_methods.findfirst(hsc >= (lr * ngs)) - 1) * 0.0001 * (mms(2) - mms(1));
        rs2 = mms(1) + ne_methods.findfirst(hsc >= (ur * ngs)) * 0.0001 * (mms(2) - mms(1)) + seps;
        rt1 = mmt(1) + (ne_methods.findfirst(htc >= (lr * ngt)) - 1) * 0.0001 * (mmt(2) - mmt(1));
        rt2 = mmt(1) + ne_methods.findfirst(htc >= (ur * ngt)) * 0.0001 * (mmt(2) - mmt(1)) + seps;

        % z-transform according to values
        mms = ne_methods.minmaxmean(vs(gs & vs >= rs1 & vs <= rs2), 1);
        vs = (1 / sqrt(mms(6))) .* (vs - mms(3));
        vs(~gs) = 0;
        mmt = ne_methods.minmaxmean(vt(gt & vt >= rt1 & vt <= rt2), 1);
        vt = (1 / sqrt(mmt(6))) .* (vt - mmt(3));
        vt(~gt) = 0;
    end

    % run rbalign
    trfv = eye(4);
    trfv([1, 6, 11]) = opts.grid(cc);
    trfv(1:3, 4) = -(127 + opts.grid(cc));
    opts.smooth = [0, 0, 0];
    opts.smpl = opts.grid(cc) .* ones(1, 3);
    opts.trfv1 = trfv;
    opts.trfv2 = trfv;
    tfm = ne_methods.rbalign(vt, vs, opts) * pinv(trfv);

    % multiply
    m44 = tfm * m44;
    bc1.RunTimeVars.TrfPlus = tfm * bc1.RunTimeVars.TrfPlus;

    % update handle
    xo.C = bc1;
end

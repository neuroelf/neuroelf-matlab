function [xo, par, trf] = fmr_Realign(xo, opts)
% FMR::Realign  - perform spatial realignment
%
% FORMAT:       [fmr, par, trf] = fmr.Realign([opts])
%
% Input fields:
%
%       opts        optional settings
%        .interpe   estim interpolation, {'linear'}, 'cubic', 'lanczos3'
%        .interpr   reslice interpolation, 'linear', 'cubic', {'lanczos3'}
%        .mask      1x1 double [0...2], relative masking threshold (0.5)
%        .robust    use robust regression to detect motion (default: false)
%        .rtplot    real-time plot of params (default: true)
%        .savemean  boolean flag, save mean as one-vol FMR (def: false)
%        .savepar   boolean flag, save parameters (default: true)
%        .smooth    smoothing kernel in mm (default: twice the voxelsize)
%        .smpl      sampling width in mm (either 1x1 or 1x3)
%        .tomean    two-pass realignment (if value > 2 perform N passes)
%        .totarget  if given, must be either a volume of data (double) or a
%                   1x2 cell with FMR object (filename) and volume number
%
% Output fields:
%
%       fmr         FMR with realigned data
%       par         Vx6 parameters of realignment (same as from .savepar)
%       trf         4x4xV transformation parameters used from resampling
%                   (see resampling code in rbalign.m)
%
% Using: clustercoordsc, dilate3d, erode3d, maxpos, rbalign, spmitrf, ztrans.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:24 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library
global ne_methods;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'interpe') || ~ischar(opts.interpe) || ...
   ~any(strcmpi(opts.interpe(:)', {'linear', 'cubic', 'lanczos3', 'lanczos5'}))
    opts.interpe = 'linear';
end
if ~isfield(opts, 'interpr') || ~ischar(opts.interpr) || ...
   ~any(strcmpi(opts.interpr(:)', {'linear', 'cubic', 'lanczos3', 'lanczos5'}))
    opts.interpr = 'lanczos3';
end
if ~isfield(opts, 'mask') || ~isa(opts.mask, 'double') || numel(opts.mask) ~= 1 || ...
    isinf(opts.mask) || isnan(opts.mask)
    opts.mask = 0.5;
else
    opts.mask = min(2, max(0, opts.mask));
end
if ~isfield(opts, 'robust') || ~islogical(opts.robust) || numel(opts.robust) ~= 1
    opts.robust = false;
end
if ~isfield(opts, 'rtplot') || ~islogical(opts.rtplot) || numel(opts.rtplot) ~= 1
    opts.rtplot = true;
end
if ~isfield(opts, 'savemean') || ~islogical(opts.savemean) || numel(opts.savemean) ~= 1
    opts.savemean = false;
end
if ~isfield(opts, 'savepar') || ~islogical(opts.savepar) || numel(opts.savepar) ~= 1
    opts.savepar = true;
end
if ~isfield(opts, 'smooth') || ~isa(opts.smooth, 'double') || ~any(numel(opts.smooth) == [1, 3]) || ...
    any(isinf(opts.smooth) | isnan(opts.smooth) | opts.smooth <= 1)
    opts.smooth = [];
elseif numel(opts.smooth) == 1
    opts.smooth = opts.smooth([1, 1, 1]);
end
if ~isfield(opts, 'smpl') || ~isa(opts.smpl, 'double') || ~any(numel(opts.smpl) == [1, 3]) || ...
    any(isinf(opts.smpl) | isnan(opts.smpl) | opts.smpl <= 1)
    opts.smpl = [];
elseif numel(opts.smpl) == 1
    opts.smpl = opts.smpl([1, 1, 1]);
end
if ~isfield(opts, 'tomean') || (~isa(opts.tomean, 'double') && ...
    ~islogical(opts.tomean)) || numel(opts.tomean) ~= 1
    opts.tomean = 0;
else
    if islogical(opts.tomean)
        opts.tomean = double(opts.tomean);
    end
    if isinf(opts.tomean) || isnan(opts.tomean) || opts.tomean < 0  || opts.tomean > 8
        opts.tomean = 0;
    end
end
if opts.tomean < 2
    opts.tomean = opts.tomean + 1;
end
opts.tomean = round(opts.tomean);
if ~isfield(opts, 'totarget')
    opts.totarget = [];
elseif isnumeric(opts.totarget)
    opts.totarget = double(opts.totarget);
end
if ~isempty(opts.totarget) && iscell(opts.totarget) && numel(opts.totarget) == 2 && ...
    isa(opts.totarget{2}, 'double') && numel(opts.totarget{2}) == 1 && ...
   ~isinf(opts.totarget{2}) && ~isnan(opts.totarget{2}) && opts.totarget{2} >= 1 && ...
    opts.totarget{2} == fix(opts.totarget{2})
    if ischar(opts.totarget{1})
        refo = {[]};
        try
            refo{1} = xff(opts.totarget{1}(:)');
            if ~xffisobject(refo{1}, true, 'fmr')
                error('NO_FMR');
            end
            opts.totarget = double(aft_GetVolume(refo{1}, opts.totarget{2}));
        catch xfferror
            neuroelf_lasterr(xfferror);
            warning('neuroelf:xff:FMRNotReadable', 'Target FMR not readable: %s.', xfferror.message);
            opts.totarget = [];
        end
        clearxffobjects(refo);
    elseif numel(opts.totarget{1}) == 1 && xffisobject(opts.totarget{1}, true, 'fmr')
        try
            opts.totarget = double(aft_GetVolume(opts.totarget{1}, opts.totarget{2}));
        catch xfferror
            neuroelf_lasterr(xfferror);
            warning('neuroelf:xff:FMRNotReadable', 'Target FMR not readable: %s.', xfferror.message);
            opts.totarget = [];
        end
    end
end
if ~isempty(opts.totarget) && ~isa(opts.totarget, 'double')
    opts.totarget = [];
end

% get data
if bc.FileVersion < 5 || bc.DataStorageFormat < 2
    stcd = bc.Slice(1).STCData(:, :, :);
    sstc = size(stcd);
    sstc = [sstc(1:2), 1, sstc(3)];
    stcd = reshape(stcd, sstc);
    stcd(1, 1, numel(bc.Slice), 1) = 0;
    for sc = 2:numel(bc.Slice)
        stcd(:, :, sc, :) = reshape(bc.Slice(sc).STCData(:, :, :), sstc);
    end
else
    stcd = permute(bc.Slice.STCData(:, :, :, :), [1, 2, 4, 3]);
end
dtype = class(stcd);
stcd = single(stcd);

% compile options
m = [ ...
    bc.InplaneResolutionX, 0, 0, -0.5 * bc.InplaneResolutionX * (bc.ResolutionX + 1); ...
    0, bc.InplaneResolutionY, 0, -0.5 * bc.InplaneResolutionY * (bc.ResolutionY + 1); ...
    0, 0, bc.SliceThickness + bc.GapThickness, -0.5 * (bc.SliceThickness + bc.GapThickness) * (bc.NrOfSlices + 1); ...
    0, 0, 0, 1];
rbopts = struct('interpe', lower(opts.interpe(:)'), 'interpr', lower(opts.interpr(:)'), ...
    'mask', [], 'robust', opts.robust, 'rtplot', opts.rtplot, 'smooth', opts.smooth, ...
    'smpl', opts.smpl, 'trfv1', m, 'trfv2', m, 'tsmooth', opts.smooth);

% masking
if opts.mask > 0
    stcm = (1 / size(stcd, 4)) .* sum(stcd, 4);
    stcm = ne_methods.erode3d(stcm > (opts.mask .* mean(stcm(:))));
    [stcs, stcm] = ne_methods.clustercoordsc(stcm, 1, floor(sqrt(numel(stcm))));
    stcm = (stcm == ne_methods.maxpos(stcs));
    stcm = (ne_methods.dilate3d(ne_methods.dilate3d(stcm)) & ~ne_methods.erode3d(ne_methods.erode3d(stcm)));
    rbopts.mask = stcm;
end

% perform realignment
tgt = stcd(:, :, :, 1);

% other target
if ~isempty(opts.totarget) && isequal(size(tgt), size(opts.totarget))
    tgt = opts.totarget;
end

% perform realignment
for pc = 1:opts.tomean
    [trf, tri, ra] = ne_methods.rbalign(tgt, stcd, rbopts);
    tgt = sum(ra, 4) ./ sum(ra ~= 0, 4);
    tgt(isinf(tgt) | isnan(tgt)) = 0;
end

% put data back into format
if bc.FileVersion < 6 || ~strcmpi(dtype, 'single')
    bc.DataType = 1;
    ra = uint16(round(max(ra, 0)));
end
if bc.FileVersion < 5 || bc.DataStorageFormat < 2
    sstc = size(bc.Slice(1).STCData);
    for sc = 1:numel(bc.Slice)
        bc.Slice(sc).STCData = reshape(ra(:, :, sc, :), sstc);
    end
else
    bc.Slice.STCData = permute(single(ra), [1, 2, 4, 3]);
end

% but content into array
xo.C = bc;

% parameters
if nargout > 1 || opts.savepar

    % get filename
    f = xo.F;

    % compile paramters
    par = zeros(size(tri, 3), 6);
    for vc = 1:size(par, 1)
        p = ne_methods.spmitrf(tri(:, :, vc));
        par(vc, :) = [p{1}, (180 / pi) * p{2}];
    end

    % save params?
    if opts.savepar

        % create SDM
        sdm = xff('new:sdm');
        sdmc = sdm.C;

        % put into SDM
        sdmc.NrOfPredictors = 6;
        sdmc.NrOfDataPoints = size(par, 1);
        sdmc.IncludesConstant = 0;
        sdmc.FirstConfoundPredictor = 1;
        sdmc.PredictorColors = [255, 0, 0; 0, 255, 0; 0, 0, 255; ...
            255, 255, 0; 0, 255, 255; 255, 0, 255];
        sdmc.PredictorNames = {'trans_x_mm', 'trans_y_mm', 'trans_z_mm', ...
            'rot_x_deg', 'rot_y_deg', 'rot_z_deg'};
        sdmc.SDMMatrix = par;
        sdmc.RTCMatrix = zeros(size(par, 1), 0);
        sdm.C = sdmc;
        aft_SaveAs(sdm, [f(1:end-4) '_MCparams.sdm']);

        % also create z-transformed and squared version
        sdmc.NrOfPredictors = 12;
        sdmc.PredictorNames = {'ztrans_x_mm', 'ztrans_y_mm', 'ztrans_z_mm', ...
            'zrot_x_deg', 'zrot_y_deg', 'zrot_z_deg', ...
            'SQztrans_x_mm', 'SQztrans_y_mm', 'SQztrans_z_mm', ...
            'SQzrot_x_deg', 'SQzrot_y_deg', 'SQzrot_z_deg'};
        sdmc.PredictorColors = [255, 0, 0; 0, 255, 0; 0, 0, 255; ...
            255, 255, 0; 0, 255, 255; 255, 0, 255; ...
            192, 128, 128; 128, 192, 128; 128, 128, 192; ...
            192, 192, 128; 128, 192, 192; 192, 128, 192];
        sdmc.SDMMatrix = [ne_methods.ztrans(par), ne_methods.ztrans(par .* par)];
        sdm.C = sdmc;
        aft_SaveAs(sdm, [f(1:end-4) '_MCzparams.sdm']);
        delete(sdm);
    end
end

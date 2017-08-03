function cvmp = vtc_Correlate(xo, regr, opts)
% VTC::Correlate  - create correlation map(s) for regressor(s)
%
% FORMAT:       cvmp = vtc.Correlate(regr [, opts])
%
% Input fields:
%
%       regr        Tx1 or TxR regressor(s)
%       opts        options settings
%        .nuisreg   nuisance regressors (variance removed from regr/data)
%        .regnames  regressor names
%        .stat      either of {'r'} or 't'
%        .tfiltfrq  number of temporal filtering frequencies (default: 0)
%        .tfilttyp  filtering type, either of {'DCT'}, 'Fourier'
%        .trobust   perform regressions (incl. filtering) robustly
%
% Output fields:
%
%       cvmp        correlation r-VMP (with 1/R maps)
%
% Note: the toolbox internal cov_nd function is used which gives
%       slightly different r values than corrcoef.
%
% Using: applyfdr, correlinvtstat, correltstat, cov_nd,
%        fitrobustbisquare_img, newnatresvmp, robustt, tempfilter, varc,
%        ztrans.

% Version:  v1.1
% Build:    16021412
% Date:     Feb-14 2016, 12:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% importing from neuroelf library
using(neuroelf, {'applyfdr', 'correlinvtstat', 'correltstat', 'cov_nd', ...
    'fitrobustbisquare_img', 'newnatresvmp', 'robustt', 'tempfilter', ...
    'varc', 'ztrans'});

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc') || ...
   ~isa(regr, 'double') || isempty(regr) || ndims(regr) > 2 || ...
    any(isinf(regr(:)) | isnan(regr(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if size(regr, 1) ~= size(bc.VTCData, 1)
    error('neuroelf:xff:badArgument', 'NrOfVolumes must match size(regr, 1).');
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'nuisreg') || isempty(opts.nuisreg)
    opts.nuisreg = [];
end
if ~isfield(opts, 'regnames') || ~iscell(opts.regnames) || numel(opts.regnames) ~= size(regr, 2)
    opts.regnames = cell(1, size(regr, 2));
else
    opts.regnames = opts.regnames(:)';
end
for rc = 1:numel(opts.regnames)
    if ~ischar(opts.regnames{rc}) || isempty(opts.regnames{rc})
        opts.regnames{rc} = sprintf('regressor %d', rc);
    end
end
if ~isfield(opts, 'stat') || ~ischar(opts.stat) || isempty(opts.stat) || ...
   ~any('rt' == lower(opts.stat(1)))
    opts.stat = 'r';
else
    opts.stat = lower(opts.stat(1));
end
if ~isfield(opts, 'tfiltfrq') || numel(opts.tfiltfrq) ~= 1 || ~isa(opts.tfiltfrq, 'double') || ...
    isinf(opts.tfiltfrq) || opts.tfiltfrq < 0 || opts.tfiltfrq > 12
    opts.tfiltfrq = 0;
end
if ~isfield(opts, 'tfilttyp') || ~ischar(opts.tfilttyp) || ...
   ~any(strcmpi(opts.tfilttyp(:)', {'dct', 'fourier'}))
    opts.tfilttyp = 'dct';
else
    opts.tfilttyp = lower(opts.tfilttyp(:)');
end
if ~isfield(opts, 'trobust') || ~islogical(opts.trobust) || numel(opts.trobust) ~= 1
    opts.trobust = false;
end
fdrl = [0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001];

% get data
vd = bc.VTCData;
if istransio(vd)
    vd = resolve(vd);
end

% filter content
if opts.tfiltfrq > 0 || ~isempty(opts.nuisreg)

    % prepare tempfilter options
    topts = opts;
    topts.nuisreg = opts.nuisreg;
    topts.spat = false;
    topts.tdim = 1;
    topts.temp = true;
    topts.tempdt = false;
    topts.robust = opts.trobust;
    if opts.tfilttyp(1) == 'd'
        if opts.tfiltfrq > 0
            topts.tempdct = ceil(0.001 * bc.TR * size(bc.VTCData, 1) / opts.tfiltfrq);
        else
            topts.tempdct = Inf;
        end
        topts.tempsc = 0;
    else
        topts.tempdct = Inf;
        topts.tempsc = opts.tfiltfrq;
    end

    % temp filter data of first object
    vd = tempfilter(vd, topts);

    % and temp filter regressors
    regr = tempfilter(regr, topts);
end

% for correlation, use z-trans on both data and regressors
regr = ztrans(regr);

% compute voxels to use
vvd = varc(vd);
usev = (vvd > 0 & ~isinf(vvd) & ~isnan(vvd));
usev = find(usev(:));
nusev = numel(usev);

% permute data
vds = size(vd);
nvol = vds(1);
vds(1) = [];
if nusev ~= prod(vds)
    vd = vd(:, usev)';
else
    vd = reshape(vd, nvol, nusev)';
end

% create vmp
bb = aft_BoundingBox(xo);
cvmp = newnatresvmp(bb.BBox, bc.Resolution, 2);
bcvmp = cvmp.C;
bcvmp.OriginatingVTC = xo.F;
if ischar(bc.NameOfLinkedPRT)
    bcvmp.LinkedPRT = bc.NameOfLinkedPRT(:)';
elseif iscell(bc.NameOfLinkedPRT) && ~isempty(bc.NameOfLinkedPRT) && ischar(bc.NameOfLinkedPRT{1})
    bcvmp.LinkedPRT = bc.NameOfLinkedPRT{1}(:)';
end
[p, f1] = fileparts(xo.F);
if opts.stat == 't'
    bcvmp.Map.Type = 1;
else
    bcvmp.Map.Type = 2;
end
bcvmp.Map.Name = sprintf('Corr %s <-> %s', f1, opts.regnames{1});
bcvmp.Map.DF1 = nvol - 2;
bcvmp.Map.BonferroniValue = nusev;
bcvmp.Map.NrOfFDRThresholds = numel(fdrl);
bcvmp.Map = bcvmp.Map(1, ones(1, size(regr, 2)));

% iterate over regressors
for rc = 1:numel(opts.regnames)

    % build tmap
    tmap = zeros(vds);

    % non-robust
    if ~opts.trobust

        % cov_nd
        [cnd, rt] = cov_nd(vd, repmat(regr(:, rc)', nusev, 1));
        rt(isinf(rt) | isnan(rt)) = 0;

        % convert ?
        if opts.stat == 't'
            rt = correltstat(rt, nvol);
        end

    % robust
    else

        % build model
        md = [regr(:, rc), ones(nvol, 1)];

        % use fitrobustbisquare_img and robustt
        [b, w] = fitrobustbisquare_img(md, double(vd));
        rt = robustt(md, vd(dc, :, :, :), b, w, [1, 0]);
        rt = reshape(rt, [1, vds(2:3)]);
        rt(isinf(rt) | isnan(rt)) = 0;

        % convert ?
        if opts.stat == 'r'
            rt = correlinvtstat(rt, nvol);
        end
    end

    % put into tmap
    tmap(usev) = rt;

    % set to map
    bcvmp.Map(rc).Name = sprintf('Corr %s <-> %s', f1, opts.regnames{rc});
    bcvmp.Map(rc).VMPData = single(tmap);
    fdrt2 = applyfdr(double(bcvmp.Map(rc).VMPData), opts.stat, fdrl(:), ...
        bcvmp.Map(rc).DF1, 0, 3);
    bcvmp.Map(rc).FDRThresholds = [fdrl(:), fdrt2];
    bcvmp.Map(rc).LowerThreshold = fdrt2(end);
    bcvmp.Map(rc).UpperThreshold = 2 * fdrt2(end);
end

% set content
cvmp.C = bcvmp;

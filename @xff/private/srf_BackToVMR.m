function xo = srf_BackToVMR(xo, opts)
% SRF::BackToVMR  - back-project vertices to VMR space
%
% FORMAT:       [srf = ] srf.BackToVMR([opts])
%
% Input fields:
%
%       opts        optional settings
%        .fillmode  either of {'nearest'} or 'linear'
%        .fillvmax  maximal filling value (combination, default: 200)
%        .nfrom     Vertex + from * Normal (default: -0.25)
%        .nstep     step along normal(default: 0.5)
%        .nto       Vertex + to * Normal (default: nfrom + nstep)
%        .res       VMR resolution, either of {1} or 0.5
%        .smp       SMP map object (if different marking is requested)
%        .smpmap    sub-map of SMP (default: 1)
%        .tcode     target color code (default: 235)
%        .triovsmp  triangular oversampling factor (default: 2)
%        .vmr       insert data into VMR (after sampling)
%
% Output fields:
%
%       srf         SRF object with VertexVMRData set accordingly
%
% Note: for BV-objects, normals point inward!

% Version:  v1.1
% Build:    16031614
% Date:     Mar-16 2016, 2:03 PM EST
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsbf', 'srf'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'fillmode') || ~ischar(opts.fillmode) || ...
   ~any(strcmpi(opts.fillmode(:)', {'linear', 'nearest'}))
    opts.fillmode = 'nearest';
else
    opts.fillmode = 'linear';
end
if ~isfield(opts, 'fillvmax') || numel(opts.fillvmax) ~= 1 || ~isnumeric(opts.fillvmax) || ...
    isinf(opts.fillmax) || isnan(opts.fillmax)
    opts.fillvmax = 200;
else
    opts.fillvmax = min(255, max(0, (double(opts.fillvmax))));
end
if ~isfield(opts, 'nfrom') || ~isa(opts.nfrom, 'double') || numel(opts.nfrom) ~= 1 || ...
    isinf(opts.nfrom) || isnan(opts.nfrom)
    opts.nfrom = -0.25;
else
    opts.nfrom = min(15, max(-15, opts.nfrom));
end
if ~isfield(opts, 'nstep') || ~isa(opts.nstep, 'double') || numel(opts.nstep) ~= 1 || ...
    isinf(opts.nstep) || isnan(opts.nstep) || opts.nstep < 0.1 || opts.nstep > 4
    opts.nstep = 0.5;
end
if ~isfield(opts, 'nto') || ~isa(opts.nto, 'double') || numel(opts.nto) ~= 1 || ...
    isinf(opts.nto) || isnan(opts.nto) || opts.nto < opts.nfrom
    opts.nto = opts.nfrom + opts.nstep;
else
    opts.nto = min(15, max(opts.nfrom + sqrt(eps), opts.nto));
end
if ~isfield(opts, 'res') || ~isa(opts.res, 'double') || numel(opts.res) ~= 1 || ...
    isinf(opts.res) || isnan(opts.res) || opts.res ~= 0.5
    opts.res = 1;
else
    opts.res = 2;
end
if ~isfield(opts, 'smp') || numel(opts.smp) ~= 1 || ~xffisobject(opts.smp, true, 'smp')
    opts.smp = [];
else
    opts.smp = opts.smp.C;
    if isempty(opts.smp.Map)
        opts.smp = [];
    end
end
if ~isempty(opts.smp) && (~isfield(opts, 'smpmap') || ~isa(opts.smpmap, 'double') || ...
    numel(opts.smpmap) ~= 1 || isinf(opts.smpmap) || isnan(opts.smpmap) || ...
    opts.smpmap < 1 || opts.smpmap > numel(opts.smp.Map))
    opts.smpmap = 1;
elseif ~isempty(opts.smp)
    opts.smpmap = round(opts.smpmap);
end
if ~isfield(opts, 'tcode') || ~isa(opts.tcode, 'double') || numel(opts.tcode) ~= 1 || ...
    isinf(opts.tcode) || isnan(opts.tcode) || opts.tcode < 0 || opts.tcode > 255
    opts.tcode = 235;
else
    opts.tcode = round(opts.tcode);
end
if ~isfield(opts, 'triovsmp') || ~isa(opts.triovsmp, 'double') || numel(opts.triovsmp) ~= 1 || ...
    isinf(opts.triovsmp) || isnan(opts.triovsmp)
    opts.triovsmp = 2;
else
    opts.triovsmp = round(min(10, max(1, opts.triovsmp)));
end
if ~isfield(opts, 'vmr') || numel(opts.vmr) ~= 1 || ~xffisobject(opts.vmr, true, 'vmr')
    opts.vmr = [];
end
o = 0:1/opts.triovsmp:1.01;
ores = opts.res;

% get content
bc = xo.C;
c = bc.VertexCoordinate;
n = bc.VertexNormal;
if lower(xo.S.Extensions{1}(1)) == 'f'
    c = 128 - c(:, [2, 3, 1]);
    n = -n(:, [2, 3, 1]);
end
t = bc.TriangleVertex;

% check SMP
tcode = opts.tcode;
if ~isempty(opts.smp)
    if numel(opts.smp.Map(opts.smpmap).SMPData) ~= size(c, 1)
        opts.smp = [];
    else
        opts.smp = opts.smp.Map(opts.smpmap);
        lt = opts.smp.LowerThreshold;
        ut = opts.smp.UpperThreshold;
        if ut <= lt
            ut = lt + sqrt(eps);
        end
        opts.smp = max(0, min(9, (9 / (ut - lt)) .* (opts.smp.SMPData - lt)));
        tcode = 240;
    end
end

% create VMRData
if strcmpi(opts.fillmode, 'nearest')
    vmrd = uint8(0);
else
    vmrd = single(0);
end
if opts.res == 1
    vmrd(1:256, 1:256, 1:256) = vmrd;
    mxs = 256;
else
    vmrd(1:512, 1:512, 1:512) = vmrd;
    mxs = 512;
end

% keep the interface responsive
nupi = 1 / 172800;
nup = now + nupi;

% iterate over range
for sc = opts.nfrom:opts.nstep:opts.nto

    % for nearest
    if strcmpi(opts.fillmode, 'nearest')

        % get coordinate (+ sc * n)
        vc = c + sc .* n;

        % alter out of bounds to bounds (max patch)
        vc(vc < 0.5) = 0.5;
        vc(vc > 256.5) = 256.5;

        % triangle vectors
        tv1  = vc(t(:, 1), :);
        tv12 = vc(t(:, 2), :) - tv1;
        tv13 = vc(t(:, 3), :) - tv1;

        % loop oversampling
        for oc2 = 1:numel(o)
            for oc3 = 1:(numel(o) + 1 - oc2)

                % compute coordinate within triangle
                tc = 1 + round(ores .* (tv1 + o(oc2) .* tv12 + o(oc3) .* tv13));

                % remove invalid indices
                tc(any(tc < 1 | tc > mxs, 2), :) = [];

                % get index
                tx = sub2ind(size(vmrd), tc(:, 1), tc(:, 2), tc(:, 3));

                % and set to 235 (or 240)
                vmrd(tx) = tcode;

                % update interface
                if now >= nup
                    drawnow;
                    nup = now + nupi;
                end
            end
        end

    % for linear
    else

        % compute fill value
        fvm = 2 * opts.nstep / (numel(o) * (numel(o) - 1)) * ...
            (1 + min(0, sc / abs(opts.nfrom - opts.nstep))) * opts.fillvmax;

        % get coordinate (+ sc * n)
        vc = c + sc .* n;

        % alter out of bounds to bounds (max patch)
        vc(vc < 0.5) = 0.5;
        vc(vc > 256.5) = 256.5;

        % triangle vectors
        tv1  = vc(t(:, 1), :);
        tv12 = vc(t(:, 2), :) - tv1;
        tv13 = vc(t(:, 3), :) - tv1;
        tv23 = vc(t(:, 2), :) - vc(t(:, 3), :);
        tsa = sqrt(sum(tv12 .* tv12, 2));
        tsb = sqrt(sum(tv13 .* tv13, 2));
        tsc = sqrt(sum(tv23 .* tv23, 2));
        ta = 0.5 * (tsa + tsb + tsc);
        ta = max(0.0001, real(2 .* sqrt(ta .* (ta - tsa) .* (ta - tsb) .* (ta - tsc))));

        % loop oversampling
        for oc2 = 1:numel(o)
            for oc3 = 1:(numel(o) + 1 - oc2)

                % adapt sampling value
                fvms = fvm;
                if any(oc2 == [1, numel(o)])
                    fvms = 0.5 * fvms;
                end

                if oc3 == 1
                    fvms = 0.5 * fvms;
                end

                % compute coordinate within triangle
                tc = 1 + ores .* (tv1 + o(oc2) .* tv12 + o(oc3) .* tv13);

                % get lower bound and diff
                tcl = floor(tc);
                tc = tc - tcl;

                % remove all that cannot make it
                tcx = any(tcl < 1 | tcl >= mxs, 2);
                tcl(tcx, :) = [];
                tc(tcx, :) = [];
                tac = fvms .* ta(~tcx);

                % get indices (for sorting)
                tcx = sub2ind(size(vmrd), tcl(:, 1), tcl(:, 2), tcl(:, 3));

                % sort indices
                [tcx, tcs] = sort(tcx);
                tc = tc(tcs, :);
                tcl = tcl(tcs, :);
                tac = tac(tcs);

                % iterate for x, y, z
                for xd = 0:1
                    for yd = 0:1
                        for zd = 0:1

                            % get indices
                            tcx = sub2ind(size(vmrd), tcl(:, 1) + xd, tcl(:, 2) + yd, tcl(:, 3) + zd);

                            % create copies
                            tcc = tc;
                            tab = tac;

                            % take care of all values
                            while ~isempty(tcx)

                                % unique values
                                tcld = [true; diff(tcx) > 0];

                                % set to higher value
                                vmrd(tcx(tcld)) = min(opts.fillvmax, vmrd(tcx(tcld)) + tab(tcld) .* prod(...
                                    [abs(tcc(tcld, 1) - xd), abs(tcc(tcld, 2) - yd), abs(tcc(tcld, 3) - zd)], 2));

                                 % set as done
                                 tcx(tcld) = [];
                                 tcc(tcld, :) = [];
                                 tab(tcld) = [];
                            end
                        end
                    end
                end

                % update interface
                if now >= nup
                    drawnow;
                    nup = now + nupi;
                end
            end
        end
    end
end

% convert back to uint8?
if ~isa(vmrd, 'uint8')

    % max reached ?
    maxv = max(vmrd(:));

    % if not
    if maxv < opts.fillvmax

        % scale to max
        vmrd = uint8(round(single(opts.fillvmax ./ maxv) .* vmrd));

    % reached
    else

        % simply round off
        vmrd = uint8(round(vmrd));
    end
end

% add SMP data?
if ~isempty(opts.smp)

    % patch coordinates
    c = 1 + c + mean(opts.nfrom:opts.nstep:opts.nto) .* n;

    % thresholded from 0 through 9
    t = [eps, 1:9, Inf];
    for tc = 1:numel(t)-1

        % find vertices that match
        tv = round(ores .* c(opts.smp >= t(tc) & opts.smp < t(tc+1), :));

        % get index
        tx = sub2ind(size(vmrd), tv(:, 1), tv(:, 2), tv(:, 3));

        % and set to 225 + tc
        vmrd(tx) = 225 + tc;
    end
end

% put back into VOI and return!
bc.VertexVMRData = vmrd;
xo.C = bc;

% mark in VMR?
if ~isempty(opts.vmr)

    % get content
    vc = opts.vmr.C;

    % resolution mismatch
    if (ores == 1 && any([vc.VoxResX, vc.VoxResY, vc.VoxResZ] ~= 1)) || ...
       (ores == 2 && any([vc.VoxResX, vc.VoxResY, vc.VoxResZ] ~= 0.5))

        % don't do anything
        return;
    end

    % overlapping bounding box
    xfrom = 1 + vc.OffsetX;
    xto = xfrom + vc.DimX - 1;
    yfrom = 1 + vc.OffsetY;
    yto = yfrom + vc.DimY - 1;
    zfrom = 1 + vc.OffsetZ;
    zto = zfrom + vc.DimZ - 1;

    % set in data
    if any([xfrom, yfrom, zfrom, xto, yto, zto] ~= [1, 1, 1, size(vc.VMRData)])
        vmrd = vmrd(xfrom:xto, yfrom:yto, zfrom:zto);
    end
    if istransio(vc.VMRData)
        vc.VMRData = resolve(vc.VMRData);
    end
    if ~isfield(vc.RunTimeVars, 'UndoBuffer')
        vc.RunTimeVars.UndoBuffer = vc.VMRData;
    end
    if strcmpi(opts.fillmode, 'nearest')
        vc.VMRData(vmrd == tcode) = tcode;
    else
        vc.VMRData = max(vc.VMRData, vmrd);
    end
    opts.vmr.C = vc;
end

function [xo2, xo3, xo4, xo5] = vmr_GradientVMR(xo, gdir, gcomp)
% VMR::GradientVMR  - compute gradient VMR(s)
%
% FORMAT:       [gvmr, gvmrx, gvmry, gvmrz] = vmr.GradientVMR([gdir]);
%
% Input fields:
%
%       gdir        flag whether to create directional VMRs too
%
% Output fields:
%
%       gvmr        gradient VMR (intensity of gradient)
%       gvmrx       X-gradient VMR (TAL X axis!)
%       gvmry       Y-gradient VMR (TAL Y axis!)
%       gvmrz       Z-gradient VMR (TAL Z axis!)

% Version:  v1.1
% Build:    16012712
% Date:     Jan-27 2016, 12:23 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bcs = xo.C;
if nargin < 2 || ~islogical(gdir) || isempty(gdir)
    gdir = false;
else
    gdir = gdir(1);
end
grad = ones([3, 3, 3]);
if nargin < 3 || ~ischar(gcomp)
    gcomp = '';
else
    gcomp = gcomp(:)';
    try
        gce = eval(gcomp);
    catch xfferror
        neuroelf_lasterr(xfferror);
        gce = [];
    end
    if ~isequal(size(gce), size(grad))
        gcomp = '';
    end
end

% first make a copy
vsz = size(bcs.VMRData);
xo2 = aft_CopyObject(xo);
xo2.F = '';
xo2.C.VMRData = uint8([]);
xo2.C.VMRData(1:vsz(1), 1:vsz(2), 1:vsz(3)) = uint8(0);
xo2.C.VMRData16 = [];
if ~gdir
    xo3 = [];
    xo4 = [];
    xo5 = [];
else
    xo3 = aft_CopyObject(xo2);
    xo4 = aft_CopyObject(xo2);
    xo5 = aft_CopyObject(xo2);
end
clobs = {xo2, xo3, xo4, xo5};

% take resolution into account
try
    xres = single(bcs.VoxResZ);
    yres = single(bcs.VoxResX);
    zres = single(bcs.VoxResY);
catch xfferror
    neuroelf_lasterr(xfferror);
    xres = single(1);
    yres = single(1);
    zres = single(1);
end

% calculus
onestep = true;
usev16 = false;
if prod(vsz) < 2e7
    try
        if isempty(bcs.VMRData16) || ...
            numel(size(bcs.VMRData16)) ~= numel(size(bcs.VMRData)) || ...
            any(size(bcs.VMRData16) ~= size(bcs.VMRData))
            usev16 = true;
            [grz, gry, grx] = gradient(single(bcs.VMRData), zres, yres, xres);
        else
            [grz, gry, grx] = gradient(single(bcs.VMRData16), zres, yres, xres);
        end
        grl = sqrt(grx .* grx + gry .* gry + grz .* grz);
    catch xfferror
        neuroelf_lasterr(xfferror);
        onestep = false;
    end
else
    onestep = false;
end
if ~onestep
    try
        if gdir
            grx = single([]);
            gry = single([]);
            grz = single([]);
            grx(1:vsz(1), 1:vsz(2), 1:vsz(3)) = 0;
            gry(1:vsz(1), 1:vsz(2), 1:vsz(3)) = 0;
            grz(1:vsz(1), 1:vsz(2), 1:vsz(3)) = 0;
        end
        grl = single([]);
        grl(1:vsz(1), 1:vsz(2), 1:vsz(3)) = 0;

        % make in packs
        for zc = 2:(vsz(3) - 1)
            if ~usev16
                [pgz, pgy, pgx] = gradient(single( ...
                    bcs.VMRData(:, :, (zc-1:zc+1))), zres, yres, xres);
            else
                [pgz, pgy, pgx] = gradient(single( ...
                    bcs.VMRData16(:, :, (zc-1:zc+1))), zres, yres, xres);
            end
            if zc == 2
                pgx = pgx(:, :, 1:2);
                pgy = pgy(:, :, 1:2);
                pgz = pgz(:, :, 1:2);
                pgl = sqrt(pgx .* pgx + pgy .* pgy + pgz .* pgz);
                if gdir
                    grx(:, :, 1:2) = pgx;
                    gry(:, :, 1:2) = pgy;
                    grz(:, :, 1:2) = pgz;
                end
                grl(:, :, 1:2) = pgl;
            elseif zc == (vsz(3) - 1)
                pgx = pgx(:, :, 2:3);
                pgy = pgy(:, :, 2:3);
                pgz = pgz(:, :, 2:3);
                pgl = sqrt(pgx .* pgx + pgy .* pgy + pgz .* pgz);
                if gdir
                    grx(:, :, end-1:end) = pgx;
                    gry(:, :, end-1:end) = pgy;
                    grz(:, :, end-1:end) = pgz;
                end
                grl(:, :, end-1:end) = pgl;
            else
                pgx = pgx(:, :, 2);
                pgy = pgy(:, :, 2);
                pgz = pgz(:, :, 2);
                pgl = sqrt(pgx .* pgx + pgy .* pgy + pgz .* pgz);
                if gdir
                    grx(:, :, zc) = pgx;
                    gry(:, :, zc) = pgy;
                    grz(:, :, zc) = pgz;
                end
                grl(:, :, zc) = pgl;
            end
        end
        clear pg*;
    catch xfferror
        clearxffobjects(clobs);
        error('neuroelf:xff:outOfMemory', 'Out of memory (%s).', xfferror.message);
    end
end

% computation
if ~isempty(gcomp)
    grad = grl;
    try
        grl = eval(gcomp);
    catch xfferror
        fprintf('Error computing on gradient: %s.\n', xfferror.message);
        grl = grad;
    end
end

% useful limitting depends on VMR type
mxg = single(ceil(max(grl(:))));
mxt = single(220);
[hn, hx] = hist(grl(grl > 0), 10 * mxg);
hnc = cumsum(hn);
ofx = find(hnc >= (0.999 * numel(grl)));
if isempty(ofx)
    ofx = floor(0.999 * numel(hn));
end
mxg = hx(ofx(1));
grl = min(mxt, grl ./ (mxg / mxt));
xo2.C.VMRData(1:vsz(1), 1:vsz(2), 1:vsz(3)) = round(grl);

% directions
if gdir

    % thresholding
    mxt = mxt / single(2);
    mxx = max(abs(grx(:)));
    mxy = max(abs(gry(:)));
    mxz = max(abs(grz(:)));
    mxg = 2 * max([mxx, mxy, mxz]);
    if mxg > mxt
        grx = (-grx - abs(min(grx(:)))) ./ (mxg ./ mxt);
        gry = (-gry - abs(min(gry(:)))) ./ (mxg ./ mxt);
        grz = (-grz - abs(min(grz(:)))) ./ (mxg ./ mxt);
    end

    % absolute values, in TAL notation !
    xo3.C.VMRData(:, :, :) = grz;
    xo4.C.VMRData(:, :, :) = grx;
    xo5.C.VMRData(:, :, :) = gry;
end

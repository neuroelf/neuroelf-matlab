function [aplot, splot, means] = avg_Average(xo, tpfile, tcd, setype, bas)
% AVG::Average  - average time course according to AVG info
%
% FORMAT:       [aplot, splot, means] = avg.Average(tpfile, tcd [, set, bas])
%
% Input fields:
%
%       tpfile      use trigger points from file # (0-based!)
%       tcd         TxV time course data
%       set         standard error type ({'BV'}, 'biased', or 'unbiased')
%       bas         optional baseline specification, can be either
%                   a Bx1 (file based), a BxC (condition-per-file based),
%                   or a BxCxE (epoch based) double array with TC indices
%
% Output fields:
%
%       aplot       average plots PxVxC
%       splot       std error plots PxVxC
%       means       mean values, 1xV or 1xVxC (depending on baseline)
%
% Note: this function only works for Volume-based AVG files
%       since without the FMR/VTC/MTC the TR is unknown

% Version:  v1.1
% Build:    16020218
% Date:     Feb-02 2016, 6:18 PM EST
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

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'avg') || ...
   ~isa(tpfile, 'double') || numel(tpfile) ~= 1 || isinf(tpfile) || isnan(tpfile) || ...
    tpfile < 0 || ~isa(tcd, 'double') || isempty(tcd)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if (tpfile + 1) > numel(bc.FileNames)
    error('neuroelf:xff:badArgument', 'Invalid tpfile argument.');
end
if ~strcmpi(bc.ProtocolTimeResolution, 'volumes') || ...
   ~ischar(bc.ResolutionOfDataPoints) || ~strcmpi(bc.ResolutionOfDataPoints, 'volumes')
    error('neuroelf:xff:notImplemented', ...
        'Averages can only be computed on Volume-based AVG files.');
end
tpfile = 1 + fix(real(tpfile));
if nargin < 4 || ~ischar(setype) || isempty(setype) || strcmpi(setype(:)', 'bv')
    setype = 0;
elseif ~strcmpi(setype(:)', 'unbiased')
    setype = 1;
else
    setype = 2;
end

% get sizes
voxtp = size(tcd, 1);
numvox = numel(tcd) / voxtp;

% get avg options
numcurves = numel(bc.Curve);
numtp = bc.NrOfTimePoints;
pretp = bc.PreInterval;
postp = numtp - (1 + pretp);
blmode = bc.BaselineMode;
if nargin > 4 && isa(bas, 'double') && ~isempty(bas)
    basz = [size(bas), 1];
    if blmode > 1 && basz(2) ~= numcurves
        bas = bas(:,1);
        blmode = 1;
        basz(2:end) = 1;
    end
    blfrom = 1;
    blupto = basz(1);
else
    bas = [];
    basz = [];
    blfrom = bc.AverageBaselineFrom;
    blupto = bc.AverageBaselineTo;
end
numbas = 1 + blupto - blfrom;
curves = bc.Curve;

% get number of curve points
numcps = zeros(1, numcurves);
for cc = 1:numcurves
    curve = curves(cc);
    if numel(curve.File) < tpfile
        error('neuroelf:xff:invalidObject', 'Invalid object, too few Files in Curve %d.', cc);
    end
    cfpts = curves(cc).File(tpfile).Points;
    numcps(cc) = numel(cfpts);
    if max(cfpts) >= voxtp
        error('neuroelf:xff:invalidObject', ...
            'Last trigger point goes beyond time course data.');
    end
end
maxcps = max(numcps);
if blmode > 2 && ~isempty(basz) && basz(3) < maxcps
    bas = bas(:, :, 1);
    basz(3:end) = 1;
    blmode = 2;
end
if blmode > 2 && ~isempty(basz) && basz(3) > 1 && basz(3) ~= maxcps
    bas = bas(:, :, 1:maxcps);
end
if ~isempty(basz)
    bas(bas(:) < 1 | bas(:) > voxtp) = 0;
    bas = round(bas);
    if basz(3) == 1
        bas = repmat(bas, [1, 1, maxcps]);
    end
    if basz(2) == 1
        bas = repmat(bas, [1, numcurves]);
    end
end

% create sample data arrays
pdata = zeros(numtp, numvox, numcurves, maxcps);
pgood = zeros(numtp, 1, numcurves, maxcps);
bdata = zeros(numbas, numvox, numcurves, maxcps);
bgood = zeros(numbas, 1, numcurves, maxcps);

% iterate over curves
for cc = 1:numcurves

    % fill with indices
    cfpts = curves(cc).File(tpfile).Points;
    for pc = 1:numel(cfpts)
        cpnt = cfpts(pc);
        sidx = (cpnt - pretp):(cpnt + postp);
        gidx = find((sidx > 0) & (sidx <= voxtp));
        pgood(gidx, 1, cc, pc) = 1;
        pdata(gidx, :, cc, pc) = tcd(sidx(gidx), :);
        if isempty(basz)
            sidx = (cpnt + blfrom):(cpnt + blupto);
        else
            sidx = bas(:, cc, pc)';
        end
        gidx = find((sidx > 0) & (sidx <= voxtp));
        bgood(gidx, 1, cc, pc) = 1;
        bdata(gidx, :, cc, pc) = tcd(sidx(gidx), :);
    end
end

% first-epoch based, replicate over points
if blmode > 3 && maxcps > 1
    bgood(:, :, :, 2:maxcps) = bgood(:, :, :, ones(1, maxcps - 1));
    bdata(:, :, :, 2:maxcps) = bdata(:, :, :, ones(1, maxcps - 1));
end

% if not epoch based, sum between curve points
if blmode < 3
    bgood = repmat(sum(bgood, 4), [1, 1, 1, maxcps]);
    bdata = repmat(sum(bdata, 4), [1, 1, 1, maxcps]);
end

% if file based, sum between conditions (curves)
if blmode < 2
    bgood = repmat(sum(bgood, 3), [1, 1, numcurves]);
    bdata = repmat(sum(bdata, 3), [1, 1, numcurves]);
end

% sum baseline between points within range
bgood = sum(bgood, 1);
bdata = sum(bdata, 1);
bdata = bdata ./ repmat(bgood, [1, numvox]);
bdata(~bgood) = 0;
bdata = repmat(bdata, [numtp, 1]);
means = bdata(1, :, :, 1);

% make transformation
if blmode > 0 || ~isempty(basz)
    pdata = (pdata - bdata) ./ (bdata / 100);
    pdata(isnan(pdata)) = 0;
end

% standard error
pgood = repmat(pgood, [1, numvox]);
pmean = repmat(sum(pdata .* pgood, 4) ./ sum(pgood, 4), [1, 1, 1, maxcps]);
pdata2 = pdata - pmean;
pdata2 = sum(pdata2 .* pdata2 .* pgood, 4);
if setype == 0
    splot = sqrt(pdata2) ./ sum(pgood, 4);
elseif setype == 1
    splot = sqrt(pdata2 ./ sum(pgood, 4));
else
    splot = sqrt(pdata2 ./ (sum(pgood, 4) - 1));
end
splot(isinf(splot(:)) | isnan(splot(:))) = 0;

% sum data where good
aplot = sum(pdata .* pgood, 4) ./ sum(pgood, 4);

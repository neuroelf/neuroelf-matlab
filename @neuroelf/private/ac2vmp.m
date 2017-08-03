function vmp = ac2vmp(acfile, opts)
% ac2vmp  - convert clusters file (TW tools) to VMP format
%
% FORMAT:       vmp = ac2vmp(acfile [, opts])
%
% Input fields:
%
%       acfile      clusters filename
%       opts        optional settings
%        .colors    if given must be 2x3xS with RGB values (positive only)
%        .output    either of 'clusters', or {'thresh'}
%        .outset    create a map of all clusters in a set (default: true)
%        .outsetf   factor for the set map (default: 1)
%        .sdf       stats df, must match type, default: [1000]
%        .stype     stats type, either 'F', 't', or {'Z'}
%        .zthresh   z-threshold to initialize VMP maps with (default: 0.05)
%
% Output fields:
%
%       vmp         VMP object with cluster maps
%
% Note: given that the Z values in the clusters are unsigned, only
%       positive-tail stats can be produced
%       several filenames can be given as cell array to create a common
%       VMP file (only if all clusters have the same dimensions!)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:12 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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
if nargin < 1 || ...
   (~ischar(acfile) && ...
    ~iscell(acfile)) || ...
    isempty(acfile)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'colors') || ...
   ~isa(opts.colors, 'double') || ...
    size(opts.colors, 1) ~= 2 || ...
    size(opts.colors, 2) ~= 3 || ...
    any(isinf(opts.colors(:)) | isnan(opts.colors(:)) | ...
        opts.colors(:) < 0 | opts.colors(:) > 255)
    opts.colors = zeros(2, 3, 0);
else
    opts.colors = fix(opts.colors);
end
numcols = size(opts.colors, 3);
if ~isfield(opts, 'output') || ...
   ~ischar(opts.output) || ...
   ~any(strcmpi(opts.output(:)', {'c', 'clusters', 't', 'thresh'}))
    opts.output = 't';
else
    opts.output = lower(opts.output(1));
end
if ~isfield(opts, 'outset') || ...
   ~islogical(opts.outset) || ...
    numel(opts.outset) ~= 1
    opts.outset = true;
end
if ~isfield(opts, 'outsetf') || ...
   ~isa(opts.outsetf, 'double') || ...
    numel(opts.outsetf) ~= 1 || ...
    isinf(opts.outsetf) || ...
    isnan(opts.outsetf) || ...
    opts.outsetf == 0
    opts.outsetf = 1;
end
if ~isfield(opts, 'sdf') || ...
   ~isa(opts.sdf, 'double') || ...
   ~any(numel(opts.sdf) == [1, 2]) || ...
    any(isinf(opts.sdf) | isnan(opts.sdf) | opts.sdf < 1 | opts.sdf ~= fix(opts.sdf))
    opts.sdf = [1, 1000];
end
if ~isfield(opts, 'stype') || ...
   ~ischar(opts.stype) || ...
    numel(opts.stype) ~= 1 || ...
   ~any(lower(opts.stype) == 'ftz')
    opts.stype = 'z';
else
    opts.stype = lower(opts.stype);
end
if opts.stype ~= 'f'
    opts.sdf = max(opts.sdf);
elseif numel(opts.sdf) ~= 2
    opts.sdf = [1, opts.sdf];
end
if ~isfield(opts, 'zthresh') || ...
   ~isa(opts.zthresh, 'double') || ...
    numel(opts.zthresh) ~= 1 || ...
    isinf(opts.zthresh) || ...
    isnan(opts.zthresh) || ...
    opts.zthresh <= 0 || ...
    opts.zthresh >= 0.5
    opts.zthresh = 0.05;
end
if opts.stype == 't'
    lowt = -sdist('tinv', 0.5 * opts.zthresh, opts.sdf);
    uppt = -sdist('tinv', 0.005 * opts.zthresh, opts.sdf);
elseif opts.stype == 'f'
    lowt = sdist('finv', opts.zthresh, opts.sdf(1), opts.sdf(2), true);
    uppt = sdist('finv', 0.01 * opts.zthresh, opts.sdf(1), opts.sdf(2), true);
else
    lowt = opts.zthresh;
    uppt = 2 * lowt;
end
if ~iscell(acfile)
    acfile = {acfile(:)'};
else
    acfile = acfile(:)';
end
cdim = [];
cM = [];
nummaps = 0;
for c = 1:numel(acfile)
    try
        filename = acfile{c}(:)';
        acfile{c} = load(filename);
        if ~isfield(acfile{c}, 'cl') || ...
           ~iscell(acfile{c}.cl) || ...
            isempty(acfile{c}.cl) || ...
            all(cellfun('isempty', acfile{c}.cl))
            error( ...
                'neuroelf:BadFile', ...
                'MAT file %s is not a clusters file.', ...
                filename ...
            );
        end
        nec = findfirst(~cellfun('isempty', acfile{c}.cl));
        if ~isstruct(acfile{c}.cl{nec}) || ...
            isempty(acfile{c}.cl{nec}) || ...
           ~isfield(acfile{c}.cl{nec}, 'M') || ...
           ~isfield(acfile{c}.cl{nec}, 'dim') || ...
           ~isequal(size(acfile{c}.cl{nec}(1).M), [4, 4]) || ...
           (~isequal(size(acfile{c}.cl{nec}(1).dim), [1, 3]) && ...
            ~isequal(size(acfile{c}.cl{nec}(1).dim), [1, 4]))
            error( ...
                'neuroelf:BadFile', ...
                'MAT file %s is not a clusters file.', ...
                filename ...
            );
        end
        for cc = 1:numel(acfile{c}.cl)
            if isempty(acfile{c}.cl{cc}) || ...
               ~isstruct(acfile{c}.cl{cc})
                acfile{c}.cl{cc} = acfile{c}.cl{nec}(1);
                acfile{c}.cl{cc}.numVox = 0;
                acfile{c}.cl{cc}.XYZ = zeros(3, 0);
                acfile{c}.cl{cc}.XYZmm = zeros(3, 0);
                acfile{c}.cl{cc}.Z = zeros(1, 0);
            end
        end
        if isempty(cdim)
            cdim = acfile{c}.cl{1}(1).dim(1:3);
        elseif ~isequal(cdim, acfile{c}.cl{1}(1).dim(1:3))
            error( ...
                'neuroelf:BadFile', ...
                'Dimensions mismatch between multiple clusters files.' ...
            );
        end
        if isempty(cM)
            cM = acfile{c}.cl{1}(1).M;
        elseif ~isequal(cM, acfile{c}.cl{1}(1).M)
            error( ...
                'neuroelf:BadFile', ...
                'Spatial layout mismatches between multiple clusters files.' ...
            );
        end
        if opts.output == 'c'
            for cc = 1:numel(acfile{c}.cl)
                nummaps = nummaps + numel(acfile{c}.cl{cc});
            end
        else
            nummaps = nummaps + numel(acfile{c}.cl);
        end
        [filepath, acfile{c}.filename] = fileparts(filename);
    catch ne_eo;
        rethrow(ne_eo);
    end
end
if opts.outset
    nummaps = nummaps + numel(acfile);
end

% check settings
if numel(cdim) ~= 3 || ...
   ~isa(cdim, 'double') || ...
    any(isinf(cdim) | isnan(cdim) | cdim < 1 | cdim ~= fix(cdim))
    error( ...
        'neuroelf:BadFile', ...
        'Invalid data dimension.' ...
    );
end
if ~isequal(size(cM), [4, 4]) || ...
   ~isa(cM, 'double') || ...
    any(isinf(cM(:)) | isnan(cM(:))) || ...
    any(cM(4, :) ~= [0, 0, 0, 1])
    error( ...
        'neuroelf:BadFile', ...
        'Invalid transformation matrix.' ...
    );
end
ncrd = cM * [1, cdim(1), 1, cdim(1), 1, cdim(1), 1, cdim(1); ...
    1, 1, cdim(2), cdim(2), 1, 1, cdim(2), cdim(2); ...
    ones(1, 4), cdim(3) .* ones(1, 4); ones(1, 8)];
bcrd = 128 - max(ncrd(1:3, :), [], 2);
ecrd = 128 - min(ncrd(1:3, :), [], 2);
res = round(min(max(cM(1:3, 1:3))));

% create new VMP
vmp = xff('new:vmp');

% make initial settings
vmp.NrOfMaps = nummaps;
vmp.XStart = bcrd(2);
vmp.XEnd = bcrd(2) + res * ceil(1 + (ecrd(2) - bcrd(2)) / res);
vmp.YStart = bcrd(3);
vmp.YEnd = bcrd(3) + res * ceil(1 + (ecrd(3) - bcrd(3)) / res);
vmp.ZStart = bcrd(1);
vmp.ZEnd = bcrd(1) + res * ceil(1 + (ecrd(1) - bcrd(1)) / res);
vmp.Resolution = res;
if opts.stype == 'f'
    vmp.Map.Type = 2;
end
vmp.Map.LowerThreshold = lowt;
vmp.Map.UpperThreshold = uppt;
vmp.Map.UseRGBColor = 1;
vmp.Map.DF1 = opts.sdf(1);
if opts.stype == 'f'
    vmp.Map.DF2 = opts.sdf(2);
end
cdim = ([vmp.XEnd - vmp.XStart, vmp.YEnd - vmp.YStart, vmp.ZEnd - vmp.ZStart]) ./ res;
vmp.Map.VMPData = single(zeros(cdim));
if nummaps > 1
    vmp.Map = vmp.Map(ones(1, nummaps));
end
bb = vmp.BoundingBox;

% iterate over clusters and structures
fm = 1;
tm = 1;
for c = 1:numel(acfile)
    for cc = 1:numel(acfile{c}.cl)
        for ccc = 1:numel(acfile{c}.cl{cc})
            mapc = bvcoordconv(acfile{c}.cl{cc}(ccc).XYZmm', 'tal2bvx', bb);
            mapz = acfile{c}.cl{cc}(ccc).Z(:);
            if isempty(mapz)
                continue;
            end
            if opts.stype == 't'
                mapz = -sdist('tinv', 0.5 .* mapz, opts.sdf);
            elseif opts.stype == 'f'
                mapz = sdist('finv', mapz, opts.sdf(1), opts.sdf(2), true);
            end
            vmp.Map(tm).VMPData(mapc) = mapz;
            if opts.output == 'c'
                vmp.Map(tm).Name = sprintf('Cluster %d from set %d in file %s', ...
                    ccc, cc, acfile{c}.filename);
                if cc <= numcols
                    vmp.Map(tm).RGBLowerThreshPos = opts.colors(1, :, cc);
                    vmp.Map(tm).RGBUpperThreshPos = opts.colors(2, :, cc);
                end
                tm = tm + 1;
            end
        end
        if opts.output == 't'
            vmp.Map(tm).Name = sprintf('Set %d in file %s', cc, acfile{c}.filename);
            if cc <= numcols
                vmp.Map(tm).RGBLowerThreshPos = opts.colors(1, :, cc);
                vmp.Map(tm).RGBUpperThreshPos = opts.colors(2, :, cc);
            end
            tm = tm + 1;
        end
    end
    if opts.outset
        vmp.Map(tm).Name = sprintf('All clusters in file %s', cc, acfile{c}.filename);
        mapz = vmp.Map(fm).VMPData;
        for mc = (fm + 1):(tm - 1)
            mapz = max(mapz, vmp.Map(mc).VMPData);
        end
        vmp.Map(tm).LowerThreshold = opts.outsetf * lowt;
        vmp.Map(tm).UpperThreshold = opts.outsetf * uppt;
        vmp.Map(tm).VMPData = single(opts.outsetf .* mapz);
        tm = tm + 1;
        fm = tm;
    end
end

% set extended colors (for display)
vmp.SetColors([], 'xauto');

function msk = vmp_CreateMSK(xo, mapsel, opts)
% VMP::CreateMSK  - create MSK from the selected maps in a VMP
%
% FORMAT:       msk = vmp.CreateMSK([mapsel [, opts]])
%
% Input fields:
%
%       mapsel      map selection, default: all maps in VMP
%       opts        optional settings
%        .clconn    cluster connectivity
%        .combtype  combination type, either of {'and'}, 'mean', 'or'
%        .mthresh   mean-threshold (only if combtype = 'mean', default: 0.5)
%
% Output fields:
%
%       msk         MSK object with voxels that survive thresholding := 1
%
% Note: the default combination type is 'and' which means conjunction of
%       the selected maps; to obtain voxels present in any of the maps
%       set combination type to 'or'; 'mean' is if the mean of the maps
%       reach mthresh (e.g. for single-subject t-maps)
%
% Using: clustercoordsc.

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:15 PM EST
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

% neuroelf library
global ne_methods;
clustercoordsc = ne_methods.clustercoordsc;

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isa(mapsel, 'double') || isempty(mapsel) || ...
    any(isinf(mapsel(:)) | isnan(mapsel(:)) | mapsel(:) ~= fix(mapsel(:)) | mapsel(:) < 1)
    mapsel = 1:numel(bc.Map);
else
    mapsel = unique(min(numel(bc.Map), floor(mapsel(:)')));
end
nummaps = numel(mapsel);
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'clconn') || ~ischar(opts.clconn) || isempty(opts.clconn) || ...
   ~any(strcmpi(opts.clconn(:)', {'edge', 'face', 'vertex'}))
    opts.clconn = 'edge';
else
    opts.clconn = lower(opts.clconn(:)');
end
opts.clconn = find(strcmp(opts.clconn, {'face', 'edge', 'vertex'}));
if ~isfield(opts, 'combtype') || ~ischar(opts.combtype) || isempty(opts.combtype) || ...
   ~any(strcmpi(opts.combtype(:)', {'and', 'mean', 'or'}))
    opts.combtype = 'and';
else
    opts.combtype = lower(opts.combtype(:)');
end
if ~isfield(opts, 'mthresh') || ~isa(opts.mthresh, 'double') || numel(opts.mthresh) ~= 1 || ...
    isinf(opts.mthresh) || isnan(opts.mthresh) || opts.mthresh <= 0 || opts.mthresh > 1
    opts.mthresh = 0.5;
end

% generate MSK
msk = xff('new:msk');
mskc = msk.C;

% copy properties
mskc.Resolution = bc.Resolution;
mskc.XStart = bc.XStart;
mskc.XEnd = bc.XEnd;
mskc.YStart = bc.YStart;
mskc.YEnd = bc.YEnd;
mskc.ZStart = bc.ZStart;
mskc.ZEnd = bc.ZEnd;

% start with usefully sized mask
masksize = floor([bc.XEnd - bc.XStart, bc.YEnd - bc.YStart, bc.ZEnd - bc.ZStart] ./ bc.Resolution);
switch (opts.combtype)
    case 'and'
        mask = true(masksize);
    case 'mean'
        mask = zeros(masksize);
    case 'or'
        mask = false(masksize);
end

% iterate over maps
for mc = 1:nummaps

    % get source map and values
    srcmap = bc.Map(mapsel(mc));
    mapval = double(srcmap.VMPData);
    srcmask = false(size(mask));

    % cluster if necessary
    if srcmap.EnableClusterCheck > 0

        % cluster > 0 tail
        if mod(srcmap.ShowPositiveNegativeFlag, 2) > 0
            [srccls, srcmaskp] = clustercoordsc(mapval >= srcmap.LowerThreshold, ...
                opts.clconn, srcmap.ClusterSize);
            srcmask = (srcmaskp > 0);
        end

        % cluster < 0 tail
        if srcmap.ShowPositiveNegativeFlag >= 2
            [srccls, srcmaskn] = clustercoordsc(mapval <= -srcmap.LowerThreshold, ...
                opts.clconn, srcmap.ClusterSize);
            srcmask = srcmask | (srcmaskn > 0);
        end

    % no clustering
    else
        % use > 0 tail
        if mod(srcmap.ShowPositiveNegativeFlag, 2) > 0
            srcmask = (mapval >= srcmap.LowerThreshold);
        end

        % use < 0 tail
        if srcmap.ShowPositiveNegativeFlag >= 2
            srcmask = srcmask | (mapval <= -srcmap.LowerThreshold);
        end
    end

    % combine with mask
    switch (opts.combtype)
        case 'and'
            mask = mask & srcmask;
        case 'mean'
            mask = mask + double(srcmask);
        case 'or'
            mask = mask | srcmask;
    end
end

% mean threshold
if strcmp(opts.combtype, 'mean')
    mask = (mask >= (opts.mthresh * nummaps));
end

% put content in new object
mskc.Mask = uint8(mask);
msk.C = mskc;

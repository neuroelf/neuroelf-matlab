function xo = vmp_PerformFDR(xo, opts)
% VMP::PerformFDR  - perform FDR correction (set properties)
%
% FORMAT:       [vmp = ] vmp.PerformFDR([opts])
%
% Input fields:
%
%       opts        1x1 struct with optional fields
%        .bonfval   alternative BonferroniValue (override map value)
%        .mapsel    map selection, default: all supported maps in VMP
%        .useske    use smoothing kernel estimate (default: true)
%
% Output fields:
%
%       vmp         VMP object
%
% Using: applyfdr, mapestsmooth, smoothdata3, smoothkern.

% Version:  v1.1
% Build:    16021110
% Date:     Feb-11 2016, 10:39 AM EST
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

% neuroelf library and settings
global ne_methods xffsngl;
applyfdr     = ne_methods.applyfdr;
mapestsmooth = ne_methods.mapestsmooth;
smoothdata3  = ne_methods.smoothdata3;
smoothkern   = ne_methods.smoothkern;
fdrc = xffsngl.CONF.settings.Statistics.FDR.Thresholds;

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if bc.NativeResolutionFile == 0
    return;
end
res = bc.Resolution;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bonfval') || ~isa(opts.bonfval, 'double') || numel(opts.bonfval) ~= 1 || ...
    isinf(opts.bonfval) || isnan(opts.bonfval) || opts.bonfval < 2
    opts.bonfval = [];
else
    opts.bonfval = round(opts.bonfval);
end
if ~isfield(opts, 'mapsel') || ~isa(opts.mapsel, 'double') || ...
    any(isinf(opts.mapsel(:)) | isnan(opts.mapsel(:)) | ...
        opts.mapsel(:) ~= fix(opts.mapsel(:)) | opts.mapsel(:) < 1)
    opts.mapsel = 1:numel(bc.Map);
else
    opts.mapsel = unique(min(opts.mapsel(:), numel(bc.Map)));
end
if ~isfield(opts, 'useske') || ~islogical(opts.useske) || numel(opts.useske) ~= 1
    opts.useske = true;
end

% iterate over maps
for mc = 1:numel(opts.mapsel)

    % get map handle
    map = bc.Map(opts.mapsel(mc));
    mmsk = (map.VMPData ~= 0 & ~isinf(map.VMPData) & ~isnan(map.VMPData));

    % continue if map type is not supported
    if ~any([1, 2, 4] == map.Type)
        continue;
    end

    % get map values we care about
    if ~isempty(opts.bonfval)
        bonfval = opts.bonfval;
    else
        if opts.useske
            if ~isfield(map, 'RunTimeVars') || ~isstruct(map.RunTimeVars) || numel(map.RunTimeVars) ~= 1
                map.RunTimeVars = struct;
            end
            rtv = map.RunTimeVars;
            if isfield(rtv, 'FWHMResEst') && isa(rtv.FWHMResEst, 'double') && numel(rtv.FWHMResEst) == 3
                smest = rtv.FWHMResEst;
            elseif isfield(rtv, 'FWHMMapEst') && isa(rtv.FWHMMapEst, 'double') && numel(rtv.FWHMMapEst) == 3
                smest = rtv.FWHMMapEst;
            else
                rtv.FWHMMapEst = mapestsmooth(double(map.VMPData), res);
                smest = rtv.FWHMMapEst;
            end
            bc.Map(opts.mapsel(mc)).RunTimeVars = rtv;

            % compute the map we care about
            smkm = max(lsqz(smoothkern(smest ./ res, 0, false, 'lanczos8')));
            bonfval = min(sum(mmsk(:)), sum(lsqz((smoothdata3(double(mmsk), ...
                smest ./ res) >= smkm))) ./ prod(smest ./ res));
            bc.Map(opts.mapsel(mc)).BonferroniValue = round(bonfval);
        end
        bonfval = round(bc.Map(opts.mapsel(mc)).BonferroniValue);
    end
    if mod(map.ShowPositiveNegativeFlag, 2) == 0
        mmsk(map.VMPData > 0) = false;
    end
    if map.ShowPositiveNegativeFlag < 2
        mmsk(map.VMPData < 0) = false;
    end
    mpv = map.VMPData(mmsk);
    if isempty(mpv)
        bc.Map(opts.mapsel(mc)).NrOfFDRThresholds = 0;
        bc.Map(opts.mapsel(mc)).FDRThresholds = zeros(0, 3);
        continue;
    end

    % computation depends on map type
    switch (map.Type)
        case {1} % t-map
            fdrv = [fdrc(:), applyfdr(abs(double(mpv)), 't', fdrc(:), map.DF1, [], 3, bonfval)];
        case {2} % r-map
            fdrv = [fdrc(:), applyfdr(abs(double(mpv)), 'r', fdrc(:), map.DF1, [], 3, bonfval)];
        case {4} % F-map
            fdrv = [fdrc(:), applyfdr(abs(double(mpv)), 'F', fdrc(:), map.DF1, map.DF2, 3, bonfval)];
    end

    % store
    bc.Map(opts.mapsel(mc)).NrOfFDRThresholds = size(fdrv, 1);
    bc.Map(opts.mapsel(mc)).FDRThresholds = fdrv;
end

% put content back in object
xo.C = bc;


% local implementation of ne_methods.lsqueeze
function v = lsqz(v)
v = v(:);

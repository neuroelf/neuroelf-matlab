function xo = vtc_Filter(xo, opts)
% VTC::Filter  - filter a VTC
%
% FORMAT:       [vtc =] vtc.Filter(opts)
%
% Input fields:
%
%       opts        mandatory struct but with optional fields
%        .nuisreg   either a VxN double matrix or single/list of SDM/s
%        .orthpoly  orthogonalize polynomials (faster computation later)
%        .spat      enable spatial filtering (default: false)
%        .spkern    smoothing kernel in mm (default: [6, 6, 6])
%        .temp      enable temporal filtering (default: false)
%        .tempdct   DCT-based (min. wavelength in seconds, default: Inf)
%        .tempdt    detrend (default: true, is overriden by dct/sc)
%        .temphp    temporal highpass (inv. smoothing) in units (def: 0)
%        .templp    temporal lowpass (smoothing) kernel in secs (def: 0)
%        .temppoly  set of orthogonal polynomials
%        .tempsc    sin/cos set of frequencies (number of pairs, def: 0)
%
% Output fields:
%
%       vtc         filtered VTC
%
% Using: tempfilter.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:04 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    return;
end
if ~isfield(opts, 'nuisreg') || ((numel(opts.nuisreg) ~= 1 || ...
    ~xffisobject(opts.nuisreg, true, 'sdm')) && ~iscell(opts.nuisreg))
    opts.nuisreg = [];
end
if ~isempty(opts.nuisreg)
    if ~iscell(opts.nuisreg)
        opts.nuisreg = {opts.nuisreg};
    end
    nr = zeros([size(bc.VTCData, 1), 0]);
    for nc = 1:numel(opts.nuisreg)
        if numel(opts.nuisreg{nc}) == 1 && xffisobject(opts.nuisreg{nc}, true, 'sdm')
            sdmc = opts.nuisreg{nc}.C;
            if size(sdmc.SDMMatrix, 1) == size(nr, 1)
                nr = [nr, sdmc.SDMMatrix];
            end
        end
    end
    nr(:, var(nr) <= sqrt(eps) | any(isinf(nr)) | any(isnan(nr))) = [];
    opts.nuisreg = nr;
else
    opts.nuisreg = [];
end
if ~isfield(opts, 'spat') || (~islogical(opts.spat) && ~isnumeric(opts.spat)) || ...
    numel(opts.spat) ~= 1
    opts.spat = false;
else
    opts.spat = (opts.spat ~= 0);
end
if ~isfield(opts, 'temp') || (~islogical(opts.temp) && ~isnumeric(opts.temp)) || ...
    numel(opts.temp) ~= 1
    opts.temp = false;
else
    opts.temp = (opts.temp ~= 0);
end
if ~isfield(opts, 'spkern') || ~isa(opts.spkern, 'double') || numel(opts.spkern) ~= 3 || ...
    any(isinf(opts.spkern) | isnan(opts.spkern) | opts.spkern < 0 | opts.spkern > 20)
    opts.spkern = [6, 6, 6];
else
    opts.spkern = opts.spkern(:)';
end
opts.spkern = opts.spkern ./ bc.Resolution;
if all(opts.spkern < 0.5)
    opts.spat = false;
end
if opts.temp
    if ~isfield(opts, 'orthpoly') || ~islogical(opts.orthpoly) || numel(opts.orthpoly) ~= 1
        opts.orthpoly = false;
    end
    if ~isfield(opts, 'tempdct') || ~isa(opts.tempdct, 'double') || numel(opts.tempdct) ~= 1 || ...
        isinf(opts.tempdct) || isnan(opts.tempdct) || opts.tempdct < (0.004 * bc.TR)
        opts.tempdct = Inf;
    else
        opts.tempdct = opts.tempdct ./ (0.001 * bc.TR);
    end
    if ~isfield(opts, 'tempdt') || ~islogical(opts.tempdt) || numel(opts.tempdt) ~= 1
        opts.tempdt = true;
    end
    if ~isfield(opts, 'temphp') || ~isa(opts.temphp, 'double') || numel(opts.temphp) ~= 1 || ...
        isinf(opts.temphp) || isnan(opts.temphp) || opts.temphp < 0
        opts.temphp = 0;
    else
        opts.temphp = 1000 * opts.temphp / bc.TR;
    end
    if ~isfield(opts, 'templp') || ~isa(opts.templp, 'double') || numel(opts.templp) ~= 1 || ...
        isinf(opts.templp) || isnan(opts.templp) || opts.templp < 0
        opts.templp = 0;
    else
        opts.templp = opts.templp ./ (0.001 * bc.TR);
    end
    if ~isfield(opts, 'temppoly') || ~isa(opts.temppoly, 'double') || numel(opts.temppoly) ~= 1 || ...
        isinf(opts.temppoly) || isnan(opts.temppoly) || opts.temppoly < 0
        opts.temppoly = 0;
    else
        opts.temppoly = min(size(bc.VTCData, 1) - 2, round(opts.temppoly));
    end
    if ~isfield(opts, 'tempsc') || ~isa(opts.tempsc, 'double') || numel(opts.tempsc) ~= 1 || ...
       ~any(opts.tempsc == (1:12))
        opts.tempsc = 0;
    end
    if isinf(opts.tempdct) && ~opts.tempdt && opts.temphp == 0 && ...
        opts.templp == 0 && opts.tempsc == 0
        opts.temp = false;
    elseif opts.temppoly > 0
        opts.tempsc = 0;
        opts.tempdct = Inf;
        opts.tempdt = false;
    elseif opts.tempsc > 0
        opts.tempdct = Inf;
        opts.tempdt = false;
    end
end
if ~opts.spat && ~opts.temp && isempty(opts.nuisreg)
    return;
end

% patch opts
opts.tdim = 1;

% and do work
bc.DataType = 2;
bc.VTCData = single(ne_methods.tempfilter(bc.VTCData, opts));

% but content into array
xo.C = bc;

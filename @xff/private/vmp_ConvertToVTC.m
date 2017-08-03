function vo = vmp_ConvertToVTC(xo, opts)
% VMP::ConvertToVTC  - converts a VMP (series of maps) to a VTC
%
% FORMAT:       [vtc =] vmp.ConvertToVTC([opts])
%
% Input fields:
%
%       opts        1x1 struct with optional fields
%       .dtype      either 1 (uint16) or {2} (single/float)
%       .mapsel     either 1xN list or regexp for map names
%       .scaling    if given scaling factor other than auto-detect
%       .vtcfile    output VTC filename
%
% Output fields:
%
%       vtc         VTC object containing VMP data
%
% Note: the scaling will always place VMP data 0.0 at 16384!

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:24 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get content of VMP
bc = xo.C;
nm = numel(bc.Map);

% check options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dtype') || ~isa(opts.dtype, 'double') || numel(opts.dtype) ~= 1 || ...
    isinf(opts.dtype) || isnan(opts.dtype) || ~any((1:2) == opts.dtype)
    opts.dtype = 2;
end
if ~isfield(opts, 'mapsel') || isempty(opts.mapsel) || (~ischar(opts.mapsel) && ~isa(opts.mapsel, 'double'))
    opts.mapsel = 1:nm;
else
    if ischar(opts.mapsel)
        msel = false(1, nm);
        for mc = 1:numel(msel)
            if ~isempty(regexpi(bc.Map(mc).Name, opts.mapsel(:)'))
                msel(mc) = true;
            end
        end
        opts.mapsel = find(msel);
    else
        opts.mapsel = opts.mapsel(:)';
        if any(isinf(opts.mapsel) | isnan(opts.mapsel) | opts.mapsel < 1 | opts.mapsel > nm)
            opts.mapsel = 1:nm;
        else
            opts.mapsel = floor(opts.mapsel);
        end
    end
end
ms = opts.mapsel;
if isempty(ms)
    error('neuroelf:xff:badArgument', 'No maps selected.');
end
if ~isfield(opts, 'scaling') || ~isa(opts.scaling, 'double') || numel(opts.scaling) ~= 1 || ...
    isinf(opts.scaling) || isnan(opts.scaling) || opts.scaling > 16384
    opts.scaling = Inf;
end
ds = opts.scaling;
if isfield(opts, 'vtcfile') && ischar(opts.vtcfile) && ~isempty(opts.vtcfile)
    vtcfile = opts.vtcfile(:)';
else
    vtcfile = '';
end

% check map types
tm = bc.Map(ms(1));
for mc = 2:numel(ms)
    if bc.Map(ms(mc)).Type ~= tm.Type
        error('neuroelf:xff:badArgument', 'Maps of different types cannot be extracted together.');
    end
end
if tm.Type == 3
    error('neuroelf:xff:notYetImplemented', 'Conversion of CC-maps not yet implemented.');
end

% auto scaling
if opts.dtype == 1
    mnv = Inf;
    mxv = -Inf;
    icp = 16384;
    for mc = 1:numel(ms)
        mnv = min([mnv, bc.Map(ms(mc)).VMPData(:)']);
        mxv = max([mxv, bc.Map(ms(mc)).VMPData(:)']);
    end
    if isinf(ds)
        switch (tm.Type)
            case {1, 11, 12, 15} % t-Map, %-Map, z-Map, (raw) beta-Map
                mxv = max(mxv, -mnv);
                ds = 16383 / mxv;
            case 2 % r-Map
                ds = 16383;
            case 3 % CC-Map (no code yet, see above)
            case {4, 13} % F-Map / TH-Map
                ds = 32767 / mxv;
                icp = 0;
            otherwise
                error('neuroelf:xff:unsupported', ...
                    'Unsupported map type for conversion: %d.', tm.Type);
        end
    end
end

% create VTC
vo = xff('new:vtc');
vbc = vo.C;

% make settings
vbc.FileVersion = 3;
vbc.DataType = opts.dtype;
vbc.NrOfVolumes = numel(ms);
vbc.Resolution = bc.Resolution;
vbc.XStart = bc.XStart;
vbc.YStart = bc.YStart;
vbc.ZStart = bc.ZStart;
if bc.Resolution > 1 || bc.FileVersion > 4
    vbc.XEnd = bc.XEnd;
    vbc.YEnd = bc.YEnd;
    vbc.ZEnd = bc.ZEnd;
else
    vbc.XEnd = bc.XEnd + 1;
    vbc.YEnd = bc.YEnd + 1;
    vbc.ZEnd = bc.ZEnd + 1;
end
tms = size(tm.VMPData);
vbc.VTCData(numel(ms), tms(1), tms(2), tms(3)) = 0;

% copy map data
if opts.dtype == 1
    for mc = 1:numel(ms)
        vbc.VTCData(mc, :, :, :) = reshape(uint16(icp + ds .* bc.Map(ms(mc)).VMPData), [1, tms]);
    end
else
    for mc = 1:numel(ms)
        vbc.VTCData(mc, :, :, :) = reshape(bc.Map(ms(mc)).VMPData(:, :, :), [1, tms]);
    end
end
% store back
vo.C = vbc;

% save file ?
if ~isempty(vtcfile)
    try
        aft_SaveAs(vo, vtcfile);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

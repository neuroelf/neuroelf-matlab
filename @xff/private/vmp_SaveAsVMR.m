function vo = vmp_SaveAsVMR(xo, vmrfile, mapno)
% VMP::SaveAsVMR  - saves a VMP as VMR
%
% FORMAT:       [vmr] = vmp.SaveAsVMR(vmrfile [, mapno])
%
% Input fields:
%
%       vmrfile     filename of VMR to write
%       mapno       number of map to write (default: 1)
%
% Output fields:
%
%       vmr         VMR xff object
%
% Note: currently only up to 256^3 VMR's are supported

% Version:  v1.1
% Build:    16021315
% Date:     Feb-13 2016, 3:10 PM EST
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp')
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~ischar(vmrfile) || isempty(vmrfile)
    vmrfile = '';
else
    vmrfile = vmrfile(:)';
end
if nargin < 3 || ~isa(mapno, 'double') || numel(mapno) ~= 1 || ...
    isinf(mapno) || isnan(mapno) || mapno < 1 || mapno > numel(bc.Map) || mapno ~= fix(mapno)
    mapno = 1;
end

% get position and resolution for data
xpos = bc.XStart + round((256 - bc.VMRDimX) / 2) + 1;
xend = bc.XEnd + round((256 - bc.VMRDimX) / 2) + 1;
ypos = bc.YStart + round((256 - bc.VMRDimY) / 2) + 1;
yend = bc.YEnd + round((256 - bc.VMRDimY) / 2) + 1;
zpos = bc.ZStart + round((256 - bc.VMRDimZ) / 2) + 1;
zend = bc.ZEnd + round((256 - bc.VMRDimZ) / 2) + 1;
pr = bc.Resolution;

% get map range
mapmin = min(bc.Map(mapno).VMPData(:));
mapmax = max(bc.Map(mapno).VMPData(:));

% check range
if isinf(mapmin) || isnan(mapmin) || isinf(mapmax) || isnan(mapmax)
    error('neuroelf:xff:badFileContent', 'No Inf or NaN in Map supported.');
elseif mapmin == mapmax
    mapmax = mapmax + eps;
end
maprng = mapmax - mapmin;

% create new VMR
vo = xff('new:vmr');
vc = vo.C;

% fill VMR
try
    map = uint8(25 + round(200 * (bc.Map(mapno).VMPData - mapmin) / maprng));
    for xc = 1:pr
        for yc = 1:pr
            for zc = 1:pr
                vc.VMRData(xpos:pr:xend, ypos:pr:yend, zpos:pr:zend) = map;
            end
        end
    end
    if size(vc.VMRData, 1) > 256
        vc.VMRData(257:end, :, :) = [];
    end
    if size(vc.VMRData, 2) > 256
        vc.VMRData(:, 257:end, :) = [];
    end
    if size(vc.VMRData, 3) > 256
        vc.VMRData(:, :, 257:end) = [];
    end
catch xfferror
    neuroelf_lasterr(xfferror);
    delete(vo);
    error('neuroelf:xff:badFileContent', 'VMP inner dimensions mismatch.');
end

% save VMR?
vo.C = vc;
if ~isempty(vmrfile)
    try
        aft_SaveAs(vo, vmrfile);
    catch xfferror
        neuroelf_lasterr(xfferror);
        warning('neuroelf:xff:errorWritingFile', ...
            'Couldn''t write VMR file to disk: ''%s''.', vmrfile);
    end
end

function xo = vmr_SaveV16(xo, v16fname)
% VMR::SaveV16  - save matching V16 file into 16-bit VMR file
%
% FORMAT:       [vmr] = vmr.SaveV16([v16fname])
%
% Input fields:
%
%       v16fname    alternative filename, otherwise use VMR's filename
%
% Output fields:
%
%       vmr         VMR object (with possibly updated transio reference!)

% Version:  v1.1
% Build:    16021321
% Date:     Feb-13 2016, 9:20 PM EST
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
bc = xo.C;
if strcmpi(class(bc.VMRData), 'uint16')
    error('neuroelf:xff:invalidObject', 'Method only valid for 8-bit VMRs.');
end
if nargin < 2
    [vmrfname{1:3}] = fileparts(xo.F);
    if isempty(vmrfname{2})
        error('neuroelf:xff:invalidObject', ...
            'This method only works without arguments on loaded VMRs.');
    end
    if strcmp(vmrfname{3}, '.VMR')
        v16fname = [vmrfname{1} '/' vmrfname{2} '.V16'];
    else
        v16fname = [vmrfname{1} '/' vmrfname{2} '.v16'];
    end
elseif ~ischar(v16fname) || numel(v16fname) < 4 || ~strcmpi(v16fname(end-3:end), '.v16')
    error('neuroelf:xff:badArgument', 'Invalid V16 filename argument.');
end
v16fname = v16fname(:)';

% check dimensions
if isempty(bc.VMRData16) || numel(size(bc.VMRData)) ~= numel(size(bc.VMRData16)) || ...
    any(size(bc.VMRData) ~= size(bc.VMRData16))
    error('neuroelf:xff:invalidObject', ...
        'Dimensions between VMRData and VMRData16 must match.');
end

% check transio case, in which the file IS saved
if istransio(bc.VMRData16) && ((ispc && strcmpi(filename(bc.VMRData16), v16fname)) || ...
    strcmp(filename(bc.VMRData16), v16fname))
    return;
end

% try to save V16
try
    vmr16 = [];
    vmr16 = xff('new:v16');
    vmr16c = vmr16.C;
    vmr16c.FileVersion = 1;
    vmr16c.DimX = size(bc.VMRData16, 1);
    vmr16c.DimY = size(bc.VMRData16, 1);
    vmr16c.DimZ = size(bc.VMRData16, 1);
    vmr16c.VMR8bit = false;
    vmr16c.VMRData = bc.VMRData16;
    vmr16.C = vmr16c;
    aft_SaveAs(vmr16, v16fname);

    % update transio reference if needed
    if istransio(vmr16c.VMRData)
        vmr16c = vmr16.C;
        bc.VMRData16 = vmr16c.VMRData;
        xo.C = bc;
    end
catch xfferror
    if ~isempty(vmr16)
        delete(vmr16);
    end
    error('neuroelf:xff:badObject', 'V16 file not writable: %s.', xfferror.message);
end
delete(vmr16);

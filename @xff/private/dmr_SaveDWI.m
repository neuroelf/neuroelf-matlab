function xo = dmr_SaveDWI(xo, totio)
% DMR::SaveDWI  - save memory-bound DWI data
%
% FORMAT:       [dmr] = dmr.SaveDWI([totio]);
%
% Input fields:
%
%       totio       if given and true convert data to transio object
%
% No output fields.
%
% Note: if the DWIData field is a transio, the function does nothing.

% Version:  v1.1
% Build:    16020110
% Date:     Feb-01 2016, 10:59 AM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'dmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~islogical(totio) || isempty(totio)
    totio = false;
else
    totio = totio(1);
end

% try saving the DWI
bc = xo.C;
if istransio(bc.DWIData)
    return;
end
[dfpn{1:3}] = fileparts(xo.F);
dsf = bc.DataStorageFormat;
nov = bc.NrOfVolumes;
nos = bc.NrOfSlices;
rsx = bc.ResolutionX;
rsy = bc.ResolutionY;
switch (dsf)
    case 2
        dds = [rsx, rsy, nov, nos];
    case 3
        dds = [rsx, rsy, nos, nov];
    case 4
        dds = [nov, rsx, rsy, nos];
    otherwise
        error('neuroelf:xff:invalidObject', 'Unsupported DataStorageFormat field.');
end
if ~isequal(size(bc.DWIData), dds) || ...
   (~isa(bc.DWIData, 'uint16') && ~isa(bc.DWIData, 'single') && ~istransio(bc.DWIData))
    error('neuroelf:xff:invalidObject', ...
        'Invalid DataStorageFormat, dimensions and/or datatype.');
end

% create filename
if any(dfpn{3}(2:end) == upper(dfpn{3}(2:end)))
    dfpn{3} = '.DWI';
else
    dfpn{3} = '.dwi';
end
dwif = [dfpn{1} '/' dfpn{2} dfpn{3}];
try
    dwii = fopen(dwif, 'w', 'ieee-le');
    if dwii < 1
        error('FILEOPEN_FAILED');
    end
    if bc.DataType == 1
        fwrite(dwii, bc.DWIData(:), 'uint16');
    else
        fwrite(dwii, bc.DWIData(:), 'single');
    end
    fclose(dwii);
catch xfferror
    neuroelf_lasterr(xfferror);
    error('neuroelf:xff:fileNotWritable', 'File not writable: ''%s''.', dwif);
end

% reload as transio?
if totio
    try
        if bc.DataType == 1
            bc.DWIData = transio(dwif, 'ieee-le', 'uint16', 0, dds);
        else
            bc.DWIData = transio(dwif, 'ieee-le', 'single', 0, dds);
        end
    catch xfferror
        error('neuroelf:xff:transioError', ...
            'Error reopening the written file as transio: %s.', xfferror.message);
    end
    xo.C = bc;
end

function xo = fmr_SaveSTC(xo, totio)
% FMR::SaveSTC  - save memory-bound STC data
%
% FORMAT:       [fmr] = fmr.SaveSTC([totio]);
%
% Input fields:
%
%       totio       if given and true convert data to transio object
%
% No output fields.
%
% Note: for STCData fields being transio, the function does nothing.

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:17 AM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~islogical(totio) || isempty(totio)
    totio = false;
else
    totio = totio(1);
end

% try saving the STC
bc = xo.C;
[dfpn{1:3}] = fileparts(xo.F);
dsf = bc.DataStorageFormat;
if ~any([1, 2] == dsf)
    error('neuroelf:xff:invalidObject', 'Unsupported DataStorageFormat field.');
end
nov = bc.NrOfVolumes;
nos = bc.NrOfSlices;
rsx = bc.ResolutionX;
rsy = bc.ResolutionY;
dds = [rsx, rsy, nov, nos];

% create filename template
if any(dfpn{3}(2:end) == upper(dfpn{3}(2:end)))
    stcx = '.STC';
else
    stcx = '.stc';
end
stctf = [dfpn{1} '/' bc.Prefix];

% for old format
if dsf == 1
    for tsc = 1:numel(bc.Slice)
        if istransio(bc.Slice(tsc).STCData)
            continue;
        end
        try
            stcf = sprintf('%s%d%s', stctf, tsc, stcx);
            stci = fopen(stcf, 'w', 'ieee-le');
            if stci < 1
                error('FILEOPEN_FAILED');
            end
            fwrite(stci, dds(1:2), 'uint16');
            if bc.DataType == 1
                fwrite(stci, bc.Slice(tsc).STCData(:), 'uint16');
            else
                fwrite(stci, bc.Slice(tsc).STCData(:), 'single');
            end
            fclose(dwii);
        catch xfferror
            error('neuroelf:xff:fileNotWritable', ...
                'File not writable: ''%s'' (%s).', dwif, xfferror.message);
        end

        % reload as transio?
        if totio
            try
                if bc.DataType == 1
                    bc.Slice(bc).STCData = transio(stctf, 'ieee-le', 'uint16', 4, dds(1:3));
                else
                    bc.Slice(bc).STCData = transio(stctf, 'ieee-le', 'single', 4, dds(1:3));
                end
            catch xfferror
                error('neuroelf:xff:transIOError', ...
                    'Error reopening the written file as transio: %s.', xfferror.message);
            end
        end
    end

    % set content back
    xo.C = bc;

% for new format
elseif ~istransio(bc.Slice.STCData)
    try
        stci = fopen([stctf stcx], 'w', 'ieee-le');
        if stci < 1
            error('FILEOPEN_FAILED');
        end
        if bc.DataType == 1
            fwrite(stci, bc.Slice.STCData(:), 'uint16');
        else
            fwrite(stci, bc.Slice.STCData(:), 'single');
        end
        fclose(stci);
    catch xfferror
        error('neuroelf:xff:fileNotWritable', ...
            'File not writable: ''%s'' (%s).', dwif, xfferror.message);
    end

    % reload as transio?
    if totio
        try
            if bc.DataType == 1
                bc.Slice.STCData = transio([stctf stcx], 'ieee-le', 'uint16', 0, dds);
            else
                bc.Slice.STCData = transio([stctf stcx], 'ieee-le', 'single', 0, dds);
            end
        catch xfferror
            error('neuroelf:xff:transIOError', ...
                'Error reopening the written file as transio: %s.', xfferror.message);
        end
    end

    % set content back
    xo.C = bc;
end

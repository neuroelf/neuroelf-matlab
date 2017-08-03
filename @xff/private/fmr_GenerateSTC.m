function xo = fmr_GenerateSTC(xo, stcprefix)
% FMR::GenerateSTC  - create empty STC files
%
% FORMAT:       [fmr] = fmr.GenerateSTC([stcprefix])
%
% Input fields:
%
%       stcprefix   if given, use this as new STC prefix
%
% Output fields:
%
%       fmr         altered FMR object
%
% Note: the FMR must be saved prior to calling this method!

% Version:  v1.1
% Build:    16021412
% Date:     Feb-14 2016, 12:58 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~ischar(stcprefix) || isempty(stcprefix) || numel(stcprefix) > 128 || ...
    any(stcprefix(:)' == '.' | stcprefix(:)' == '/' | stcprefix(:)' == '\')
    stcprefix = bc.Prefix(:)';
else
    stcprefix = stcprefix(:)';
end

% is FMR saved?
fmrfile = xo.F;
if isempty(fmrfile) || exist(fmrfile(:)', 'file') ~= 2
    error('neuroelf:xff:badArgument', 'The FMR object must be saved first.');
end
fmrpath = fileparts(fmrfile(:)');

% unload current slices
bc.Slice = [];
xo.C = bc;

% get some settings
dt = bc.NrOfVolumes;
dx = bc.ResolutionX;
dy = bc.ResolutionY;
dz = bc.NrOfSlices;

% depending on FileVersion -> < 4 multiple files
fv = bc.FileVersion;
if fv < 5 || bc.DataStorageFormat == 1

    % loop over slices
    for sc = 1:dz

        % generate STC file
        stcfile = sprintf('%s/%s%d.stc', fmrpath, stcprefix, sc);
        try
            stci = -1;
            stci = fopen(stcfile, 'w', 'ieee-le');
            if stci < 1
                error('FILE_OPEN_ERROR');
            end
            fwrite(stci, [dx, dy], 'uint16');
            fwrite(stci, zeros(dx * dy * dt, 1), 'uint16');
            fclose(stci);
            bc.Slice(sc).STCData = transio(stcfile, 'ieee-le', 'uint16', 4, [dx, dy, dt]);
        catch xfferror
            if stci > 1
                fclose(stci);
            end
            error('neuroelf:xff:errorSavingFile', ...
                'Error saving STC ''%s'': %s.', stcfile, xfferror.message);
        end
    end

% otherwise one file only
else

    % generate STC filename
    stcfile = sprintf('%s/%s.stc', fmrpath, stcprefix);

    % open and write content
    stcfid = fopen(stcfile, 'wb');
    if stcfid < 1
        error('neuroelf:xff:errorSavingFile', 'Error saving STC ''%s''.', stcfile);
    end

    % write content
    ne = dx * dy * dz * dt;
    fz = uint16(zeros(1048576,1));
    while ne > 1048576
        fwrite(stcfid, fz, 'uint16');
        ne = ne - 1048576;
    end
    if ne > 0
        fwrite(stcfid, fz(1:ne), 'uint16');
    end

    % close file
    fclose(stcfid);

    % check file size
    stcdir = dir(stcfile);
    if isempty(stcdir) || stcdir(1).bytes ~= (2 * dx * dy * dz * dt)
        error('neuroelf:xff:errorSavingFile', 'Error saving STC ''%s''.', stcfile);
    end

    % set in array
    switch (bc.DataStorageFormat)
        case 2
            tio = transio(stcfile, 'ieee-le', 'uint16', 0, [dx, dy, dt, dz]);
        case 3
            tio = transio(stcfile, 'ieee-le', 'uint16', 0, [dx, dy, dz, dt]);
        case 4
            tio = transio(stcfile, 'ieee-le', 'uint16', 0, [dt, dx, dy, dz]);
        otherwise
            error('neuroelf:xff:fieldUnsupported', 'Unsupported DataStorageFormat.');
    end
    bc.Slice.STCData = tio;
end

% set new prefix
bc.Prefix = stcprefix;

% set back
xo.C = bc;

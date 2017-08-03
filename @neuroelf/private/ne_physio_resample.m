% FUNCTION ne_physio_resample: resample selected channels
function ne_physio_resample(varargin)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-10 2011, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% try to get handle
try
    fPhysio = xfigure(get(varargin{1}, 'Parent'));
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    return;
end
ud = fPhysio.UserData;
hTag = ud.hTag;

% get data
try
    coi = ddeblank(hTag.ED_physio_coi.String);
    if ~strcmpi(coi, 'all')
        coi = evalin('base', coi);
    else
        coi = 0;
    end
    sfreq = str2double(hTag.ED_physio_sfreq.String);
    rfreq = str2double(hTag.ED_physio_rfreq.String);
    if numel(sfreq) ~= 1 || ...
        isinf(sfreq) || ...
        isnan(sfreq) || ...
        sfreq < 1 || ...
        sfreq > 1e6 || ...
        numel(rfreq) ~= 1 || ...
        isinf(rfreq) || ...
        isnan(rfreq) || ...
        rfreq < 1 || ...
        rfreq > 1e3 || ...
        isempty(coi) || ...
       ~isa(coi, 'double') || ...
        numel(coi) ~= max(size(coi)) || ...
        any(isinf(coi) | isnan(coi) | coi < 0 | coi ~= fix(coi))
        error('BAD_INPUT');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    uiwait(warndlg('Invalid setting; please correct.', ...
        'NeuroElf - input error', 'modal'));
    return;
end
ud.obj.RunTimeVars.Frequency = sfreq;

% channel selection
try
    if coi(1) ~= 0
        switch lower(ud.obj.Filetype)
            case {'acq'}
                coi(coi > ud.obj.NrOfChannels) = [];
                if ~isempty(ud.obj.RawData)
                    ud.obj.RawData = ud.obj.RawData(coi, :);
                end
                ud.obj.Channel = ud.obj.Channel(coi);
                ud.NrOfChannels = numel(coi);
            case {'mat', 'ntt'}
                if diff(size(ud.obj.Data)) <= 0
                    coi(coi > size(ud.obj.Data, 2)) = [];
                    ud.obj.Data = ud.obj.Data(:, coi);
                else
                    coi(coi > size(ud.obj.Data, 1)) = [];
                    ud.obj.Data = ud.obj.Data(coi, :);
                end
        end
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    uiwait(warndlg('Error selecting channels, please re-load and re-try!', ...
        'NeuroElf - Error', 'modal'));
end

% resample data
try
    if strcmpi(ud.obj.Filetype, 'acq')
        ud.obj.Resample(rfreq);
    else
        ud.obj.Resample(sfreq, rfreq);
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
ud.obj.RunTimeVars.Resampled = rfreq;

% disable resampled group
hTag.ED_physio_sfreq.String = sprintf(' %g', rfreq);
fPhysio.SetGroupEnabled('Resampled', 'off');

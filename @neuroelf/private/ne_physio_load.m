% FUNCTION ne_physio_load: load physio file into buffer object variable
function ne_physio_load(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:39 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% release any loaded file
ud = fPhysio.UserData;
if ~isempty(ud.obj) && ...
    isxff(ud.obj, true)
    ud.obj.ClearObject;
    ud.obj = [];
    fPhysio.UserData = ud;
    ud.hTag.ED_physio_sfreq.String = ' <from file>';
end

% disable groups
fPhysio.SetGroupEnabled('DLoaded', 'off');
fPhysio.SetGroupEnabled('DECG', 'off');
fPhysio.SetGroupEnabled('DGSR', 'off');
ud.hTag.TX_physio_ecgrtvinfo.String = 'RTV:';
drawnow;

% try to load object
try
    % get filename
    pf = ud.hTag.ED_physio_sf.String;
    pft = {};

    % if not an ACQ or MAT file
    if numel(pf) < 4 || ...
       ~any(strcmpi(pf(end-3:end), {'.acq', '.mat'}))

        % force load as NTT
        pft{1} = 'ntt';
    end

    % try loading
    ud.obj = xff(pf, pft{:});

% bail out on error
catch ne_eo;
    uiwait(warndlg(sprintf('Error loading requested file: %s', ne_eo.message), ...
        'NeuroElf - Error message', 'modal'));
    return;
end
ud.obj.RunTimeVars.Resampled = [];
rtv = ud.obj.RunTimeVars;

% store in UserData
fPhysio.UserData = ud;

% update info from file
if isfield(rtv, 'heartbeats') && ...
    isstruct(rtv.heartbeats)
    info = 'RTV: ';
    if isfield(rtv.heartbeats, 'bppre')
        info = sprintf('%s %d beats ', info(1:end), ...
            numel(rtv.heartbeats.bppre));
    end
    if isfield(rtv.heartbeats, 'freq')
        info = sprintf('%s f=%.1fHz ', info(1:end), ...
            rtv.heartbeats.freq);
    end
    ud.hTag.TX_physio_ecgrtvinfo.String = info(1:end-1);
end

% enable correct groups
fPhysio.SetGroupEnabled('DLoaded', 'on');
if ud.hTag.DD_physio_datakind.Value == 1
    fPhysio.SetGroupEnabled('DECG', 'on');
    ud.hTag.CB_physio_ecgrwchan.Value = 0;
    fPhysio.SetGroupEnabled('DECGRWC', 'off');
    if ud.hTag.CB_physio_ecgusertv.Value > 0
        if isfield(rtv, 'ECGChannel') && ...
            numel(rtv.ECGChannel) == 1
            ud.hTag.ED_physio_ecgschan.String = ...
                sprintf('%d', rtv.ECGChannel);
        end
        if isfield(rtv, 'LimitRange') && ...
            numel(rtv.LimitRange) == 2
            ud.hTag.ED_physio_ecglrmin.String = ...
                ddeblank(sprintf('%10.4g', rtv.LimitRange(1)));
            ud.hTag.ED_physio_ecglrmax.String = ...
                ddeblank(sprintf('%10.4g', rtv.LimitRange(2)));
        end
        if isfield(rtv, 'PosHelp') && ...
            ischar(rtv.PosHelp) && ...
            any(strcmp(rtv.PosHelp, ud.hTag.DD_physio_ecgmposhlp.String))
            ud.hTag.DD_physio_ecgmposhlp.Value = find(strcmp( ...
                rtv.PosHelp, ud.hTag.DD_physio_ecgmposhlp.String));
        end
        if isfield(rtv, 'RWChannel') && ...
            numel(rtv.RWChannel) == 1
            ud.hTag.ED_physio_ecgrwchan.String = ...
                sprintf('%d', rtv.RWChannel);
        end
        if isfield(rtv, 'RWCode') && ...
            numel(rtv.RWCode) == 1
            ud.hTag.ED_physio_ecgrwcode.String = ...
                sprintf('%d', rtv.RWCode);
        end
    end
else
    fPhysio.SetGroupEnabled('DGSR', 'on');
end

% get frequency
freq = [];
if isfield(rtv, 'Frequency') && ...
    numel(rtv.Frequency) == 1 && ...
    isa(rtv.Frequency, 'double')
    freq = rtv.Frequency;
else
    switch lower(ud.obj.Filetype)
        case {'acq'}
            freq = 1000 / acq.MilliSecsPerSample;
        case {'mat'}
            if isfield(ud.obj.Info, 'isi')
                if isfield(ud.obj.Info, 'isi_units') && ...
                    strcmpi(ud.obj.Info.isi_units, 'ms')
                    freq = 1000 / ud.obj.Info.isi;
                else
                    freq = 1 / ud.obj.Info.isi;
                end
            elseif isfield(ud.obj.Info, 'freq')
                freq = ud.obj.Info.freq;
            end
        case {'ntt'}
    end
    ud.obj.RunTimeVars.Frequency = freq;
    ud.obj.RunTimeVars.Resampled = freq;
end
if numel(freq) == 1 && ...
    isnumeric(freq)
    ud.hTag.ED_physio_sfreq.String = sprintf(' %g', freq);
else
    ud.hTag.ED_physio_sfreq.String = ' <specify!>';
end

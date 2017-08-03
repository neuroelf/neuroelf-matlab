% FUNCTION ne_physio_process: process physio data
function ne_physio_process(varargin)

% Version:  v0.9b
% Build:    11052514
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

% global variable
global ne_gcfg;

% try to get handle
try
    fPhysio = xfigure(get(varargin{1}, 'Parent'));
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    return;
end
ud = fPhysio.UserData;
hTag = ud.hTag;
obj = ud.obj;
if ~isempty(obj.RunTimeVars.Resampled)
    freq = obj.RunTimeVars.Resampled;
else
    freq = obj.RunTimeVars.Frequency;
end

% build title
[null, objtitle, objext] = fileparts(obj.FilenameOnDisk);
objtitle = [objtitle objext];


% depending on type of computation
switch (hTag.DD_physio_datakind.Value)

    % ECG
    case {1}

        % prepare options for first and second pass
        try
            ecgchan = str2double(hTag.ED_physio_ecgschan.String);
            ecgcalc = 'none';
            switch hTag.DD_physio_ecgptrans.Value
                case {2}
                    ecgcalc = 'squarez';
                case {3}
                    ecgcalc = 'thirdz';
                case {4}
                    ecgcalc = 'fourthz';
                case {5}
                    ecgcalc = 'absdiff';
                case {6}
                    ecgcalc = 'diffsq';
                case {7}
                    ecgcalc = {'none', 'squarez', 'thirdz', 'fourthz', ...
                               'absdiff', 'diffsq'};
            end
            ecgpfreps = hTag.DD_physio_ecgpfiltp.Value - 1;
            ecgpflength = str2double(hTag.ED_physio_ecgpfiltk.String);
            lrmin = str2double(hTag.ED_physio_ecglrmin.String);
            lrmax = str2double(hTag.ED_physio_ecglrmax.String);
            if numel(lrmin) ~= 1 || ...
                isnan(lrmin) || ...
                numel(lrmax) ~= 1 || ...
                isnan(lrmax) || ...
                lrmin >= lrmax
                error( ...
                    'neuroelf:BadInput', ...
                    'Invalid limit range input.' ...
                );
            end
            limrange = [lrmin, lrmax];
            ecgcleanup = (hTag.CB_physio_ecgmcleaup.Value > 0);
            if numel(ecgchan) ~= 1 || ...
               ~isa(ecgchan, 'double') || ...
                isinf(ecgchan) || ...
                isnan(ecgchan) || ...
                ecgchan < 1 || ...
                ecgchan ~= fix(ecgchan)
                error('BAD_CHANNEL');
            end
            if numel(ecgpflength) ~= 1 || ...
               ~isa(ecgpflength, 'double') || ...
                isinf(ecgpflength) || ...
                isnan(ecgpflength) || ...
                ecgpflength < 0 || ...
                ecgpflength > 1
                error('BAD_PREFILTER_KERNEL');
            end
            if hTag.CB_physio_ecgsavertv.Value > 0
                ud.obj.RunTimeVars.ECGChannel = ecgchan;
            end
            try
                ecgchan = ud.obj.ChannelData(ecgchan);
            catch ne_eo;
                rethrow(ne_eo);
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            uiwait(warndlg( ...
                'Error in ECG configuration. Please correct and try again.', ...
                'NeuroElf - input error', 'modal'));
            return;
        end

        % bppre available from previous session
        if hTag.CB_physio_ecgusertv.Value > 0 && ...
            isfield(obj.RunTimeVars, 'heartbeats') && ...
            isstruct(obj.RunTimeVars.heartbeats) && ...
            isfield(obj.RunTimeVars.heartbeats, 'bppre') && ...
           ~isempty(obj.RunTimeVars.heartbeats.bppre) && ...
            isfield(obj.RunTimeVars.heartbeats, 'freq') && ...
            obj.RunTimeVars.heartbeats.freq == freq

            % get beats from there and plot
            bp = heartbeats(ecgchan, struct( ...
                'bppre',    obj.RunTimeVars.heartbeats.bppre, ...
                'correct',  false, ...
                'freq',     freq, ...
                'limrange', limrange, ...
                'plot',     ecgcleanup, ...
                'title',    objtitle));

        % R-wave code not in file
        elseif hTag.CB_physio_ecgrwchan.Value <= 0

            % first pass
            try
                ne_gcfg.h.MainFig.Pointer = 'watch';
                fPhysio.Pointer = 'watch';
                drawnow;
                if iscell(ecgcalc)
                    calco = zeros(numel(ecgcalc), 2);
                    for calcc = 1:numel(ecgcalc)
                        bp = heartbeats(ecgchan, struct( ...
                            'calc',     ecgcalc{calcc}, ...
                            'freq',     freq, ...
                            'limrange', limrange, ...
                            'pflength', ecgpflength, ...
                            'pfreps',   ecgpfreps, ...
                            'plot',     false, ...
                            'title',    objtitle));
                        bpd = diff(bp);
                        if numel(bpd) > 3
                            calco(calcc, 1) = mean(bpd) / std(bpd);
                            bph = computehrv(bp, struct('freq', freq));
                            calco(calcc, 2) = bph.rmssd;
                        else
                            calco(calcc, 2) = 1;
                        end
                    end
                    calco = calco(:, 1) ./ calco(:, 2);
                    ecgcalc = ecgcalc{maxpos(calco)};
                    disp(sprintf('Auto-detected calculation type for %s: %s', ...
                        ud.obj.FilenameOnDisk, ecgcalc));
                end
                bp = heartbeats(ecgchan, struct( ...
                    'calc',     ecgcalc, ...
                    'freq',     freq, ...
                    'limrange', limrange, ...
                    'pflength', ecgpflength, ...
                    'pfreps',   ecgpfreps, ...
                    'plot',     ecgcleanup, ...
                    'title',    objtitle));
            catch ne_eo;
                ne_gcfg.h.MainFig.Pointer = 'arrow';
                fPhysio.Pointer = 'arrow';
                drawnow;
                uiwait(warndlg(['An error occurred while using heartbeats: ' ...
                    ne_eo.message], 'NeuroElf - error message', 'modal'));
                return;
            end

        % otherwise detect from file
        else
            try
                rwc = str2double(hTag.ED_physio_ecgrwchan.String);
                rwv = str2double(hTag.ED_physio_ecgrwcode.String);
                if numel(rwc) ~= 1 || ...
                    isinf(rwc) || ...
                    isnan(rwc) || ...
                    rwc < 1 || ...
                    rwc ~= fix(rwc)
                    error('Invalid R-wave channel.');
                end
                if numel(rwv) ~= 1
                    error('Invalid R-wave code.');
                end
                if isinf(rwv) || ...
                    isnan(rwv)
                    rwvu = unique(ud.obj.ChannelData(rwc));
                    if numel(rwvu) > 16
                        error('Non-unique data in R-wave channel.');
                    end
                    rwvn = zeros(numel(rwvu), 2);
                    for rwvc = 1:numel(rwvu)
                        rwvp = find(ud.obj.ChannelData(rwc) == rwvu(rwvc));
                        rwvn(rwvc, :) = ...
                            [numel(rwvp), mean(diff(rwvp)) / std(diff(rwvp))];
                    end
                    rwvn = rwvn(:, 1) .* rwvn(:, 2);
                    rwv = rwvu(maxpos(rwvn));
                end
                if hTag.CB_physio_ecgsavertv.Value > 0
                    ud.obj.RunTimeVars.RWChannel = rwc;
                    ud.obj.RunTimeVars.RWCode = rwv;
                end
                bp = find(ud.obj.ChannelData(rwc) == rwv);
                bp = heartbeats(ecgchan, struct( ...
                    'bppre', bp, ...
                    'freq',  freq, ...
                    'limrange', limrange, ...
                    'plot',  ecgcleanup, ...
                    'title', objtitle));
            catch ne_eo;
                ne_gcfg.h.MainFig.Pointer = 'arrow';
                fPhysio.Pointer = 'arrow';
                drawnow;
                uiwait(warndlg(['An error occurred while extracting the R-wave markers: ' ...
                    ne_eo.message], 'NeuroElf - error message', 'modal'));
                return;
            end
        end

        % second pass
        ne_physio_outnames = {'bp', 'bs', 'bf', 'bv', 'cp', 'wgd', 'wd'};
        if hTag.CB_physio_ecghrv.Value > 0
            ne_physio_out = cell(1, 8);
            ne_physio_outnames{end+1} = 'hrv';
        else
            ne_physio_out = cell(1, 7);
        end
        [ne_physio_out{1:numel(ne_physio_out)}] = heartbeats(ecgchan, struct( ...
            'bppre',   bp, ...
            'correct', false, ...
            'cleanup', ecgcleanup, ...
            'limrange', limrange, ...
            'freq',    freq, ...
            'plot',    ecgcleanup, ...
            'poshlp',  hTag.DD_physio_ecgmposhlp.String{hTag.DD_physio_ecgmposhlp.Value}, ...
            'title',   objtitle));

        % store information in object
        if hTag.CB_physio_ecgsavertv.Value > 0
            obj.RunTimeVars.PosHelp = ...
                hTag.DD_physio_ecgmposhlp.String{hTag.DD_physio_ecgmposhlp.Value};
            obj.RunTimeVars.heartbeats.bppre = ne_physio_out{1};
            obj.RunTimeVars.heartbeats.freq = freq;
            try
                obj.SaveRunTimeVars;
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                uiwait(warndlg(['Error saving RunTimeVars of ' ...
                    obj.FilenameOnDisk '.'], 'Neuroelf - user information', 'modal'));
            end
        end

    % GSR
    case {2}
end
ne_gcfg.h.MainFig.Pointer = 'arrow';
fPhysio.Pointer = 'arrow';
drawnow;

% pass on output arguments
for outc = 1:numel(ne_physio_out)
    ud.obj.RunTimeVars.(ne_physio_outnames{outc}) = ne_physio_out{outc};
    try
        assignin('base', ne_physio_outnames{outc}, ne_physio_out{outc});
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        warning( ...
            'neuroelf:AssignIn', ...
            'Couldn''t assign ''%s'' in BASE workspace.', ...
            ne_physio_outnames{outc} ...
        );
    end
end

% close dialog (using the waitfor in the main function)
fPhysio.Visible = 'off';

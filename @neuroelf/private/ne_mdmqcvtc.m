% FUNCTION ne_mdmqcvtc: create QC-VTC from files referenced in MDM
function varargout = ne_mdmqcvtc(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:37 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2016, Jochen Weber
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
mainfig = ne_gcfg.h.MainFig;
mfp = mainfig.Pointer;
pbar = ne_gcfg.h.Progress;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% load MDM
mdm = {[]};
try
    mainfig.Pointer = 'watch';
    drawnow;
    mdm{1} = xff('*.mdm', 'Please select MDM containing time-course references...');
    if ~isxff(mdm{1}, 'mdm') || ...
        isempty(mdm{1}.XTC_RTC) || ...
        size(mdm{1}.XTC_RTC, 2) ~= 2 || ...
       ~strcmpi(mdm{1}.TypeOfFunctionalData, 'vtc')
        error('BAD_MDM');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    if ~isempty(mdm{1})
        uiwait(warndlg('Invalid input file selected.', ...
            'NeuroElf - Error', 'modal'));
    end
    mainfig.Pointer = mfp;
    drawnow;
    clearxffobjects(mdm);
    return;
end
mdmf = mdm{1}.FilenameOnDisk;
xtc = mdm{1}.XTC_RTC(:, 1);
mdm{1}.ClearObject;

% load first object
tc = {[]};
try
    tc{1} = xff(xtc{1});
    if ~isxff(tc{1}, {'hdr', 'vtc'})
        error('BAD_TC');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    clearxffobjects(tc);
    uiwait(warndlg('Invalid or unfound time-course file in MDM.', ...
        'NeuroElf - Error', 'modal'));
    mainfig.Pointer = mfp;
    drawnow;
    return;
end

% if no VTC, sample as VTC
if ~isxff(tc{1}, 'vtc')
    vtc = importvtcfromanalyze(tc{1});
    tc{1}.ClearObject;
    tc{1} = vtc;
else
    if istransio(tc{1}.VTCData)
        tc{1}.VTCData = resolve(tc{1}.VTCData);
    end
    vtc = tc{1}.CopyObject;
    tc{1}.ClearObject;
    tc{1} = vtc;
end

% ask for type of averaging
atype = questdlg('Which volumes to store?', 'NeuroElf - request', ...
    'Mean only', '1st, mean, & last', 'Cancel', 'Mean only');
if ~ischar(atype) || ...
    isempty(atype) || ...
    strcmpi(atype, 'cancel')
    tc{1}.ClearObject;
    mainfig.Pointer = mfp;
    drawnow;
    return;
end
atype = ~(lower(atype(1)) == 'm');

% progress initialization
[tcp, tcf] = fileparts(xtc{1});
cprog = ne_progress(0, 0, {true, 0, sprintf('Averaging TC #1 (%s)...', tcf)});

% build averages
mvtc = mean(tc{1}.VTCData);
if atype
    fvtc = tc{1}.VTCData(1, :, :, :);
    lvtc = tc{1}.VTCData(end, :, :, :);
end

% initialize VTC
vtc = tc{1};
tc = {[]};
vsz = size(vtc.VTCData);
vtc.FileVersion = 3;
vtc.NameOfSourceFMR = mdmf;
vtc.NrOfLinkedPRTs = 1;
vtc.NameOfLinkedPRT = '';
vtc.DataType = 2;
vtc.NrOfVolumes = (1 + 2 * double(atype)) * numel(xtc);
vtc.VTCData = single(0);
vtc.VTCData(vtc.NrOfVolumes, vsz(2), vsz(3), vsz(4)) = 0;

% start with first
if atype
    vtc.VTCData(1, :, :, :) = fvtc;
    vtc.VTCData(2, :, :, :) = mvtc;
    vtc.VTCData(3, :, :, :) = lvtc;
else
    vtc.VTCData(1, :, :, :) = mvtc;
end

% loop over objects
for tcc = 2:numel(xtc)
    try
        tc{1} = xff(xtc{tcc});
        if ~isxff(tc{1}, {'hdr', 'vtc'})
            error('BAD_TC');
        end
        if isxff(tc{1}, 'hdr')
            tvtc = importvtcfromanalyze(tc{1});
            tc{1}.ClearObject;
            tc{1} = tvtc;
        elseif istransio(tc{1}.VTCData)
            tc{1}.VTCData = resolve(tc{1}.VTCData);
        end
        [tcp, tcf] = fileparts(xtc{tcc});
        pbar.Progress(tcc / numel(xtc), sprintf('Averaging TC #%d (%s)...', tcc, tcf));
        mvtc = mean(tc{1}.VTCData);
        if atype
            vtc.VTCData(3 * tcc - 2, :, :, :) = tc{1}.VTCData(1, :, :, :);
            vtc.VTCData(3 * tcc - 1, :, :, :) = mvtc;
            vtc.VTCData(3 * tcc, :, :, :) = tc{1}.VTCData(end, :, :, :);
        else
            vtc.VTCData(tcc, :, :, :) = mvtc;
        end
        tc{1}.ClearObject;
        tc{1} = [];
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        vtc.ClearObject;
        clearxffobjects(tc);
        uiwait(warndlg('Invalid or unfound time-course file in MDM.', ...
            'NeuroElf - Error', 'modal'));
        mainfig.Pointer = mfp;
        ne_progress(0, 0, cprog);
        return;
    end
end

% open average VTC
ne_openfile(0, 0, vtc, true);

% reset pointer and progress bar
mainfig.Pointer = mfp;
ne_progress(0, 0, cprog);

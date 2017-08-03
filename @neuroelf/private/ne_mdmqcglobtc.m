% FUNCTION ne_mdmqcvtc: create QC-global TC from files referenced in MDM
function varargout = ne_mdmqcglobtc(varargin)

% Version:  v1.0
% Build:    16010821
% Date:     Jan-08 2016, 9:29 PM EST
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
gtc = cell(numel(xtc), 1);
mdm{1}.ClearObject;

% ask for type of transformation
ttype = questdlg('Apply transformation?', 'NeuroElf - request', ...
    'none', '% change (PSC)', 'z-transform', 'none');
if ~ischar(ttype) || ...
    isempty(ttype)
    mainfig.Pointer = mfp;
    drawnow;
    return;
end
ttype = lower(ttype(1));

% initialize progress bar
cprog = ne_progress(0, 0, {true, 0, 'Averaging TCs'});

% loop over objects
tc = {[]};
for tcc = 1:numel(xtc)
    try
        [tcp, tcf] = fileparts(xtc{tcc});
        pbar.Progress((tcc - 1)/ numel(xtc), ...
            sprintf('Averaging TC #%d (%s)...', tcc, tcf));
        tc{1} = xff(xtc{tcc});
        if ~isxff(tc{1}, {'hdr', 'vtc'})
            error('BAD_TC');
        end
        if isxff(tc{1}, 'hdr')
            mtc = lsqueeze(mean(mean(mean(tc{1}.VoxelData(:, :, :, :), 3), 2), 1));
        else
            mtc = mean(mean(mean(tc{1}.VTCData(:, :, :, :), 4), 3), 2);
        end
        tc{1}.ClearObject;
        tc{1} = [];
        if ttype == '%'
            mtc = psctrans(mtc);
        elseif ttype == 'z'
            mtc = ztrans(mtc);
        end
        gtc{tcc} = mtc;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        clearxffobjects(tc);
        uiwait(warndlg('Invalid or unfound time-course file in MDM.', ...
            'NeuroElf - Error', 'modal'));
        mainfig.Pointer = mfp;
        ne_progress(0, 0, cprog);
        return;
    end
end

% concatenate time course
ctc = cat(1, gtc{:});
if ttype == '%'
    ctc = ctc - 100;
end
mctc = minmaxmean(ctc);

% assign in base workspace
assignin('base', 'mtc_global_tc', ctc);
assignin('base', 'mtc_global_tcs', gtc);

% plot time course
f = figure;
set(f, 'NumberTitle', 'off', ...
    'Name', sprintf('NeuroElf - global time course - %s', mdmf));
x = axes;
plot(1:numel(ctc), ctc);

% plot subject markers
hold(x, 'on');
tcc = 0;
lastid = '::INVALID::';
for cc = 1:numel(gtc)
    [tcp, tcf] = fileparts(xtc{cc});
    newid = regexprep(tcf, '^([^_]+)_.*$', '$1');
    if ~strcmpi(newid, lastid)
        lastid = newid;
        plot(tcc .* ones(2, 1), mctc(1:2)', 'r-');
    else
        plot(tcc .* ones(2, 1), mctc(1:2)', 'r--');
    end
    tcc = tcc + numel(gtc{cc});
end
plot(tcc .* ones(2, 1), mctc(1:2)');

% reset pointer and progress bar
mainfig.Pointer = mfp;
ne_progress(0, 0, cprog);

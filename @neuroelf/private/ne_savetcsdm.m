% FUNCTION ne_savetcsdm: saves the displayed timecourse as a SDM (or RTC)
function ne_savetcsdm(varargin)
%
% FORMAT:       ne_savetcsdm(src, evt [, filename])
%
% Input fields:
%
%       src         handle of UI object issuing the call (0 for GUI)
%       evt         event data (currently unused)
%       filename    optional filename (otherwise requested)

% Version:  v0.9c
% Build:    11051414
% Date:     May-14 2011, 2:32 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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

% time course data
tcd = ne_gcfg.fcfg.tcplotdata(:);
if isempty(tcd)
    return;
end

% transform accordingly
tcn = ne_gcfg.c.ini.Statistics.TCnorm;
if ~isempty(tcn) && ...
    lower(tcn(1)) == 'p'
    tcd = psctrans(tcd);
    if mean(tcd) > 50
        tcd = tcd - 100;
    end
elseif ~isempty(tcn) && ...
    lower(tcn(1)) == 'z'
    tcd = ztrans(tcd);
end

% output filename
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
    isempty(regexpi(varargin{3}(:)', '\.(rtc|sdm|txt)$'))
    tlist = {'*.sdm', 'Single-study design matrix (*.sdm)'; ...
         '*.rtc', 'Reference time-course (*.rtc)'; ...
         '*.txt', 'Plain text file (*.txt)'};
    [savefile, savepath, saveidx] = uiputfile(tlist, ...
        'Save time course data as...');
    if isequal(savefile, 0) || ...
        isequal(savepath, 0) || ...
        isempty(savefile)
        return;
    end
    if isempty(savepath)
        savepath = pwd;
    end
    [nullpath, savefile, saveext] = fileparts(savefile);
    if isempty(saveext) || ...
       ~any(strcmpi(saveext, {'.rtc', '.sdm', '.txt'}))
        saveext = {'.sdm', '.rtc', '.txt'};
        saveext = saveext{saveidx};
    end
    savefile = [savepath, '/', savefile, saveext];
else
    savefile = varargin{4}(:)';
end

% try to save text
sdm = [];
try
    if ~strcmpi(savefile(end-2:end), 'txt')
        slvar = ne_gcfg.fcfg.SliceVar;
        if numel(slvar) == 1 && ...
            isxff(slvar, {'dmr', 'fmr', 'hdr', 'head', 'vtc'})
            slfile = slvar.FilenameOnDisk;
        else
            slfile = 'extracted';
        end
        [slpath, slfile, slext] = fileparts(slfile);
        slfile = [slfile slext];
        sdm = xff('new:sdm');
        sdm.NrOfPredictors = 1;
        sdm.NrOfDataPoints = numel(tcd);
        sdm.IncludesConstant = 0;
        sdm.FirstConfoundPredictor = 2;
        sdm.PredictorNames = {slfile};
        sdm.PredictorColors = floor(255.999 .* rand(1, 3));
        sdm.SDMMatrix = tcd;
        sdm.RTCMatrix = tcd;
        if strcmpi(savefile(end-2:end), 'rtc')
            sdm.FileVersion = 2;
        end
        sdm.SaveAs(savefile);
        sdm.ClearObject;
        sdm = [];
    else
        save(savefile, tcd, '-ascii');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg(['Error writing text output file: ', ne_eo.message], ...
        'NeuroElf - error', 'modal'));
end
if numel(sdm) == 1 && ...
    isxff(sdm, true)
    sdm.ClearObject;
end

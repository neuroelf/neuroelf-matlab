% FUNCTION ne_mdm_ffiles: find MDM relevant files
function ne_mdm_ffiles(varargin)

% Version:  v0.9d
% Build:    14061710
% Date:     Jun-17 2014, 10:04 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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
cc = ne_gcfg.fcfg.MDM;
cf = ne_gcfg.h.MDM.MDMFig;
ch = ne_gcfg.h.MDM.h;

% get base folder and patterns
bf = ch.Basefolder.String;
fpat = ddeblank(ch.FuncPattern.String);
dpat = ddeblank(ch.DsgnPattern.String);
mpat = ddeblank(ch.MParPattern.String);
userps = (ch.UseMotParms.Value > 0);
[ffld, fpat, fext] = fileparts(fpat);
if ~strcmpi(fext, '.vtc')
    uiwait(warndlg('Currently only VTC files are supported.', ...
        'NeuroElf - user information', 'modal'));
    return;
end
[dfld, dpat, dext] = fileparts(dpat);
if ~any(strcmpi(dext, {'.prt', '.rtc', '.sdm'}))
    uiwait(warndlg('Currently only PRT/RTC/SDM files are supported.', ...
        'NeuroElf - user information', 'modal'));
    return;
end
[mfld, mpat, mext] = fileparts(mpat);
if userps && ...
   ~any(strcmpi(mext, {'.rtc', '.sdm', '.txt'}))
    uiwait(warndlg('Currently only RTC/SDM/TXT files are supported.', ...
        'NeuroElf - user information', 'modal'));
    return;
end

% clear current selection
cc.mdm.XTC_RTC = cell(0, 2);
cc.mdm.NrOfStudies = 0;
cc.mdm.RunTimeVars.MotionParameters = {};
cf.SetGroupEnabled('FFound', 'off');
cf.SetGroupEnabled('PRTDsgn', 'off');
ch.FuncFiles.String = {'<no files selected>'};
ch.FuncFiles.Value = [];
ch.DsgnFiles.String = {'<no files selected>'};
ch.DsgnFiles.Value = [];
ch.MParFiles.String = {'<no files selected>'};
ch.MParFiles.Value = [];

% set pointer
cf.Pointer = 'watch';
drawnow;

% with error handling
try

    % locate files
    fpat = [fpat fext];
    if isempty(ffld)
        vtcs = findfiles(bf, fpat);
    else
        vtcs = findfiles([bf '/' ffld], fpat, 'depth=1');
    end
    numstudy = numel(vtcs);
    if numstudy == 0
        uiwait(warndlg('No VTC files found using this folder/pattern.', ...
            'NeuroElf - user information', 'modal'));
        cf.Pointer = 'arrow';
        return;
    end
    dpat = [dpat dext];
    if isempty(dfld)
        dess = findfiles(bf, dpat);
    else
        dess = findfiles([bf '/' dfld], dpat, 'depth=1');
    end
    if isempty(dess)
        uiwait(warndlg('No design files found using this folder/pattern.', ...
            'NeuroElf - user information', 'modal'));
        cf.Pointer = 'arrow';
        return;
    end
    if numstudy ~= numel(dess)
        uiwait(warndlg(sprintf('Number of VTC and design files mismatch: %d ~= %d.', ...
            numel(vtcs), numel(dess)), 'NeuroElf - user information', 'modal'));
        cf.Pointer = 'arrow';
        return;
    end
    sids = regexprep(vtcs, '.*[\\\/]([^\\\/]+)$', '$1');
    sids = regexprep(sids, '_.*$', '');
    sidb = false(1, numstudy);
    for sc = 1:numstudy
        if isempty(strfind(dess{sc}, sids{sc}))
            sidb(sc) = true;
        end
    end
    if any(sidb)
        bsid = sprintf('\n%s', dess{find(sidb)});
        uiwait(warndlg(sprintf( ...
            'Some design files seem to mismatch the subject:\n%s\nPlease confirm!', bsid), ...
            'NeuroElf - user information', 'modal'));
    end
    if userps
        mpat = [mpat mext];
        if isempty(mfld)
            rpss = findfiles(bf, mpat);
        else
            rpss = findfiles([bf '/' mfld], mpat, 'depth=1');
        end
        if isempty(rpss)
            uiwait(warndlg('No motion parameters found using this folder/pattern.', ...
                'NeuroElf - user information', 'modal'));
            userps = false;
        elseif numel(rpss) ~= numstudy
            uiwait(warndlg(sprintf( ...
                'Number of VTC and motion parameters files mismatch: %d ~= %d.', ...
                numstudy, numel(rpss)), 'NeuroElf - user information', 'modal'));
            cf.Pointer = 'arrow';
            return;
        end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    cf.Pointer = 'arrow';
    return;
end

% re-set pointer and groups
cf.Pointer = 'arrow';
cf.SetGroupEnabled('FFound', 'on');
if ~any(cellfun('isempty', regexpi(dess, '\.prt$')))
    cf.SetGroupEnabled('PRTDsgn', 'on');
end

% put files into listboxes and MDM
ch.FuncFiles.String = vtcs;
ch.DsgnFiles.String = dess;
if userps
    ch.MParFiles.String = rpss;
else
    ch.UseMotParms.Value = 0;
    cf.SetGroupEnabled('MotParm', 'off');
end

% FUNCTION ne_mdm_showopt: show advanced options
function ne_mdm_showopt(varargin)

% Version:  v1.1
% Build:    16053116
% Date:     May-31 2016, 4:09 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, 2014, 2016, Jochen Weber
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

% check preliminaries
if ~isfield(ne_gcfg.h, 'MDM') || ~isstruct(ne_gcfg.h.MDM) || ~isfield(ne_gcfg.h.MDM, 'MDMFig') || ...
   ~isxfigure(ne_gcfg.h.MDM.MDMFig, true) || ~isfield(ne_gcfg.fcfg, 'MDM') || ...
   ~isstruct(ne_gcfg.fcfg.MDM) || ~isfield(ne_gcfg.fcfg.MDM, 'mdm') || ~isxff(ne_gcfg.fcfg.MDM.mdm, 'mdm')
    return;
end
hFig = ne_gcfg.h.MDM.MDMFig;
cc = ne_gcfg.fcfg.MDM;
ch = ne_gcfg.h.MDM.h;
mdm = cc.mdm;

% main page again
if nargin > 2 && islogical(varargin{3}) && numel(varargin{3}) == 1 && ~varargin{3}
    hFig.ShowPage(1);
    return;
end

% make sure all design files are PRTs
prts = mdm.XTC_RTC(:, end);
if isempty(prts) || any(cellfun('isempty', regexpi(prts, '\.prt$')))
    uiwait(warndlg('Extended options requires all design files to be PRTs.', ...
        'NeuroElf - information', 'modal'));
    return;
end

% show second page
hFig.ShowPage(2);

% no more work to be done
if ~isempty(ch.PRTConds.UserData)
    return;
end

% set cursor to wait
hfp = hFig.Pointer;
hFig.Pointer = 'watch';
drawnow;

% get conditions
conds = cell(numel(prts, 1));
condcols = conds;
for pc = 1:numel(prts)

    % read PRT
    try
        prt = [];
        prt = xff(prts{pc});
        if ~isxff(prt, 'prt')
            error('neuroelf:xffError:invalidPRT', 'Invalid PRT file for study %d.', pc);
        end
        [conds{pc}, condcols{pc}] = prt.ConditionNames;
        prt.ClearObject;
    catch ne_eo;
        if ~isempty(prt) && isxff(prt, true)
            prt.ClearObject;
        end
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(errordlg(ne_eo.message, 'NeuroElf - error', 'modal'));
        hFig.ShowPage(1);
        hFig.Pointer = hfp;
        drawnow;
        return;
    end
end

% consolidate conds
conds = cat(1, conds{:});
condcols = cat(1, condcols{:});
dconds = false(numel(conds), 1);
for cc = 2:numel(conds)
    dconds(cc:end) = dconds(cc:end) | strcmpi(conds(cc:end), conds{cc-1});
end
conds(dconds) = [];
condcols(dconds, :) = [];
ch.PRTConds.String = conds;
ch.PRTConds.UserData = conds;
ccol = condcols;
dcol = round(0.333 .* double(ccol));
ne_gcfg.fcfg.MDM.condcols = [dcol; ccol; 255 - dcol; 255 - ccol];

% set control states
hFig.RadioGroupSetOne('Deriv', 1);
hFig.SetGroupEnabled('Derivs', 'off');
%ch.VTCAvgCondCols.Callback = {@ne_mdm_acond, 'setcolors'};
ch.VTCAvgCondCols.Enable = 'on';

% re-set pointer
hFig.Pointer = hfp;
drawnow;

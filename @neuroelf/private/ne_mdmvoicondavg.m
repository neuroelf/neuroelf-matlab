function varargout = ne_mdmvoicondavg(varargin)
% ne_mdmvoicondavg  - initialize MDM-based VOI cond average plots
%
% FORMAT:       [x, f, fid, data] = ne_mdmvoicondavg(SRC, EVT [, reload [, file [, conds]]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       reload      if set to 'reload', use MAT file instead of MDM
%       file        either an MDM or MAT file (see reload flag)
%       conds       for MDM-mode, set conditions (otherwise per UI)
%
% Output fields:
%
%       x           axes handle (MATLAB, not xff)
%       f           figure handle (MATLAB, not xff)
%       fid         figure ID (used as fieldname for ne_gcfg.cc)
%       data        mvc_data structure (saved in MAT file)
%
% Examples:
%
%   [h_ax, h_fig, fid, data] = ne_mdmvoicondavg(0, 0, 'mdm', 'study_prts.mdm', ...
%       {'reapp_neg', 'look_neg', 'look_neu'});
%   neuroelf_gui('screenshot', h_fig, 'high-q');
%   delete(h_fig);

% Version:  v1.1
% Build:    16021715
% Date:     Feb-17 2016, 3:17 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2015, 2016, Jochen Weber
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
cini = ne_gcfg.c.ini;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% reloading
voi = {[]};
if nargin > 2 && ...
    ischar(varargin{3}) && ...
    strcmpi(varargin{3}(:)', 'reload')
    reload = true;

% creating
else
    reload = false;

    % only valid if VOI active
    voi{1} = ne_gcfg.voi;
    if numel(voi{1}) ~= 1 || ...
       ~isxff(voi{1}, 'voi') || ...
        isempty(voi{1}.VOI)
        uiwait(warndlg('No VOI/clusters defined.', 'Neuroelf GUI - info', 'modal'));
        return;
    end
end

% load MDM
mdm = {[]};
try
    if reload
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{4}) || ...
            isempty(regexpi(varargin{4}(:)', '\.mat$')) || ...
            exist(varargin{4}(:)', 'file') ~= 2
            [mvc_data, mvc_path] = uigetfile( ...
                {'*.mat', 'MDM voi-extracted MAT-files (*.mat)'}, ...
                'Please select the datafile containing the extracts...');
            if isequal(mvc_path, 0) || ...
                isequal(mvc_data, 0) || ...
                isempty(mvc_data)
                return;
            end
        else
            [mvc_path, mvc_data, mvc_dext] = fileparts(varargin{4}(:)');
            mvc_data = [mvc_data, mvc_dext];
        end
        if isempty(mvc_path)
            mvc_path = pwd;
        end
        ne_gcfg.h.MainFig.Pointer = 'watch';
        drawnow;
        mvc_data = load([mvc_path '/' mvc_data], '-mat');
        mvc_data = mvc_data.mvc_data;
        if ~isfield(mvc_data, 'btc') || ...
           ~isa(mvc_data.btc, 'single') || ...
            ndims(mvc_data.btc) ~= 5 || ...
           ~isfield(mvc_data, 'mtc') || ...
           ~isa(mvc_data.mtc, 'single') || ...
            ndims(mvc_data.mtc) ~= 5 || ...
           ~isfield(mvc_data, 'rtc') || ...
           ~iscell(mvc_data.rtc) || ...
            isempty(mvc_data.rtc) || ...
           ~isfield(mvc_data, 'avgopt') || ...
           ~isstruct(mvc_data.avgopt) || ...
            numel(mvc_data.avgopt) ~= 1 || ...
            ((numel(fieldnames(mvc_data.avgopt)) ~= 12 || ...
             ~all(strcmp(fieldnames(mvc_data.avgopt), ...
                  {'avgwin'; 'baseline'; 'basewin'; 'motpars'; 'group'; ...
                   'remnuis'; 'robust'; 'samptr'; 'sdse'; 'tfilter'; ...
                   'tfilttype'; 'xconfound'}))) && ...
             (numel(fieldnames(mvc_data.avgopt)) ~= 13 || ...
             ~all(strcmp(fieldnames(mvc_data.avgopt), ...
                  {'avgwin'; 'baseline'; 'basewin'; 'motpars'; 'group'; ...
                   'remnuis'; 'robust'; 'samptr'; 'sdse'; 'tfilter'; ...
                   'tfilttype'; 'xconfound'; 'avgbegin'})))) || ...
           ~isfield(mvc_data, 'colls') || ...
           ~iscell(mvc_data.colls) || ...
            size(mvc_data.colls, 2) ~= 3 || ...
           ~isfield(mvc_data, 'conds') || ...
           ~iscell(mvc_data.conds) || ...
            isempty(mvc_data.conds) || ...
            size(mvc_data.conds, 2) ~= 1 || ...
           ~isfield(mvc_data, 'condcol') || ...
           ~isa(mvc_data.condcol, 'double') || ...
           ~isequal(size(mvc_data.condcol), [size(mvc_data.conds, 1), 3]) || ...
            any(isinf(mvc_data.condcol(:)) | isnan(mvc_data.condcol(:))) || ...
           ~isfield(mvc_data, 'mdm') || ...
           ~isstruct(mvc_data.mdm) || ...
            numel(mvc_data.mdm) ~= 1 || ...
           ~isfield(mvc_data, 'onsets') || ...
           ~iscell(mvc_data.onsets) || ...
           ~isequal(size(mvc_data.onsets), [numel(mvc_data.rtc), numel(mvc_data.conds)]) || ...
           ~isfield(mvc_data, 'sdms') || ...
           ~iscell(mvc_data.sdms) || ...
            numel(mvc_data.sdms) ~= numel(mvc_data.rtc) || ...
           ~isfield(mvc_data, 'tcf') || ...
           ~iscell(mvc_data.tcf) || ...
            numel(mvc_data.tcf) ~= numel(mvc_data.rtc) || ...
           ~isfield(mvc_data, 'tr') || ...
           ~isa(mvc_data.tr, 'double') || ...
            numel(mvc_data.tr) ~= numel(mvc_data.rtc) || ...
           ~isfield(mvc_data, 'voi') || ...
           ~isstruct(mvc_data.voi) || ...
            numel(mvc_data.voi) ~= 1
            error('BAD_EXTRACT');
        end
        if ~isfield(mvc_data.avgopt, 'avgbegin')
            mvc_data.avgopt.avgbegin = 0;
        end
        mdm{1} = xff('new:mdm');
        setcont(mdm{1}, mvc_data.mdm);
        mdm{1}.SetHandle('FilesChecked', true);
        voi{1} = xff('new:voi');
        setcont(voi{1}, mvc_data.voi);
        if ~isxff(ne_gcfg.voi, 'voi') || ...
            isempty(ne_gcfg.voi.VOI)
            ne_loadcluster(0, 0, voi{1});
        end
    else
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{4}) || ...
            isempty(regexpi(varargin{4}(:)', '\.mdm$')) || ...
            exist(varargin{4}(:)', 'file') ~= 2
            mdm{1} = xff('*.mdm', 'Please select PRT+VTC-based MDM for averaging plots...');
        else
            mdm{1} = xff(varargin{4}(:)');
        end
        if ~isxff(mdm{1}, 'mdm') || ...
            isempty(mdm{1}.XTC_RTC) || ...
            size(mdm{1}.XTC_RTC, 2) ~= 2 || ...
            any(cellfun('isempty', regexpi(mdm{1}.XTC_RTC(:, 2), '\.prt$')))
            error('BAD_MDM');
        end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    if ~isempty(mdm{1})
        uiwait(warndlg('Invalid input file selected.', ...
            'NeuroElf - Error', 'modal'));
    end
    ne_gcfg.h.MainFig.Pointer = 'arrow';
    drawnow;
    clearxffobjects(mdm);
    if reload
        clearxffobjects(voi);
    end
    return;
end
rtv = mdm{1}.RunTimeVars;
xtc = mdm{1}.XTC_RTC;

% set pointer to watch for main figure
ne_gcfg.h.MainFig.Pointer = 'watch';
drawnow;

% remove nuisance variance
motpars = {};
remnuis = false;
tfilter = Inf;
tfilttype = 'dct';
xconfound = {};
if isfield(rtv, 'MotionParameters') && ...
    iscell(rtv.MotionParameters) && ...
    numel(rtv.MotionParameters) == size(xtc, 1)
    remnuis = true;
    motpars = rtv.MotionParameters;
end
if isfield(rtv, 'Confounds') && ...
    iscell(rtv.Confounds) && ...
    numel(rtv.Confounds) == size(xtc, 1)
    remnuis = true;
    xconfound = rtv.Confounds;
end
if isfield(rtv, 'TempFilterType') && ...
    ischar(rtv.TempFilterType) && ...
    ~isempty(rtv.TempFilterType) && ...
    isfield(rtv, 'TempFilterCutoff') && ...
    isa(rtv.TempFilterCutoff, 'double') && ...
    numel(rtv.TempFilterCutoff) == 1 && ...
    ~isinf(rtv.TempFilterCutoff) && ...
    ~isnan(rtv.TempFilterCutoff)
    remnuis = true;
    tfilter = rtv.TempFilterCutoff;
    tfilttype = rtv.TempFilterType;
end

% try to load the plotting figure
try
    hFig = neuroelf_file('f', 'ne_mdmcondavg');
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    mdm{1}.ClearObject;
    if reload
        clearxffobjects(voi);
    end
    uiwait(warndlg('Error loading MDM condition averaging figure.', ...
        'NeuroElf - Error', 'modal'));
    ne_gcfg.h.MainFig.Pointer = 'arrow';
    drawnow;
    return;
end
hTag = hFig.TagStruct;
tstr = hFig.Tag(1:8);

% get preliminary data
try
    if reload
        avgopt = mvc_data.avgopt;
        mtc = mvc_data.mtc;
        btc = mvc_data.btc;
        rtc = mvc_data.rtc;
        preds = mvc_data.conds;
        predcols = mvc_data.condcol;
        onsets = mvc_data.onsets;
        sdms = mvc_data.sdms;
        tcf = mvc_data.tcf;
        tr = mvc_data.tr;
    else
        cprog = ne_progress(0, 0, {true, 0, 'Reading condition list...'});
        mdm{1}.ConditionOnsets;
        ne_progress(0, 0, cprog);
        mdmh = handles(mdm{1});
        preds = mdmh.OOstlist;
        predcols = mdmh.OOstcols;
        predons = mdmh.OnOffsets;
        if nargin < 5 || ...
           ~iscell(varargin{5}) || ...
            isempty(varargin{5}) || ...
            any(cellfun('isempty', varargin{5}(:))) || ...
           ~all(cellfun(@ischar, varargin{5}(:)))
            [predi, ok] = listdlg( ...
                'PromptString', 'Please select the conditions to extract...', ...
                'ListString',   preds(:), ...
                'ListSize',     [min(300, max(600, 10 .* size(char(preds(:)), 2))), 200], ...
                'InitialValue', (1:numel(preds))', ...
                'Name',         'NeuroElf - condition selection');
        else
            predi = multimatch(varargin{5}(:), preds(:));
            if all(predi > 0) && ...
                numel(predi) == numel(unique(predi))
                ok = 1;
            else
                ok = 0;
            end
        end
        if ~isequal(ok, 1) || ...
            isempty(predi)
            mdm{1}.ClearObject;
            ne_gcfg.h.MainFig.Pointer = 'arrow';
            drawnow;
            return;
        end
        preds = preds(predi(:));
        predcols = predcols(predi, :);
        mdm{1}.SetHandle('OOstlist', preds(:));
        mdm{1}.SetHandle('OOstcols', predcols);
        mdm{1}.SetHandle('OnOffsets', predons(:, predi));
        avgopt = struct( ...
            'avgwin',    20000, ...
            'baseline',  'none', ...
            'basewin',   -2000:1000:0, ...
            'motpars',   {motpars}, ...
            'group',     'off', ...
            'pbar',      ne_gcfg.h.Progress, ...
            'remnuis',   remnuis, ...
            'robust',    false, ...
            'samptr',    100, ...
            'sdse',      'sd', ...
            'tfilter',   tfilter, ...
            'tfilttype', tfilttype, ...
            'xconfound', {xconfound}, ...
            'avgbegin',  0);
        cprog = ne_progress(0, 0, {true, 0, 'Extracting VOI time courses...'});
        [mtc, mtcse, btc, rtc, tcf, tr, onsets, sdms] = mdm{1}.VOICondAverage(voi{1}, preds, avgopt);
        ne_progress(0, 0, cprog);
    end
    btc = meannoinfnan(btc, 1);
    btc(btc == 0) = NaN;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    hFig.Delete;
    mdm{1}.ClearObject;
    if reload
        clearxffobjects(voi);
    end
    uiwait(warndlg(sprintf('Error extracting timecourses or reading PRTs (%s)', ...
        ne_eo.message), 'NeuroElf - Error', 'modal'));
    ne_gcfg.h.MainFig.Pointer = 'arrow';
    drawnow;
    return;
end

% tag shortcuts
hTag.AvgTR = hTag.(['ED_' tstr '_AVTR']);
hTag.AvgWin = hTag.(['ED_' tstr '_AWin']);
hTag.Axes = hTag.(['AX_' tstr '_Plot']);
hTag.BackColor = hTag.(['BT_' tstr '_BGrC']);
hTag.BaseWin = hTag.(['ED_' tstr '_BWin']);
hTag.CData = hTag.(['UIM_' tstr '_CData']);
hTag.CPlot = hTag.(['UIM_' tstr '_CPlot']);
hTag.CRaw = hTag.(['UIM_' tstr '_CRaw']);
hTag.CollAdd = hTag.(['BT_' tstr '_CPls']);
hTag.CollColors = hTag.(['BT_' tstr '_ColC']);
hTag.CollDelete = hTag.(['BT_' tstr '_CMns']);
hTag.Collapsed = hTag.(['LB_' tstr '_Coll']);
hTag.CondColors = hTag.(['BT_' tstr '_CndC']);
hTag.Conditions = hTag.(['LB_' tstr '_Cond']);
hTag.Diff = hTag.(['UIM_' tstr '_OptD']);
hTag.Diff1 = hTag.(['UIM_' tstr '_OptD1']);
hTag.Diff2 = hTag.(['UIM_' tstr '_OptD2']);
hTag.EBands = hTag.(['CB_' tstr '_EBnd']);
hTag.EBandAlpha = hTag.(['ED_' tstr '_EAlp']);
hTag.ECurves = hTag.(['CB_' tstr '_ECrv']);
hTag.Groups = hTag.(['LB_' tstr '_Grps']);
hTag.LLNE = hTag.(['UIM_' tstr '_OptLNE']);
hTag.LLNW = hTag.(['UIM_' tstr '_OptLNW']);
hTag.LLSE = hTag.(['UIM_' tstr '_OptLSE']);
hTag.LLSW = hTag.(['UIM_' tstr '_OptLSW']);
hTag.OptVis = hTag.(['CB_' tstr '_OVis']);
hTag.P05 = hTag.(['UIM_' tstr '_OptSP']);
hTag.Robust = hTag.(['CB_' tstr '_Rob']);
hTag.SAxes = hTag.(['UIM_' tstr '_SAxes']);
hTag.SData = hTag.(['UIM_' tstr '_SData']);
hTag.SDSE = hTag.(['CB_' tstr '_SDSE']);
hTag.SFigure = hTag.(['UIM_' tstr '_SFig']);
hTag.SubSel = hTag.(['UIM_' tstr '_OptSS']);
hTag.TSmooth = hTag.(['UIM_' tstr '_OptTS']);
hTag.Type = hTag.(['DD_' tstr '_Type']);
hTag.VOI = hTag.(['DD_' tstr '_VOI']);
hTag.VOIUpdate = hTag.(['CB_' tstr '_UVOI']);
hTag.XFrom = hTag.(['ED_' tstr '_AxX1']);
hTag.XTo = hTag.(['ED_' tstr '_AxX2']);
hTag.YFrom = hTag.(['ED_' tstr '_AxY1']);
hTag.YTo = hTag.(['ED_' tstr '_AxY2']);

% set update callbacks
hTag.AvgTR.Callback = {@ne_mdmvoicondavggui, tstr, 'AvgTR'};
hTag.AvgWin.Callback = {@ne_mdmvoicondavggui, tstr, 'AvgWin'};
hTag.BackColor.Callback = {@ne_mdmvoicondavggui, tstr, 'BackColor'};
hTag.BaseWin.Callback = {@ne_mdmvoicondavggui, tstr, 'BaseWin'};
hTag.CData.Callback = {@ne_mdmvoicondavggui, tstr, 'CopyData'};
hTag.CPlot.Callback = {@ne_mdmvoicondavggui, tstr, 'CopyPlot'};
hTag.CRaw.Callback = {@ne_mdmvoicondavggui, tstr, 'CopyRawData'};
hTag.CollAdd.Callback = {@ne_mdmvoicondavggui, tstr, 'CollAdd'};
hTag.CollColors.Callback = {@ne_mdmvoicondavggui, tstr, 'CollColors'};
hTag.CollDelete.Callback = {@ne_mdmvoicondavggui, tstr, 'CollDelete'};
hTag.Collapsed.Callback = {@ne_mdmvoicondavgup, tstr};
hTag.CondColors.Callback = {@ne_mdmvoicondavggui, tstr, 'CondColors'};
hTag.Conditions.Callback = {@ne_mdmvoicondavgup, tstr};
hTag.Diff1.Callback = {@ne_mdmvoicondavggui, tstr, 'Diff1'};
hTag.Diff2.Callback = {@ne_mdmvoicondavggui, tstr, 'Diff2'};
hTag.EBandAlpha.Callback = {@ne_mdmvoicondavggui, tstr, 'EBandAlpha'};
hTag.EBands.Callback = {@ne_mdmvoicondavggui, tstr, 'EBands'};
hTag.ECurves.Callback = {@ne_mdmvoicondavgup, tstr};
hTag.Groups.Callback = {@ne_mdmvoicondavggui, tstr, 'Groups'};
hTag.LLNE.Callback = {@ne_mdmvoicondavggui, tstr', 'LegLoc', 'NorthEast'};
hTag.LLNW.Callback = {@ne_mdmvoicondavggui, tstr', 'LegLoc', 'NorthWest'};
hTag.LLSE.Callback = {@ne_mdmvoicondavggui, tstr', 'LegLoc', 'SouthEast'};
hTag.LLSW.Callback = {@ne_mdmvoicondavggui, tstr', 'LegLoc', 'SouthWest'};
hTag.OptVis.Callback = {@ne_mdmvoicondavgrsz, tstr, 'OptVis'};
hTag.P05.Callback = {@ne_mdmvoicondavggui, tstr, 'p05'};
hTag.Robust.Callback = {@ne_mdmvoicondavggui, tstr, 'Robust'};
hTag.SAxes.Callback = {@ne_mdmvoicondavggui, tstr, 'SaveAxes'};
hTag.SDSE.Callback = {@ne_mdmvoicondavgup, tstr};
hTag.SData.Callback = {@ne_mdmvoicondavggui, tstr, 'SaveData'};
hTag.SFigure.Callback = {@ne_screenshot, hFig.MLHandle, '', 'high-q'};
hTag.SubSel.Callback = {@ne_mdmvoicondavggui, tstr, 'SubSel'};
hTag.TSmooth.Callback = {@ne_mdmvoicondavggui, tstr, 'TSmooth'};
hTag.Type.Callback = {@ne_mdmvoicondavgup, tstr};
hTag.VOI.Callback = {@ne_mdmvoicondavggui, tstr, 'VOI'};
hTag.XFrom.Callback = {@ne_mdmvoicondavggui, tstr, 'XRange'};
hTag.XTo.Callback = {@ne_mdmvoicondavggui, tstr, 'XRange'};
hTag.YFrom.Callback = {@ne_mdmvoicondavggui, tstr, 'YRange'};
hTag.YTo.Callback = {@ne_mdmvoicondavggui, tstr, 'YRange'};

% make sure MDM has minimal options
if ~isfield(mdm{1}.RunTimeVars, 'Groups') || ...
   ~iscell(mdm{1}.RunTimeVars.Groups) || ...
    size(mdm{1}.RunTimeVars.Groups, 2) ~= 2
    mdm{1}.RunTimeVars.Groups = cell(0, 2);
end
groups = mdm{1}.RunTimeVars.Groups;
for gc = 1:size(groups, 1)
    groups{gc, 2} = groups{gc, 2}(:);
end

% fill in figure options
hTag.Axes.Color = [1, 1, 1];
hTag.Axes.FontSize = cini.Satellites.FontSize;
hTag.VOI.String = strrep(voi{1}.VOINames, '_', ' ');
voiidx = ne_gcfg.h.Clusters.Value;
if numel(voiidx) == 1
    hTag.VOI.Value = voiidx;
else
    hTag.VOI.Value = 1;
end
hTag.Collapsed.String = {'none defined'};
hTag.Collapsed.Value = [];
hFig.SetGroupEnabled('CGrp', 'off');
hTag.Conditions.String = preds(:);
hTag.Conditions.Value = (1:numel(preds))';
if ~isempty(groups)
    hTag.Groups.String = groups(:, 1);
    hTag.Groups.Value = (1:size(groups, 1))';
else
    hTag.Groups.String = {'none defined'};
    hTag.Groups.Value = [];
    hFig.SetGroupEnabled('UGrp', 'off');
end
hTag.Robust.Value = double(avgopt.robust);
if mdm{1}.PSCTransformation ~= 0
    hTag.Trans.Value = 2;
elseif mdm{1}.zTransformation ~= 0
    hTag.Trans.Value = 3;
end

% create options
avgopts = struct( ...
    'ax',        hTag.Axes.MLHandle, ...
    'avgopt',    avgopt, ...
    'axcol',     [255, 255, 255], ...
    'axpos',     hTag.Axes.Position, ...
    'bands',     [], ...
    'btc',       btc, ...
    'colls',     {cell(0, 3)}, ...
    'condcol',   max(0, min(255, predcols)), ...
    'conds',     {preds}, ...
    'curves',    [], ...
    'diffmode',  0, ...
    'fig',       hFig, ...
    'groups',    hTag.Groups.Value, ...
    'legloc',    'NorthWest', ...
    'legends',   [], ...
    'lines',     [], ...
    'mdm',       mdm{1}, ...
    'mtc',       mtc, ...
    'p05',       false, ...
    'rtc',       {rtc}, ...
    'onsets',    {onsets}, ...
    'optvis',    true, ...
    'optvispos', hTag.OptVis.Position(1:2), ...
    'robust',    avgopt.robust, ...
    'sattype',   'mdmcondavg', ...
    'sdms',      {sdms}, ...
    'sealpha',   0.25, ...
    'subsel',    (1:numel(mdm{1}.Subjects))', ...
    'tcf',       {tcf}, ...
    'texts',     [], ...
    'tsmooth',   0, ...
    'tr',        tr, ...
    'typeuic',   hTag.Type, ...
    'voi',       voi{1}.CopyObject);
if reload
    voi{1}.ClearObject;
    if isfield(mvc_data, 'subsel') && ...
        isa(mvc_data.subsel, 'double') && ...
       ~any(isinf(mvc_data.subsel(:)) | isnan(mvc_data.subsel(:)) | ...
            mvc_data.subsel(:) < 1 | mvc_data.subsel(:) > avgopts.subsel(end)) && ...
        all(mvc_data.subsel(:) == fix(mvc_data.subsel(:)))
        avgopts.subsel = unique(mvc_data.subsel(:));
    end
    avgopts.colls = mvc_data.colls;
    if ~isempty(avgopts.colls)
        hTag.Collapsed.String = avgopts.colls(:, 1);
        hFig.SetGroupEnabled('CGrp', 'on');
    end
end
voih = handles(avgopts.voi);
if isfield(voih, 'RGBImage')
    avgopts.voi.DeleteHandle('RGBImage');
end

% set in children array
ne_gcfg.cc.(tstr) = struct( ...
    'Config',       avgopts, ...
    'Satellite',    hFig, ...
    'SatelliteMLH', hFig.MLHandle, ...
    'Tags',         hTag);

% give the figure a more appropriate name
[mdmp, mdmf] = fileparts(mdm{1}.FilenameOnDisk);
hFig.Name = sprintf('NeuroElf - MDM VOI-based condition average plot - %s', mdmf);

% and set a resize function
hFig.ResizeFcn = {@ne_mdmvoicondavgrsz, tstr};

% robust?
if cini.Statistics.RobustVOICondAvg
    hTag.Robust.Value = 1;
    ne_mdmvoicondavggui(0, 0, tstr, 'Robust');
end

% initialize plot
ne_mdmvoicondavgup(0, 0, tstr);

% set position and make visible
if any(cini.Children.MDMCondAvgPosition > 0)
    lastpos = hFig.Position;
    lastpos(1:2) = cini.Children.MDMCondAvgPosition;
    hFig.Position = lastpos;
end
ne_gcfg.h.MainFig.Pointer = 'arrow';
hFig.Visible = 'on';
drawnow;

% return values
if nargout > 0
    varargout{1} = hTag.Axes.MLHandle;
    if nargout > 1
        varargout{2} = hFig.MLHandle;
        if nargout > 2
            varargout{3} = tstr;
            if nargout > 3
                if ~reload
                    mvc_data = struct;
                    mvc_data.btc = btc;
                    mvc_data.mtc = mtc;
                    mvc_data.rtc = rtc;
                    mvc_data.avgopt = rmfield(avgopt, 'pbar');
                    mvc_data.colls = cell(0, 3);
                    mvc_data.condcol = avgopts.condcol;
                    mvc_data.conds = preds;
                    mvc_data.mdm = mdm{1};
                    mvc_data.onsets = onsets;
                    mvc_data.sdms = sdms;
                    mvc_data.subsel = avgopts.subsel;
                    mvc_data.tcf = tcf;
                    mvc_data.tr = tr;
                    mvc_data.voi = getcont(avgopts.voi);
                end
                varargout{4} = mvc_data;
            end
        end
    end
end

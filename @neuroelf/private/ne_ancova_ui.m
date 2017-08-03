% FUNCTION ne_ancova_ui: UI functions for ANCOVA dialog/UI
function ne_ancova_ui(varargin)

% Version:  v1.1
% Build:    17061220
% Date:     Jun-12 2017, 8:58 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, 2017, Jochen Weber
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

% handles
cc = ne_gcfg.fcfg.AC;
ch = ne_gcfg.h.AC.h;
hFig = ne_gcfg.h.AC.ACFig;

% argument?
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3})
    return;
end

% what to do
switch (lower(varargin{3}(:)'))

    % cell content
    case {'cellcont'}

    % cell type
    case {'celltype'}

        % required text
        if ch.CellType.Value == 1
            cnames = [{'<predictors>'}; cc.preds(:)];
        else
            if ch.Layout1.Value > 10
                ch.Layout1.Value = 2;
                ne_ancova_ui(0, 0, 'layout');
            end
            cnames = [{'<contrasts>'}; cc.cons(:, 1)];
        end
        for clc = 1:numel(ch.CellCont)
            ch.CellCont{clc}.Value = 1;
            ch.CellCont{clc}.String = cnames;
        end
        
    % close UI
    case {'closeui'}

        % save RunTimeVars on last GLM
        if numel(hFig.UserData.lastglm) == 1 && ...
            isxff(hFig.UserData.lastglm, 'glm')
            ne_ancova_ui(0, 0, 'updatertv', hFig.UserData.lastglm, true);
        end

        % update last known position
        ne_gcfg.c.ini.Children.ACPosition = hFig.Position(1:2);

        % delete figure and remove from global struct
        hFig.Delete;
        ne_gcfg.fcfg.AC = [];
        ne_gcfg.h.AC = [];

        % release contrast manager call
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'ancova_open')) = [];

    % layout
    case {'layout'}

        % parse enable/disable information
        fl1 = ch.Layout1.Value;
        fl2 = ch.LayoutX.Value;
        if fl2 == 4
            fl2 = [2, 2];
        elseif fl2 == 5
            fl2 = [3, 2];
        else
            fl2(2) = 1;
        end
        predonly = false;
        if fl1 > 10
            fls = [fl2, 1];
            predonly = true;
        else
            fls = [fl1, fl2];
        end
        for c1 = 1:10
            for c2 = 1:3
                for c3 = 1:2
                    if any([c1, c2, c3] > fls)
                        ch.CellCont{c1, c2, c3}.Enable = 'off';
                        ch.CellCont{c1, c2, c3}.Value = 1;
                    else
                        if strcmpi(ch.CellCont{c1, c2, c3}.Enable, 'off')
                            ch.CellCont{c1, c2, c3}.Value = 1;
                        end
                        ch.CellCont{c1, c2, c3}.Enable = 'on';
                    end
                end
            end
        end
        if fls(3) < 2
            hFig.Position(3) = 776;
            ch.HDivider.Position(3) = 524;
            ch.Compute.Position(1) = 656;
        else
            hFig.Position(3) = 1288;
            ch.HDivider.Position(3) = 1036;
            ch.Compute.Position(1) = 1168;
        end
        if predonly && ...
            ch.CellType.Value == 2
            ch.CellType.Value = 1;
            ne_ancova_ui(0, 0, 'celltype');
        end
        if ch.Layout1.Value == 11
            uiwait(warndlg('Please select the last lag (_Dx) you wish to include', ...
                'NeuroElf - info', 'modal'));
        end

    % load config
    case {'loadcfg'}

        % request load filename
        [loadname, loadpath] = uigetfile( ...
            {'*.ini', 'ANCOVA configuration files (*.ini)'}, ...
            'Please select the ANCOVA configuration to load (*.ini)');
        if isequal(loadname, 0) || ...
            isequal(loadpath, 0) || ...
            isempty(loadname)
            return;
        end
        cfgfile = [loadpath filesep loadname];
        try
            cfg = [];
            cfg = xini(cfgfile, 'convert');
            glmfile = cfg.GLM.Filename;
            glmxffid = cfg.GLM.xffID;
            groups = cfg.Groups;
            subjects = cfg.Subjects.Selection(:);
            if isfield(cfg.Covariates, 'Interactions')
                coviact = cfg.Covariates.Interactions;
            else
                coviact = false;
            end
            if isfield(cfg.Covariates, 'RemoveMean') && ...
                ischar(cfg.Covariates.RemoveMean) && ...
                any(strcmpi(cfg.Covariates.RemoveMean, {'none', 'group'}))
                if lower(cfg.Covariates.RemoveMean(1)) == 'n'
                    mrcov = 1;
                else
                    mrcov = 3;
                end
            else
                mrcov = 2;
            end
            covs = cfg.Covariates.Selection(:);
            factors = cfg.Factors;
            cells = cfg.Cell;
            globcfg = cfg.Global;
            cfg.Release;
            cfg = [];

            % make settings
            if ~isempty(glmfile) && ...
               ~strcmpi(glmfile, hFig.UserData.lastglm.FilenameOnDisk) && ...
               ~strcmpi(glmxffid, hFig.UserData.lastglm.RunTimeVars.xffID)
                vans = questdlg( ...
                    'Configuration was created for different GLM. Use it anyway?', ...
                    'NeuroElf - request', 'Yes', 'No', 'No');
                if ~ischar(vans) || ...
                    isempty(vans) || ...
                    lower(vans(1)) ~= 'y'
                    return;
                end
            end
            glm = hFig.UserData.lastglm;
            ch.UseGroups.Value = 0;
            if groups.UseGroups && ...
               ~isempty(glm.RunTimeVars.Groups)
                ch.UseGroups.Value = 1;
                ch.Groups.Value = groups.Groups(:);
            end
            if ~isempty(glm.RunTimeVars.Groups)
                ne_ancova_ui(0, 0, 'usegroups');
            end
            gsubjects = glm.Subjects;
            submatch = multimatch(subjects, gsubjects);
            if any(submatch == 0)
                misssubs = subjects(submatch == 0);
                misssubs = sprintf('%s, ', misssubs{:});
                uiwait(warndlg(sprintf('Missing subjects from configuration: %s.', ...
                    misssubs(1:end-2)), 'NeuroElf - warning', 'modal'));
                submatch = submatch(submatch > 0);
            end
            ch.Subjects.Value = unique(submatch(:));
            gcovs = glm.RunTimeVars.CovariatesNames(:);
            covmatch = multimatch(covs, gcovs);
            if any(covmatch == 0)
                misscovs = covs(covmatch == 0);
                misscovs = sprintf('%s, ', misscovs{:});
                uiwait(warndlg(sprintf('Missing covariates from configuration: %s.', ...
                    misscovs(1:end-2)), 'NeuroElf - warning', 'modal'));
                covmatch = covmatch(covmatch > 0);
            end
            ch.Covs.Value = unique(covmatch(:));
            ch.CovInteractions.Value = double(coviact);
            hFig.RadioGroupSetOne('MRCov', mrcov);
            f1s = [];
            if ~ischar(factors.Layout1)
                f1s = factors.Layout1;
                ch.Layout1.Value = f1s;
            elseif lower(factors.Layout1(1)) == 'f'
                ch.Layout1.Value = 11;
            else
                ch.Layout1.Value = 12;
            end
            switch sum([10, 1] .* factors.LayoutX)
                case {11}
                    ch.LayoutX.Value = 1;
                case {21}
                    ch.LayoutX.Value = 2;
                case {31}
                    ch.LayoutX.Value = 3;
                case {22}
                    ch.LayoutX.Value = 4;
                case {32}
                    ch.LayoutX.Value = 5;
            end
            fsel = [f1s, factors.LayoutX];
            if numel(fsel) < 3
                fsel(3) = 1;
            end
            ne_ancova_ui(0, 0, 'layout');
            if lower(cells.Type(1)) == 'p'
                ch.CellType.Value = 1;
            else
                ch.CellType.Value = 2;
            end
            ne_ancova_ui(0, 0, 'celltype');
            gcellconts = ch.CellCont{1}.String;
            if ~iscell(gcellconts)
                gcellconts = cellstr(gcellconts);
            end
            cvals = zeros(fsel);
            for c3 = 1:fsel(3)
                for c2 = 1:fsel(2)
                    for c1 = 1:fsel(1)
                        cvals(c1, c2, c3) = multimatch({ ...
                            cells.(sprintf('Cont%02d%d%d', c1, c2, c3))}, ...
                            gcellconts);
                    end
                end
            end
            if any(cvals(:) < 1)
                uiwait(warndlg('Invalid predictor/contrast specification.', ...
                    'NeuroElf - error', 'modal'));
                return;
            end
            if any(cvals(:) < 2) || ...
                numel(unique(cvals(:))) ~= numel(cvals)
                uiwait(warndlg('Non-unique predictor/contrast specification.', ...
                    'NeuroElf - info', 'modal'));
            end
            for c3 = 1:fsel(3)
                for c2 = 1:fsel(2)
                    for c1 = 1:fsel(1)
                        ch.CellCont{c1, c2, c3}.Value = cvals(c1, c2, c3);
                    end
                end
            end
            ch.AddGlobalMean.Value = double(globcfg.AddMean);
            ch.SmoothData.Value = double(globcfg.SmoothMaps);
            ch.SmoothDataKernel.String = sprintf('%.1f', globcfg.SmoothingKernel);
            ch.BRange.Value = double(globcfg.LimitBetas);
            ch.BRangeMin.String = sprintf('%g', globcfg.LimitRange(1));
            ch.BRangeMax.String = sprintf('%g', globcfg.LimitRange(1));
            ch.InterpMethod.Value = 1;
            switch lower(globcfg.Interpolation(1))
                case {'c'}
                    ch.InterpMethod.Value = 2;
                case {'s'}
                    ch.InterpMethod.Value = 3;
            end
        catch ne_eo;
            if ~isempty(cfg)
                cfg.Release;
            end
            uiwait(warndlg(['Error loading configuration: ' ne_eo.message], ...
                'NeuroElf - error', 'modal'));
            return;
        end

    % save config
    case {'savecfg'}

        % request save filename
        [savename, savepath] = uiputfile( ...
            {'*.ini', 'ANCOVA configuration files (*.ini)'}, ...
            'Save ANCOVA configuration as (*.ini)');
        if isequal(savename, 0) || ...
            isequal(savepath, 0) || ...
            isempty(savename)
            return;
        end
        if isempty(savepath)
            savepath = pwd;
        end
        cfgfile = [savepath filesep savename];

        % load template
        cfg = xini([neuroelf_path filesep '_core' filesep 'config' filesep 'ancova.ini'], 'convert');
        
        % try saving as
        try
            cfg.SaveAs(cfgfile);
        catch ne_eo;
            cfg.Release;
            uiwait(warndlg(['Error saving configuration: ' ne_eo.message], ...
                'NeuroElf - error', 'modal'));
            return;
        end

        % make settings
        glm = hFig.UserData.lastglm;
        cfg.GLM.Filename = glm.FilenameOnDisk;
        cfg.GLM.xffID = glm.RunTimeVars.xffID;
        cfg.Groups.UseGroups = (ch.UseGroups.Value > 0);
        cfg.Groups.Groups = ch.Groups.Value(:)';
        subsel = ch.Subjects.Value;
        if ~isempty(subsel)
            cfg.Subjects.Selection = lsqueeze(ch.Subjects.String(subsel));
        else
            cfg.Subjects.Selection = {};
        end
        covsel = ch.Covs.Value;
        if ~isempty(covsel)
            cfg.Covariates.Selection = lsqueeze(ch.Covs.String(covsel));
        else
            cfg.Covariates.Selection = {};
        end
        cfg.Covariates.Interactions = (ch.CovInteractions.Value > 0);
        if ch.MRCovNone.Value > 0
            cfg.Covariates.RemoveMean = 'none';
        elseif ch.MRCovGroup.Value > 0
            cfg.Covariates.RemoveMean = 'group';
        else
            cfg.Covariates.RemoveMean = 'global';
        end
        f1sel = ch.Layout1.Value;
        if f1sel <= 10
            f1seln = f1sel;
        elseif f1sel == 10
            f1sel = [];
            f1seln = 'FIR';
        else
            f1sel = [];
            f1seln = 'ST';
        end
        cfg.Factors.Layout1 = f1seln;
        f2sel = ch.LayoutX.Value;
        switch f2sel
            case {1}
                f2sel = [1, 1];
            case {2}
                f2sel = [2, 1];
            case {3}
                f2sel = [3, 1];
            case {4}
                f2sel = [2, 2];
            case {5}
                f2sel = [3, 2];
        end
        cfg.Factors.LayoutX = f2sel;
        if ch.CellType.Value == 1
            cfg.Cell.Type = 'predictors';
        else
            cfg.Cell.Type = 'contrasts';
        end
        fsel = [f1sel, f2sel];
        if numel(fsel) < 3
            fsel(3) = 1;
        end
        ccont = ch.CellCont{1}.String;
        if ~iscell(ccont)
            ccont = cellstr(ccont);
        end
        for c3 = 1:fsel(3)
            for c2 = 1:fsel(2)
                for c1 = 1:fsel(1)
                    cval = ch.CellCont{c1, c2, c3}.Value;
                    if cval == 1
                        cfg.Release;
                        delete(cfgfile);
                        uiwait(warndlg('Not all cells filled in.', ...
                            'NeuroElf - info', 'modal'));
                        return;
                    end
                    cfg.Cell.(sprintf('Cont%02d%d%d', c1, c2, c3)) = ccont{cval};
                end
            end
        end
        cfg.Global.AddMean = (ch.AddGlobalMean.Value > 0);
        interps = {'linear', 'cubic', 'sinc3'};
        cfg.Global.Interpolation = interps{ch.InterpMethod.Value};
        cfg.Global.LimitBetas = (ch.BRange.Value > 0);
        brangemin = str2double(ch.BRangeMin.String);
        brangemax = str2double(ch.BRangeMin.String);
        if numel(brangemin) ~= 1 || ...
            isnan(brangemin) || ...
            brangemin > 0
            brangemin = -Inf;
        end
        if numel(brangemax) ~= 1 || ...
            isnan(brangemax) || ...
            brangemax < 0
            brangemax = Inf;
        end
        cfg.Global.LimitRange = [brangemin, brangemax];
        cfg.Global.SmoothMaps = (ch.SmoothData.Value > 0);
        smkern = str2double(ch.SmoothDataKernel.String);
        if numel(smkern) ~= 1 || ...
            isnan(smkern) || ...
            isinf(smkern) || ...
            smkern < 0
            smkern = 0;
        end
        cfg.Global.SmoothingKernel = 0.01 * round(100 * smkern);
        
        % save over
        cfg.Save;
        cfg.Release;
        
    % set GLM
    case {'setglm'}

        % update RunTimeVars if required
        if numel(hFig.UserData.lastglm) == 1 && ...
            isxff(hFig.UserData.lastglm, 'glm')
            ne_ancova_ui(0, 0, 'updatertv', hFig.UserData.lastglm, true);
        end

        % reset content of entire UI!
        ne_gcfg.fcfg.AC = struct( ...
            'cons',      {cell(0, 2)}, ...
            'covs',      {cell(0, 2)}, ...
            'glm',       [], ...
            'GLMs',      {ne_gcfg.fcfg.AC.GLMs}, ...
            'groups',    {cell(0, 2)}, ...
            'preds',     {cell(0, 1)}, ...
            'subsel',    {cell(0, 1)}, ...
            'usegroups', false, ...
            'vmp',       []);
        ch.Layout1.Value = 2;
        ch.LayoutX.Value = 1;
        ch.CellType.Value = 1;
        for clc = 1:numel(ch.CellCont)
            ch.CellCont{clc}.Value = 1;
        end
        ch.StoreInVMP.Value = 1;

        % error-handling
        try

            % get currently selected GLM
            glm = cc.GLMs{ch.GLMs.Value};

            % disable groups by default
            cc.usegroups = false;
            ch.UseGroups.Value = 0;
            ch.Groups.String = {'<no groups specified>'};
            ch.Groups.Value = 1;
            hFig.SetGroupEnabled('Groups', 'off');

            % get list of subjects
            subjects = glm.Subjects;

            % make sure there is a valid Contrasts field
            if ~isfield(glm.RunTimeVars, 'Contrasts') || ...
               ~iscell(glm.RunTimeVars.Contrasts) || ...
                size(glm.RunTimeVars.Contrasts, 2) ~= 2
                glm.RunTimeVars.Contrasts = cell(0, 2);
            end

            % make sure there is a valid CovariatesData and CovariatesNames field
            if ~isfield(glm.RunTimeVars, 'CovariatesData') || ...
               ~isfield(glm.RunTimeVars, 'CovariatesNames') || ...
               ~isa(glm.RunTimeVars.CovariatesData, 'double') || ...
               ~iscell(glm.RunTimeVars.CovariatesNames) || ...
                size(glm.RunTimeVars.CovariatesData, 1) ~= numel(subjects) || ...
                size(glm.RunTimeVars.CovariatesData, 2) ~= numel(glm.RunTimeVars.CovariatesNames)
                glm.RunTimeVars.CovariatesData = zeros(numel(subjects), 0);
                glm.RunTimeVars.CovariatesNames = cell(1, 0);
            end

            % make sure there is a valid Groups field
            if ~isfield(glm.RunTimeVars, 'Groups') || ...
               ~iscell(glm.RunTimeVars.Groups) || ...
                size(glm.RunTimeVars.Groups, 2) ~= 2
                glm.RunTimeVars.Groups = cell(0, 2);
            end

            % make sure there is a valid SubSel field
            if ~isfield(glm.RunTimeVars, 'SubSel') || ...
               ~iscell(glm.RunTimeVars.SubSel) || ...
                size(glm.RunTimeVars.SubSel, 2) ~= 1
                glm.RunTimeVars.SubSel = glm.Subjects;
            end

            % make sure there is a valid UseGroups field
            if ~isfield(glm.RunTimeVars, 'UseGroups') || ...
               ~islogical(glm.RunTimeVars.UseGroups) || ...
                numel(glm.RunTimeVars.UseGroups) ~= 1
                glm.RunTimeVars.UseGroups = false;
            end

            % set subjects
            cc.glm = glm;
            subjects = glm.Subjects;
            cc.nsubs = numel(subjects);
            ch.Subjects.ListboxTop = 1;
            ch.Subjects.Value = 1;
            ch.Subjects.String = subjects(:);
            ch.Subjects.Value = 1:numel(subjects);

            % and get SubjectPredictors (without the mean predictor)
            preds = glm.SubjectPredictors;
            if glm.NrOfSubjectPredictors == numel(preds)
                preds(end) = [];
            end

            % set in config
            cc.cons = glm.RunTimeVars.Contrasts;
            cc.covs = glm.RunTimeVars.CovariatesNames(:);
            cc.groups = glm.RunTimeVars.Groups;
            for cvc = 1:size(cc.covs, 1)
                cc.covs{cvc, 2} = glm.RunTimeVars.CovariatesData(:, cvc);
            end
            cc.preds = preds(:);
            cc.subsel = glm.RunTimeVars.SubSel;
            cc.usegroups = glm.RunTimeVars.UseGroups;
            ne_gcfg.fcfg.AC = cc;
            ne_ancova_ui(0, 0, 'celltype');
            ne_ancova_ui(0, 0, 'layout');

            % re-set groups controls
            if ~isempty(cc.groups) && ...
                cc.usegroups
                ch.UseGroups.Value = 1;
                ch.Groups.String = cc.groups(:, 1);
                ch.Groups.Value = (1:size(cc.groups, 1))';
                groupsubs = lsqueezecells(cc.groups(:, 2));
                ch.Subjects.Value = unique(lsqueeze(cat(1, groupsubs{:})));
                hFig.SetGroupEnabled('Groups', 'on');
                hFig.SetGroupEnabled('NGroups', 'off');
            else
                ch.Groups.Value = [];
                if isempty(cc.groups)
                    hFig.SetGroupEnabled('Groups', 'off');
                end
                hFig.SetGroupEnabled('NGroups', 'on');
                ch.Subjects.Value = lsqueeze(find(multimatch(glm.Subjects, cc.subsel) > 0));
            end
            ch.Subjects.ListboxTop = 1;

            % re-set covariates control
            ch.Covs.ListboxTop = 1;
            ch.Covs.Value = [];
            ch.Covs.String = cc.covs(:, 1);
            ch.CovInteractions.Value = 1;
            hFig.RadioGroupSetOne('MRCov', 2);

            % store in (default)
            if glm.ProjectType < 2
                ch.StoreInVMP.String = {'<new.vmp>'};
            else
                ch.StoreInVMP.String = {'<new.smp>'};
            end

            % disable beta-range setting (default!)
            ch.BRange.Value = 0;

            % update
            drawnow;
            ch.Subjects.ListboxTop = 1;
            drawnow;

            % set reference to last GLM
            hFig.UserData.lastglm = cc.glm;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            ne_ancova_ui(0, 0, 'closeui');
            uiwait(warndlg('Error using the selected GLM; closing UI.', 'NeuroElf GUI - error', 'modal'));
            return;
        end

        % which VMPs are available to store data
        glay = glm.Layout;
        grtv = glm.RunTimeVars;
        gtrf = grtv.TrfPlus;
        broot = xff;
        if glm.ProjectType < 2
            tvmp = broot.Documents('vmp');
        else
            tvmp = broot.Documents('smp');
        end
        if ~isfield(grtv, 'SubjectSPMsn') || ...
           ~isstruct(grtv.SubjectSPMsn) || ...
            numel(grtv.SubjectSPMsn) ~= 1 || ...
            isempty(fieldnames(grtv.SubjectSPMsn))
            for vmpc = numel(tvmp):-1:1
                vhnd = broot.Document(tvmp{vmpc});
                vlay = vhnd.Layout;
                vhnr = vhnd.RunTimeVars;
                if any(glay([1:3, 5:10]) ~= vlay([1:3, 5:10])) || ...
                   (isfield(vhnr, 'TrfPlus') && ~isequal(vhnr.TrfPlus, gtrf))
                    tvmp(vmpc) = [];
                end
            end
        end

        % any good
        if ~isempty(tvmp)

            % compile name of VMPs
            vmpn = cell(numel(tvmp) + 1, 1);
            if glm.ProjectType < 2
                mext = 'vmp';
            else
                mext = 'smp';
            end
            vmpn{1} = ['<new.' mext '>'];
            for vmpc = numel(tvmp):-1:1
                if ischar(tvmp{vmpc}) && ~isempty(regexpi(tvmp{vmpc}, ['\.' mext '$']))
                    [vmppath, vmpfile, vmpfext] = fileparts(tvmp{vmpc});
                    vmpn{vmpc + 1} = [vmpfile, vmpfext];
                else
                    vmpn{vmpc + 1} = sprintf('<xffID=%s>.%s', tvmp{vmpc}, mext);
                end
                tvmp{vmpc + 1} = broot.Document(tvmp{vmpc});
            end
            tvmp{1} = [];
            ch.StoreInVMP.String = vmpn;
            ch.StoreInVMP.Value = 2;
        end
        cc.vmp = tvmp;
        ne_gcfg.fcfg.AC = cc;

    % smoothing settings
    case {'smoothdata'}

        % check flag
        if ch.SmoothData.Value > 0
            hFig.SetGroupEnabled('SMData', 'on');
            try
                smk = str2double(ch.SmoothDataKernel.String);
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                smk = 6;
            end
            if numel(smk) ~= 1 || ...
                isnan(smk) || ...
                isinf(smk) || ...
                smk < 0
                smk = 6;
            else
                smk = 0.1 * round(10 * smk);
            end
            ch.SmoothDataKernel.String = sprintf('%.1f', smk);
        else
            hFig.SetGroupEnabled('SMData', 'off');
        end

    % use-groups
    case {'usegroups'}
        
        % enable?
        if nargin > 3 && ...
            numel(varargin{4}) == 1 && ...
            islogical(varargin{4}) && ...
            varargin{4} && ...
           ~isempty(hFig.UserData.lastglm.RunTimeVars.Groups)
            ch.UseGroups.Value = 1;
        end

        % if groups enabled
        if ch.UseGroups.Value > 0
            hFig.SetGroupEnabled('NGroups', 'off');
        else
            hFig.SetGroupEnabled('NGroups', 'on');
        end
end

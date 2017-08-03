function ne_neurosynth(varargin)
% ne_neurosynth  - access (older) neurosynth database
%
% FORMAT:       ne_neurosynth(SRC, EVT, action, [opts{:}])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       action      action, supported: {'load'}
%       opts        options for action
%
%    loading a map
%       opts{1}     term name
%       opts{2}     map type, one of 'fi', 'pp', 'ri (must be given!)
%
% Example:
%
%   ne_neurosynth(0, 0, 'load', 'fear', 'fi');


% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:39 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, 2013, 2014, 2016, Jochen Weber
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
ch = ne_gcfg.h;

% arguments
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3})
    return;
end
action = lower(varargin{3}(:)');

% neurosynth path
npath = neuroelf_path('nsynth');

% what to do
switch (action)

    % load map
    case {'load'}

        % request word if not given
        if nargin > 4 && ...
            ischar(varargin{4}) && ...
           ~isempty(varargin{4}) && ...
            ischar(varargin{5}) && ...
            any(strcmpi(varargin{5}(:)', {'fi', 'pp', 'ri'}))
            term = lower(varargin{4}(:)');
            type = lower(varargin{5}(:)');
        else
            term = inputdlg({'Please enter a term to download a map for...', ...
                'Forward inference (fi), reverse inferece (ri), or posterior prob (pp)'}, ...
                'NeuroElf - user input', 1, {'search term', 'fi'});
            if ~iscell(term) || ...
                numel(term) ~= 2 || ...
               ~ischar(term{1}) || ...
                isempty(term{1}) || ...
               ~ischar(term{2}) || ...
                isempty(term{2}) || ...
                any(term{1} == ' ') || ...
               ~any('fpr' == lower(term{2}(1)))
                return;
            end
            switch (lower(term{2}(1)))
                case {'f'}
                    type = 'fi';
                case {'p'}
                    type = 'pp';
                case {'r'}
                    type = 'ri';
            end
            term = lower(term{1}(:)');
        end

        % term must be available for download
        if ~any(strcmp(ne_gcfg.fcfg.nsynth.terms, term))
            uiwait(errordlg(['Term ' term ' not in the database.'], ...
                'NeuroElf - error', 'modal'));
            return;
        end

        % download term map
        tdate = datestr(now, 'mm.dd.yy');
        tname = sprintf('%s/terms/%s_%s_%s_term.nii.gz', ...
            npath, type, term, tdate);
        try
            ne_httpget(0, 0, sprintf( ...
                'http://old.neurosynth.org/downloads/terms/%s/%s', ...
                term, type), tname);
        catch ne_eo;
            uiwait(errordlg(ne_eo.message, 'NeuroElf - error', 'modal'));
            return;
        end

        % load map
        try
            ne_neurosynth(0, 0, 'reload', sprintf('%s_%s_%s_term.nii.gz', ...
                type, term, tdate));
            ne_neurosynth(0, 0, 'menu');
        catch ne_eo;
            uiwait(errordlg(ne_eo.message, 'NeuroElf - error', 'modal'));
        end

    % update menu
    case {'menu'}

        % update list of maps
        termm = findfiles([npath '/terms'], '*.nii.gz', 'depth=1', 'relative=');
        ne_gcfg.fcfg.nsynth.termm = termm;

        % get main menu item
        mmh = ch.NeuroSynth.MLHandle;

        % delete children (if any)
        mmc = get(mmh, 'Children');
        if ~isempty(mmc)
            try
                for ccc = 1:numel(mmc)
                    mmcc = get(mmc(ccc), 'Children');
                    try
                        if ~isempty(mmcc)
                            delete(mmcc);
                        end
                    catch ne_eo;
                        ne_gcfg.c.lasterr = ne_eo;
                    end
                end
                delete(mmc);
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end

        % add "load term" menu entry
        uimenu(mmh, 'Label', 'Load term map from NeuroSynth...',  ...
            'Callback', {@ne_neurosynth, 'load'});

        % if any maps
        if ~isempty(termm)
            if numel(termm) <= 20
                for tc = 1:numel(termm)
                    t = splittocellc(termm{tc}, '_');
                    switch (lower(t{1}(1)))
                        case {'f'}
                            t{1} = 'Forward-inference:';
                        case {'p'}
                            t{1} = 'Posterior-probability:';
                        case {'r'}
                            t{1} = 'Reverse-inference:';
                    end
                    if tc < 2 || ...
                        termm{tc}(1) ~= termm{tc-1}(1)
                        sep = 'on';
                    else
                        sep = 'off';
                    end
                    uimenu(mmh, 'Label', ['    ' t{1} ' ' t{2} ' (' t{3} ')'], ...
                        'Callback', {@ne_neurosynth, 'reload', termm{tc}}, ...
                        'Separator', sep);
                end
            elseif numel(termm)
                % do inference types separately
                for ttype = {'fi', 'ri', 'pp'; ...
                        'Forward inference', 'Reverse inferece', 'Posterior probability'}
                    termm = findfiles([npath '/terms'], ...
                        [ttype{1} '*.nii.gz'], 'depth=1', 'relative=');
                    if isempty(termm)
                        continue;
                    end
                    if numel(termm) <= 20
                        mmh1 = uimenu(mmh, 'Label', [ttype{2} ' maps'], ...
                            'Separator', 'on');
                        for tc = 1:numel(termm)
                            t = splittocellc(termm{tc}, '_');
                            uimenu(mmh1, 'Label', ['    ' t{2} ' (' t{3} ')'], ...
                                'Callback', {@ne_neurosynth, 'reload', termm{tc}}, ...
                                'Separator', 'off');
                        end
                    else
                        spp = 'on';
                        while ~isempty(termm)
                            ntc = min(20, numel(termm));
                            t1 = splittocellc(termm{1}, '_');
                            t2 = splittocellc(termm{ntc}, '_');
                            mmh1 = uimenu(mmh, 'Separator', spp, 'Label', ...
                                sprintf('%s maps (%s through %s)', ttype{2}, t1{2}, t2{2}));
                            spp = 'off';
                            for tc = 1:ntc
                                t = splittocellc(termm{tc}, '_');
                                uimenu(mmh1, 'Label', ['    ' t{2} ' (' t{3} ')'], ...
                                    'Callback', {@ne_neurosynth, 'reload', termm{tc}}, ...
                                    'Separator', 'off');
                            end
                            termm(1:ntc) = [];
                        end
                    end
                end
            end
        end

        % add menu to load all regions
        uimenu(mmh, 'Label', 'Download all ROIs from NeuroSynth...',  ...
            'Callback', {@ne_neurosynth, 'regions'}, 'Separator', 'on');

    % posterior probability map
    case {'ppmap'}

        % load all PP maps
        ppmaps = findfiles([npath '/terms'], 'pp_*.nii.gz', 'depth=1');

        % empty, return
        if isempty(ppmaps)
            return;
        end

        % create term names (same files without dates, etc.)
        ppterms = regexprep(regexprep(findfiles([npath '/terms'], ...
            'pp_*.nii.gz', 'depth=1', 'relative='), '^pp_', ''), '_.*$', '');

        % load first
        pppart = xff(ppmaps{1});
        ppmap = pppart.CopyObject;
        pppart.ClearObject;
        ppmap.RunTimeVars.AutoSave = true;
        ppmap.RunTimeVars.IntensityLabels = ppterms;

        % initialize volume and data
        ppmap.VoxelData = zeros(size(ppmap.VoxelData));
        ppmap.VoxelData(:, :, :, 2) = 0;
        ppmap.ImgDim.DataType = 64;
        ppmap.ImgDim.Dim(1) = 4;
        ppmap.ImgDim.Dim(5) = 2;
        ppmap.ImgDim.BitsPerPixel = 64;
        ppmap.ImgDim.ScalingSlope = 1;
        ppmap.ImgDim.ScalingIntercept = 0;
        ppplab = ppmap.VoxelData(:, :, :, 1);
        pppmax = ppmap.VoxelData(:, :, :, 2);

        % extend Map
        ppmap.RunTimeVars.Map(2) = ppmap.RunTimeVars.Map(1);
        ppmap.RunTimeVars.Map(1).Name = ...
            sprintf('Region labels by intensity (%d terms)', numel(ppterms));
        ppmap.RunTimeVars.Map(2).Name = 'Max. posterior probability';

        % load terms and set
        for tc = 1:numel(ppterms)
            try
                pppart = xff(ppmaps{tc});
                pppart.LoadVoxelData;
                pppmap = pppart.ImgDim.ScalingIntercept + ...
                    (pppart.ImgDim.ScalingSlope .* double(pppart.VoxelData));
                pppart.ClearObject;
                pppmap(isinf(pppmap) | isnan(pppmap)) = 0;
                ppmask = (pppmap >= 0.5 & pppmap > pppmax);
                if any(ppmask(:))
                    ppplab(ppmask) = tc;
                    pppmax(ppmask) = pppmap(ppmask);
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end

        % set result
        ppmap.VoxelData = cat(4, ppplab, pppmax);

        % show
        ppmap.Browse;

    % regions
    case {'regions'}

        % load region terms
        regions = splittocellc(asciiread([npath '/rterms.txt']), ',');

        % progress bar
        pbar = ch.Progress;
        cprog = ne_progress(0, 0, {true, 0, 'Loading terms...'});

        % load
        for rc = 1:numel(regions)
            for tc = {'fi', 'ri', 'pp'}
                try
                    ne_neurosynth(0, 0, 'load', regions{rc}, tc{1});
                    drawnow;
                    ne_gcfg.fcfg.StatsVar.ClearObject;
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end
            end
            pbar.Progress(rc / numel(regions), sprintf('Loaded %s (%d/%d)', ...
                regions{rc}, rc, numel(regions)));
        end

        % reset progress bar
        ne_progress(0, 0, cprog);

    % reload
    case {'reload'}

        % requires filename
        if nargin < 4 || ...
           ~ischar(varargin{4}) || ...
            isempty(varargin{4}) || ...
            exist([npath '/terms/' varargin{4}(:)'], 'file') ~= 2
            return;
        end

        % try to load file
        try
            h = xff([npath '/terms/' varargin{4}(:)']);
            if ~isxff(h, 'hdr')
                if isxff(h, true)
                    h.ClearObject;
                end
                error( ...
                    'neuroelf:BadFilename', ...
                    'Not a HDR/NII file.' ...
                );
            end
            h.LoadVoxelData;
            if ~isfield(h.RunTimeVars, 'AutoSave') || ...
               ~h.RunTimeVars.AutoSave
                h.RunTimeVars.AutoSave = true;
                h.RunTimeVars.StatsObject = true;
                h.RunTimeVars.Map.DF1 = 1000;
                h.RunTimeVars.Map.UpperThreshold = double(max(h.VoxelData(:))) + sqrt(eps);
                if lower(varargin{4}(1)) ~= 'p'
                    h.RunTimeVars.Map.Type = 1;
                    h.RunTimeVars.Map.LowerThreshold = 0.25 .* h.RunTimeVars.Map.UpperThreshold;
                else
                    h.RunTimeVars.Map.Type = 11;
                    h.RunTimeVars.Map.LowerThreshold = ...
                        double(min(h.VoxelData(h.VoxelData > 0.5)));
                    if isempty(h.RunTimeVars.Map.LowerThreshold)
                        h.RunTimeVars.Map.LowerThreshold = 0.25 .* h.RunTimeVars.Map.UpperThreshold;
                    end
                end
                h.RunTimeVars.Map.ClusterSize = 10;
                h.RunTimeVars.Map.EnableClusterCheck = 1;
            end

            % ensure MapSelection
            if ~isfield(h.RunTimeVars, 'MapSelection')
                h.RunTimeVars.MapSelection = {{}, []};
            end

            ne_openfile(0, 0, h, true);
        catch ne_eo;
            rethrow(ne_eo);
        end
end

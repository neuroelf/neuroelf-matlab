% FUNCTION ne_averagenii: average NII files (from preprocessing)
function varargout = ne_averagenii(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:32 AM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% action?
if nargin > 2 && ...
    ischar(varargin{3}) && ...
   ~isempty(varargin{3})
    action = lower(varargin{3}(:)');
else
    action = 'opengui';
end

% all actions but opengui require figure to be open!
if ~strcmp(action, 'opengui') && ...
   (~isfield(ne_gcfg.cc, 'AverageNII') || ...
    ~isxfigure(ne_gcfg.cc.AverageNII.Satellite, true))
    return;
elseif ~strcmp(action, 'opengui')
    tags = ne_gcfg.cc.AverageNII.Tags;
end

% depending on action
switch (action)

    % edit anatomical folder/file pattern
    case {'anatedit'}

        % disable search button
        ne_gcfg.cc.AverageNII.Satellite.SetGroupEnabled('FoldOK', 'off');
        drawnow;

        % check pattern for validity
        [aff, afp, afe] = fileparts(tags.ED_averagenii_anat.String);
        if isempty(aff)
            aff = pwd;
        end
        afpe = [afp afe];
        try
            afiles = findfiles(aff, afpe, 'depth=1');
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            afiles = {};
        end
        tags.ED_averagenii_norm.String = '';
        if ~isempty(afiles)
            ne_gcfg.cc.AverageNII.Satellite.SetGroupEnabled('FoldOK', 'on');
            try
                nfiles = findfiles(aff, [afp '_seg_sn.mat'], 'depth=1');
                if numel(nfiles) == numel(afiles)
                    tags.ED_averagenii_norm.String = [aff filesep afp '_seg_sn.mat'];
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end

    % search for files
    case {'afsearch'}

        % disable search button (no double click)
        ne_gcfg.cc.AverageNII.Satellite.SetGroupEnabled('FoldOK', 'off');
        drawnow;

        % find files
        [aff, afp, afe] = fileparts(tags.ED_averagenii_anat.String);
        if isempty(aff)
            aff = pwd;
        end
        afp = [afp afe];
        try
            afiles = findfiles(aff, afp, 'depth=1');
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            afiles = {};
        end
        if tags.CB_averagenii_usenorm.Value > 0
            [aff, afp, afe] = fileparts(tags.ED_averagenii_norm.String);
            if isempty(aff)
                aff = pwd;
            end
            afp = [afp afe];
            try
                nfiles = findfiles(aff, afp, 'depth=1');
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
                nfiles = repmat({''}, numel(afiles), 1);
            end
            if numel(afiles) ~= numel(nfiles)
                afiles = {};
                nfiles = {};
            end
        else
            nfiles = repmat({''}, numel(afiles), 1);
        end
        if ~isempty(afiles)
            ne_gcfg.cc.AverageNII.Config.afiles = ...
                [ne_gcfg.cc.AverageNII.Config.afiles(:); afiles(:)];
        end
        if ~isempty(nfiles)
            ne_gcfg.cc.AverageNII.Config.nfiles = ...
                [ne_gcfg.cc.AverageNII.Config.nfiles(:); nfiles(:)];
        end
        afiles = ne_gcfg.cc.AverageNII.Config.afiles;
        nfiles = ne_gcfg.cc.AverageNII.Config.nfiles;
        sfiles = cell(numel(afiles), 1);
        for fc = 1:numel(afiles)
            if numel(afiles{fc}) > 80
                afiles{fc} = [afiles{fc}(1:37) '...' afiles{fc}(end-36:end)];
            end
            if ~isempty(nfiles{fc})
                [afp, nfiles{fc}] = fileparts(nfiles{fc});
                sfiles{fc} = sprintf('%s (norm: %s)', afiles{fc}, nfiles{fc});
            else
                sfiles{fc} = afiles{fc};
            end
        end
        if isempty(nfiles)
            sfiles = afiles;
        end
        tags.LB_averagenii_af_found.String = sfiles;
        tags.LB_aveargenii_af_found.Value = [];
        tags.LB_averagenii_af_found.ListboxTop = 1;
        if ~isempty(afiles)
            ne_gcfg.cc.AverageNII.Satellite.SetGroupEnabled('FilesOK', 'on');
        end

    % do the work
    case {'average'}

        % get configuration
        acfg = ne_gcfg.cc.AverageNII.Config;

        % parse some options
        avgtypes = {'mean', 'mean50', 'robmean', 'median'};
        bboxval = tags.DD_averagenii_wbb.Value;
        switch (bboxval)
            case {1}
                bbox = [-127, -127, -127; 128, 128, 128];
            case {2}
                bbox = [-78, -111, -71; 78, 84, 84];
            case {3}
                bbox = [-89, -119, -71; 90, 96, 96];
        end
        imeths = strrep(tags.DD_averagenii_imeth.String, 'sinc', 'lanczos');
        if ~iscell(imeths)
            imeths = cellstr(imeths);
        end
        imeths = strrep(imeths, 'sinc', 'lanczos');
        indivvmr = '';
        if tags.CB_averagenii_writesn.Value > 0
            indivvmr = '_SN.vmr';
        end
        maskc3 = false;
        if tags.DD_averagenii_mask.Value > 1
            maskc1c2 = true;
            if tags.DD_averagenii_mask.Value > 2
                maskc3 = true;
            end
        else
            maskc1c2 = false;
        end
        if tags.DD_averagenii_wvox.Value == 1
            res = 0.5;
        else
            res = 1;
        end

        % create configuration for function
        cprog = ne_progress(0, 0, {true, 0, 'averagenii'});
        aopt = struct( ...
            'avgtype',  avgtypes{tags.DD_averagenii_mmeth.Value}, ...
            'bbox',     bbox, ...
            'clean',    false, ...
            'cnorm',    true, ...
            'ihc',      (tags.CB_averagenii_inhcorr.Value > 0), ...
            'imeth',    imeths{tags.DD_averagenii_imeth.Value}, ...
            'indivvmr', indivvmr, ...
            'maskc1c2', maskc1c2, ...
            'maskc3',   maskc3, ...
            'outtype',  'VMR', ...
            'pbar',     ne_gcfg.h.Progress, ...
            'res',      res, ...
            'snorm',    (tags.CB_averagenii_usenorm.Value > 0));

        % make figure invisible
        ne_gcfg.cc.AverageNII.Satellite.Visible = 'off';
        drawnow;

        % attempt processing
        try
            vmr = averagestruct(acfg.afiles, aopt);

            % if OK, delete dialog
            ne_averagenii(0, 0, 'closegui');

            % then browse result
            vmr.Browse;

        % error occurred...
        catch ne_eo;
            ne_gcfg.cc.AverageNII.Satellite.Visible = 'on';
            warndlg(ne_eo.message, 'NeuroElf - error', 'modal');
        end
        ne_progress(0, 0, cprog);

    % browse for (single) file
    case {'browse'}

        % request file
        [affile, afpath] = uigetfile( ...
            {'*.hdr;*.nii', 'Analyze files (*.hdr, *.nii)'; ...
             '*.vmr',       'BrainVoyager VMR files (*.vmr)'}, ...
            'Please select one of the anatomical files...', ...
            'MultiSelect', 'off');
        if isequal(affile, 0) || ...
            isequal(afpath, 0) || ...
            isempty(affile)
            return;
        end
        if isempty(afpath)
            afpath = pwd;
        elseif any(afpath(end) == '/\')
            afpath(end) = [];
        end
        affile = [afpath '/' affile];
        tags.ED_averagenii_anat.String = affile;
        [aff, afp, afe] = fileparts(affile);
        if any(strcmpi(afe, {'.hdr', '.img', '.nii'}))
            tags.ED_averagenii_norm.String = [aff filesep afp '_seg_sn.mat'];
        else
            tags.ED_averagenii_norm.String = '';
        end
        ne_gcfg.cc.AverageNII.Satellite.SetGroupEnabled('FoldOK', 'on');

    % close GUI
    case {'closegui'}

        % close GUI
        ne_gcfg.cc.AverageNII.Satellite.Delete;

        % remove field in global struct
        ne_gcfg.cc = rmfield(ne_gcfg.cc, 'AverageNII');

    % delete from list
    case {'delete'}

        % remove from lists
        didx = tags.LB_averagenii_af_found.Value;
        if isempty(didx)
            return;
        end
        ne_gcfg.cc.AverageNII.Config.afiles(didx) = [];
        ne_gcfg.cc.AverageNII.Config.nfiles(didx) = [];
        tags.LB_averagenii_af_found.Value = [];
        dstr = tags.LB_averagenii_af_found.String;
        if ~iscell(dstr)
            dstr = cellstr(dstr);
        end
        dstr(didx) = [];
        tags.LB_averagenii_af_found.String = dstr;
        tags.LB_averagenii_af_found.Value = [];
        tags.LB_averagenii_af_found.ListboxTop = ...
            max(1, min(tags.LB_averagenii_af_found.ListboxTop, numel(dstr)));
        if isempty(dstr)
            ne_gcfg.cc.AverageNII.Satellite.SetGroupEnabled('FilesOK', 'off');
        end

    % open GUI
    case {'opengui'}

        % only allow one instance to run
        if any(strcmpi(ne_gcfg.c.blockcb, 'averagenii'))
            return;
        end
        ne_gcfg.c.blockcb{end+1} = 'averagenii';

        % safety net
        try

            % open GUI
            hGUI = [];
            hGUI = xfigure([neuroelf_path('tfg') '/ne_averagenii.tfg']);
            hGUI.WindowStyle = 'modal';

            % set the close-request callback
            hGUI.CloseRequestFcn = {@ne_averagenii, 'closegui'};

            % get tags
            tags = hGUI.TagStruct;

            % set callbacks with function handle
            tags.ED_averagenii_anat.Callback = {@ne_averagenii, 'anatedit'};
            tags.BT_averagenii_af_browse.Callback = {@ne_averagenii, 'browse'};
            tags.BT_averagenii_af_search.Callback = {@ne_averagenii, 'afsearch'};
            tags.CB_averagenii_usenorm.Callback = {@ne_averagenii, 'usenorm'};
            tags.ED_averagenii_norm.Callback = {@ne_averagenii, 'normedit'};
            tags.BT_averagenii_nf_browse.Callback = {@ne_averagenii, 'nfbrowse'};
            tags.BT_averagenii_delsubj.Callback = {@ne_averagenii, 'delete'};
            tags.ED_averagenii_vmrfpat.Callback = {@ne_averagenii, 'checkout'};
            tags.BT_averagenii_cancel.Callback = {@ne_averagenii, 'closegui'};
            tags.BT_averagenii_average.Callback = {@ne_averagenii, 'average'};

            % correct settings for list of files
            tags.LB_averagenii_af_found.ListboxTop = 1;
            tags.LB_averagenii_af_found.Value = [];
            tags.LB_averagenii_af_found.String = {};

            % store in global struct
            ne_gcfg.cc.AverageNII = struct( ...
                'Config',       struct( ...
                    'afiles',   {{}}, ...
                    'nfiles',   {{}}), ...
                'Satellite',    hGUI, ...
                'SatelliteMLH', hGUI.MLHandle, ...
                'Tags',         tags);

            % wait until GUI is closed
            hGUI.HandleVisibility = 'callback';
            hGUI.Visible = 'on';
            uiwait(hGUI.MLHandle);

        % safety net (part two)
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            if isxfigure(hGUI, true)
                try
                    hGUI.Delete;
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end
            end
        end

        % unblock calls
        ne_gcfg.c.blockcb(strcmpi(ne_gcfg.c.blockcb, 'averagenii')) = [];

	% switch use-norm
    case {'usenorm'}

        % clear list of files
        ne_gcfg.cc.AverageNII.Config.afiles = {};
        ne_gcfg.cc.AverageNII.Config.nfiles = {};
        tags.LB_averagenii_af_found.ListboxTop = 1;
        tags.LB_averagenii_af_found.Value = [];
        tags.LB_averagenii_af_found.String = {};

        % disable/enable edit box
        if tags.CB_averagenii_usenorm.Value > 0
            une = 'on';
        else
            une = 'off';
        end
        ne_gcfg.cc.AverageNII.Satellite.SetGroupEnabled('UseNorm', une);
end

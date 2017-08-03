function [fig, rpf] = showspmrparams(folder, opts)
% showspmrparams  - display the rp_*.txt files of an SPM preprocessing
%
% FORMAT:       [fig, rp] = showspmrparams([folder [, opts]]);
%
% Input fields:
%
%       folder      folder containing the aligned volumes / rp_*.txt files
%       opts        optional settings
%        .color     1x6 cell array with valid color tokens for plot(...)
%        .legend    boolean flag, add legend to the axes (default: true)
%        .markruns  boolean flag, add vertical lines at runs (default: true)
%        .source    filename pattern, default 'rp_*.txt'
%        .split     either off 'complete', 'none', or {'spm'}
%                   whereas 'complete' shows 6 plots, 'none' shows 1 plot,
%                   and 'spm' shows 2 plots, translation and rotation
%
% Output fields:
%
%       fig         Matlab figure handle to output figure
%       rp          read parameters

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
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

% variable for UI stuff
global ne_ui;

% argument check
if nargin < 1 || ...
   ~ischar(folder) || ...
    isempty(folder) || ...
   (~any(folder(:)' == '*') && ...
    exist(folder(:)', 'dir') ~= 7)

    % be friendly
    folder = uigetdir(pwd, ...
        'Please select a subject''s (functional) folder...');
    if isempty(folder) || ...
        isequal(folder, 0) || ...
        exist(folder, 'dir') ~= 7
        fig = [];
        return;
    end
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'color') || ...
   ~iscell(opts.color) || ...
    numel(opts.color) ~= 6
    opts.color = ...
        {[0, 0, 0.75], [0, 0.75, 0], [0.75, 0, 0], ...
         [0.2, 0.2, 1], [0.2, 1, 0.2], [1, 0.2, 0.2]};
end
if ~isfield(opts, 'legend') || ...
   ~islogical(opts.legend) || ...
    numel(opts.legend) ~= 1
    opts.legend = true;
end
if ~isfield(opts, 'markruns') || ...
   ~islogical(opts.markruns) || ...
    numel(opts.markruns) ~= 1
    opts.markruns = true;
end
if ~isfield(opts, 'source') || ...
   ~ischar(opts.source) || ...
    numel(opts.source) < 5 || ...
   ~any(strcmpi(lsqueeze(opts.source(numel(opts.source)-3:end))', ...
        {'.hdr', '.img', '.mat', '.nii', '.txt'}))
    opts.source = 'rp_*.txt';
else
    opts.source = opts.source(:)';
end
ext = lower(opts.source(end-2:end));
if ~isfield(opts, 'split') || ...
   ~ischar(opts.split) || ...
   ~any(strcmpi(opts.split, {'f', 'full', 'n', 'none', 's', 'spm'}))
    opts.split = 's';
else
    opts.split = lower(opts.split(1));
end

% create or re-use figure
if isfield(ne_ui, 'showspmrparams') && ...
    isstruct(ne_ui.showspmrparams) && ...
    numel(ne_ui.showspmrparams) == 1 && ...
    isfield(ne_ui.showspmrparams, 'fig') && ...
    ishandle(ne_ui.showspmrparams.fig)
    fig = ne_ui.showspmrparams.fig;
else
    fig = figure;
    ne_ui.showspmrparams = struct('fig', fig);
end

% delete all figure children
delete(get(fig, 'Children'));

% find required files
rpf = findfiles(folder, opts.source);

% on empty, just return
if isempty(rpf)
    return;
end

% set figure title
set(fig, 'NumberTitle', 'off', 'Name', ...
    sprintf('Motion parameters of %s', folder));

% initialize some variables
rpl = zeros(numel(rpf), 1);
nruns = [];
iquat = eye(4);

% depending on filetype
switch (ext)

    % images
    case {'hdr', 'img', 'nii'}

        % if img, replace with hdr
        if strcmp(ext, 'img')
            rpf = regexprep(rpf, 'img$', 'hdr', 'preservecase');
        end

        % if number <= 16, think of this as runs, not volumes!
        if numel(rpf) <= 16

            % for each run get the time of last write access (order)
            rpt = zeros(numel(rpf), 1);
            for rpc = 1:numel(rpf)
                rpd = dir(rpf{rpc});
                rpt(rpc) = rpd.datenum;
            end

            % sort according to date
            [rpt, rpts] = sort(rpt, 'ascend');
            rpf = rpf(rpts);

        % otherwise, we need to detect the run lengths from filenames
        else
            rpd = rpf;
            for rpc = 1:numel(rpd)
                rpd{rpc} = fileparts(rpd{rpc});
            end
            nruns = numel(unique(rpd));
            nrun = 1;
            nvol = 1;
            rstr = rpd{1};
            for rpc = 2:numel(rpd)
                if ~strcmp(rpd{rpc}, rstr)
                    nruns(nrun) = nvol;
                    nrun = nrun + 1;
                    nvol = 1;
                    rstr = rpd{rpc};
                else
                    nvol = nvol + 1;
                end
            end
            nruns(nrun) = nvol;
        end

        % load images, one by one
        ih = cell(1, 1);
        for ic = 1:numel(rpf)
            try
                rpfile = rpf{ic};
                ih{1} = xff(rpfile);
                if ~isxff(ih{1}, 'hdr')
                    delete(fig);
                    error( ...
                        'neuroelf:BadFileContent', ...
                        'Not an Analyze/NIftI file.' ...
                    );
                end
            catch ne_eo;
                clearxffobjects(ih);
                rethrow(ne_eo);
            end

            % get coordinate frame and number of volumes
            rpf{ic} = ih{1}.CoordinateFrame.Trf;
            fds = size(ih{1}.VoxelData, 4);

            % clear object
            clearxffobjects(ih);

            % for 4D datasets
            if fds > 1

                % try to load accompanying mat file with trf matrices
                try
                    rpmatf = regexprep(rpfile, '(hdr|nii)$', 'mat', ...
                        'preservecase');
                    rpmat = load(rpmatf);

                    % check matfile
                    if ~isstruct(rpmat) || ...
                       ~isfield(rpmat, 'mat') || ...
                       ~isa(rpmat.mat, 'double') || ...
                       ~isequal(size(rpmat.mat), [4, 4, fds]) || ...
                        any(isinf(rpmat.mat(:)) | isnan(rpmat.mat(:))) || ...
                        any(lsqueeze(rpmat.mat(4, 1:3, :)) ~= 0) || ...
                        any(rpmat.mat(4, 4, :) ~= 1)
                        error('BAD_MATFILE');
                    end
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    delete(fig);
                    error( ...
                        'neuroelf:BadFileContent', ...
                        '4D files require a mat file with matching info.' ...
                    );
                end
                rpf{ic} = rpmat.mat;
            end

            % for first file get inverse
            if ic == 1
                iquat = inv(rpf{1}(:, :, 1));
            end
        end

    % MAT/RP text files
    case {'mat', 'txt'}

        % if number <= 16, think of this as runs, not volumes!
        if numel(rpf) <= 16

            % for each plot get the time of last write access (order)
            rpt = zeros(numel(rpf), 1);
            for rpc = 1:numel(rpf)
                rpd = dir(rpf{rpc});
                rpt(rpc) = rpd.datenum;
            end

            % sort according to date
            [rpt, rpts] = sort(rpt, 'ascend');
            rpf = rpf(rpts);
        end

        % load files
        for rpc = 1:numel(rpf)
            try
                rpf{rpc} = load(rpf{rpc});
                if (strcmp(ext, 'mat') && ...
                    (~isstruct(rpf{rpc}) || ...
                     ~isfield(rpf{rpc}, 'mat') || ...
                     ~isa(rpf{rpc}.mat, 'double') || ...
                      size(rpf{rpc}.mat, 1) ~= 4 || ...
                      size(rpf{rpc}.mat, 2) ~= 4 || ...
                      any(isinf(rpf{rpc}.mat(:)) | isnan(rpf{rpc}.mat(:))) || ...
                      any(lsqueeze(rpf{rpc}.mat(4, 1:3, :)) ~= 0) || ...
                      any(rpf{rpc}.mat(4, 4, :) ~= 1))) || ...
                   (strcmp(ext, 'txt') && ...
                    (isempty(rpf{rpc}) || ...
                     ~isa(rpf{rpc}, 'double') || ...
                      ndims(rpf{rpc}) > 2 || ...
                      size(rpf{rpc}, 2) ~= 6 || ...
                      any(isinf(rpf{rpc}(:)) | isnan(rpf{rpc}(:))) || ...
                      any(abs(rpf{rpc}(1, :)) > sqrt(eps))))
                    error('BAD_RP_FILE');
                end
                if isstruct(rpf{rpc})
                    rpf{rpc} = rpf{rpc}.mat;
                end
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                delete(fig);
                error( ...
                    'neuroelf:FileIOError', ...
                    'Realignment parameter file does not match specifications.' ...
                );
            end

            % for first file/run get inverse (for mat files)
            if rpc == 1 && ...
                strcmp(ext, 'mat')
                iquat = inv(rpf{1}(:, :, 1));
            end
        end
end

% now, we parse this information
for rpc = 1:numel(rpf)

    % this information is a 4-by-4-by-volumes mat
    if size(rpf{rpc}, 2) == 4

        % number of volumes
        rpl(rpc) = size(rpf{rpc}, 3);

        % then we compute a table with as many rows as volumes
        rft = zeros(rpl(rpc), 6);
        for volc = 1:rpl(rpc)

            % compute elements
            rfvol = spmitrf(rpf{rpc}(:, :, volc) * iquat);

            % and put into table
            rft(volc, :) = [rfvol{1}, rfvol{2}];
        end

        % then replace in rpf
        rpf{rpc} = rft;

    % this is a table from the txt files
    elseif strcmp(ext, 'txt')

        % number of volumes
        rpl(rpc) = size(rpf{rpc}, 1);

        % for the first run, we don't need anything
        if rpc == 1
            continue;
        end

        % for runs 2..N, create SPM-like (!) output (the actual true parameters
        % cannot be retrieved fully, as the misplacement between the last image
        % of the preceding run and the first image of the next run is LOST!)
        % get TRF of last image of previous run
        trf = spmtrf(rpf{rpc-1}(end, 1:3), rpf{rpc-1}(end, 4:6));

        % for this run, iterate over volumes
        rpv = rpf{rpc};
        for vc = 1:size(rpv, 1)

            % compute new TRF as product of old two, then decompose again
            itrf = spmitrf(spmtrf(rpv(vc, 1:3), rpv(vc, 4:6)) * trf);

            % then store in table
            rpv(vc, :) = [itrf{1}, itrf{2}];
        end

        % store table
        rpf{rpc} = rpv;
    end
end

% concatenate tables
rpf = cat(1, rpf{:});

% if necessary override number of volumes per run
if ~isempty(nruns)
    rpl = nruns;
end

% number of total volumes
nvol = (1:sum(rpl))';

% check against concatenated table
if opts.markruns && ...
    size(rpf, 1) ~= numel(nvol)
    opts.markruns = false;
    nvol = (1:size(rpf, 1))';
    warning( ...
        'neuroelf:InternalError', ...
        'Error detecting run lengths.' ...
    );
end

% do the radiens to degrees transform
rpf(:, 4:6) = (180 / pi) .* rpf(:, 4:6);

% create axes accordingly
switch (opts.split)

    % full split
    case {'f'}

        % create axes for 3 translation + 3 rotation parameters
        tx = subplot(6, 1, 1);
        hold(tx, 'on');
        ty = subplot(6, 1, 2);
        hold(ty, 'on');
        tz = subplot(6, 1, 3);
        hold(tz, 'on');
        rx = subplot(6, 1, 4);
        hold(rx, 'on');
        ry = subplot(6, 1, 5);
        hold(ry, 'on');
        rz = subplot(6, 1, 6);
        hold(rz, 'on');
        ax = [tx, ty, tz, rx, ry, rz];

        % titles
        title(tx, 'Translation-X');
        ylabel(tx, 'mm');
        title(ty, 'Translation-Y');
        ylabel(ty, 'mm');
        title(tz, 'Translation-Z');
        ylabel(tz, 'mm');
        title(rx, 'Rotation-X (pitch)');
        ylabel(rx, 'degrees');
        title(ry, 'Rotation-Y (roll)');
        ylabel(ry, 'degrees');
        title(rz, 'Rotation-Z (yaw)');
        ylabel(rz, 'degrees');

        % plot
        plot(tx, nvol, rpf(:, 1), 'Color', opts.color{1});
        plot(ty, nvol, rpf(:, 2), 'Color', opts.color{2});
        plot(tz, nvol, rpf(:, 3), 'Color', opts.color{3});
        plot(rx, nvol, rpf(:, 4), 'Color', opts.color{4});
        plot(ry, nvol, rpf(:, 5), 'Color', opts.color{5});
        plot(rz, nvol, rpf(:, 6), 'Color', opts.color{6});

        % legends
        if opts.legend
            legend(tx, 'trans-X');
            legend(ty, 'trans-Y');
            legend(tz, 'trans-Z');
            legend(rx, 'pitch');
            legend(ry, 'roll');
            legend(rz, 'yaw');
        end

    % no split
    case {'n'}
        tx = subplot(1, 1, 1);
        hold(tx, 'on');
        ax = tx;
        title(tx, 'Translation (mm) + Rotation (degrees)');
        trl = plot(tx, nvol, rpf);
        set(trl(1), 'Color', opts.color{1});
        set(trl(2), 'Color', opts.color{2});
        set(trl(3), 'Color', opts.color{3});
        set(trl(4), 'Color', opts.color{4});
        set(trl(5), 'Color', opts.color{5});
        set(trl(6), 'Color', opts.color{6});
        if opts.legend
            legend(tx, 'trans-X', 'trans-Y', 'trans-Z', 'pitch', 'roll', 'yaw');
        end

    % SPM-typical split
    case {'s'}
        tx = subplot(2, 1, 1);
        hold(tx, 'on');
        rx = subplot(2, 1, 2);
        hold(rx, 'on');
        ax = [tx, rx];
        title(tx, 'Translation');
        ylabel(tx, 'mm');
        title(rx, 'Rotation');
        ylabel(rx, 'degrees');
        tl = plot(tx, nvol, rpf(:, 1:3));
        rl = plot(rx, nvol, rpf(:, 4:6));
        set(tl(1), 'Color', opts.color{1});
        set(tl(2), 'Color', opts.color{2});
        set(tl(3), 'Color', opts.color{3});
        set(rl(1), 'Color', opts.color{4});
        set(rl(2), 'Color', opts.color{5});
        set(rl(3), 'Color', opts.color{6});
        if opts.legend
            legend(tx, 'trans-X', 'trans-Y', 'trans-Z');
            legend(rx, 'pitch', 'roll', 'yaw');
        end
end

% mark runs?
if opts.markruns

    % break volumes
    bvol = cumsum(rpl(1:end-1));

    % for each axes
    for ac = 1:numel(ax)

        % get y limits
        al = get(ax(ac), 'YLim');

        % for each break volume
        for bc = 1:numel(bvol)

            % plot vertical line
            plot(ax(ac), [bvol(bc); bvol(bc)], al(:), ...
                'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
        end
    end
end

% label X axes of bottom-most plot
rpls = sprintf('%d, ', rpl);
xlabel(ax(end), sprintf('%d runs with %d total volumes (%s volumes/run)', ...
    numel(rpl), size(rpf, 1), rpls(1:end-2)));

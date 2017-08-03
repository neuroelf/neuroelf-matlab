function xo = voi_ImportClusters(xo, file, opts)
% VOI::ImportClusters  - import clusters from a separate file
%
% FORMAT:       [voi = ] voi.ImportClusters(file [, opts]);
%
% Input fields:
%
%       file        filename containing additional clusters
%       opts        optional settings struct
%        .color     specific color for all added VOIs (default: [])
%        .ithresh   image value threshold (default: 1)
%        .kthresh   size threshold (only for image formats, default: 1)
%        .labels    filename with label information (e.g. from FreeSurfer)
%        .radius    sphere radius (for text files with < 4 columns, 10)
%        .sepclus   create separate VOIs for each cluster (default: true)
%        .sort      sorting criteria, either {'size'} or 'value
%
% Output fields:
%
%       voi         VOI with added VOI/s (integer coordinates only!)
%
% Note: if .color is empty (or not a valid RGB color), each VOI will have
%       a different, random color
%
%       the .sepclus flag is ignored when a VOI file is imported
%
% Using: asciiread, clustercoords, emptystruct, joinstructs, splittocellc,
%        u8str2double.

% Version:  v1.1
% Build:    16121617
% Date:     Dec-16 2016, 5:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011 - 2016, Jochen Weber
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

% neuroelf library
global ne_methods;

% argument check
if nargin ~= 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'voi') || ...
   ~ischar(file) || numel(file) < 5 || size(file, 2) ~= numel(file) || ...
   ~any(strcmpi(file(end-3:end), {'.hdr', 'i.gz', '.img', '.msk', '.nii', '.txt', '.voi'})) || ...
    exist(file, 'file') ~= 2
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
fext = lower(file(end-2:end));
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'color') || ~isa(opts.color, 'double') || ~isequal(size(opts.color), [1, 3]) || ...
    any(isinf(opts.color) | isnan(opts.color) | opts.color < 0 | opts.color > 255)
    opts.color = [];
elseif all(opts.color <= 1)
    opts.color = round(255 .* opts.color);
else
    opts.color = round(opts.color);
end
if ~isfield(opts, 'ithresh') || ~isa(opts.ithresh, 'double') || ...
    numel(opts.ithresh) ~= 1 || isinf(opts.ithresh) || isnan(opts.ithresh)
    opts.ithresh = 1 - sqrt(eps) * 16;
end
if ~isfield(opts, 'kthresh') || ~isa(opts.kthresh, 'double') || numel(opts.kthresh) ~= 1 || ...
    isinf(opts.kthresh) || isnan(opts.kthresh) || opts.kthresh < 1
    opts.kthresh = 1;
else
    opts.kthresh = floor(opts.kthresh);
end
if ~isfield(opts, 'labels') || ~ischar(opts.labels) || isempty(opts.labels) || ...
    exist(opts.labels(:)', 'file') ~= 2
    opts.labels = '';
end
if ~isfield(opts, 'radius') || ~isa(opts.radius, 'double') || numel(opts.radius) ~= 1 || ...
    isinf(opts.radius) || isnan(opts.radius) || opts.radius < 0
    opts.radius = 10;
else
    opts.radius = min(30, opts.radius);
end
if ~isfield(opts, 'sepclus') || ~islogical(opts.sepclus) || numel(opts.sepclus) ~= 1
    opts.sepclus = true;
end
if ~isfield(opts, 'sort') || ~ischar(opts.sort) || isempty(opts.sort) || ...
   ~any(lower(opts.sort(1)) == 'sv')
    opts.sort = 's';
else
    opts.sort = lower(opts.sort(1));
end

% get content
bc = xo.C;

% coordinate valid in current system
if ~strcmpi(bc.ReferenceSpace, 'tal')
    error('neuroelf:xff:badArgument', 'Import requires VOI to be in TAL space.');
end
bb = [-128, 128];

% what to do
switch (fext)

    % for image files
    case {'.gz', 'hdr', 'img', 'msk', 'nii'}

        % replace img with hdr
        if strcmp(fext, 'img')
            file = regexprep(file, '\.img$', '.hdr', 'preservecase');
        end

        % try to load MSK file
        try
            msk = [];
            msk = xff(file);
            if strcmp(fext, 'msk')
                if ~xffisobject(msk, true, 'msk')
                    error('neuroelf:xff:badFileContent', 'Not a MSK file: ''%s''.', file);
                end
                mskc = msk.C;
                mskb = aft_BoundingBox(msk);
                trf = mskb.QuatB2T';
                mskc = (mskc.Mask ~= 0);
                msku = 1;
            else
                if ~xffisobject(msk, true, 'hdr')
                    error('neuroelf:xff:badFileContent', 'Not an Analyze file: ''%s''.', file);
                end
                mskv = aft_GetVolume(msk, 1);
                mskb = hdr_CoordinateFrame(msk);
                trf = mskb.Trf';
                mskc = (mskv >= opts.ithresh);
                msku = unique(mskv(mskc));

                % patch extension
                if numel(unique(mskv(:))) < 4
                    fext = 'msk';
                end
            end
            delete(msk);
        catch xfferror
            if numel(msk) == 1 && xffisobject(msk, true)
                delete(msk);
            end
            rethrow(xfferror);
        end
        [fpath, file, filee] = fileparts(file);
        file = [file filee];

        % cluster as single cluster
        if ~opts.sepclus || numel(msku) == 1 || ~all(msku == round(msku))
            [mskcrd, mskcrdc, msks] = ne_methods.clustercoords(mskc, 2, opts.kthresh);
        else
            mskcrdc = cell(1, numel(msku));
            msks = cell(1, numel(msku));
            for cc = 1:numel(msku)
                [mskcrd, mskcrdc{cc}, msks{cc}] = ...
                    ne_methods.clustercoords(mskv == msku(cc), 2, opts.kthresh);
                if numel(mskcrdc{cc}) > 1
                    mskcrdc{cc} = {cat(1, mskcrdc{cc}{:})};
                    msks{cc} = size(mskcrdc{cc}{1}, 1);
                end
            end
            mskcrdc = cat(2, mskcrdc{:});
            msks = cat(2, msks{:});
        end

        % and sort by size
        if opts.sort == 's' && numel(msks) > 1
            [msks, msksi] = sort(msks(:), 'descend');
            mskcrdc = mskcrdc(msksi);
            msku = msku(msksi);
        end
        if isempty(mskcrdc)
            return;
        end

        % as one cluster
        if ~opts.sepclus

            msks = sum(msks);
            mskcrdc = {cat(1, mskcrdc{:})};
        end

        % create space
        voi = ne_methods.emptystruct({'Name', 'Color', 'NrOfVoxels', 'Voxels'}, ...
            [numel(mskcrdc), 1]);

        % load labels file
        if ~isempty(opts.labels)
            lfcont = ne_methods.splittocellc(ne_methods.asciiread(opts.labels), ...
                char([10, 13]), true, true);
            lfcont(cellfun('isempty', lfcont) | ~cellfun('isempty', regexpi(lfcont, '^\s*\#'))) = [];
            lfnum = regexpi(lfcont, '^\s*(\d+)\s+', 'match');
            lfnum(cellfun('isempty', lfnum)) = {{'NaN'}};
            lfnum = str2double(cat(1, lfnum{:}));
            lfcont = regexprep(lfcont, '^\s*\d+\s+', '');
            lflab = regexpi(lfcont, '^(\S+)', 'match');
            lflab(cellfun('isempty', lflab)) = {{'Unknown'}};
            lflab = cat(1, lflab{:});
            cnames = cell(numel(mskcrdc), 1);
            for cc = 1:numel(cnames)
                cnamel = find(lfnum == msku(cc));
                if ~isempty(cnamel)
                    cnames{cc} = lflab{cnamel(1)};
                else
                    cnames{cc} = sprintf('Cluster %04d (value %d) from %s', ...
                        cc, lfnum(cc), file);
                end
            end
        else
            lfnum = msku;
            cnames = cell(numel(mskcrdc), 1);
            for cc = 1:numel(cnames)
                cnames{cc} = sprintf('Cluster %04d (value %d) from %s', cc, lfnum(cc), file);
            end
        end

        % set color if given
        if ~isempty(opts.color)
            color = opts.color;
        end

        % transform to TAL space
        for cc = 1:numel(mskcrdc)

            % for msk
            if strcmp(fext, 'msk')

                % sort by distance to mean
                mc = mean(mskcrdc{cc}(:, 1:3), 1);
                dv = sqrt(sum((mskcrdc{cc}(:, 1:3) - mc(ones(1, msks(cc)), :)) .^ 2, 2));
                [dv, dvi] = sort(dv);
                mskcrd = mskcrdc{cc}(dvi, :);

            % for Analyze/NIFTI
            else

                % sort by descending value
                [mv, mvi] = sort(mskv(sub2ind(size(mskv), mskcrdc{cc}(:, 1), ...
                    mskcrdc{cc}(:, 2), mskcrdc{cc}(:, 3))), 'descend');
                mskcrd = mskcrdc{cc}(mvi, :);
            end

            % set fourth column to 1 for multiplication
            mskcrd(:, 4) = 1;

            % multiply (and round off)
            mskcrd = round(mskcrd * trf);

            % then fill VOI struct
            voi(cc).Name = cnames{cc};
            if isempty(opts.color)
                color = floor(255.999 .* rand(1, 3));
            end
            voi(cc).Color = color;
            voi(cc).NrOfVoxels = msks(cc);
            voi(cc).Voxels = mskcrd(:, 1:3);
        end

    % for VOI files
    case 'voi'

        % try to load file
        try
            newvoi = [];
            newvoi = xff(file);
            if ~xffisobject(newvoi, true, 'voi');
                error('neuroelf:xff:badFileContent', 'Not a VOI file: ''%s''.', file);
            end
            voi = newvoi.C;
            voi = voi.VOI(:);
            delete(newvoi);
        catch xfferror
            if numel(newvoi) == 1 && xffisobject(newvoi, true)
                delete(newvoi);
            end
            rethrow(xfferror);
        end

    % for text files
    case 'txt'

        % load file
        filecont = ne_methods.asciiread(file);
        if any(uint16(filecont) > 127)
            error('neuroelf:xff:badArgument', 'File contains invalid characters.');
        end

        % split to lines
        filecont = ne_methods.splittocellc(filecont, char([10, 13]), true, true);
        if numel(filecont) > 1000
            error('neuroelf:xff:badArgument', ...
                'Too many lines in file. Unlikely to be a cluster file!');
        end
        cluscont = repmat({zeros(0, 3)}, numel(filecont), 1);
        clussize = opts.radius(ones(numel(filecont), 1), 1);

        % for each line
        for lc = 1:numel(filecont)

            % contains a cluster?
            cfound = regexp(filecont{lc}, '([\+\-]?\d+[\,\s]+[\+\-]?\d+[\,\s]+[\+\-]?\d+)(\,?\s+\d+)?', 'tokens');
            if ~isempty(cfound)

                % center
                cc = ne_methods.u8str2double(cfound{1}{1});

                % fill cluscont accordingly
                if isempty(cfound{2})
                    cluscont{lc} = sphroi(cc, opts.radius, bb);
                else
                    clussize(lc) = min(30, max(0, ne_methods.u8str2double(cfound{1}{2})));
                    cluscont{lc} = sphroi(cc, clussize(lc), bb);
                end
            else
                clussize(lc) = -1;
            end
        end

        % don't join
        if opts.sepclus

            % remove empty clusters
            cluscont(clussize < 0) = [];
            clussize(clussize < 0) = [];

            % create space
            voi = ne_methods.emptystruct({'Name', 'Color', 'NrOfVoxels', 'Voxels'}, ...
                [numel(cluscont), 1]);

            % create each cluster separately
            if ~isempty(opts.color)
                color = opts.color;
            end
            for lc = 1:numel(voi)
                if isempty(opts.color)
                    color = floor(255.999 .* rand(1, 3));
                end
                voi(lc).Name = sprintf('%d-mm sphere around %d, %d, %d', ...
                    clussize(lc), cluscont{lc}(1, :));
                voi(lc).Color = color;
                voi(lc).NrOfVoxels = size(cluscont{lc}, 1);
                voi(lc).Voxels = cluscont{lc};
            end

        % join
        else
            if isempty(opts.color)
                opts.color = floor(255.999 .* rand(1, 3));
            end
            cluscont = unique(cat(1, cluscont{:}), 'rows');
            voi = struct('Name', sprintf('%d spheres', numel(clussize)), ...
                'Color', opts.color, 'NrOfVoxels', size(cluscont, 1), 'Voxels', cluscont);
        end

end

% update VOI and NrOfVOIs
if numel(fieldnames(voi)) == numel(fieldnames(bc.VOI)) && ...
    all(strcmp(fieldnames(voi), fieldnames(bc.VOI)))
    bc.VOI(end+1:end+numel(voi)) = voi;
else
    bc.VOI = ne_methods.joinstructs(bc.VOI(:), voi(:));
end
bc.NrOfVOIs = numel(bc.VOI);

% set back into storage
xo.C = bc;


% sub-functions


function xg = sphroi(c, r, bb)

% create grid
rc = round(c);
cr = ceil(r + 0.5);
[xg, yg, zg] = ndgrid( ...
    rc(1)-cr:rc(1)+cr, rc(2)-cr:rc(2)+cr, rc(3)-cr:rc(3)+cr);

% fill with grid
xg = [xg(:), yg(:), zg(:)];

% remove voxels beyond bounding box
xg(any(xg < bb(1), 2) | any(xg > bb(2), 2), :) = [];

% compute distance
yg = sqrt(sum((xg - c(ones(1, size(xg, 1)), :)) .^ 2, 2));

% sort by distance
[yg, ygi] = sort(yg);
xg = xg(ygi, :);

% remove further away stuff
xg(yg > r, :) = [];

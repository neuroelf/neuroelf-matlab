function [rstring, xlabel, xlabelc] = tdlocal2(rtype, xc, yc, zc, cs)
% tdlocal2  - lookup local TD NiFTi image file for TalLabel
%
%       FORMAT:     tdlabel = tdlocal2(rtype, x, y, z, size)
%  or
%       FORMAT:     tdvoi   = tdlocal2(rtype, voiname)
%
% Input fields:
%       rtype       numeric type selection: 2: TalLabel, 3: TalCube, 5: NGM
%                   6: TalVOI
%       x,y,z       Talairach coordinates of request
%       size        only required for types 3 and 5; where with type 3
%                   this is the width, height, and depth of the searched
%                   cube, while with type 5 its a maximum distance;
%                   for type 3, its maximum is 13, for type 5 it's 8.5
%
%       voiname     string that is matched against the list of labels
%
% as with the original network database, tdlocal2 returns its reply in
% one continuous row. hence, it's preferable to use tdclient.m instead for
% requests.
%
% for an extension, NGM search does **not** search in a cubed pattern, but
% rather spherically, hence it is far more accurate than the same function
% in some other localized TD client versions.
%
% See also tdclient, tdlocal.

% Version:  v1.1
% Build:    16060811
% Date:     Jun-08 2016, 11:19 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% enough arguments ?
if nargin < 2 || ...
   ~isa(rtype, 'double') || ...
    numel(rtype) ~= 1 || ...
    isinf(rtype) || ...
    isnan(rtype) || ...
    rtype < -6 || ...
    rtype > 8 || ...
   (~ischar(xc) && nargin < 4)
    error( ...
        'neuroelf:TooFewArguments',...
        'Bad or too few arguments.' ...
    );
end

% persistent variable for internal database representation
persistent tnii_opts;
if isempty(tnii_opts)

    % db is not loaded yet, loading (and saving) should be tried though
    tnii_opts.dbloaded = false;

    % load talairach image
    tdsrcfld = neuroelf_path('tal');
    talimage = [];
    try
        talimage = xff([tdsrcfld '/talairach.nii']);
        tnii_opts.taldb = talimage.VoxelData(:, :, :);
        xyzsize = size(tnii_opts.taldb);
        tnii_opts.xsize  = xyzsize(1);
        tnii_opts.ysize  = xyzsize(2);
        tnii_opts.zsize  = xyzsize(3);
        tnii_opts.s2i = [1, xyzsize(1), xyzsize(1) * xyzsize(2)];
        tnii_opts.xtrans = 1 - talimage.DataHist.NIftI1.QuatOffsetX;
        tnii_opts.ytrans = 1 - talimage.DataHist.NIftI1.QuatOffsetY;
        tnii_opts.ztrans = 1 - talimage.DataHist.NIftI1.QuatOffsetZ;
        labcount = double(max(tnii_opts.taldb(:)));
        tnii_opts.labels = splittocell( ...
            strrep(char(talimage.IntermedData(13:end)), '.', ','), ...
            char(10), 1);
        tnii_opts.labels(1) = [];
        tnii_opts.labels(end) = [];
        talimage.ClearObject;
        if length(tnii_opts.labels) ~= labcount
            error('WRONG_LABEL_NUMBER');
        end
    catch ne_eo;
        clearxffobjects({talimage});
        rethrow(ne_eo);
    end

    % clear no longer needed memory
    clear talimage;

    % parsing labels
    plabels = cell(labcount, 5);
    try

        % split labels into particles
        for lc = 1:labcount
            plabels(lc,:) = splittocell(lower(tnii_opts.labels{lc}), ',');
        end

        % check each column for unique labels
        labtargets = struct;
        for cc = 1:5
            ulab = unique(plabels(:, cc));
            for lc = 1:length(ulab)
                if ~strcmp(ulab{lc}, '*')
                    labtargets.(makelabel(ulab{lc})) = find( ...
                        strcmp(plabels(:, cc), ulab{lc}));
                end
            end
        end

        % store in persistent variable
        tnii_opts.labtargets = labtargets;
    catch ne_eo;
        error( ...
            'neuroelf:InternalError', ...
            'Invalid Talairach label found (%s).', ...
            ne_eo.message ...
        );
    end
end

% performing argument parsing
if ischar(rtype)
    rtype = str2double(rtype);
end
if nargin > 3 && ischar(xc) && ischar(yc) && ischar(zc)
    try
        xc = str2double(xc);
        yc = str2double(yc);
        zc = str2double(zc);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error('neuroelf:general:badArgument', 'Invalid char argument given.');
    end
end
if nargin > 4 && ischar(cs)
    try
        cs = eval(['[' strrep(cs, '#', ',') ']']);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        cs = [7, 3];
    end
end
if rtype < 0
    rtype = -rtype;
    btype = 1;
else
    btype = 0;
end

% default set label components (empty) and return string (*)
rstring = '*,*,*,*,*';
xlabel  = 'l_0';
xlabelc = 0;
if numel(xc) > 1
    rstring = repmat({rstring}, numel(xc), 1);
    if nargout > 1
        xlabel = repmat({xlabel}, numel(xc), 1);
        if nargout > 2
            xlabelc = zeros(numel(xc), 1);
        end
    end
end

% different lookup types
switch (round(rtype))

    % Talairach label
    case {2}

        % get lookup x,y,z coordinates and sizes
        xc = round(xc + tnii_opts.xtrans);
        yc = round(yc + tnii_opts.ytrans);
        zc = round(zc + tnii_opts.ztrans);
        xm = tnii_opts.xsize;
        ym = tnii_opts.ysize;
        zm = tnii_opts.zsize;

        % check coordinates
        goodc = (xc > 0 & xc <= xm & yc > 0 & yc <= ym & zc > 0 & zc <= zm);
        if ~any(goodc)
            return;
        end

        % find label
        xlabelc(goodc) = double(tnii_opts.taldb( ...
            1 + ([xc(goodc), yc(goodc), zc(goodc)] - 1) * tnii_opts.s2i(:)));

        % any content
        if numel(xc) == 1 && xlabelc > 0
            rstring = tnii_opts.labels{xlabelc};
            if nargout > 1
                xlabel = sprintf('l_%d', xlabelc);
            end
        elseif numel(xc) > 1
            rstring(xlabelc > 0) = tnii_opts.labels(xlabelc(xlabelc > 0));
            if nargout > 1
                for goodc = find(goodc(:)')
                    xlabel{goodc} = sprintf('l_%d', xlabelc(goodc));
                end
            end
        end

    % Talairach cube lookup
    case {3}

        % only valid with 5 arguments, cs being double
        if nargin < 5 || ~isa(cs, 'double') || numel(cs) ~= 1
            error( ...
                'neuroelf:TooFewArguments',...
                'Argument ''size'' missing.' ...
            );
        end

        % get good cs from [3..13]
        cs  = min(13, max(3, cs));
        cs  = floor(cs / 2);

        % create data holding structs
        lss = struct;
        lsn = struct;

        % iterate over all coordinates
        for xr = (xc-cs):(xc+cs)
            for yr = (yc-cs):(yc+cs)
                for zr = (zc-cs):(zc+cs)

                    % recursively call
                    [lc, lcf] = tdlocal2(2, xr, yr, zr);

                    % check struct
                    if ~isfield(lss,lcf)
                        lss.(lcf) = 1;
                        lsn.(lcf) = lc;
                    else
                        lss.(lcf) = lss.(lcf) + 1;
                    end
                end
            end
        end

        % get filled in, distinct labels, init to zero (to sort by number)
        lt  = fieldnames(lsn);
        lti = zeros(size(lt));

        % iterate over fields (unique labels)
        for lout = 1:length(lt)
            lti(lout) = lss.(lt{lout});
        end

        % sort and init return string
        [ltii{1:2}] = sort(lti);
        ltii = ltii{2};
        rstring = '';

        % fill return string
        for lout = length(lt):-1:1
            rstring = sprintf('%s%d:%s:', rstring, lss.(lt{ltii(lout)}), lsn.(lt{ltii(lout)}));
        end

    % Nearest gray matter search
    case {5}

        % only valid with 5 arguments, cs being double
        if nargin < 5 || ~isa(cs, 'double') || numel(cs) ~= 2 || ...
            any(isinf(cs) | isnan(cs) | cs < 1 | cs > 13);
            error('neuroelf:general:tooFewArguments',...
                'Argument size missing or invalid.');
        end

        % get good range for sphere
        xcs      = min(12.501, max(1.19, cs(1)));
        radius   = fix(xcs + 0.8);

        % check direct coordinate first
        tresult = tdlocal2(2, xc, yc, zc);
        rstring = '*:*';

        % if found
        if strfind(lower(tresult), 'gray matter')
            if btype
                tresult = tdlocal2(-2, xc, yc, zc);
                rstring = {{[xc, yc, zc], [0, 0, 0], 0, tresult}};
            else
                rstring = ['+0,+0,+0  (dist=0.000): ' tresult];
            end
            return;
        end

        % get label over several positions
        [xr, yr, zr] = ndgrid(-radius:radius, -radius:radius, -radius:radius);
        dist = sqrt(xr(:) .* xr(:) + yr(:) .* yr(:) + zr(:) .* zr(:));
        distc = (dist <= xcs);
        dist = dist(distc);
        xr = xr(distc);
        yr = yr(distc);
        zr = zr(distc);
        tresult = tdlocal2(2, xr + xc, yr + yc, zr + zc);
        rfound = cellfun('isempty', regexpi(tresult, 'gray matter'));

        % no coordinate found
        if all(rfound)
            return;
        end

        % find distances
        dist(rfound) = [];
        xr(rfound) = [];
        yr(rfound) = [];
        zr(rfound) = [];
        [rsort, rindex] = sort(dist);

        % how many distinct labels at most?
        cs(2) = fix(cs(2) + 0.1);

        % init return string
        if btype
            rstring = {};
        else
            rstring = '';
        end

        % init struct
        lss = struct;

        % iterate over number of found occurrences
        for rcount = 1:numel(rsort)

            % get coordinate for index
            rii = rindex(rcount);
            mxc  = xc + xr(rii);
            myc  = yc + yr(rii);
            mzc  = zc + zr(rii);

            % on requested output type
            if btype

                % call tdlocal2 again
                [tresult, lcf] = tdlocal2(-2, mxc, myc, mzc);

                % put into output
                if ~isfield(lss, lcf)

                    % add to cell output
                    rstring{end+1} = {[mxc, myc, mzc], ...
                        [xr(rii), yr(rii), zr(rii)], rsort(rcount), tresult};

                    % make field available and reduce cs(2)
                    lss.(lcf) = 1;
                    cs(2) = cs(2)-1;
                end

            % direct output
            else

                % call tdlocal2 again
                [tresult, lcf] = tdlocal2(2, mxc, myc, mzc);

                % put into output
                if ~isfield(lss, lcf)
                    rstring = sprintf('%s%+d,%+d,%+d  (dist=%.3f): %s (coords=%d,%d,%d):', ...
                                      rstring, xr(rii), yr(rii), zr(rii), ...
                                      dist(rii), tresult, ...
                                      mxc, myc, mzc);

                    % make field available and reduce cs(2)
                    lss.(lcf) = 1;
                    cs(2) = cs(2)-1;
                end
            end

            % leave if number of hits was found
            if cs(2) < 1
                break;
            end
        end

    % Talairach label search
    case {6}

        % prepare output
        rstring = [zeros(0,1), zeros(0,1), zeros(0,1)];

        % check second arg
        if ~ischar(xc) || ...
            isempty(xc) || ...
           ~isfield(tnii_opts.labtargets, makelabel(lower(xc(:)')))
            return;
        end

        % find matching strings
        ccm = tnii_opts.labtargets.(makelabel(lower(xc(:)')));

        % find matching voxels
        mvx = false;
        mvx(1:tnii_opts.xsize, 1:tnii_opts.ysize, 1:tnii_opts.zsize) = false;
        for lc = ccm(:)'
            mvx(tnii_opts.taldb == uint16(lc)) = true;
        end
        [mvx, mvy, mvz] = ind2sub( ...
            [tnii_opts.xsize, tnii_opts.ysize, tnii_opts.zsize], ...
            find(mvx(:)));

        % create output array
        rstring = [ ...
            mvx - tnii_opts.xtrans, ...
            mvy - tnii_opts.ytrans, ...
            mvz - tnii_opts.ztrans];

    % return available labels (for search)
    case {8}

        % return sorted, camel-cased labels
        rstring = camelcase(sort(strrep(fieldnames(tnii_opts.labtargets), '_', ' ')));

        % ensure good order of BAs
        bas = find(~cellfun('isempty', regexpi(rstring, 'brodmann\s+area')));
        sd = cellfun('isempty', regexpi(rstring(bas), '\s+\d\d$'));
        rstring(bas) = rstring([lsqueeze(bas(sd)); lsqueeze(bas(~sd))]);

    % invalid type
    otherwise
        rstring = '* (unknown rtype in request!)';
end

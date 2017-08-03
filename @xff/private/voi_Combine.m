function [xo, res] = voi_Combine(xo, vspec, cmdname, cmdspec)
% VOI::Combine  - combine regions within a VOI object
%
% FORMAT:       [voi, res] = voi.Combine(vspec, cmdname [,cmdspec])
%
% Input fields:
%
%       vspec       VOI region selection, 1xN list or regexp pattern
%       cmdname     either of
%                   - 'intersect'  (BrainVoyager QX's "a AND b")
%                   - 'overlap'    (probabilistic overlap between VOIs)
%                   - 'restrict'   (restrict VOI to shape)
%                   - 'setdiff'    (voxels in a NOT b)
%                   - 'union'      (BrainVoyager QX's "a OR b")
%       cmdspec     optional struct with settings
%        .bbox      2x3 mask/bounding box (in BVSystem coordinates)
%        .bothdirs  combine in both dirs (where useful, default false)
%        .cmbsubj   boolean flag, over-subject combination (default true)
%        .nostore   boolean flag, do not store combined VOIs into object,
%                   only return the VOI struct (*as first output*)
%        .rcenter   restriction center, if not given, use first coordinate
%        .rinplace  restrict in place, default: true
%        .rshape    restriction shape, one of 'box', {'sphere'}
%        .rsize     restriction size, default 8mm (radius)
%        .thresh    relative value, when e.g. overlap is set (default 0.5)
%
% Output fields:
%
%       voi         altered object (unless nostore is set to true)
%       res         resulting VOI structures (1xN struct)

% Version:  v1.1
% Build:    16012718
% Date:     Jan-27 2016, 6:18 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'voi') || ...
   ((~isa(vspec, 'double') || isempty(vspec) || ...
     any(isinf(vspec(:)) | isnan(vspec(:)) | vspec(:) < 0 | vspec(:) ~= fix(vspec(:)))) && ...
    (~ischar(vspec) || isempty(vspec) || ~any(vspec(:) == '+' | vspec(:) == '*'))) || ...
   ~ischar(cmdname) || isempty(cmdname) || ...
   ~any(strcmpi(cmdname(:)', {'intersect', 'overlap', 'restrict', 'setdiff', 'union'}))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
vspec = vspec(:)';
if isa(vspec, 'double') && any(vspec > numel(bc.VOI))
    error('neuroelf:xff:badArgument', 'Selected VOI(s) out of bounds.');
end
cmdname = lower(cmdname(:)');
if isempty(vspec) && strcmp(cmdname, 'restrict')
    vspec = 1:numel(bc.VOI);
end
if nargin < 4 || ~isstruct(cmdspec) || numel(cmdspec) ~= 1
    cmdspec = struct;
end
if ~isfield(cmdspec, 'bbox') || ~isa(cmdspec.bbox, 'double') || numel(size(cmdspec.bbox)) ~= 2 || ...
    any(size(cmdspec.bbox) ~= [2, 3]) || ...
    any(isnan(cmdspec.bbox(:)) | cmdspec.bbox(:) < 0 | cmdspec.bbox > 255)
    cmdspec.bbox = [];
else
    cmdspec.bbox = round(cmdspec.bbox);
end
bb = cmdspec.bbox;
if ~isfield(cmdspec, 'bothdirs') || ~islogical(cmdspec.bothdirs) || numel(cmdspec.bothdirs) ~= 1
    cmdspec.bothdirs = false;
end
if ~isfield(cmdspec, 'cmbsubj') || ~islogical(cmdspec.cmbsubj) || numel(cmdspec.cmbsubj) ~= 1
    cmdspec.cmbsubj = true;
end
if ~isfield(cmdspec, 'nostore') || ~islogical(cmdspec.nostore) || numel(cmdspec.nostore) ~= 1
    cmdspec.nostore = false;
end
if ~isfield(cmdspec, 'rcenter') || ~isa(cmdspec.rcenter, 'double') || numel(cmdspec.rcenter) ~= 3 || ...
    any(isinf(cmdspec.rcenter) | isnan(cmdspec.rcenter))
    cmdspec.rcenter = [];
end
if ~isfield(cmdspec, 'rinplace') || ~islogical(cmdspec.rinplace) || numel(cmdspec.rinplace) ~= 1
    cmdspec.rinplace = true;
end
if ~isfield(cmdspec, 'rshape') || ~ischar(cmdspec.rshape) || ...
   ~any(strcmpi(cmdspec.rshape, {'box', 'sphere'}))
    cmdspec.rshape = 'sphere';
else
    cmdspec.rshape = lower(cmdspec.rshape(:)');
end
if ~isfield(cmdspec, 'rsize') || ~isa(cmdspec.rsize, 'double') || numel(cmdspec.rsize) ~= 1 || ...
    isinf(cmdspec.rsize) || isnan(cmdspec.rsize) || cmdspec.rsize < 1 || cmdspec.rsize > 100
    cmdspec.rsize = 8;
end
if ~isfield(cmdspec, 'thresh') || ~isa(cmdspec.thresh, 'double') || numel(cmdspec.thresh) ~= 1 || ...
    isnan(cmdspec.thresh) || cmdspec.thresh <= 0 || cmdspec.thresh > 1
    cmdspec.thresh = 0.5;
end

% check vspec
vst = bc.VOI;
nvoi = numel(vst);
if ischar(vspec)
    vspecname = vspec;
    vspecname(vspec == '.' | vspec == '+' | vspec == '*' | vspec == '(' | vspec == ')' | ...
              vspec == '|' | vspec == '[' | vspec == ']' | vspec == '-' | vspec == '{' | ...
              vspec == '}' | vspec == ',' | vspec == '^' | vspec == '$' | vspec == '?') = [];
    nsv = false(1, nvoi);
    for vc = 1:nvoi
        rxr = regexpi(vst(vc).Name, vspec);
        if ~isempty(rxr)
            nsv(vc) = true;
        end
    end
    vspec = find(nsv);
else
    vspecname = '';
end
nsel = numel(vspec);
if nsel < 1 || (~strcmpi(cmdname, 'restrict') && nsel < 2)
    error('neuroelf:xff:badArgument', 'VOI selection too small');
end

% good new name pattern given
if isfield(bc, 'SubjectVOINamingConvention') && strcmpi(bc.SubjectVOINamingConvention, '<subj>_<voi>')
    subjvoi = true;
else
    subjvoi = false;
end
if isempty(vspecname)
    vnames = cell(nsel, 1);
    emptyname = false;
    if subjvoi
        vname1 = regexprep(vst(vspec(1)).Name, '^[^_]*_', '');
        if isempty(vname1)
            emptyname = true;
        end
        for vc = 1:nsel
            vnames{vc} = regexprep(vst(vspec(vc)).Name, '^[^_]*_', '');
            if isempty(vnames{vc})
                emptyname = true;
                break;
            end
        end
    else
        vname1 = regexprep(vst(vspec(1)).Name, '_[^_]*$', '');
        if isempty(vname1)
            emptyname = true;
        end
        for vc = 1:nsel
            vnames{vc} = regexprep(vst(vspec(vc)).Name, '_[^_]*$', '');
            if isempty(vnames{vc})
                emptyname = true;
                break;
            end
        end
    end
    if ~emptyname
        for cc = numel(vname1):-1:0
            if cc > 0 && numel(strmatch(vname1(1:cc), vnames)) == numel(vnames)
                break;
            end
        end
    else
        cc = 0;
    end
    if cc > 0
        vspecname = vname1(1:cc);
    else
        if ~isempty(vname1)
            vspecname = vname1;
        else
            vspecname = vst(vspec(1)).Name;
        end
    end
end

% for restriction, we follow a special line of code
if strcmp(cmdname, 'restrict')

    % get valid vspec
    vspec = intersect(1:numel(bc.VOI), vspec);
    nsel = numel(vspec);

    % not in place ?
    if ~cmdspec.rinplace
        tspec = (numel(bc.VOI) + 1):(numel(bc.VOI) + nsel);
        bc.VOI(end + nsel).Name = '';
    end

    % iterate over VOIs
    for vc = 1:nsel

        % get coordinates
        c = bc.VOI(vspec(vc)).Voxels;

        % subtract center
        if isempty(cmdspec.rcenter)
            cdiff = c - c(ones(size(c, 1), 1), :);
        else
            cdiff = c - cmdspec.rcenter(ones(size(c, 1), 1), :);
        end

        % compute distance
        if strcmp(cmdspec.rshape, 'box')
            cdist = max(abs(cdiff), [], 2);
        else
            cdist = sqrt(sum(cdiff .* cdiff, 2));
        end

        % remove bad coords
        c(cdist > cmdspec.rsize, :) = [];

        % put into target
        if cmdspec.rinplace
            if isempty(strfind(bc.VOI(vspec(vc)).Name, '_sph_'))
                bc.VOI(vspec(vc)).Name = ...
                    sprintf('%s_sph_%.1f', bc.VOI(vspec(vc)).Name, cmdspec.rsize);
            end
            bc.VOI(vspec(vc)).NrOfVoxels = size(c, 1);
            bc.VOI(vspec(vc)).Voxels = c;
        else
            bc.VOI(tspec(vc)) = bc.VOI(vspec(vc));
            bc.VOI(tspec(vc)).Name = ...
                sprintf('%s_sph_%.1f', bc.VOI(tspec(vc)).Name, cmdspec.rsize);
            bc.VOI(tspec(vc)).NrOfVoxels = size(c, 1);
            bc.VOI(tspec(vc)).Voxels = c;
        end
    end

    % set back to content array and return
    xo.C = bc;
    return;
end

% get coords
vistal = (lower(bc.ReferenceSpace(1)) == 't');
crdnum = zeros(1, nsel);
coords = cell(nsel, 1);
maxc = 1;
minc = 16777216;
for vc = 1:nsel
    c = vst(vspec(vc)).Voxels;
    if vistal
        c = 128 - c;
    end
    if ~isempty(bb)
        c = c(:, ( ...
            c(:, 1) >= bb(1, 1) & c(:, 1) <= bb(2, 1) & ...
            c(:, 2) >= bb(1, 2) & c(:, 2) <= bb(2, 2) & ...
            c(:, 3) >= bb(1, 3) & c(:, 3) <= bb(2, 3)));
    end
    c = round(c(:, 1)) + 256 .* round(c(:, 2)) + 65536 .* round(c(:, 3)) + 1;
    if any(c < 1 | c > 16777216)
        error('neuroelf:xff:invalidObject', ...
            'Invalid coordinates in VOI (out of VMR space).');
    end
    if ~isempty(c)
        maxc = max(max(c), maxc);
        minc = min(min(c), minc);
    end
    coords{vc} = c;
    crdnum(vc) = numel(c);
end
minc = minc - 1;
if minc > 0
    for vc = 1:nsel
        coords{vc} = coords{vc} - minc;
    end
end
minc = minc - 1;

% build resulting array and intermittent arrays
res = vst(vspec(1));
res.NrOfVoxels = 0;
res.Voxels = zeros(0, 3);
msk = false(1, maxc);
msk(coords{1}) = true;
vec = false(1, maxc);

% rest depends of command
switch (cmdname)

    % intersect (AND)
    case 'intersect'

        % perform intersection with boolean operators
        for vc = 2:nsel
            vec(:) = false;
            vec(coords{vc}) = true;
            msk = msk & vec;
        end
        res.Voxels = find(msk) + minc;
        res.NrOfVoxels = numel(res.Voxels);
        if subjvoi == cmdspec.cmbsubj
            res.Name = ['ISECT_' vspecname];
        else
            res.Name = [vspecname '_ISECT'];
        end

    % overlap
    case 'overlap'

        % perform intersection with boolean operators
        clear vec;
        msk = uint16(msk);
        for vc = 2:nsel
            msk(coords{vc}) = msk(coords{vc}) + uint16(1);
        end
        res.Voxels = find(msk >= uint16(ceil(nsel * cmdspec.thresh))) + minc;
        res.NrOfVoxels = numel(res.Voxels);
        if subjvoi == cmdspec.cmbsubj
            res.Name = ['OVERLAP_' vspecname];
        else
            res.Name = [vspecname '_OVERLAP'];
        end

    % setdiff (NOT)
    case 'setdiff'

        % vec not needed
        clear vec;

        % without both/all directions
        if ~cmdspec.bothdirs
            for vc = 2:nsel
                msk(coords{vc}) = false;
            end
            res.Voxels = find(msk) + minc;
            res.NrOfVoxels = numel(res.Voxels);
            if subjvoi == cmdspec.cmbsubj
                res.Name = [res.Name '_UNIQUE'];
            else
                res.Name = ['UNIQUE_' res.Name];
            end

        % otherwise perform setdiff with each as a start!
        else

            % grow result
            res(2:nsel) = res;

            % iterate over first part of combination
            for vc = 1:nsel

                % setup array
                msk(:) = false;
                msk(coords{vc}) = true;

                % iterate over others
                for svc = 1:nsel
                    if svc ~= vc
                        msk(coords{svc}) = false;
                    end
                end

                % fill output
                res(vc).Voxels = find(msk) + minc;
                res(vc).NrOfVoxels = numel(res(vc).Voxels);
                if subjvoi == cmdspec.cmbsubj
                    res(vc).Name = [vst(vspec(vc)).Name '_UNIQUE'];
                else
                    res(vc).Name = ['UNIQUE_' vst(vspec(vc)).Name];
                end
            end
        end

    % union (OR)
    case 'union'

        % perform union with boolean operators
        for vc = 2:nsel
            vec(coords{vc}) = true;
            msk = msk | vec;
        end
        res.Voxels = find(msk) + minc;
        res.NrOfVoxels = numel(res.Voxels);
        if subjvoi == cmdspec.cmbsubj
            res.Name = ['UNION_' vspecname];
        else
            res.Name = [vspecname '_UNION'];
        end

end

% recompute coordinate system
for vc = 1:numel(res)
    c = res(vc).Voxels(:);
    c = [mod(c, 256), mod(floor(c ./ 256), 256), floor(c ./ 65536)];
    if vistal
        c = 128 - c(end:-1:1, :);
    end
    res(vc).Voxels = c;
    res(vc).Name = strrep(res(vc).Name, '__', '_');
end

% only output
if cmdspec.nostore
    xo = res;
    return;
end

% add to VOI file
bc.VOI = [bc.VOI(:)', res(:)'];
bc.NrOfVOIs = numel(bc.VOI);
xo.C = bc;

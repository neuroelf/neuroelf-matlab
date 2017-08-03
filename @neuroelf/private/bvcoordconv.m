function c = bvcoordconv(c, ct, bbox)
% bvcoordconv  - converting coordinates for BV
%
% FORMAT:       c = bvcoordconv(c, ct, bbox)
%
% Input fields:
%
%       c           Cx3 coordinates (or Cx1 / 1xC 1-based indices)
%       ct          conversion type, one of
%                   'bvc2tal'  convert BV internal volume coords to TAL
%                   'bvi2tal'  convert BV internal system coords to TAL
%                   'bvs2tal'  convert BV external system coords to TAL
%                   'bvx2tal'  convert BV internal volume indices to TAL
%                   'tal2bvc'  convert TAL to BV internal volume coords
%                   'tal2bvi'  convert TAL to BV internal system coords
%                   'tal2bvs'  convert TAL to BV external system coords
%                   'tal2bvx'  convert TAL to BV internal volume indices
%       bbox        1x1 bounding box struct with (at least) fields
%        .BBox      hardcoded bounding box
%        .DimXYZ    size of volume to sample (needed for indexing)
%        .FCube     framing cube (needed to decide about center)
%        .ResXYZ    resolution of voxels (in BV internal notation)
%
% Output fields:
%
%       c           converted coordinates
%
% For further explanations see 'dbtype bvcoordconv' !

% Version:  v1.1
% Build:    16012313
% Date:     Jan-23 2016, 1:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% Note: the documentation in this file has been influenced by
%       Graham Wideman's very good notes about coordinate systems; see
%
%       http://www.grahamwideman.com/gw/brain/orientation/orientterms.htm
%  and  http://www.grahamwideman.com/gw/brain/analyze/formatdoc.htm
%
% Abstract:
%
% Coordinate systems (or their differences) commonly lead to
% misinterpretations of spatial references and/or wrongly visualized
% data (sliced images). This can only be solved by knowing about the
% "good" interpretation of data as found in the several data formats
%
% BrainVoyager files (different for FMR/STC, VMR/VTC/VDW/GLM/VMP/...)
% SPM/Analyze/NIftI
% AFNI HEAD/BRIK
% DICOM (scanner space)
%
% There are a few different ways to describe the orientation and
% layout of coordinate systems. In this file, the "from"-naming scheme
% is used, which gives the three letters (e.g. RPI for right-to-left
% in posterior-to-anterior in inferior-to-superior order) to describe
% how data is ordered in a chunk of data which then also determines the
% order and orientation of the coordinate system axes. So, each letter
% denotes the "side" with the **LOWER** coordinate indices (so that for
% RPI -X is right, +X is left, etc.). Other guides use the reverse
% nomenclature, so the same system could be described as LAS! Note
% that in this file, the "from"-scheme is used henceforth.
%
% For all datasets, a center-sampling is assumed (that is, a voxel at
% position (0;0;0) with resulution (r;r;r) represents the "average"
% contrast of the volume spanning from (-r/2;-r/2;-r/2) through
% (r/2;r/2;r/2). This is consistent with the idea of interpolation and
% holds true for at least SPM and BrainVoyager coordinate systems.
%
% Explanation of working mechanism:
%
% The TAL coordinate system (LPI) works (all in 1mm resolution)
%
% - 1st dimension:       left (-127) to right    (128)
% - 2nd dimension:  posterior (-127) to anterior (128)
% - 3rd dimension:   inferior (-127) to superior (128)
%
% with a center coordinate at (0;0;0)
%
% The BV system coordinates (RAS, external, as presented to the user!) have
% the same general axes logic, but with reversed directions. To
% slighlty complicate matters, the interpretation (center) depends on
% the framing cube of the data (in older versions fixed to 256):
%
% - 1st dimension:    right (0) to left      (FCube - 1)
% - 2nd dimension: anterior (0) to posterior (FCube - 1)
% - 3rd dimension: superior (0) to inferior  (FCube - 1)
%
% with the center at (FCube / 2;FCube / 2;FCube / 2), usually at
% (128;128;128) which is still fixed for non VMR/VMP files!
%
% Please note: in the voxel-based (VMR) space, any rotations to
% coordinates are transformed around the center coordinate - 0.5
% (so the virtual point at the center of the data!), assuming that
% the voxel center is at (0.5;0.5;0.5) !! In contrast, for the surface
% space (SRF), the center coordinate *is* at (0;0;0)!
%
% The *second* BV system coordinates (ASR, internal ordering only,
% rarely shown to users!) on top permute the axes order as [2, 3, 1]:
%
% - 1st dimension:  anterior (0) to posterior (FCube - 1)
% - 2nd dimension:  superior (0) to inferior  (FCube - 1)
% - 3rd dimension:     right (0) to left      (FCube - 1)
%
% This leads to a sagittal slicing, and all three standard sliced
% images in BrainVoyager have a growing indexing from (screen-space!)
% left-to-right and top-to-bottom (SAG: X,Y; COR: Z,Y; TRA: Z,X)
%
% The BV internal volume coordinates (0-based) then represent the
% indexing information into a 3-D array, so that from the BV internal
% system coordinates the Offset must be subtracted and the result
% devided by the resolution (finally adding 1 for ML's indexing).
%
% The BV volume indexing, eventually, represents this information
% in a simple number using sub2ind(DataDimXYZ, c).
%
% Note: this notation differs from the BrainVoyager QX coordinates in
%       one (more or less crucial) point:
%
% the toolbox always represents system coordinates in a 1x1x1 ISO-voxel
% resolution, while BrainVoyager does not do so for hi-res VMRs (all other
% filetypes assume an underlying 1x1x1 resolution of the hosting object!)
%
% NB: Up to version 1.10.x, in VMR files only BrainVoyager also allows
%     Neurological convention of the data, leading to a
%     LAS/ASL external/internal coordinate system. All other files are
%     unaware of the actual orientation! Visualizing any VTC/GLM/VMP/SRF
%     file hence depends on loading the (a) correctly oriented VMR.
%     Due to this fact, it is *strongly suggested* not to use the
%     Neurological convention to allow simple data exchange between
%     users. Newer version (2.x.x) also have a corresponding flag in
%     the other voxel-based filetypes, but the recommendation remains...
%
% Other formats:
%
% SPM/Analyze:
%
% The Analyze format has (like BrainVoyager file formats) grown over time.
% One of the changes is the "flipping" that has to be guessed for older
% files that do not contain this information in the Analyze file header.
% SPM assumes (for files where this information is not available) that
% data comes in a flipped orientation (see spm_defaults.m).
%
% The general logic for data storage follows a (from) RPI scheme, while the
% display is in "neurological convention" (so, in the GUI, left IS left!)
% hence, data that is stored in the radiological convention must be
% accessed in a "flipped" manner. This is denoted in the Analyze header
% either
%
% - by using a negatively oriented direction vector in the Quaternion
%   and/or AffineTransX matrix (header with nii/ni1/n+1 magic token)
% - by setting the first (otherwise unused) PixSpacing to a negative
%   value (intermediate solution)
% - by using the corresponding DataHist.Orientation value (3) for old
%   Analyze format (header without nii/ni1/n+1 magic token)
%
% If neither of these conditions is met, the toolbox's default is yet
% to assume flipped data ! To disable this behavior, use
%
% x = xff();
% x.Config('hdr', 'assumeflipped', false);

% argument check
if nargin < 2 || ...
   ~isa(c, 'double') || ...
    ndims(c) > 2 || ...
   ~ischar(ct) || ...
    numel(ct) ~= 7 || ...
    size(ct, 2) ~= 7 || ...
   ~any(strcmpi(ct(1:3), {'bvc', 'bvi', 'bvs', 'bvx', 'tal'})) || ...
   ~any(strcmpi(ct(5:7), {'bvc', 'bvi', 'bvs', 'bvx', 'tal'})) || ...
    ct(4) ~= '2' || ...
   (~strcmpi(ct(1:3), 'bvx') && ...
    ~any(size(c) == 3))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
    isempty(bbox)
    bbox = struct( ...
        'BBox',   [0, 0, 0; 255, 255, 255], ...
        'DimXYZ', [256, 256, 256], ...
        'FCube',  [256, 256, 256], ...
        'ResXYZ', [1, 1, 1]);
elseif isa(bbox, 'double') && ...
    isequal(size(bbox), [2, 3]) && ...
   ~any(isinf(bbox(:)) | isnan(bbox(:)) | bbox(:) < 0)
    bbox = struct( ...
        'BBox',   bbox, ...
        'DimXYZ', 1 + diff(bbox), ...
        'FCube',  [256, 256, 256], ...
        'ResXYZ', [1, 1, 1]);
elseif ~isstruct(bbox) || ...
    numel(bbox) ~= 1 || ...
   ~isfield(bbox, 'BBox') || ...
   ~isfield(bbox, 'DimXYZ') || ...
   ~isfield(bbox, 'FCube') || ...
   ~isfield(bbox, 'ResXYZ')
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid bbox argument.' ...
    );
end
ct = lower(ct);
cf = ct(1:3);
ci = ct(5:7);
if size(c, 1) == 3 && ...
    size(c, 2) ~= 3
    c = c';
end
bbx = bbox.BBox(1, 1:3);
xyz = bbox.DimXYZ;
fcb = 0.5 .* bbox.FCube;
res = bbox.ResXYZ;

% which processing direction
switch (ct)
    case {'bvc2tal', 'bvi2tal', 'bvs2tal', 'bvx2tal', ...
          'bvc2bvs', 'bvi2bvs', 'bvx2bvs', 'bvc2bvi', 'bvx2bvi', 'bvx2bvc'}
        total = true;
    case {'tal2bvc', 'tal2bvi', 'tal2bvs', 'tal2bvx', ...
          'bvs2bvc', 'bvs2bvi', 'bvs2bvx', 'bvi2bvc', 'bvi2bvx', 'bvc2bvx'}
        total = false;
end

% for empty coordinates, return transformation matrix
if isempty(c)

    % set to unity
    c = eye(4);

    % in TAL direction
    if total

        % starting at indexing
        if strcmp(cf, 'bvx')
            cf = 'bvc';
        end

        % coordinates in internal volume coords but further conversion
        if strcmp(cf, 'bvc') && ...
           ~strcmp(ci, 'bvc')

            % convert to internal system coords
            c = [res(1), 0, 0, bbx(1) - res(1); ...
                 0, res(2), 0, bbx(2) - res(2); ...
                 0, 0, res(3), bbx(3) - res(3); ...
                 0, 0, 0, 1] * c;
            cf = 'bvi';
        end

        % coordinates in internal system coords but further conversion
        if strcmp(cf, 'bvi') && ...
           ~strcmp(ci, 'bvi')

            % convert to external system coords
            c = [0, 0, 1, 0; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 0, 1] * c;
            cf = 'bvs';
        end

        % coordinates in external system coords but further conversion
        if strcmp(cf, 'bvs') && ...
           ~strcmp(ci, 'bvs')

            % convert to TAL space
            c = [-1,  0,  0, fcb(3); ...
                  0, -1,  0, fcb(1); ...
                  0,  0, -1, fcb(2); ...
                  0,  0,  0,    1  ] * c;
        end

    % coming from TAL direction
    else

        % starting at TAL
        if strcmp(cf, 'tal')

            % convert to external system coords
            c = [-1,  0,  0, fcb(3); ...
                  0, -1,  0, fcb(1); ...
                  0,  0, -1, fcb(2); ...
                  0,  0,  0,    1  ] * c;
            cf = 'bvs';
        end

        % coordinates in external system coords but further conversion
        if strcmp(cf, 'bvs') && ...
           ~strcmp(ci, 'bvs')

            % convert to internal system coords
            c = [0, 1, 0, 0; 0, 0, 1, 0; 1, 0, 0, 0; 0, 0, 0, 1] * c;
            cf = 'bvi';
        end

        % coordinates in internal system coords but further conversion
        if strcmp(cf, 'bvi') && ...
           ~strcmp(ci, 'bvi')

            % convert to internal volume coords
            ires = 1 ./ res;
            ibbx = 1 - bbx ./ res;
            c = [ires(1), 0, 0, ibbx(1); ...
                 0, ires(2), 0, ibbx(2); ...
                 0, 0, ires(3), ibbx(3); ...
                 0, 0, 0, 1] * c;
        end
    end

    % return early
    return;
end

% in TAL direction
if total

    % starting at indexing
    if strcmp(cf, 'bvx')

        % convert to internal volume
        [cc{1:3}] = ind2sub(xyz, c(:));
        c = [cc{1}, cc{2}, cc{3}];
        cf = 'bvc';
    end

    % coordinates in internal volume coords but further conversion
    if strcmp(cf, 'bvc') && ...
       ~strcmp(ci, 'bvc')

        % convert to internal system coords
        c = [res(1) * (c(:, 1) - 1) + bbx(1), ...
             res(2) * (c(:, 2) - 1) + bbx(2), ...
             res(3) * (c(:, 3) - 1) + bbx(3)];
        cf = 'bvi';
    end

    % coordinates in internal system coords but further conversion
    if strcmp(cf, 'bvi') && ...
       ~strcmp(ci, 'bvi')

        % convert to external system coords
        c = [c(:, 3), c(:, 1), c(:, 2)];
        cf = 'bvs';
    end

    % coordinates in external system coords but further conversion
    if strcmp(cf, 'bvs') && ...
       ~strcmp(ci, 'bvs')

        % convert to TAL space
        if all(fcb == fcb(1))
            c = fcb(1) - c;
        else
            c = fcb(ones(1, size(c, 1)), [3, 1, 2]) - c;
        end
    end

% coming from TAL direction
else

    % starting at TAL
    if strcmp(cf, 'tal')

        % convert to external system coords
        if all(fcb == fcb(1))
            c = fcb(1) - c;
        else
            c = fcb(ones(1, size(c, 1)), [3, 1, 2]) - c;
        end
        cf = 'bvs';
    end

    % coordinates in external system coords but further conversion
    if strcmp(cf, 'bvs') && ...
       ~strcmp(ci, 'bvs')

        % convert to internal system coords
        c = [c(:, 2), c(:, 3), c(:, 1)];
        cf = 'bvi';
    end

    % coordinates in internal system coords but further conversion
    if strcmp(cf, 'bvi') && ...
       ~strcmp(ci, 'bvi')

        % convert to internal volume coords
        c = [(c(:, 1) - bbx(1)) ./ res(1) + 1, ...
             (c(:, 2) - bbx(2)) ./ res(2) + 1, ...
             (c(:, 3) - bbx(3)) ./ res(3) + 1];
        cf = 'bvc';
    end

    % coordinates in internal volume coords but further conversion
    if strcmp(cf, 'bvc') && ...
       ~strcmp(ci, 'bvc')

        % convert to indexes
        c = round(c);
        gc = (all(c > 0, 2) & ...
            c(:, 1) <= xyz(1) & c(:, 2) <= xyz(2) & c(:, 3) <= xyz(3));
        c(gc, 1) = sub2ind(xyz, c(gc, 1), c(gc, 2), c(gc, 3));
        c = c(:, 1);
        c(~gc) = NaN;
    end
end

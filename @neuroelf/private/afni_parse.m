function fields = afni_parse(info, filename)
% afni_parse  - internal parser for AFNI HEADer file structures
%
% FORMAT:       fields = afni_parse(Info, filename);
%
% Input fields:
%
%       Info        Nx1 struct with fields
%        .Name      name of field
%        .Type      type of field content
%        .Size      number of elements
%        .Data      actual data of field
%       filename    filename of HEADer file (to get to the BRIK file)
%
% Output field:
%
%       fields      1x1 struct with mandatory and optional settings
%                   extracted from Info (as well as a copy of Info!)
%                   as well as a Brick field that allows access to the
%                   binary data in the BRIK file

% Note: the information incorporated here-in comes from
%       http://afni.nimh.nih.gov/pub/dist/src/README.attributes

% Version:  v0.9b
% Build:    11051100
% Date:     Apr-08 2011, 10:18 PM EST
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

% persistent memory
persistent afni_pm;
if isempty(afni_pm)
    afni_pm = struct( ...
        'TypeStrings', {{'3DIM_HEAD_ANAT', '3DIM_HEAD_FUNC', ...
                '3DIM_GEN_ANAT', '3DIM_GEN_FUNC'}});
end

% argument check
if nargin ~= 2 || ...
   ~isstruct(info) || ...
    numel(fieldnames(info)) ~= 4 || ...
   ~all(strcmp(fieldnames(info), {'Name'; 'Type'; 'Size'; 'Data'})) || ...
    isempty(info) || ...
   ~ischar(filename)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid call to function.' ...
    );
end

% initiate output
fields = struct(...
    'Info', info, ...
    'Endianness', 'ieee-le', ...
    'IDCode', sprintf('%c', 65 + floor(25 * rand(1, 15))), ...
    'IDDate', datestr(now), ...
    'NrOfVolumes', 1, ...
    'NrOfTimePoints', 0, ...
    'DataDimensions', [0, 0, 0], ...
    'AxesOrientation', [1, 3, 5], ...
    'Origin', [0, 0, 0], ...
    'Resolution', [1, 1, 1], ...
    'TimeUnits', 'msec', ...
    'TypeOfVolumes', '3DIM_GEN_FUNC', ...
    'SubtypeOfVolumes', 0, ...
    'ScalingFactor', 0, ...
    'SceneView', 'Original', ...
    'SceneFuncType', 'FIM', ...
    'SceneVolumeType', 3, ...
    'NrOfSliceTimes', 0, ...
    'TimeOrigin', 0, ...
    'TimeStep', 1, ...
    'TimeAcquisition', 0, ...
    'TimeZOffset', 0, ...
    'TimeZStep', 1, ...
    'TimeSliceOffsets', [], ...
    'WarpParent', '', ...
    'WarpType', -1, ...
    'WarpValues', emptystruct({ ...
        'ForwardMat', 'BackwardMat', 'ForwardVec', 'BackwardVec', 'From', 'To'}), ...
    'MarkersType', 1, ...
    'Markers', emptystruct({'Label', 'Coordinate', 'Help'}), ...
    'Brick', struct( ...
        'Label', '<none>', 'DataType', 1, 'FuncType', 0, 'FuncParams', [], ...
        'MinMaxValue', [0, 225], 'ScalingFactor', 0, 'Data', [], 'DataCT', logical([])));

% set intermediate data fields
nvals = 1;
bstats = [0, 225];
btypes = 1;
bffacs = 0;
blabls = {'<none>'};
bstaux = [0, 0, 0];
b2saux = [0, 0];
wptype = 0;

% split filename
[fp, fn] = fileparts(filename);
if isempty(fp)
    fp = pwd;
end

% now parse each Info element
for ic = 1:numel(info)

    % get element handle
    ie = info(ic);
    ien = upper(ie.Name);
    ied = ie.Data;

    % switch over field name
    switch (ien)

        % brick float factors
        case {'BRICK_FLOAT_FACS'}
            bffacs = ied;

        % brick labels
        case {'BRICK_LABS'}
            blabls = splittocellc(ied, '~');

        % brick stats aux field
        case {'BRICK_STATAUX'}
            bstaux = ied;

        % brick stats
        case {'BRICK_STATS'}
            bstats = ied;

        % brick types
        case {'BRICK_TYPES'}
            btypes = ied;

        % byte order
        case {'BYTEORDER_STRING'}
            if numel(ied) < 3
                continue;
            end
            switch (lower(ied(1:3)))
                case {'lsb'}
                    fields.Endianness = 'ieee-le';
                case {'msb'}
                    fields.Endianness = 'ieee-be';
                otherwise
                    warning( ...
                        'neuroelf:AFNIWarning', ...
                        'Invalid Endianness specified: %s.', ...
                        ied ...
                    );
            end

        % dimensions
        case {'DATASET_DIMENSIONS'}
            if numel(ied) < 3
                continue;
            end
            if any(ied == 0)
                ied(findfirst(ied == 0):end) = [];
            end
            fields.DataDimensions = ied;

        % rank
        case {'DATASET_RANK'}
            if numel(ied) < 2 || ...
                ied(1) ~= 3
                continue;
            end
            nvals = round(ied(2));
            fields.NrOfVolumes = nvals;
            fields.Brick = fields.Brick(ones(1, nvals));

        % resolution
        case {'DELTA'}
            if numel(ied) < 3
                continue;
            end
            fields.Resolution = ied;

        % ID code
        case {'IDCODE_STRING'}
            if ied(end) == '~'
                ied(end) = [];
            end
            fields.IDCode = ied;

        % ID date
        case {'IDCODE_DATE'}
            if ied(end) == '~'
                ied(end) = [];
            end
            fields.IDDate = ied;

        % orientation specification
        case {'ORIENT_SPECIFIC'}
            if numel(ied) < 3
                continue;
            end
            fields.AxesOrientation = ied;

        % origin
        case {'ORIGIN'}
            if numel(ied) < 3
                continue;
            end
            fields.Origin = ied;

        % scene data
        case {'SCENE_DATA'}
            if numel(ied) < 3
                continue;
            end
            switch (ied(1))
                case {0}
                    fields.SceneView = 'Original';
                case {1}
                    fields.SceneView = 'ACPC';
                case {2}
                    fields.SceneView = 'Talairach';
            end
            fields.SubtypeOfVolumes = ied(2);
            if ied(3) < 0 || ...
                ied(3) > 3 || ...
               ~strcmpi(afni_pm.TypeStrings{round(ied(3)) + 1}, fields.TypeOfVolumes)
                warning( ...
                    'neuroelf:AFNIWarning', ...
                    'Invalid SCENE_DATA[2] field.' ...
                );
            else
                field.TypeOfVolumes = afni_pm.TypeStrings{round(ied(3)) + 1};
            end

        % stat aux field (brick #2)
        case {'STAT_AUX'}
            if numel(ied) < 2
                continue;
            end
            b2saux = ied;

        % time axis dims
        case {'TAXIS_NUMS'}
            if numel(ied) < 3
                continue;
            end
            fields.NrOfTimePoints = ied(1);
            fields.NrOfSliceTimes = ied(2);
            switch (ied(3))
                case {77001}
                    fields.TimeUnits = 'msec';
                case {77002}
                    fields.TimeUnits = 'msec';
                case {77003}
                    fields.TimeUnits = 'Hz';
            end

        % time axis settings
        case {'TAXIS_FLOATS'}
            if numel(ied) < 5
                continue;
            end
            fields.TimeOrigin = ied(1);
            fields.TimeStep = ied(2);
            fields.TimeAcquisition = ied(3);
            fields.TimeZOffset = ied(4);
            fields.TimeZStep = ied(5);

        % (slice) time axis offsets
        case {'TAXIS_OFFSETS'}
            fields.TimeSliceOffsets = ied;

        % type string
        case {'TYPESTRING'}
            if ied(end) == '~'
                ied(end) = [];
            end
            if ~any(strcmp(ied, afni_pm.TypeStrings))
                warning( ...
                    'neuroelf:AFNIError', ...
                    'Invalid TypeString: %s.', ...
                    ied ...
                );
            end
            fields.TypeOfVolumes = ied;

        % warping data
        case {'WARP_DATA'}

            % depending on warp type
            switch (wptype)
                case {0}
                    fields.WarpValues(1).ForwardMat = [];
                    fields.WarpValues = fields.WarpValues(1);
                    if numel(ied) == 30
                        fields.WarpValues.ForwardMat = ...
                            reshape(ied(1:9), 3, 3)';
                        fields.WarpValues.BackwardMat = ...
                            reshape(ied(10:18), 3, 3)';
                        fields.WarpValues.ForwardVec = lsqueeze(ied(19:21));
                        fields.WarpValues.BackwardVec = lsqueeze(ied(22:24));
                        fields.WarpValues.From = lsqueeze(ied(25:27));
                        fields.WarpValues.To = lsqueeze(ied(28:30));
                    else
                        fields.WarpValues.From = ied(:)';
                        warning( ...
                            'neuroelf:AFNIError', ...
                            'Invalid number of warping parameters.' ...
                        );
                    end
                case {1}
                    fields.WarpValues(12).ForwardMat = [];
                    fields.WarpValues = fields.WarpValues(1:12);
                    if numel(ied) == 360
                        for wf = 1:12
                            wfi = (wf - 1) * 30;
                            fields.WarpValues(wf).ForwardMat = ...
                                reshape(ied(wfi+1:wfi+9), 3, 3)';
                            fields.WarpValues(wf).BackwardMat = ...
                                reshape(ied(wfi+10:wfi+18), 3, 3)';
                            fields.WarpValues(wf).ForwardVec = ...
                                lsqueeze(ied(wfi+19:wfi+21));
                            fields.WarpValues(wf).BackwardVec = ...
                                lsqueeze(ied(wfi+22:wfi+24));
                            fields.WarpValues(wf).From = ...
                                lsqueeze(ied(wfi+25:wfi+27));
                            fields.WarpValues(wf).To = ...
                                lsqueeze(ied(wfi+28:wfi+30));
                        end
                    else
                        fields.WarpValues = ied(:)';
                        warning( ...
                            'neuroelf:AFNIError', ...
                            'Invalid number of warping parameters.' ...
                        );
                    end
                otherwise
                    fields.WarpValues = ied;
            end

        % warping type
        case {'WARP_TYPE'}
            if numel(ied) == 8 && ...
               ~isinf(ied(1)) && ...
               ~isnan(ied(1)) && ...
                any(ied(1) == [0, 1])
                wptype = ied(1);
                fields.WarpType = wptype;
            else
                warning( ...
                    'neuroelf:AFNIError', ...
                    'Invalid WARP_TYPE flag.' ...
                );
            end

        % warped dataset parent name
        case {'WARP_PARENTNAME'}
            if ied(end) == '~'
                ied(end) = [];
            end
            fields.WarpParent = ied;
            if ~isempty(ied)
                [wppath, wpname] = fileparts(ied);
                if isempty(wppath)
                    wppath = fileparts(filename);
                end
                if exist([fp filesep fn '.BRIK'], 'file') ~= 2 && ...
                    exist([fp filesep fn '.brik'], 'file') ~= 2 && ...
                    exist([fp filesep fn '.BRIK.gz'], 'file') ~= 2 && ...
                    exist([fp filesep fn '.brik.gz'], 'file') ~= 2 && ...
                   (exist([wppath filesep wpname '.BRIK'], 'file') == 2 || ...
                    exist([wppath filesep wpname '.brik'], 'file') == 2 || ...
                    exist([wppath filesep wpname '.BRIK.gz'], 'file') == 2 || ...
                    exist([wppath filesep wpname '.brik.gz'], 'file') == 2)
                    fp = wppath;
                    fn = wpname;
                end
            end
    end
end

% some sanity checks
if numel(bffacs) < nvals
    bffacs(end+1:nvals) = 0;
end
if numel(blabls) < nvals
    blabls(end+1:nvals) = {'<none>'};
end
if numel(bstats) < (2 * nvals)
    bstats(end+1:2*nvals) = 0;
end
if numel(btypes) < nvals
    btypes(end+1:nvals) = 1;
end
bstaux(end+1) = -1;

% does the BRIK file exist
briktmp = false;
if exist([fp filesep fn '.BRIK'], 'file') == 2
    brik = [fp filesep fn '.BRIK'];
elseif exist([fp filesep fn '.brik'], 'file') == 2
    brik = [fp filesep fn '.brik'];
elseif exist([fp filesep fn '.BRIK.gz'], 'file') == 2
    try
        tmpfile = sprintf('%s%s%06x_%s.BRIK.gz', ...
            tempdir, filesep, round(2^24 * rand(1, 1)), fn);
        cpfile([fp filesep fn '.BRIK.gz'], tmpfile);
        ne_gzip('-d', tmpfile);
        tmpfilu = tmpfile(1:end-3);
        if exist(tmpfilu, 'file') == 2
            brik = tmpfilu;
            briktmp = true;
        else
            if exist(tmpfile, 'file') == 2
                try
                    delete(tmpfile);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
            if exist(tmpfilu, 'file') == 2
                try
                    delete(tmpfilu);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
            error('UngzipError');
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        warning( ...
            'neuroelf:GzipError', ...
            ['Error decompressing BRIK.gz file. ' ...
             'The system binary ''gzip'' must be available and working.'] ...
        );
        fields.Brick(:) = [];
        return;
    end
else
    warning( ...
        'neuroelf:FileNotFound', ...
        'BRIK file not found (and transformed access not implemented).' ...
    );
    fields.Brick(:) = [];
    return;
end

% try to create BRIK access object
[dt, fac, siz] = brik_datatype(btypes(1));
try
    tio = transio(brik, fields.Endianness, dt, 0, [fac, fields.DataDimensions]);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    warning( ...
        'neuroelf:transioError', ...
        'Unable to create transio access...' ...
    );
    fields.Brick(:) = [];
    return;
end

% put object into Brick field
stap = 1;
for vc = 1:nvals
    if vc > 1
        tios = struct(tio);
        tios.IOBuffer = {{}, []};
        tios.IOOffset = tios.IOOffset + siz * prod(tios.DataDims);
        [dt, fac, siz] = brik_datatype(btypes(vc));
        tios.DataType = dt;
        tios.TypeSize = siz;
        tios.DataDims = [fac, fields.DataDimensions];
        tio = transio(0, 'makeobject', tios);
    end
    fields.Brick(vc).Label = blabls{vc};
    fields.Brick(vc).DataType = btypes(vc);
    if bstaux(stap) == (vc - 1)
        fields.Brick(vc).FuncType = bstaux(stap + 1);
        if bstaux(stap + 2) > 0
            fields.Brick(vc).FuncParams = bstaux(stap+3:stap+2+bstaux(stap + 2));
            stap = stap + 3 + bstaux(stap + 2);
        else
            stap = stap + 3;
        end
    end
    fields.Brick(vc).MinMaxValue = bstats(2*vc-1:2*vc);
    fields.Brick(vc).ScalingFactor = bffacs(vc);
    if briktmp
        try
            fields.Brick(vc).Data = resolve(tio);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    else
        fields.Brick(vc).Data = tio;
    end
end

% remove temporary object
if briktmp
    try
        delete(brik);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% brick 2 special info
if any(b2saux > 0)
    fields.Brick(2).FuncType = b2saux(1);
    fields.Brick(2).FuncParams = b2saux(3:end);
end



% subfunction to get matlab class from BRIK type
function [mcn, fac, siz] = brik_datatype(t)
fac = [];
switch (t)
    case {0}
        mcn = 'uint8';
        siz = 1;
    case {1}
        mcn = 'int16';
        siz = 2;
    case {2}
        mcn = 'int32';
        siz = 4;
    case {3}
        mcn = 'single';
        siz = 4;
    case {4}
        mcn = 'double';
        siz = 8;
    case {5}
        mcn = 'single';
        fac = 2;
        siz = 4;
    case {6}
        mcn = 'uint8';
        fac = 3;
        siz = 1;
end
% end of function brik_datatype

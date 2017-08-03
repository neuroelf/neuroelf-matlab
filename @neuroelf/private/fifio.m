function [varargout] = fifio(varargin)
% fifio  - read FIF(F) files (NeuroMag)
%
% FORMAT:       fifobject = fifio(filename)
%
% Input fields:
%
%       filename    filename of FIF(F) file to read
%
% Output fields:
%
%       fifobject   xff object
%
% See also xff

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-07 2011, 4:23 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% persistent variables
persistent fif_kinds fif_itypes fif_types;
if isempty(fif_kinds) || isempty(fif_itypes) || isempty(fif_types)

    % structure for blockkinds, datakinds, and matrix-kinds, etc.
    fif_kinds = struct;

    % block kinds
    fif_kinds.B0064 = 'Measurement';
    fif_kinds.B0065 = 'MetaInformation';
    fif_kinds.B0066 = 'RawData';
    fif_kinds.B0067 = 'CookedData';
    fif_kinds.B0068 = 'EvokedData';
    fif_kinds.B0069 = 'AspectData';
    fif_kinds.B006A = 'SubjectInfomation';
    fif_kinds.B006B = 'ISOTrak';
    fif_kinds.B006D = 'HPIResult';
    fif_kinds.B0073 = 'Comments';
    fif_kinds.B0139 = 'ProjectInformation';
    fif_kinds.B013A = 'ProjectItem';

    % coil types
    fif_kinds.C07D0 = {'PointMagnetometer'};
    fif_kinds.C0BC4 = {'VectorviewPlanarGradiometer', 'T1', 0.02639, 0.0168};
    fif_kinds.C0BC5 = {'VectorviewPlanarGradiometer', 'T2', 0.02639, 0.0168};
    fif_kinds.C0BCE = {'VectorviewMagnetometer', 'T1', 0.0258};
    fif_kinds.C0BCF = {'VectorviewMagnetometer', 'T2', 0.0258};
    fif_kinds.C0BD0 = {'VectorviewMagnetometer', 'T3', 0.021};
    fif_kinds.C0FA1 = {'MagnesWH2500Magnetometer', '', 0.0115};
    fif_kinds.C0FA2 = {'MagnesWH3600Gradiometer', '', 0.018, 0.050};
    fif_kinds.C1389 = {'CTFAxialGradiometer', '', 0.018, 0.050};

    % frame types
    fif_kinds.F0000 = 'Unknown';
    fif_kinds.F0001 = 'Device';
    fif_kinds.F0002 = 'ISOTrak';
    fif_kinds.F0003 = 'HPI';
    fif_kinds.F0004 = 'Head';
    fif_kinds.F0005 = 'MRI';
    fif_kinds.F0006 = 'MRISlice';
    fif_kinds.F0007 = 'MRIDisplay';
    fif_kinds.F0008 = 'DICOMDevice';
    fif_kinds.F0009 = 'ImagingDevice';

    % datakinds
    fif_kinds.K0064 = 'FileID';
    fif_kinds.K0065 = 'DirPointer';
    fif_kinds.K0066 = 'Directory';
    fif_kinds.K0067 = 'BlockID';
    fif_kinds.K0068 = 'BlockStart';
    fif_kinds.K0069 = 'BlockEnd';
    fif_kinds.K006A = 'FreeList';
    fif_kinds.K006C = 'NoOPerationElement';
    fif_kinds.K006D = 'ParentFileID';
    fif_kinds.K006E = 'ParentBlockID';
    fif_kinds.K0071 = 'HPIComment';
    fif_kinds.K0096 = 'Settings';
    fif_kinds.K00C8 = 'NrOfChannels';
    fif_kinds.K00C9 = 'SamplingFrequency';
    fif_kinds.K00CA = 'SomeValue00CA'; % 4          ???
    fif_kinds.K00CB = 'ChannelInfo';
    fif_kinds.K00CC = 'UTCuSec';
    fif_kinds.K00CE = 'Comment';
    fif_kinds.K00CF = 'Nave';          % ???
    fif_kinds.K00D0 = 'FirstSample';
    fif_kinds.K00D1 = 'LastSample';
    fif_kinds.K00D2 = 'AspectKind';
    fif_kinds.K00D4 = 'SubjectName';   % dubious
    fif_kinds.K00D5 = 'DigPoint';
    fif_kinds.K00D7 = 'SomeList00D7';  % +-1e-12    ???
    fif_kinds.K00D8 = 'SomeValue00D8'; % 4          ???
    fif_kinds.K00DB = 'LowPass';
    fif_kinds.K00DC = 'SomeList00DC';  % -1 terminated, in averaged data (Hanna's data)
    fif_kinds.K00DE = 'TransformationMatrix';
    fif_kinds.K00DF = 'HighPass';
    fif_kinds.K00E2 = 'SomeList00E2';  % 0. (...)   ???
    fif_kinds.K00F0 = 'SomeMatrix00F0';% 3x4 float  ??? Hanna's data
    fif_kinds.K00F1 = 'SomeList00F1';  % 0. (...)   ???
    fif_kinds.K00F2 = 'SomeValue00F2'; % 0          ???
    fif_kinds.K00F3 = 'SomeValue00F3'; % 0. (...)   ???
    fif_kinds.K00F4 = 'SomeValue00F4'; % 0.00 (...) ???
    fif_kinds.K00F5 = 'SomeValue00F5'; % 1          ???
    fif_kinds.K00F6 = 'SomeList00F6';  % 1, 2, 3, 4 ???
    fif_kinds.K00F7 = 'SomeList00F7';  % 2, 3, 1, 4 ???
    fif_kinds.K0107 = 'SomeValue0107'; % 4          ???
    fif_kinds.K0108 = 'SomeValue0108'; % 2          ???
    fif_kinds.K0109 = 'SomeList0109';  % [0,0,0.04] ???
    fif_kinds.K010A = 'NrOfStimChannels'; %         ???
    fif_kinds.K010B = 'NrOfExtraChannels'; %        ???
    fif_kinds.K010C = 'NrOfSensorChannels'; %       ???
    fif_kinds.K010D = 'SomeList010D';  % [0 | 1]    ???
    fif_kinds.K010E = 'SomeMatrix010E';%            ???
    fif_kinds.K010F = 'SomeMatrix010F';%            ???
    fif_kinds.K0110 = 'SomeValue0110'; % 0.9        ???
    fif_kinds.K012C = 'DataBuffer';
    fif_kinds.K012D = 'DataSkip';
    fif_kinds.K012E = 'Epoch';
    fif_kinds.K012F = 'DataSkipInSamples';
    fif_kinds.K0190 = 'SomeValue0190'; % 1310       ???
    fif_kinds.K0191 = 'InvestigatorName'; % dubious
    fif_kinds.K0192 = 'InvestigatorMiddleName'; % very dubious, Hanna's data
    fif_kinds.K0193 = 'InvestigatorGivenName'; % dubious
    fif_kinds.K0194 = 'SomeValue0194'; % 2444011    ???
    fif_kinds.K0195 = 'SomeValue0195'; % 2          ???
    fif_kinds.K0196 = 'SomeValue0196'; % 1          ???
    fif_kinds.K0198 = 'SomeValue0198'; % 1.67       ???
    fif_kinds.K0199 = 'SomeValue0199'; % "Ambidextrous", Hanna's data
    fif_kinds.K01F4 = 'SomeValue01F4'; % 2.12       ???
    fif_kinds.K01F5 = 'ExperimentName'; %           ???
    fif_kinds.K0258 = 'StimulusChannels';
    fif_kinds.K0259 = 'EventList';     % dubious
    fif_kinds.K0320 = 'SignalSquareMatrix'; %       ???
    fif_kinds.K0D53 = 'NrOfMEGChannels'; % dubious, Hanna's data
    fif_kinds.K0D54 = 'SomeValue0D54'; % 0.0        ??? Hanna's data
    fif_kinds.K0D56 = 'SomeValue0D56'; %            ??? Hanna's data
    fif_kinds.K0D57 = 'SomeValue0D57'; %            ??? Hanna's data
    fif_kinds.K0D59 = 'MEGChannelIDs';
    fif_kinds.KRAWD = 'RawData';       % used to read raw data

    % matrix encoding schemes
    fif_kinds.M0000 = 'Normal';
    fif_kinds.M4000 = 'Dense';
    fif_kinds.M4010 = 'Sparse';

    % point kinds
    fif_kinds.P0001 = 'HeadFiducial';
    fif_kinds.P0002 = 'HPIPoint';
    fif_kinds.P0004 = 'ExtraPoint';

    % build lookup
    lup = fieldnames(fif_kinds);
    blup = struct;
    klup = struct;
    mlup = struct;
    for lc = 1:numel(lup)
        lpk = lup{lc};
        if numel(lpk) < 10 && ...
            ~isempty(regexpi(lpk(2:end), '^[0-9a-f]+$'))
            lkv = hex2dec(lpk(2:end));
        else
            lkv = -1;
        end
        switch (lpk(1))
            case {'B'}
                blup.(fif_kinds.(lpk)) = lkv;
            case {'K'}
                klup.(fif_kinds.(lpk)) = lkv;
            case {'M'}
                mlup.(fif_kinds.(lpk)) = lkv;
        end
    end
    fif_kinds.Lookup = struct( ...
        'Block', blup, ...
        'Datakind', klup, ...
        'MatrixEnconding', mlup);

    % structure for internal datatypes
    fif_itypes = struct;
    fif_itypes.double = 8;
    fif_itypes.int8   = 1;
    fif_itypes.int16  = 2;
    fif_itypes.int32  = 4;
    fif_itypes.single = 4;
    fif_itypes.uint8  = 1;
    fif_itypes.uint16 = 2;
    fif_itypes.uint32 = 4;

    % structure for datatypes
    fif_types = cell2struct({ ...
         'uint8',  'uint8',  1, {}; ...
         'int16',  'int16',  2, {}; ...
         'int32',  'int32',  4, {}; ...
        'single', 'single',  4, {}; ...
        'double', 'double',  8, {}; ...
        'uint32', 'uint32',  4, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
        'string',  'uint8',  1, {@uint8_to_string, @string_to_uint8}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'raw',  'uint8',  1, {}; ...  % created to support raw reading
        'pack16',  'int16',  2, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
        'oldpak',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
        'cpxsng', 'single',  8, {@single_to_complex, @complex_to_single}; ...
        'cpxdbl', 'double', 16, {@double_to_complex, @complex_to_double}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
           'n/a',  'uint8',  1, {}; ...
        'chinfo',  'int32', 96, {@int32_to_chinfo, @chinfo_to_int32}; ...
            'id',  'uint8', 20, {@uint8_to_id, @id_to_uint8}; ...
           'dir', 'uint32', 16, {@uint32_to_dir, @dir_to_uint32}; ...
         'point',  'int32', 20, {@int32_to_point, @point_to_int32}; ...
           'n/a',  'uint8',  1, {}; ...
         'trans',  'int32',104, {@int32_to_trans, @trans_to_int32}  ...
        }, {'Name', 'Disktype', 'Factor', 'Conversion'}, 2);

end

% at least one input must be given!
if nargin == 0
    error( ...
        'xff:BadArgument', ...
        'Invalid call to fifio.' ...
    );
end

% what to do...
action = 0;

% special tag
if nargin > 1 && ...
    isa(varargin{1}, 'double') && ...
    numel(varargin{1}) == 1 && ...
    varargin{1} == 0 && ...
    ischar(varargin{2})

    % what command
    switch (lower(varargin{2}(:)'))
        case {'fif_kinds'}
            varargout{1} = fif_kinds;
            return;
        case {'fif_types'}
            varargout{1} = fif_types;
            return;
        otherwise
            error( ...
                'xff:BadArgument', ...
                'Unknown command for fifio.' ...
            );
    end

% filename given as first argument
elseif nargin == 1 && ...
    ischar(varargin{1}) && ...
   ~isempty(varargin{1}) && ...
    exist(varargin{1}(:)', 'file') == 2
    action = 1;

% FIF struct given, plus command
elseif nargin > 1 && ...
    numel(varargin{1}) == 1 && ...
    isstruct(varargin{1}) && ...
    isfield(varargin{1}, 'Directory') && ...
    ischar(varargin{2}) && ...
    any(strcmpi(varargin{2}(:)', ...
        {'elemblock', 'freeelem', 'readelem'}))

    % get FIF structure
    FIF = varargin{1};
    cmd = lower(varargin{2}(:)');

    % what command
    switch (cmd)

        % get element block(s)
        case {'elemblock'}
            action = 21;

        % read elements
        case {'readelem'}
            action = 5;

        % free elements
        case {'freeelem'}
            action = 11;

    end

end
if action == 0
    error( ...
        'xff:BadArgument', ...
        'Invalid call to fifio.' ...
    );
end

% what action
switch (action)

    % open file
    case {1}
        filename = varargin{1}(:)';
        fid = fopen(filename, 'r', 'ieee-be');
        if fid < 1
            error( ...
                'xff:FileNotReadable', ...
                'File ''%s'' is not readable.', ...
                filename ...
            );
        end
        fseek(fid, 0, 1);
        filesz = ftell(fid);
        fseek(fid, 0, -1);

        % check first few bytes
        try
            fcheck = fread(fid, [1, 4], 'uint32=>double');
            if ~all(fcheck == [100, 31, 20, 0])
                error( ...
                    'xff:BadFile', ...
                    'Not a FIF file: ''%s''.', ...
                    filename ...
                );
            end
            fseek(fid, 0, -1);
        catch ne_eo;
            fclose(fid);
            rethrow(ne_eo);
        end

        % read header
        try
            FIF = readtag(fid, fif_kinds, fif_itypes, fif_types);
        catch ne_eo;
            fclose(fid);
            error( ...
                'xff:ErrorReadingTag', ...
                'Error reading FileID from ''%s'' (%s).', ...
                filename, ne_eo.message ...
            );
        end

        % read second tag
        try
            firstleaf = readtag(fid, fif_kinds, fif_itypes, fif_types);
        catch ne_eo;
            fclose(fid);
            error( ...
                'xff:ErrorReadingTag', ...
                'Error reading first data tag from ''%s'' (%s).', ...
                filename, ne_eo.message ...
            );
        end

        % file has a directory built-in ? then read dir first !
        if strcmp(firstleaf.Datakind, 'DirPointer') && ...
            numel(firstleaf.Value) == 1 && ...
            firstleaf.Value > 0

            try
                % store current file position and go to directory
                cfpos = ftell(fid);
                fseek(fid, firstleaf.Value, -1);

                % read DIR
                try
                    DIR = readtag(fid, fif_kinds, fif_itypes, fif_types);
                catch ne_eo;
                    error(rethrow(ne_eo));
                end

                % go back to position
                fseek(fid, cfpos, -1);

            % on error bail out !
            catch ne_eo;
                fclose(fid);
                error( ...
                    'xff:ErrorReadingDir', ...
                    'Error reading Directory of FIF file at 0x%08X (%s).', ...
                    firstleaf.Value, ne_eo.message ...
                );
            end

            % get true DIR and remove stop tags
            DIR = DIR.Value;
            while all(sprintf('%08X', DIR(end).IOPos) == 'F')
                DIR(end) = [];
            end

        % otherwise
        else

            % build data structure for directory (up to 65536 entries)
            ds = uint32(0);
            ds(1:262144) = ds;

            % set position back to 0
            fseek(fid, 0, -1);
            readpos = ftell(fid);
            maxpos = filesz - 16;

            % go until end of file is hit
            rc = 1;
            while readpos < maxpos

                % read and set information
                ds(rc:rc + 2) = fread(fid, [1, 3], 'uint32');
                ds(rc + 3) = uint32(readpos);
                rawsize = double(ds(rc + 2));

                % read additional offset
                offs = fread(fid, [1, 1], 'int32=>double');
                if offs < 0
                    offs = 0;
                end

                % set new filepos
                try
                    fseek(fid, rawsize + double(offs), 0);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    warning( ...
                        'xff:FileTooShort', ...
                        'File too short for Element at position 0x%08X.', ...
                        readpos ...
                    );
                end

                % increase read counter
                readpos = readpos + 16 + rawsize;
                rc = rc + 4;
            end

            % convert DIR
            DIR = uint32_to_dir(ds, floor(rc / 4));
        end

        % put FIF filename and directory into top leaf
        FIF.Filename = filename;
        FIF.Directory = DIR;

        % parse TREE into FIF as well
        [tinfo{1:4}] = parsedir(DIR, 1, fid, zeros(3, 0), fif_kinds);
        FIF.Tree = tinfo{1};
        FIF.Tree.BlockID = FIF.Value;

        % close file
        fclose(fid);

        % get short lookup list
        lup = zeros(1, numel(DIR));
        for dc = 1:numel(DIR)
            lup(dc) = DIR(dc).Datakind;
        end
        FIF.Lookup = lup;
        FIF.BlockLookup = tinfo{4};

        % keep track of read elements
        FIF.ReadFlag = false(1, numel(DIR));


    % read elements
    case {5}

        % get element specification
        DIR = FIF.Directory;
        nel = numel(DIR);
        if nargin < 3 || ...
           ~isa(varargin{3}, 'double')
            espec = 1:nel;
        else
            espec = fix(varargin{3}(:)');
            espec(espec < 1 | espec > nel) = [];
        end

        % grow children if necessary
        if numel(FIF.Children) < nel
            FIF.Children(nel).Value = [];
        end

        % get all elements that are still to be read
        toread = espec(~FIF.ReadFlag(espec));
        if isempty(toread)
            varargout{1} = FIF;
            return;
        end

        % sort elements
        toread = unique(toread);

        % open file and go to last position
        try
            fid = fopen(FIF.Filename, 'r', 'ieee-be');
            if fid < 1
                error('FILE_OPEN_ERROR');
            end
            etest = toread(end);
            if etest < nel
                etest = etest + 1;
            end
            fseek(fid, DIR(etest).IOPos, -1);
            if ftell(fid) ~= DIR(etest).IOPos
                fclose(fid);
                error('FILE_SEEK_ERROR');
            end
        catch ne_eo;
            error( ...
                'xff:FileError', ...
                'Error opening/seeking in ''%s'': %s.', ...
                FIF.Filename, ne_eo.message ...
            );
        end

        % loop over elements to read
        CHL = FIF.Children;
        for trc = toread

            % go to position
            try
                fseek(fid, DIR(trc).IOPos, -1);
                CHL(trc) = readtag(fid, fif_kinds, fif_itypes, fif_types);
            catch ne_eo;
                fclose(fid);
                error( ...
                    'xff:FileReadError', ...
                    'Error reading element %d in ''%s'': %s.', ...
                    trc, FIF.Filename, ne_eo.message ...
                );
            end
        end

        % put Children back into FIF and keep track of read elements
    	fclose(fid);
        FIF.Children = CHL;
        FIF.ReadFlag(toread) = true;


    % free elements
    case {11}

        % get element specification
        DIR = FIF.Directory;
        nel = numel(DIR);
        if nargin < 3 || ...
           ~isa(varargin{3}, 'double')
            espec = 1:nel;
        else
            espec = fix(varargin{3}(:)');
            espec(espec < 1 | espec > nel) = [];
        end

        % grow children if necessary
        if numel(FIF.Children) < nel
            FIF.Children(nel).Value = [];
        end

        % clear value
        for ec = espec
            FIF.Children(ec).Value = [];
        end

        % keep track of changes
        FIF.ReadFlag(espec) = false;


    % element block(s) lookup
    case {21}

        % check argument
        isinb = zeros(3, 0);
        varargout{1} = isinb;
        if nargin < 3 || ...
           ~isa(varargin{3}, 'double') || ...
            numel(varargin{3}) ~= 1 || ...
            varargin{3} ~= fix(varargin{3}) || ...
            varargin{3} < 1
            return;
        else
            espec = varargin{3};
        end

        % get dir
        DIR = FIF.Directory;
        nel = numel(DIR);
        if espec > nel
            return;
        end

        % get block lookup and traverse...
        bup = FIF.BlockLookup;
        maxb = size(bup, 2);
        for bc = maxb:-1:1
            if espec >= bup(2, bc) && ...
                espec <= bup(3, bc)
                isinb = [isinb, bup(:, bc)];
            end
        end

        % return looked up value
        varargout{1} = isinb;
        return;
end

% by default, return the FIF object
varargout{1} = FIF;


% % % % % %


%
% internal functions (conversion, reading)
%


%function nv = chinfo_to_uint8(v, ne, k)
%    nv = uint8(0);
%    nv(1:64 * ne) = nv;
% end of function v = chinfo_to_uint8(v, ne, k)
%

%function v = complex_to_double(v, ne, k)
%
% end of function v = complex_to_double(v, ne, k)
%

%function v = complex_to_single(v, ne, k)
%
% end of function v = complex_to_single(v, ne, k)
%

%function v = dir_to_uint32(v, ne, k)
%
% end of function v = dir_to_uint32(v, ne, k)
%

%function v = double_to_complex(v, ne, k)
%
% end of function v = double_to_complex(v, ne, k)
%

%function v = id_to_uint8(v, ne, k)
%
% end of function v = id_to_uint8(v, ne, k)
%

function nv = int32_to_chinfo(v, ne, k)
    v = double(v);
    nv = newchinfo;
    if ne > 1
        nv(1:ne) = nv;
    end
    for ec = 1:ne
        ecc = 24 * (ec - 1) + 1;
        nv(ec).ScanNumber   = v(ecc);
        nv(ec).LogNumber    = v(ecc + 1);
        nv(ec).Kind         = v(ecc + 2);
        nv(ec).Range        = x_int32_to_float32(v(ecc + 3));
        nv(ec).Calibration  = x_int32_to_float32(v(ecc + 4));
        nv(ec).CoilType     = v(ecc + 5);
        ctstr = sprintf('C%04X', double(nv(ec).CoilType));
        if isfield(k, ctstr)
            nv(ec).CoilType = k.(ctstr);
        end
        nv(ec).Location     = x_int32_to_float32(v(ecc+6:ecc+17));
        nv(ec).Unit         = v(ecc + 18);
        nv(ec).UnitExponent = v(ecc + 19);
        nv(ec).ChannelName  = x_int32_to_string(v(ecc+20:ecc+23));
        switch (nv(ec).Kind)
            case {1, 301}
                loc = nv(ec).Location;
                nv(ec).CoilTransformation = ...
                    [[loc(4:6); loc(7:9); loc(10:12); loc(1:3)]';[0, 0, 0, 1]];
                nv(ec).CoordinateFrame = 'Device';
            case {2}
                loc = nv(ec).Location;
                if any(loc(4:6) ~= 0)
                    nv(ec).EEGLocation = [loc(1:3); loc(4:6)]';
                else
                    nv(ec).EEGLocation = loc(1:3)';
                end
        end
    end
% end of function nv = int32_to_chinfo(v, ne, k)


function nv = int32_to_point(v, ne, k)
    v = double(v);
    vk = (v < 0);
    vk(1:5:end) = false;
    vk(2:5:end) = false;
    v(vk) = v(vk) + 2 ^ 32;
    nv = newpoint;
    if ne > 1
        nv(1:ne) = nv;
        for ec = 1:ne
            ecc = 5 * (ec - 1) + 1;
            nv(ec).Kind  = v(ecc);
            pkind = sprintf('P%04X', double(nv(ec).Kind));
            if isfield(k, pkind)
                nv(ec).Kind = k.(pkind);
            end
            nv(ec).Ident = v(ecc + 1);
            nv(ec).Coord = x_int32_to_float32(v(ecc+2:ecc+4));
        end
    else
        nv.Kind  = v(1);
        pkind = sprintf('P%04X', double(nv.Kind));
        if isfield(k, pkind)
            nv.Kind = k.(pkind);
        end
        nv.Ident = v(2);
        nv.Coord = x_int32_to_float32(v(3:5));
    end
% end of function nv = int32_to_point(v, ne, k)


function nv = int32_to_trans(v, ne, k)
    v = double(v);
    vk = (v < 0);
    vk(1:26:end) = false;
    vk(2:26:end) = false;
    v(vk) = v(vk) + 2 ^ 32;
    nv = newtrans;
    if ne > 1
        nv(1:ne) = nv;
        for ec = 1:ne
            ecc = 26 * (ec - 1) + 1;
            nv(ec).From = v(ecc);
            tystr = sprintf('F%04X', double(nv(ec).From));
            if isfield(k, tystr)
                nv(ec).From = k.(tystr);
            end
            nv(ec).To   = v(ecc + 1);
            tystr = sprintf('F%04X', double(nv(ec).To));
            if isfield(k, tystr)
                nv(ec).To = k.(tystr);
            end
            rotate = reshape(x_int32_to_float32(v(ecc+2:ecc+10)), [3, 3])';
            transl = reshape(x_int32_to_float32(v(ecc+11:ecc+13)), [3, 1]);
            nv(ec).Trans = [[rotate, transl]; 0, 0, 0, 1];
            rotate = reshape(x_int32_to_float32(v(ecc+14:ecc+22)), [3, 3])';
            transl = reshape(x_int32_to_float32(v(ecc+23:ecc+25)), [3, 1]);
            nv(ec).InvTrans = [[rotate, transl]; 0, 0, 0, 1];
        end
    else
        nv.From  = v(1);
        tystr = sprintf('F%04X', double(nv.From));
        if isfield(k, tystr)
            nv.From = k.(tystr);
        end
        nv.To    = v(2);
        tystr = sprintf('F%04X', double(nv.To));
        if isfield(k, tystr)
            nv.To = k.(tystr);
        end
        rotate = reshape(x_int32_to_float32(v(3:11)), [3, 3])';
        transl = reshape(x_int32_to_float32(v(12:14)), [3, 1]);
        nv.Trans = [[rotate, transl]; 0, 0, 0, 1];
        rotate = reshape(x_int32_to_float32(v(15:23)), [3, 3])';
        transl = reshape(x_int32_to_float32(v(24:26)), [3, 1]);
        nv.InvTrans = [[rotate, transl]; 0, 0, 0, 1];
    end
% end of function nv = int32_to_trans(v, ne, k)


function t = newchinfo
    t = emptystruct({ ...
        'ScanNumber', 'LogNumber', 'Kind', 'Range', 'Calibration', ...
        'CoilType', 'CoilTransformation', 'Location', 'EEGLocation', ...
        'CoordinateFrame', 'Unit', 'UnitExponent', 'ChannelName'});
    t(1).ChannelName = '';
% end of function t = newchinfo


function t = newdir
    t = struct( ...
        'IOPos',    0, ...
        'Datakind', 0, ...
        'Datatype', 0, ...
        'Rawsize',  0, ...
        'DKindTxt', 'n/a' ...
    );
% end of function t = newdir


function t = newid
    t = struct( ...
        'Version',   [1, 1], ...
        'MachineID', '00:00:00:00:00:00', ...
        'UTCDate',    0, ...
        'UTCuSecs',   0 ...
    );
% end of function t = newid


function t = newleaf(iopos, ikind, iname, isize, ioffs)

    % create leaf
    t = emptystruct({ ...
        'IOPos', 'Datakind', 'Datatype', 'RawSize', 'Value', ...
        'IOPosOffset', 'Children', 'BlockID', 'Meta'});

    % start to fill leaf
    t(1).IOPos = iopos;
    t.Datakind = ikind;
    t.Datatype = iname;
    t.RawSize = isize;
    t.IOPosOffset = ioffs;

    % generate empty Children list
    t.Children = t;
    t.Children(:) = [];

% end of function t = newleaf


function t = newpoint
    t = struct( ...
        'Kind',  0, ...
        'Ident', 0, ...
        'Coord', [0; 0; 0]);
% end of function t = newpoint


function t = newtrans
    t = struct( ...
        'From',     1, ...
        'To',       2, ...
        'Trans',    eye(4, 4), ...
        'InvTrans', eye(4, 4));
% end of function t = newtrans


function [TREE, dc, en, bup] = parsedir(DIR, dc, fid, bup, fif_kinds)
    % prepare TREE structure
    TREE = newleaf(20, 'FileID', 'int32', 4, 0);
    TREE.IOPos = DIR(dc).IOPos;
    TREE.BlockID = newid;
    dirsize = numel(DIR);
    en = 0;
    cdc = dc;
    bp = size(bup, 2) + 1;
    bup(:, bp) = [-1; 0; dirsize];

    % set start dir, maxdir and bogus block num
    TREE.Value = [dc, dirsize, -1];

    % parse DIR structure to TREE
    while dc <= dirsize

        % get next element's datakind
        dt = DIR(dc).Datakind;

        % look for ...
        switch (dt)

            % block ID
            case {103}

                % go to position
                try
                    fseek(fid, DIR(dc).IOPos, -1);
                    iblockid = fread(fid, [1, 36], 'uint8=>double');
                catch ne_eo;
                    error( ...
                        'xff:BadFilePosition', ...
                        'Could not read block ID at position 0x%08X: %s.', ...
                        DIR(dc).IOPos, ne_eo.message ...
                    );
                end
                if any(iblockid(4:4:12) ~= [103, 31, 20])
                    error( ...
                        'xff:BadFileContent', ...
                        'Bad block ID at position 0x%08X.', ...
                        DIR(dc).IOPos ...
                    );
                end
                TREE.BlockID = uint8_to_id(iblockid(17:end), 1);

            % block start
            case {104}

                % get block ID
                try
                    fseek(fid, DIR(dc).IOPos, -1);
                    iblocknum = fread(fid, [1, 5], 'int32=>double');
                catch ne_eo;
                    error( ...
                        'xff:BadFilePosition', ...
                        'Could not read block start at position 0x%08X: %s.', ...
                        DIR(dc).IOPos, ne_eo.message ...
                    );
                end
                if any(iblocknum(1:3) ~= [104, 3, 4])
                    error( ...
                        'xff:BadFileContent', ...
                        'Bad block number at position 0x%08X.', ...
                        DIR(dc).IOPos ...
                    );
                end

                % get children list
                [newchild, dc, endnum, bup] = ...
                    parsedir(DIR, dc + 1, fid, bup, fif_kinds);

                % check block number
                if endnum ~= iblocknum(end)
                    error( ...
                        'xff:NestingError', ...
                        'Block nesting error for block ending at 0x%08X.', ...
                        DIR(dc).IOPos ...
                    );
                end

                % fill in values
                TREE.Children(end + 1) = newchild;

            % block end
            case {105}

                % get block ID
                try
                    fseek(fid, DIR(dc).IOPos, -1);
                    iblocknum = fread(fid, [1, 5], 'int32=>double');
                catch ne_eo;
                    error( ...
                        'xff:BadFilePosition', ...
                        'Could not read block start at position 0x%08X: %s.', ...
                        DIR(dc).IOPos, ne_eo.message ...
                    );
                end
                if any(iblocknum(1:3) ~= [105, 3, 4])
                    error( ...
                        'xff:BadFileContent', ...
                        'Bad block number at position 0x%08X.', ...
                        DIR(dc).IOPos ...
                    );
                end
                en = iblocknum(end);

                % set end dir marker
                TREE.Value(2:3) = [dc - 1, en];
                ikind = sprintf('B%04X', en);
                if isfield(fif_kinds, ikind)
                    TREE.Datakind = fif_kinds.(ikind);
                else
                    TREE.Datakind = sprintf('Unknown_%d', en);
                end
                bup(:, bp) = [en; cdc; dc - 1];
                return;
        end

        % increase counter
        dc = dc + 1;
    end

    % keep track of block lookup
    bup(:, bp) = [TREE.Value(3); cdc; dc - 1];

% end of function TREE = parsedir(DIR)


%function v = point_to_uint8(v, ne, k)
%
% end of function v = point_to_uint8(v, ne, k)
%

function t = readtag(fid, fif_kinds, fif_itypes, fif_types, skip)

    % get current IO position
    iopos = ftell(fid);

    % read Datakind and build lookup string
    ikind = fread(fid, [1, 2], 'uint16=>double');
    tkind = sprintf('K%04X', ikind(2));

    % unknown Datakind ?
    if ~isfield(fif_kinds, tkind)
        warning( ...
            'xff:UnknownTagInFile', ...
            'Unknown datakind 0x%04x at position 0x%08X.', ...
            ikind(2), iopos ...
        );

        % then read raw
        tkind = 'KRAWD';
    end

    % read Datatype
    itype = fread(fid, [1, 2], 'uint16=>double');
    imenc = itype(1);
    itype = itype(2);

    % unknown matrix format ?
    if ~any([0, 16384, 16400] == imenc)
        warning( ...
            'xff:UnknownMatrixEncoding', ...
            'Unknown matrix encoding scheme 0x%04x at position 0x%08X.', ...
            ikind(1), iopos ...
        );

        % equally read raw
        tkind = 'KRAWD';
    end

    % for NOP use double
    if itype == 0 && ...
        ikind(2) == 108
        itype = 5;
    end

    % invalid or unknown Datatype
    if itype < 1 || ...
        itype > 35 || ...
        strcmp(fif_types(itype).Name, 'n/a')
        warning( ...
            'xff:UnknownDatatype', ...
            'Unknown datatype 0x%04x at position 0x%08X.', ...
            itype, iopos ...
        );

        % read raw then
        itype = 15;

    end

    % get Datatype structure
    itype = fif_types(itype);

    % read size and IOPositionOffset
    isize = fread(fid, [1, 1], 'uint32=>double');
    ioffs = fread(fid, [1, 1], 'int32=>double');

    % create new leaf object
    t = newleaf(iopos, fif_kinds.(tkind), itype.Name, isize, ioffs);

    % calculate number of target elements
    ielem = floor(isize / itype.Factor);

    % calculate final IOPosition
    tpos = iopos + 16 + isize(1);
    if ioffs ~= 0 && ioffs ~= -1
        tpos = tpos + ioffs;
    end

    % on skip reading, do so
    if nargin > 4 && skip
        fseek(fid, tpos, -1);
        return;
    end

    % only read if elements to read
    if ielem > 0

        % matrix encoding
        switch (imenc)

            % none
            case {0}
                ivalue = fread(fid, ...
                    [1, floor(isize / fif_itypes.(itype.Disktype))], ...
                    ['*' itype.Disktype]);

            % dense (dimensional)
            case {16384}
                t.Datakind = ['DENSE_' t.Datakind];

                % keep current position
                diopos = ftell(fid);

                % read number of dimensions and dimensions
                fseek(fid, isize - 4, 0);
                inumd = fread(fid, [1, 1], 'int32=>double');
                fseek(fid, -4 * (inumd + 1), 0);
                idims = fread(fid, [1, inumd], 'int32=>double');

                % go back to original data position
                fseek(fid, diopos, -1);

                % now do the read
                if inumd < 2
                    idims(2) = 1;
                end
                ivalue = fread(fid, idims, ['*' itype.Disktype]);

                % make sure to go to good position !
                fseek(fid, tpos, -1);

            % sparse (dimensional)
            case {16400}
                t.Datakind = ['SPARSE_' t.Datakind];

                % keep current position
                diopos = ftell(fid);

                % read number of dimensions and dimensions
                fseek(fid, isize - 4, 0);
                inumd = fread(fid, [1, 1], 'int32=>double');
                fseek(fid, -4 * (inumd + 1), 0);
                idims = fread(fid, [1, inumd], 'int32=>double');
                fseek(fid, -4 * (inumd + 2), 0);
                inset = fread(fid, [1, 1], 'int32=>double');

                % go back to original data position
                fseek(fid, diopos, -1);

                % read the data to store
                tvalue = fread(fid, [1, inset], ['*' itype.Disktype]);

                % read where to store it
                target = fread(fid, [1, inset], 'int32=>double') + 1;

                % read where lines start
                nlines = prod(idims(1:end-1));
                tlines = fread(fid, [1, nlines + 1], 'int32=>double') + 1;

                % iterate over lines
                ilinel = idims(1);
                for lc = 1:nlines
                    ton  = tlines(lc);
                    toff = tlines(lc+1) - 1;
                    target(ton:toff) = target(ton:toff) + ilinel * (lc - 1);
                end

                % now create the target matrix
                if inumd < 2
                    idims(2) = 1;
                end
                ivalue = zeros(idims);
                eval(['ivalue=' itype.Disktype '(ivalue);'], '');

                % set elements
                try
                    ivalue(target) = tvalue;
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    error( ...
                        'xff:AssignmentFailed', ...
                        'Error assigning sparse matrix values.' ...
                    );
                end

                % make sure to go to good position !
                fseek(fid, tpos, -1);
        end
    else
        eval(['ivalue = ' itype.Disktype '([]);'], 'ivalue=[];');
        fseek(fid, tpos, -1);
        return;
    end

    % convert with true conversion function ?
    if ~isempty(itype.Conversion)
        ivalue = feval(itype.Conversion{1}, ivalue, ielem, fif_kinds);

    % all single-element value go to double for simplicity
    elseif ielem == 1
        ivalue = double(ivalue);
    end
    t.Value = ivalue;
% end of function t = readtag(...)


%function v = single_to_complex(v, ne, k)
%
% end of function v = single_to_complex(v, ne, k)
%

%function v = string_to_uint8(v, ne, k)
%
% end of function v = string_to_uint8(v, ne, k)
%

%function v = trans_to_uint8(v, ne, k)
%
% end of function v = trans_to_uint8(v, ne, k)
%

function nv = uint32_to_dir(v, ne, k)
    v = double(v);
    nv = newdir;
    if ne > 1
        nv(1:ne) = nv;
    end
    for ec = 1:ne
        ecc = 4 * (ec - 1) + 1;
        nv(ec).Datakind = v(ecc);
        nv(ec).Datatype = v(ecc + 1);
        nv(ec).Rawsize  = v(ecc + 2);
        nv(ec).IOPos    = v(ecc + 3);
        dktxt = sprintf('K%04X', double(nv(ec).Datakind));
        if isfield(k, dktxt)
            nv(ec).DKindTxt = k.(dktxt);
        end
    end
% end of function v = uint32_to_dir(v, ne, k)


function v = uint8_to_id(v, ne, k)
    v = double(v);
    nv = newid;
    if ne > 1
        nv(1:ne) = nv;
        for ec = 1:ne
            ecc = 20 * (ec - 1) + 1;
            nv(ec).Version   = ...
                x_uint8_to_version(v(ecc:ecc+3));
            nv(ec).MachineID = ...
                x_uint8_to_machineid(v(ecc+4:ecc+11));
            nv(ec).UTCDate   = ...
                x_uint8_to_utcdate(v(ecc+12:ecc+15));
            nv(ec).UTCuSecs  = ...
                x_uint8_to_long(v(ecc+16:ecc+19));
        end
    else
        nv.Version   = x_uint8_to_version(v(1:4));
        nv.MachineID = x_uint8_to_machineid(v(5:12));
        nv.UTCDate   = x_uint8_to_utcdate(v(13:16));
        nv.UTCuSecs  = x_uint8_to_long(v(17:20));
    end
    v = nv;
% end of function v = uint8_to_id(v, ne, k)


function v = uint8_to_string(v, ne, k)
    v = char(double(v));
    v = v(1:min(numel(v), ne));
% end of function v = uint8_to_string(v, ne, k)


function nv = x_int32_to_float32(v)
    v = double(v);
    vk = (v < 0);
    v(vk) = v(vk) + 2 ^ 32;
    nv = double(hxsingle(sprintf('%08X', v)));
% end of function nv = x_uint8_to_long(v)


function nv = x_int32_to_string(v)
    v = double(v(:)');
    vk = (v < 0);
    v(vk) = v(vk) + 2 ^ 32;
    v1 = floor(v / 16777216);
    v2 = mod(floor(v / 65536), 256);
    v3 = mod(floor(v / 256), 256);
    v4 = mod(v, 256);
    nv = reshape(char([v1;v2;v3;v4]), [1, 4 * numel(v)]);
    pnull = find(nv == char(0));
    if ~isempty(pnull)
        nv(pnull(1):end) = [];
    end
    if isempty(nv)
        nv = '';
    end
% end of function nv = x_uint8_to_long(v)


function nv = x_uint8_to_long(v)
    nv = 16777216 * v(1) + 65536 * v(2) + 256 * v(3) + v(4);
% end of function nv = x_uint8_to_long(v)


function nv = x_uint8_to_machineid(v)
    if any(v(7:8) ~= 0) && all(v(1:2) == 0)
        machid = v(3:8);
    elseif all(v(7:8) == 0)
        machid = v(1:6);
    else
        machid = v(5:8);
    end
    machid = sprintf('%02X:', machid);
    nv = machid(1:end-1);
% end of function nv = x_uint8_to_machineid(v)


function nv = x_uint8_to_version(v)
    nv = [256 * v(1) + v(2), 256 * v(3) + v(4)];
% end of function nv = x_uint8_to_version(v)


function nv = x_uint8_to_utcdate(v)
    nv = datestr(719529 + ...
        (16777216 * v(1) + 65536 * v(2) + 256 * v(3) + v(4)) / 86400, 31);
% end of function nv = x_uint8_to_utc(v)

function [varargout] = dcmio(varargin)
% dcmio  - read DICOM files
%
% FORMAT:       dcmobject = dcmio(filename)
%
% Input fields:
%
%       filename    filename of DICOM file to read
%
% Output fields:
%
%       dcmobject   xff object
%
% See also xff

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:16 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% persistent VR dict and empty file
persistent my_dicom_vr my_dicom_dic my_newdcm;
if isempty(my_dicom_vr)
    my_dicom_vr = dicom_vr;
    my_dicom_dic = dicom_dic(4);
    newdcm = xff('new:dcm');
    my_newdcm = getcont(newdcm);
    delete(newdcm);
end
vrSQ = my_dicom_vr.SQ;
vrUN = my_dicom_vr.UN;

% argument check
if nargin < 1 || ...
    ~ischar(varargin{1}) || ...
    isempty(varargin{1})
    error( ...
        'xff:BadArgument', ...
        'Bad or missing argument for dcmio.' ...
    );
end
filename = varargin{1}(:)';

% default options
writemode = false;

% dcm content structure given
if nargin > 1 && ...
    xffisobject(varargin{2}, true, 'dcm') && ...
    numel(varargin{2}) == 1
    try
        dcmcont = getcont(varargin{2});
    catch ne_eo;
        rethrow(ne_eo);
    end
    writemode = true;
end

% reading DCM content
if ~writemode

    % create new object
    dcmcont = my_newdcm;

    % check file (opening in little endian syntax by default)
    fid = fopen(filename, 'rb', 'ieee-le');
    if fid < 1
        error( ...
            'xff:FileNotReadable', ...
            'File not readable: ''%s''.', ...
            filename ...
        );
    end

    % get file size
    fseek(fid, 0, 1);
    flen = ftell(fid);
    flenmin = flen - 8;
    fseek(fid, 0, -1);
    if flen < 132
        fclose(fid);
        error( ...
            'xff:InvalidFilesize', ...
            'File too short: ''%s''.', ...
            filename ...
        );
    end

    % get file content (little endian, up to 16MB at most)
    frmost = min(flen, 16777216);
    fc_uint8 = fread(fid, [1, frmost], '*uint8');
    frmost = floor(frmost / 2);
    fseek(fid, 0, -1);
    fc_uint16_le = fread(fid, [1, frmost], '*uint16');
    fseek(fid, 0, -1);

    % get file content (big endian)
    fidbe = fopen(filename, 'rb', 'ieee-be');
    fc_uint16_be = fread(fidbe, [1, frmost], '*uint16');
    fseek(fidbe, 0, -1);

% writing BFF content
else

    % check open writable file
    if dcmcont.MetaLittleEndian
        fid = fopen(filename, 'wb', 'ieee-le');
    else
        fid = fopen(filename, 'wb', 'ieee-be');
    end
    if fid < 1
        error( ...
            'xff:FileNotWritable', ...
            'File not writable: ''%s''.', ...
            filename ...
        );
    end
    fclose(fid);

    % however, still not yet implemented
    error( ...
        'xff:NotYetImplemented', ...
        'Writing of DICOM files not yet implemented.' ...
    );

end

% counters
MetaKeyCount = 1;
DataKeyCount = 0;

% try used for forcemode (on indentation !!!)
try % forcemode

    % reading
    if ~writemode

        % detect magic token
        dcmcont.Magic = 'DICM';
        dcmcont.Preamble = true;

        % if magic token present
        if all(fc_uint8(129:132) == [68, 73, 67, 77])

            % start reading after token
            fpos = 133;

        % otherwise
        else

            % reset defaults and start at char 1!
            dcmcont.Magic = '';
            dcmcont.Preamble = false;
            fpos = 1;
        end

        % defaults, detection
        dcmcont.MetaTSExplicit = false;
        dcmcont.MetaLittleEndian = true;

        % if Magic token found
        if dcmcont.Preamble

            % initialize meta struct
            MetaStr = dcmcont.Meta;

            % check implicit explicit
            vrcheck = fc_uint8(fpos+4:fpos+5);
            dcmcont.MetaTSExplicit = all(vrcheck > 64 & vrcheck < 91);

            % check endian only on explicit VR syntax
            MetaExplicit = dcmcont.MetaTSExplicit;
            if MetaExplicit
                lecheck = fc_uint8(fpos:fpos+4);
                if lecheck(1) < lecheck(2)
                    dcmcont.MetaLittleEndian = false;
                elseif lecheck(1) == lecheck(2)
                    if lecheck(3) < lecheck(4)
                        dcmcont.MetaLittleEndian = false;
                    elseif lecheck(3) == lecheck(4)
                        error( ...
                            'xff:FileDetectionFailed', ...
                            'Error detecting Endian type in Meta header.' ...
                        );
                    end
                end
            end
            MetaLEndian = dcmcont.MetaLittleEndian;
            if MetaLEndian
                rfid = fid;
            else
                rfid = fidbe;
            end

            % read until key is over group 7
            while true

                % get correct key format
                if MetaLEndian
                    DicomKey = double(fc_uint16_le(round(fpos/2):round(fpos/2)+1));
                else
                    DicomKey = double(fc_uint16_be(round(fpos/2):round(fpos/2)+1));
                end

                % break on condition
                if DicomKey(1) > 7
                    break;
                end
                fpos = fpos + 4;
                MetaStr(MetaKeyCount).Key = DicomKey;
                if DicomKey(1) > 65533
                    MetaSequence = true;
                    DataVR = vrSQ;
                else
                    MetaSequence = false;
                end

                % VR given ?
                if ~MetaSequence && ...
                    MetaExplicit
                    DataVR = upper(char(fc_uint8(fpos:fpos+1)));
                    fpos = fpos + 2;
                    if ~isfield(my_dicom_vr, DataVR)
                        error( ...
                            'xff:InvalidToken', ...
                            'Invalid VR found: %s.', ...
                            DataVR ...
                        );
                    end
                    DataVR = my_dicom_vr.(DataVR);
                elseif ~MetaSequence
                    DataVR = vrUN;
                end
                MetaStr(MetaKeyCount).VR = DataVR.tag;
                if DataVR.length(2) > 32767
                    DataShortVLength = false;
                else
                    DataShortVLength = true;
                end

                % length
                if MetaLEndian
                    DataVLengthShort = double(fc_uint16_le(round(fpos/2)));
                else
                    DataVLengthShort = double(fc_uint16_be(round(fpos/2)));
                end
                if ~DataShortVLength && ...
                   ~MetaSequence && ...
                    MetaExplicit && ...
                    DataVLengthShort > 0
                    error( ...
                        'xff:InvalidFileContent', ...
                        'Invalid 16-bit VL for given VR.' ...
                    );
                elseif MetaSequence
                    fpos = fpos - 2;
                end
                if ~DataShortVLength
                    if MetaLEndian
                        DataVLengthLong = 65536 * ...
                            double(fc_uint16_le(round(fpos/2)+2)) + ...
                            double(fc_uint16_le(round(fpos/2)+1));
                    else
                        DataVLengthLong = 65536 * ...
                            double(fc_uint16_le(round(fpos/2)+1)) + ...
                            double(fc_uint16_le(round(fpos/2)+2));
                    end
                    fpos = fpos + 6;
                    dvl = DataVLengthLong;
                else
                    fpos = fpos + 2;
                    DataVLengthLong = NaN;
                    dvl = DataVLengthShort;
                end
                MetaStr(MetaKeyCount).VLShort = DataVLengthShort;
                MetaStr(MetaKeyCount).VLLong  = DataVLengthLong;
                dvl = 2 * round(dvl/2);

                % get value according to type
                ddt = lower(DataVR.datatype);
                if ~strcmp(ddt, 'sequence')

                    % go to and update position
                    fseek(rfid, fpos - 1, -1);
                    fpos = fpos + dvl;

                    % number of objects
                    switch (ddt)
                        case {'double'}
                            dvl = round(dvl / 8);
                        case {'int32', 'single', 'uint32'}
                            dvl = round(dvl / 4);
                        case {'int16', 'uint16'}
                            dvl = dvl / 2;
                    end
                    dval = fread(rfid, [1, dvl], ['*' ddt]);
                    if ischar(dval) && ...
                        isempty(dval)
                        dval = '';
                    end
                    if DataVR.deblank
                        dval = deblank(dval);
                    end
                else
                    dval = [];
                end
                MetaStr(MetaKeyCount).Value = dval;
                MetaKeyCount = MetaKeyCount + 1;

            end

            % put meta back
            dcmcont.Meta = MetaStr;
            mkey = cell(length(MetaStr), 1);
            mkys = struct;
            for c = 1:length(mkey)
                mky = sprintf('k_%04X_%04X', MetaStr(c).Key);
                mkey{c} = mky;
                mkys.(mky) = c;
            end
            dcmcont.MetaKeys = mkey;
            dcmcont.MetaKeyLookup = mkys;
        end

        % initialize data struct
        dstr = dcmcont.Data;
        DataLEndian = true;

        % check implicit/explicit TS for data part
        vrcheck = fc_uint8(fpos+4:fpos+5);
        dcmcont.DataTSExplicit = all(vrcheck > 64 & vrcheck < 91);
        dvs = dcmcont.DataTSExplicit;

        % check endian only on explicit VR syntax
        if dcmcont.DataTSExplicit
            lecheck = fc_uint8(fpos:fpos+4);
            if lecheck(1) < lecheck(2)
                dcmcont.DataLittleEndian = false;
            elseif lecheck(1) == lecheck(2)
                if lecheck(3) < lecheck(4)
                    dcmcont.DataLittleEndian = false;
                elseif lecheck(3) == lecheck(4)
                    warning( ...
                        'xff:FileDetectionFailed', ...
                        'Error detecting Endian type in Data part, assuming LE.' ...
                    );
                    dcmcont.DataLittleEndian = true;
                end
            end
            DataLEndian = dcmcont.DataLittleEndian;
        end
        if DataLEndian
            rfid = fid;
        else
            rfid = fidbe;
        end

        % data part until fpos > flenmin
        DStrKeyCount = numel(dstr);
        while fpos <= flenmin

            % get correct key format
            DataKeyCount = DataKeyCount + 1;
            fposh = round(fpos/2);
            if DataLEndian
                DicomKey = double([fc_uint16_le(fposh), fc_uint16_le(fposh+1)]);
            else
                DicomKey = double([fc_uint16_be(fposh), fc_uint16_be(fposh+1)]);
            end
            if DataKeyCount > DStrKeyCount
                if DStrKeyCount < 128
                    DStrKeyCount = 128;
                end
                dstr(ceil(1.4 * DStrKeyCount)).Key = [];
                DStrKeyCount = numel(dstr);
            end
            dstr(DataKeyCount).Key = DicomKey;
            fpos = fpos + 4;

            % sequence ?
            if DicomKey(1) > 65533
                DataSequence = true;
                DataVR = vrSQ;
            else
                DataSequence = false;
            end

            % VR ?
            if ~DataSequence && ...
                dvs
                DataVR = upper(char(fc_uint8(fpos:fpos+1)));
                fpos = fpos + 2;
                if ~isfield(my_dicom_vr, DataVR)
                    error( ...
                        'xff:InvalidToken', ...
                        'Invalid VR found: %s.', ...
                        DataVR ...
                    );
                end
                DataVR = my_dicom_vr.(DataVR);
            elseif ~dvs
                tkey = sprintf('k_%04X_%04X', DicomKey(1), DicomKey(2));
                try
                    DataVR = my_dicom_vr.(upper(my_dicom_dic.Dictionary.(tkey).vr));
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    DataVR = vrUN;
                end
            elseif ~DataSequence
                DataVR = vrUN;
            end
            dstr(DataKeyCount).VR = DataVR.tag;
            if DataVR.length(2) > 32767
                DataShortVLength = false;
            else
                DataShortVLength = true;
            end

            % length
            fposh = round(fpos/2);
            if DataLEndian
                DataVLengthShort = double(fc_uint16_le(fposh));
            else
                DataVLengthShort = double(fc_uint16_be(fposh));
            end
            if ~DataShortVLength && ...
               ~DataSequence && ...
                dvs && ...
                DataVLengthShort > 0
                error( ...
                    'xff:InvalidFileContent', ...
                    'Invalid 16-bit VL for given VR.' ...
                );
            elseif DataSequence && ...
                any([57357, 57565] == DicomKey(2)) && ...
                all(fc_uint16_le(fposh:fposh+1) == 0)
                fpos = fpos - 2;
                fposh = fposh - 1;
            end
            if ~DataShortVLength || ...
                 ~dvs
                if DataVLengthShort == 0 && ...
                    dvs
                    if DataLEndian
                        DataVLengthLong = 65536 * ...
                            double(fc_uint16_le(fposh+2)) + ...
                            double(fc_uint16_le(fposh+1));
                    else
                        DataVLengthLong = 65536 * ...
                            double(fc_uint16_le(fposh+1)) + ...
                            double(fc_uint16_le(fposh+2));
                    end
                    fpos = fpos + 6;
                else
                    if DataLEndian
                        DataVLengthLong = 65536 * ...
                            double(fc_uint16_le(fposh+1)) + ...
                            double(fc_uint16_le(fposh));
                    else
                        DataVLengthLong = 65536 * ...
                            double(fc_uint16_le(fposh)) + ...
                            double(fc_uint16_le(fposh+1));
                    end
                    fpos = fpos + 4;
                end
                dvl = DataVLengthLong;
            else
                fpos = fpos + 2;
                DataVLengthLong = NaN;
                dvl = DataVLengthShort;
            end
            dstr(DataKeyCount).VLShort = DataVLengthShort;
            dstr(DataKeyCount).VLLong  = DataVLengthLong;
            dvl = 2 * round(dvl/2);

            % get value according to type
            ddt = lower(DataVR.datatype);
            if ~strcmp(ddt, 'sequence')

                % go to and update position
                fseek(rfid, fpos - 1, -1);
                fpos = fpos + dvl;

                % number of objects
                switch (ddt)
                    case {'double'}
                        dvl = round(dvl / 8);
                    case {'int32', 'single', 'uint32'}
                        dvl = round(dvl / 4);
                    case {'int16', 'uint16'}
                        dvl = dvl / 2;
                    case {'uint816'}
                        dval = fread(rfid, [1, dvl], 'uint8=>double');
                        dval = dval(:);
                        fseek(rfid, -dvl, 0);
                        cr1 = corrcoef(dval(1:end-1), dval(2:end));
                        cr2 = 0.4 .* (corrcoef(dval(1:2:end-2), dval(3:2:end)) + ...
                            corrcoef(dval(2:2:end-2), dval(4:2:end)));
                        if cr1(2) > cr2(2)
                            ddt = 'uint8';
                        else
                            ddt = 'uint16';
                            dvl = dvl / 2;
                        end
                end
                dval = fread(rfid, [1, dvl], ['*' ddt]);
                if ischar(dval) && ...
                    isempty(dval)
                    dval = '';
                end
                if DataVR.deblank
                    dval = deblank(dval);
                end
            else
                dval = [];
            end
            dstr(DataKeyCount).Value = dval;

        end

        % put meta back
        dcmcont.Data = dstr(1:DataKeyCount);
        DicomKey = cell(numel(dcmcont.Data), 1);
        dkys = struct;
        for c = 1:length(DicomKey)
            dky = sprintf('k_%04X_%04X', dstr(c).Key);
            DicomKey{c} = dky;
            dkys.(dky) = c;
        end
        dcmcont.DataKeys = DicomKey;
        dcmcont.DataKeyLookup = dkys;

    % writing
    else

        % not yet implemented

    end

catch ne_eo;
    neuroelf_lasterr(ne_eo);
end

% close main file
fclose(fid);

% create good object
hfile = xff('new:dcm');
setcont(hfile, dcmcont);

% close big-endian file
if ~writemode
    fclose(fidbe);
end

% give correct output
if nargout > 1
    varargout = cell(1, nargout);
    varargout{1} = hfile;
else
    varargout{1} = hfile;
end

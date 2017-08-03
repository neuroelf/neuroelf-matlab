function [varargout] = rdaio(varargin)
% rdaio  - read R files
%
% FORMAT:       rdaobject = rdaio(filename)
%
% Input fields:
%
%       filename    filename of RData file to read
%
% Output fields:
%
%       rdaobject   xff object
%
% See also xff.

% See http://yetanothermathprogrammingconsultant.blogspot.com/2016/02/r-rdata-file-format.html

% Version:  v1.1
% Build:    16060209
% Date:     Jun-02 2016, 9:25 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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
persistent my_newrda;
if isempty(my_newrda)
    newrda = xff('new:rda');
    my_newrda = getcont(newrda);
    delete(newrda);
end

% argument check
if nargin < 1 || ~ischar(varargin{1}) || isempty(varargin{1})
    error('neuroelf:general:badArgument', 'Bad or missing argument for rdaio.');
end
filename = varargin{1}(:)';

% default options
writemode = false;

% dcm content structure given
if nargin > 1 && xffisobject(varargin{2}, true, 'rda') && numel(varargin{2}) == 1
    try
        rdacont = getcont(varargin{2});
    catch ne_eo;
        rethrow(ne_eo);
    end
    writemode = true;
end

% reading DCM content
tmpfile = '';
if ~writemode

    % create new object
    rdacont = my_newrda;

    % check file (opening in little endian syntax by default)
    fid = fopen(filename, 'rb', 'ieee-be');
    if fid < 1
        error('neuroelf:general:fileOpenError', 'Cannot open RDA file for reading.');
    end
    fseek(fid, 0, 1);
    fsize = ftell(fid);
    if fsize < 48
        fclose(fid);
        error('neuroelf:fileio:fileTooShort', 'RDA file too short.');
    end
    fseek(fid, 0, -1);

    % read the first 7 bytes and check content
    fhead = fread(fid, [1, 7], 'uint8=>double');
    if ~isequal(fhead, [82, 68, 65, 50, 10, 65, 10]) && ...
       ~isequal(fhead, [82, 68, 88, 50, 10, 88, 10]) && ...
       ~isequal(fhead, [88, 10, 0, 0, 0, 2, 0])
        fclose(fid);

        % test for a GZIP-ed file
        tmpfile = [tempname '.rda.gz'];
        try
            copyfile(filename, tmpfile);
            gunzip(tmpfile, fileparts(tmpfile));
            if exist(tmpfile, 'file') == 2
                delete(tmpfile);
            end
            tmpfile = tmpfile(1:end-3);
        catch ne_eo;
            if exist(tmpfile, 'file') == 2
                delete(tmpfile);
            end
            rethrow(ne_eo);
        end
        fid = fopen(tmpfile, 'rb', 'ieee-be');
        if fid < 1
            if exist(tmpfile, 'file') == 2
                delete(tmpfile);
            end
            error('neuroelf:general:fileOpenError', 'Cannot open RDA file for reading.');
        end
        fseek(fid, 0, 1);
        fsize = ftell(fid);
        if fsize < 48
            fclose(fid);
            delete(tmpfile);
            error('neuroelf:fileio:fileTooShort', 'RDA file too short.');
        end
        fseek(fid, 0, -1);
        fhead = fread(fid, [1, 7], 'uint8=>double');
        if ~isequal(fhead, [82, 68, 65, 50, 10, 65, 10]) && ...
           ~isequal(fhead, [82, 68, 88, 50, 10, 88, 10]) && ...
           ~isequal(fhead, [88, 10, 0, 0, 0, 2, 0])
            fclose(fid);
            delete(tmpfile);
            error('neuroelf:fileio:invalidHeader', 'Not an RDA file.');
        end
    end

    % compression status, file magic and type
    rdacont.FileCompressed = ~isempty(tmpfile);
    if fhead(3) ~= 0
        rdacont.FileMagic = char(fhead(1:4));
        rdacont.FileType = char(fhead(6));
    else
        rdacont.FileMagic = '';
        rdacont.FileType = char(fhead(1));
        fseek(fid, -5, 0);
    end

    % type
    ftype = rdacont.FileType;

    % header
    if ftype == 'X'

        % read version info, etc.
        rdacont.FormatVersion = fread(fid, [1, 1], 'uint32=>double');
        rdacont.RVersion = fread(fid, [1, 4], 'uint8=>double');
        rdacont.RVersion(1) = [];
        rdacont.RXVersion = fread(fid, [1, 1], 'uint32=>double');
    else
        rdacont.FormatVersion = str2double(fgetl(fid));
        rv = str2double(fgetl(fid));
        rv = [floor(rv / 65536), floor(mod(rv, 65536) / 256), mod(rv, 256)];
        rdacont.RVersion = rv;
        rdacont.RXVersion = fread(fid, [1, 1], 'uint32=>double');
    end

    % read vars
    try
        if ~isempty(rdacont.FileMagic)
            [rdacont.Vars, dok] = readvars(fid, ftype, fsize);
        else
            [v, dok] = readdata(fid, ftype, fsize);
            if ftell(fid) < fsize
                dok = false;
            end
            if isstruct(v) && numel(v) == 1 && isfield(v, 'Type') && isfield(v, 'Value') && ...
                ischar(v.Type) && strcmp(v.Type, 'CELSXP') && isstruct(v.Value) && ...
                numel(v.Value) == 1 && isfield(v.Value, 'Data') && iscell(v.Value.Data) && ...
                isfield(v.Value, 'Info') && isstruct(v.Value.Info) && numel(v.Value.Info) == 1 && ...
                isfield(v.Value.Info, 'names') && isstruct(v.Value.Info.names) && ...
                numel(v.Value.Info.names) == 1 && isfield(v.Value.Info.names, 'Value') && ...
                iscell(v.Value.Info.names.Value) && ...
                isequal(size(v.Value.Data), size(v.Value.Info.names.Value)) && ...
                all(cellfun(@ischar, v.Value.Info.names.Value(:)))
                ovn = v.Value.Info.names.Value;
                ovv = v.Value.Data;
                rdacont.Vars = cell2struct(cat(2, ovn(:), ovv(:)), {'Name', 'Data'}, 2);
            end
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        dok = false;
    end

    % something didn't go right
    if ~dok

        % read the rest of the file
        rdacont.REMAININGCONTENT = fread(fid, [1, Inf], 'uint8=>uint8');

    % special format
    elseif numel(rdacont.Vars) == 2 && isempty(rdacont.Vars(1).Name) && ...
        isstruct(rdacont.Vars(1).Data) && numel(rdacont.Vars(1).Data) == 1 && ...
        isfield(rdacont.Vars(1).Data, 'Type') && ischar(rdacont.Vars(1).Data.Type) && ...
        strcmp(rdacont.Vars(1).Data.Type(:)', 'CELSXP') && ...
        isfield(rdacont.Vars(1).Data, 'Value') && iscell(rdacont.Vars(1).Data.Value) && ...
        strcmp(rdacont.Vars(2).Name, 'names') && isstruct(rdacont.Vars(2).Data) && ...
        numel(rdacont.Vars(2).Data) == 1 && isfield(rdacont.Vars(2).Data, 'Type') && ...
        ischar(rdacont.Vars(2).Data.Type) && strcmp(rdacont.Vars(2).Data.Type, 'STRSXP') && ...
        isfield(rdacont.Vars(2).Data, 'Value') && iscell(rdacont.Vars(2).Data.Value) && ...
        isequal(size(rdacont.Vars(1).Data.Value), size(rdacont.Vars(2).Data.Value)) && ...
        all(cellfun(@ischar, rdacont.Vars(2).Data.Value(:)))

        % regenerate content
        ovv = rdacont.Vars(1).Data.Value;
        ovn = rdacont.Vars(2).Data.Value;
        rdacont.Vars = cell2struct(cat(2, ovn(:), ovv(:)), {'Name', 'Data'}, 2);
    end
end

% close main file
fclose(fid);

% remove temp file
if ~isempty(tmpfile)

    % delete (for read)
    if ~writemode
        delete(tmpfile)

    % gzip to final file otherwise
    else
    end
end

% reading
if ~writemode

    % create good object
    hfile = xff('new:rda');
    setcont(hfile, rdacont);

    % give correct output
    if nargout > 1
        varargout = cell(1, nargout);
        varargout{1} = hfile;
    else
        varargout{1} = hfile;
    end
else
end



% sub functions
function [v, dok] = readvars(fid, ftype, fsize)

% initialize vars
v = emptystruct({'Name', 'Data'});

% keep reading until end of file (or list) marker
if ftype == 'X'
    dtype = fread(fid, [1, 1], 'uint32=>double');
else
    dtype = str2double(fgetl(fid));
end
vc = 1;
while dtype ~= 254 && ftell(fid) < fsize
    try
        [d, dok] = readdata(fid, ftype, fsize, dtype);
        if dtype == 1026 && numel(d) == 1 && isstruct(d) && numel(fieldnames(d)) == 2 && ...
            isfield(d, 'Type') && isfield(d, 'Value') && ischar(d.Type) && ...
            strcmp(d.Type(:)', 'LISTSXP') && iscell(d.Value) && numel(d.Value) == 2 && ...
            numel(d.Value{1}) == 1 && isstruct(d.Value{1}) && isfield(d.Value{1}, 'Type') && ...
            ischar(d.Value{1}.Type) && strcmp(d.Value{1}.Type(:)', 'SYMSXP') && ...
            isfield(d.Value{1}, 'Value') && ischar(d.Value{1}.Value) && ...
            ~isempty(d.Value{1}.Value) && isrealvarname(d.Value{1}.Value(:)')
            v(vc).Name = d.Value{1}.Value;
            v(vc).Data = d.Value{2};
        else
            v(vc).Data = d;
        end
        if ~dok
            break;
        end
        vc = vc + 1;
        if ftype == 'X'
            dtype = fread(fid, [1, 1], 'uint32=>double');
        else
            dtype = str2double(fgetl(fid));
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        dok = false;
        break;
    end
end


% read one variable/element of data
function [d, dok] = readdata(fid, ftype, fsize, dtype, dlen)

% assume OK
dok = true;

% initialize variable data
d = struct('Type', '', 'Value', []);

% variable type
if nargin < 4 || isempty(dtype)
    if ftype == 'X'
        dtype = fread(fid, [1, 1], 'uint32=>double');
    else
        dtype = str2double(fgetl(fid));
    end
end

% we should NOT get an end marker, but just in case
if dtype == 254
    return;
end

% unknown type
if ~any(dtype == [1, 9, 10, 13, 14, 16, 19, 525, 526, 531, 781, 787, 1023, 1026, 1791, 2047, 2815, 262153])
    dok = false;
    return;
end

% type with length field upfront
if any(mod(dtype, 256) == [9, 10, 13, 14, 16, 19]) && nargin < 5
    if ftype == 'X'
        dlen = fread(fid, [1, 1], 'uint32=>double');
    else
        dlen = str2double(fgetl(fid));
    end
elseif nargin < 5
    dlen = [];
end

% depending on datatype
switch dtype

    % SYMBOL
    case 1
        d.Type = 'SYMSXP';
        try
            [dv, dok] = readdata(fid, ftype, fsize);
            d.Value = dv;
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
        end
        if ~dok
            return;
        end

    % CHRSXP
    case 9
        d.Type = 'CHRSXP';
        if ftype == 'X'
            dv = fread(fid, [1, dlen], 'uint8=>char');
        else
            dv = fgetl(fid);
        end
        if numel(dv) ~= dlen
            dok = false;
        end
        d.Value = dv;

    % TFVSXP
    case 10
        d.Type = 'TFVSXP';
        if ftype == 'X'
            dv = (fread(fid, [1, dlen], 'uint32=>double') ~= 0);
        else
            dv = false(1, dlen);
            for vc = 1:dlen
                dv(vc) = ~isequal(str2double(fgetl(fid)), 0);
            end
        end
        if numel(dv) ~= dlen
            dok = false;
        end
        d.Value = dv;

    % INTSXP
    case 13
        d.Type = 'INTSXP';
        if ftype == 'X'
            dv = fread(fid, [1, dlen], 'int32=>double');
            if numel(dv) ~= dlen
                dok = false;
            end
        else
            dv = zeros(1, dlen);
            for vc = 1:dlen
                dv(vc) = str2double(fgetl(fid));
            end
        end
        d.Value = dv;

    % DBLSXP
    case 14
        d.Type = 'DBLSXP';
        if ftype == 'X'
            dv = fread(fid, [1, dlen], 'double=>double');
            if numel(dv) ~= dlen
                dok = false;
            end
        else
            dv = zeros(1, dlen);
            for vc = 1:dlen
                dv(vc) = str2double(fgetl(fid));
            end
        end
        d.Value = dv;

    % STRSXP
    case 16
        d.Type = 'STRSXP';
        dv = cell(1, dlen);
        try
            for vc = 1:dlen
                [dv{vc}, dok] = readdata(fid, ftype, fsize);
                if ~dok
                    break;
                end
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
        end
        d.Value = dv;

    % CELSXP
    case {19, 531}
        d.Type = 'CELSXP';
        dv = cell(1, dlen);
        try
            for vc = 1:dlen
                [dv{vc}, dok] = readdata(fid, ftype, fsize);
                if ~dok
                    break;
                end
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
        end
        d.Value = dv;

        % list info
        if dok && dtype > 255
            [di, dok] = readinfo(fid, ftype, fsize);
            if dok
                d.Value = struct('Data', {dv}, 'Info', di);
            end
        end

    % DBLSXP
    case {525, 526}
        if dtype == 525
            d.Type = 'DBMSXP';
        else
            d.Type = 'DBASXP';
        end
        if ftype == 'X'
            if dtype == 525
                dv = fread(fid, [1, dlen], 'int32=>double');
            else
                dv = fread(fid, [1, dlen], 'double=>double');
            end
            if numel(dv) ~= dlen
                dok = false;
            end
        else
            dv = zeros(1, dlen);
            for vc = 1:dlen
                dv(vc) = str2double(fgetl(fid));
            end
        end
        
        % matrix info
        if dok
            [di, dok] = readinfo(fid, ftype, fsize);
        else
            di = struct;
        end
        d.Value = struct('Data', dv, 'Info', di);

    % FACSXP
    case 781
        d.Type = 'FACSXP';
        try
            [dv, dok] = readdata(fid, ftype, fsize, 13, dlen);
            dv = dv.Value;
            if iscell(dv)
                dv = {dv};
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
            d.Value = struct('Data', dv, 'Info', []);
            return;
        end

        % matrix info
        if dok
            [di, dok] = readinfo(fid, ftype, fsize);
        else
            di = struct;
        end
        d.Value = struct('Data', dv, 'Info', di);

    % DFRSXP
    case 787
        d.Type = 'DFRSXP';
        try
            [dv, dok] = readdata(fid, ftype, fsize, 19, dlen);
            dv = dv.Value;
            if iscell(dv)
                dv = {dv};
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
            d.Value = struct('Data', dv, 'Info', []);
            return;
        end

        % matrix info
        if dok
            [di, dok] = readinfo(fid, ftype, fsize);
        else
            di = struct;
        end
        d.Value = struct('Data', dv, 'Info', di);

    % class name info
    case 1023
        d.Type = 'CLISXP';
        try
            [dv, dok] = readdata(fid, ftype, fsize);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
        end
        d.Value = dv;
        
    % list (symbol + value)
    case 1026
        d.Type = 'LISTSXP';
        try
            [dv1, dok] = readdata(fid, ftype, fsize);
            if isstruct(dv1) && numel(dv1) == 1 && isfield(dv1, 'Type') && ...
                ischar(dv1.Type) && strcmp(dv1.Type, 'SYMSXP')
                d.Value = {dv1, readdata(fid, ftype, fsize)};
            else
                d.Value = {dv1};
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
        end

    % factor info (levels)
    case 1791
        d.Type = 'Levels';
        try
            [dv, dok] = readdata(fid, ftype, fsize);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
        end
        if ~dok
            dv = 'BAD';
        end
        d.Value = dv;

    % factor info (class)
    case 2047
        d.Type = 'Class';
        try
            [dv, dok] = readdata(fid, ftype, fsize);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
        end
        if ~dok
            dv = 'BAD';
        end
        d.Value = dv;

    % unnamed matrix dims
    case 2815
        d.Type = 'MDims';
        try
            [dv, dok] = readdata(fid, ftype, fsize);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            dok = false;
        end
        if ~dok
            dv = 'BAD';
        end
        d.Value = dv;

    % STRING
    case 262153
        if ftype == 'X'
            d = fread(fid, [1, dlen], 'uint8=>char');
        else
            d = fgetl(fid);
        end
        if numel(d) ~= dlen
            dok = false;
        end

    % unknown
    otherwise
        dok = false;
end

% read info fields
function [di, dok] = readinfo(fid, ftype, fsize)
dok = true;
di = struct;
try
    if ftype == 'X'
        sdtype = fread(fid, [1, 1], 'uint32=>double');
    else
        sdtype = str2double(fgetl(fid));
    end
    while dok && sdtype ~= 254 && ftell(fid) < fsize
        [div, dok] = readdata(fid, ftype, fsize, sdtype);
        if isstruct(div) && numel(div) == 1 && isfield(div, 'Type') && ...
            ischar(div.Type) && strcmp(div.Type(:)', 'LISTSXP') && ...
            isfield(div, 'Value') && iscell(div.Value) && numel(div.Value) == 2 && ...
            isstruct(div.Value{1}) && numel(div.Value{1}) == 1 && ...
            isfield(div.Value{1}, 'Type') && ischar(div.Value{1}.Type) && ...
            strcmp(div.Value{1}.Type(:)', 'SYMSXP') && isfield(div.Value{1}, 'Value') && ...
            ischar(div.Value{1}.Value) && ~isempty(div.Value{1}.Value) && ...
            isrealvarname(div.Value{1}.Value(:)')
            di.(div.Value{1}.Value(:)') = div.Value{2};
        else
            dif = sprintf('Info%04d', numel(fieldnames(di)) + 1);
            di.(dif) = div;
        end
        if dok && ftype == 'X'
            sdtype = fread(fid, [1, 1], 'uint32=>double');
        elseif dok
            sdtype = str2double(fgetl(fid));
        end
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    dok = false;
end

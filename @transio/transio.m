function [varargout] = transio(varargin)
% transio (Object Class)
%
% FORMAT:       tio_obj = transio(file, endian, class, offset, size [, enlarge]);
%
% Input fields:
%
%       file        filename where data is stored in
%       endian      endian type (e.g. 'le', or 'ieee-be')
%       class       numerical class (e.g. 'uint16', 'single')
%       offset      offset within file (in bytes)
%       size        size of array in file
%       enlarge     flag necessary to create content in file
%
% Output fields:
%
%       tio_obj     transio object that supports the methods
%                   - subsref  : tio_obj(I) or tio_obj(IX, IY, IZ)
%                   - subsasgn : same as subsref
%                   - size     : retrieving the array size
%                   - end      : for building relative indices
%                   - display  : showing information
%
% Note 1: enlarging of existing files (if the array is the last element
%         in the file) can be done by adding a (class-independent) sixth
%         parameter to the call.
%
% Note 2: both subsref and subsasgn will only work within the existing
%         limits; growing of the array as with normal MATLAB variables
%         is *NOT* supported--so tio_obj(:,:,ZI) = []; will *NOT* work!
%
% Using: makelabel.

% Version:  v1.1
% Build:    16060314
% Date:     Jun-03 2016, 2:43 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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

% This Class provides transparent IO access to files that store large
% variables at certain place.
%
% When accessed, the file is opened and closed each time! Since MATLAB
% does not offer destructors, it is virtually impossible to keep track
% of open files!
%
% The construction of a transio object can be done with a syntax of:
%
% tio = transio(file, endian, class, offset, size [, enlarge]);
%
% When the fifth argument is not given, transio will issue an error
% if the filesize is not sufficient.

% declare persistent stored variable
persistent transio_singleton;

% check class initialization
if isempty(transio_singleton) || ...
   ~isstruct(transio_singleton) || ...
   ~transio_singleton.is_initialized

    % initialize persistent struct
    transio_singleton = struct;
    transio_singleton.is_initialized = false;

    % neuroelf methods (different from ne_methods)
    makelabel = using(neuroelf, 'makelabel');
    transio_singleton.makelabel = makelabel{1};

    % perform init commands
    transio_singleton.validclasses = struct( ...
        'double',  8, ...
        'int16',   2, ...
        'int32',   4, ...
        'int64',   8, ...
        'int8',    1, ...
        'single',  4, ...
        'uint16',  2, ...
        'uint32',  4, ...
        'uint64',  8, ...
        'uint8',   1  ...
    );
    transio_singleton.validcnames = fieldnames(transio_singleton.validclasses);
    transio_singleton.validendian = struct( ...
        'be',      false, ...
        'le',      true,  ...
        'ieee_be', false, ...
        'ieee_le', true,  ...
        'ieeebe',  false, ...
        'ieeele',  true   ...
    );
    transio_singleton.validenames = fieldnames(transio_singleton.validendian);

    % say we are initialized
    transio_singleton.is_initialized = true;
end

% object creation
if nargin > 1 && isa(varargin{1}, 'double') && ischar(varargin{2}) && ...
    ~isempty(varargin{2}) && ...
    all(((double(varargin{2}(:)) - 109.5) .^ 2) < 156.5)

    % what action to perform
    switch lower(varargin{2}(:)')
        case 'makeobject'
            if isstruct(varargin{3}) && numel(varargin{3}) == 1 && ...
                length(fieldnames(varargin{3})) == 8 && ...
                all(strcmp(fieldnames(varargin{3}), { ...
                    'FileName'; ...
                    'LittleND'; ...
                    'DataType'; ...
                    'TypeSize'; ...
                    'IOOffset'; ...
                    'DataDims'; ...
                    'ReadOnly'; ...
                    'IOBuffer'}))
                varargout{1} = class(varargin{3}, 'transio');
            end
        otherwise
            if ~any(strcmpi(transio_singleton.validcnames, varargin{2}(:)'))
                error( ...
                    'transio:BadArgument', ...
                    'Unknown action string: %s.', ...
                    varargin{2}(:)' ...
                );
            end
            varargout{1} = ...
                transio_singleton.validclasses.(lower(varargin{2}(:)'));
    end

% check arguments > object creation call
elseif nargin > 4 && ...
    ischar(varargin{1}) && ...
    ~isempty(varargin{1}) && ...
    any(strcmp(transio_singleton.makelabel(lower(varargin{2}(:)')), transio_singleton.validenames)) && ...
    any(strcmpi(varargin{3}(:)', transio_singleton.validcnames)) && ...
    isa(varargin{4}, 'double') && ...
    numel(varargin{4}) == 1 && ...
   ~isinf(varargin{4}) && ...
   ~isnan(varargin{4}) && ...
    fix(varargin{4}) == varargin{4} && ...
    varargin{4} >= 0 && ...
    isa(varargin{5}, 'double') && ...
   ~isempty(varargin{5}) && ...
    numel(varargin{5}) == length(varargin{5}) && ...
    numel(varargin{5}) < 32 && ...
   ~any(isinf(varargin{5}) | isnan(varargin{5}) | fix(varargin{5}) ~= varargin{5}) && ...
   ~any(varargin{5} < 1 | varargin{5} > 2147483647)

    % make absolute filename
    [isabs, filename] = isabsolute(varargin{1}(:)');

    % open the file
    try
        rdo = false;
        if nargin < 6
            fid = fopen(filename, 'r+');
        else
            fid = fopen(filename, 'a+');
        end
        if fid < 1
            rdo = true;
            fid = fopen(filename, 'r');
            if fid < 1
                error('INVALID_FID');
            end
        end
    catch ne_eo;
        error( ...
            'transio:ErrorOpeningFile', ...
            'Error opening file ''%s'': %s.', ...
            filename, ne_eo.message ...
        );
    end

    % calculate most significant position
    lnd = transio_singleton.validendian.(transio_singleton.makelabel(lower(varargin{2}(:)')));
    cls = lower(varargin{3}(:)');
    csz = transio_singleton.validclasses.(cls);
    obs = varargin{5}(:)';
    while numel(obs) > 1 && ...
        obs(end) == 1
        obs(end) = [];
    end
    while numel(obs) < 2
        obs(2) = 1;
    end
    obs = obs(:)';
    ofs = varargin{4};
    msp = ofs + csz * prod(obs);

    % go to end position and check size
    fseek(fid, 0, 1);
    fsz = ftell(fid);

    % file too short and no enlarging
    if fsz < msp && ...
       (nargin < 6 || ...
        rdo)

        % give an error
        error( ...
            'transio:FileTooShort', ...
            'The file is too short to read/write the data.' ...
        );

    % enlarging
    elseif fsz < msp

        % create array to write
        war = uint8(0);
        war(1:1048576) = war(1);
        tfsz = fsz;
        mspf = msp - 1048576;
        try
            while tfsz < mspf
                fwrite(fid, war, 'uint8');
                tfsz = tfsz + 1048576;
            end
            fwrite(fid, war(1:(msp - tfsz)), 'uint8');
            % check again
            fseek(fid, 0, 1);
            if ftell(fid) ~= msp
                error('ENLARGE_ERROR');
            end
        catch ne_eo;
            error( ...
                'transio:ErrorEnlargingFile', ...
                'Error while enlarging the file: %s.', ...
                ne_eo.message ...
            );
        end
    end

    % close file again, we don't keep it open anyway
    fclose(fid);

    % build object
    sobj.FileName = filename;
    sobj.LittleND = lnd;
    sobj.DataType = cls;
    sobj.TypeSize = csz;
    sobj.IOOffset = ofs;
    sobj.DataDims = obs;
    sobj.ReadOnly = rdo;
    sobj.IOBuffer = {{}, []};
    obj = class(sobj, 'transio');

    % return object
    varargout{1} = obj;

% multiple filenames and offsets
elseif nargin > 4 && ...
    iscell(varargin{1}) && ...
    numel(varargin{1}) > 1 && ...
    any(strcmp(transio_singleton.makelabel(lower(varargin{2}(:)')), transio_singleton.validenames)) && ...
    any(strcmpi(varargin{3}(:)', transio_singleton.validcnames)) && ...
    isa(varargin{4}, 'double') && ...
    numel(varargin{5}) == length(varargin{5}) && ...
    numel(varargin{4}) == numel(varargin{1}) && ...
   ~any(isinf(varargin{4}) | isnan(varargin{4}) | ...
        varargin{4} ~= fix(varargin{4}) | varargin{4} < 0) && ...
    isa(varargin{5}, 'double') && ...
    numel(varargin{5}) > 1 && ...
    numel(varargin{5}) < 32 && ...
    numel(varargin{5}) == length(varargin{5}) && ...
   ~any(isinf(varargin{5}) | isnan(varargin{5}) | fix(varargin{5}) ~= varargin{5}) && ...
   ~any(varargin{5} < 1 | varargin{5} > 2147483647) && ...
    varargin{5}(end) == numel(varargin{1})

    % calculate most significant position
    lnd = transio_singleton.validendian.(transio_singleton.makelabel(lower(varargin{2}(:)')));
    cls = lower(varargin{3}(:)');
    csz = transio_singleton.validclasses.(cls);
    ofs = varargin{4}(:)';
    obs = varargin{5}(:)';
    obf = prod(obs(1:end-1));

    % make absolute filenames
    filenames = varargin{1}(:)';
    numfiles = numel(filenames);
    rdo = false(1, numfiles);
    for fc = 1:numfiles

        % not char or empty
        if ~ischar(filenames{fc}) || ...
            isempty(filenames{fc})
            error( ...
                'transio:BadArgument', ...
                'Invalid multi-file compoment (file %d).', ...
                fc ...
            );
        end
        [isabs, filenames{fc}] = isabsolute(filenames{fc}(:)');

        % open the files (one by one)
        try
            fid = fopen(filenames{fc}, 'r+');
            if fid < 1
                rdo(fc) = true;
                fid = fopen(filenames{fc}, 'r');
                if fid < 1
                    error('INVALID_FID');
                end
            end

            % go to end position and check size
            fseek(fid, 0, 1);
            fsz = ftell(fid);
            fclose(fid);
            if fsz < (ofs(fc) + csz * obf)
                error('FILE_TOO_SHORT');
            end
        catch ne_eo;
            error( ...
                'transio:ErrorOpeningFile', ...
                'Error opening file, or file too short (''%s'').', ...
                filenames{fc}, ne_eo.message ...
            );
        end
    end

    % build object
    sobj.FileName = filenames;
    sobj.LittleND = lnd;
    sobj.DataType = cls;
    sobj.TypeSize = csz;
    sobj.IOOffset = ofs;
    sobj.DataDims = obs;
    sobj.ReadOnly = any(rdo);
    sobj.IOBuffer = {{}, []};
    obj = class(sobj, 'transio');

    % return object
    varargout{1} = obj;

% default constructor
elseif nargin == 0
    varargout{1} = class(struct( ...
        'FileName', 'NONE', ...
        'LittleND', 1, ...
        'DataType', 'double', ...
        'TypeSize', 8, ...
        'IOOffset', 0, ...
        'DataDims', [0, 0], ...
        'ReadOnly', true, ...
        'IOBuffer', {{{}, []}}), 'transio');

% otherwise give an error
else
    error( ...
        'transio:BadArgument', ...
        'Invalid combination of arguments.' ...
    );
end

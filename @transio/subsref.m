function [varargout] = subsref(htio, S)
% transio::subsref  - overloaded method
%
% transio(I)       - values at position(s) I
% transio(XI, YI)  - values at position(s) XI, YI
% transio(end, YI) - values at position(s) size(transio, 1), YI

% Version:  v0.9c
% Build:    11090911
% Date:     Sep-9 2011, 11:37 AM EST
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

% argument check
if nargin < 2 || ...
   ~isstruct(S) || ...
    isempty(S) || ...
   ~isfield(S, 'type') || ...
   ~isfield(S, 'subs')
    error( ...
        'transio:BadSubsRef', ...
        'No S subsref struct given.' ...
    );
end

% allow struct subsref for certain fields
if strcmp(S(1).type, '.') && ...
   (ischar(S(1).subs) || ...
    (numel(S(1).subs) == 1 && ...
     iscell(S(1).subs) && ...
     ischar(S(1).subs{1})))

    % cell field
    if iscell(S(1).subs)
        S(1).subs = S(1).subs{1};
    end

    % valid field
    if any(strcmp(S(1).subs(:)', {'DataDims', 'DataType', 'FileName', ...
        'IOBuffer', 'IOOffset', 'LittleND', 'ReadOnly', 'TypeSize'}))

        % pass on to struct representation
        varargout{1} = subsref(struct(htio), S(1));

        % further on?
        if numel(S) > 1
            try
                varargout{1} = subsref(varargout{1}, S(2:end));
            catch ne_eo;
                rethrow(ne_eo);
            end
        end
        return;

    % unknown field
    else
        error( ...
            'transio:UnknownProperty', ...
            'Unknown property: %s.', ...
            S.subs{1}(:)' ...
        );
    end

% allow sub-file subsref
elseif strcmp(S(1).type, '{}') && ...
    numel(S(1).subs) == 1

    % try to apply to FileName
    varargout{1} = htio;
    if ~iscell(varargout{1}.FileName)
        varargout{1}.FileName = {varargout{1}.FileName};
        varargout{1}.DataDims(end+1) = 1;
    end
    varargout{1}.FileName = varargout{1}.FileName(S(1).subs{1}(:)');
    varargout{1}.FileName = varargout{1}.FileName(:)';
    varargout{1}.IOOffset = varargout{1}.IOOffset(S(1).subs{1}(:)');
    varargout{1}.IOOffset = varargout{1}.IOOffset(:)';
    if numel(varargout{1}.FileName) == 0
        error( ...
            'transio:BadArgument', ...
            'Empty transio not supported.' ...
        );
    elseif numel(varargout{1}.FileName) == 1
        varargout{1}.FileName = varargout{1}.FileName{1};
        varargout{1}.DataDims(end) = [];
        if numel(varargout{1}.DataDims) == 1
            varargout{1}.DataDims(2) = 1;
        end
    else
        varargout{1}.DataDims(end) = numel(varargout{1}.FileName);
    end

    % further on
    if numel(S) > 1
        try
            varargout{1} = subsref(varargout{1}, S(2:end));
        catch ne_eo;
            rethrow(ne_eo);
        end
    end

    % return early
    return;
end

% buffered set
if ~isempty(htio.IOBuffer{1}) && ...
    isequal(htio.IOBuffer{1}, S(1).subs)

    % return buffer
    varargout{1} = htio.IOBuffer{2};
    return;
end

% get datatype and type size for access functions
csn = ['*' htio.DataType];
csz = htio.TypeSize;
siz = htio.DataDims;

% handle multi-files
if iscell(htio.FileName)

    % disallow bad subsref
    if numel(S) > 1 || ...
       ~strcmp(S.type, '()')
        error( ...
            'transio:BadSubsRef', ...
            'Bad subsref construct.' ...
        );
    end
    idx = S.subs;
    hdd = htio.DataDims;
    nd = numel(hdd);

    % allow ':'
    if numel(idx) == 1 && ...
        ischar(idx{1}) && ...
        strcmp(idx{1}, ':')

        % use resolve
        varargout{1} = resolve(htio);

    % indexing below number of dims
    elseif numel(idx) < nd

        % pass on to first object
        varargout{1} = htio;
        varargout{1}.FileName = varargout{1}.FileName{1};
        varargout{1}.DataDims(end) = [];
        if numel(varargout{1}.DataDims) < 2
            varargout{1}.DataDims(2) = 1;
        end
        varargout{1}.IOOffset = varargout{1}.IOOffset(1);
        varargout{1} = subsref(varargout{1}, S);

    % multi-dim indexing
    else

        % remove ':' and ,1 datadims at end
        for ic = numel(idx):-1:(nd+1)
            if (ischar(idx{ic}) && ...
                strcmp(idx{ic}, ':')) || ...
                isequal(idx{ic}, 1) || ...
                isequal(idx{ic}, true)
                idx(ic) = [];
            else
                break;
            end
        end

        % too many dims?
        if numel(idx) > nd
            error( ...
                'transio:BadSubsRef', ...
                'Index exceeds matrix dimensions.' ...
            );
        end

        % how many outputs (cells)
        idnd = idx{nd};
        if ischar(idnd) && ...
            strcmp(idnd, ':')
            idnd = 1:hdd(nd);
        elseif islogical(idnd) && ...
            numel(idnd) == hdd(nd)
            idnd = find(idnd(:));
        elseif isnumeric(idnd)
            if any(isinf(idnd(:)) | isnan(idnd(:)) | ...
                idnd(:) < 1 | idnd(:) > hdd(nd) | idnd(:) ~= fix(idnd(:)))
                error( ...
                    'transio:BadSubsRef', ...
                    'Invalid indexing expression.' ...
                );
            end
            idnd = idnd(:);
        else
            error( ...
                'transio:BadSubsRef', ...
                'Invalid indexing expression.' ...
            );
        end

        % create cell array with last index
        lnd = numel(idnd);
        varargout{1} = cell(lnd, 1);

        % create copy
        cpio = htio;
        cpio.DataDims(end) = [];
        if numel(cpio.DataDims) < 2
            cpio.DataDims(2) = 1;
        end
        Ss = S;
        Ss.subs = idx(1:nd-1);

        % and then parse along
        try
            for lnc = 1:lnd
                cpio.FileName = htio.FileName{idnd(lnc)};
                cpio.IOOffset = htio.IOOffset(idnd(lnc));
                varargout{1}{lnc} = subsref(cpio, Ss);
            end
        catch ne_eo;
            rethrow(ne_eo);
        end
        varargout{1} = cat(nd, varargout{1}{:});
    end

    % return early
    return;
end

% parse S
try

    % use external function to efficiently get indices
    [ots, idx, ridx, cidx, lastcol] = ioidx(S, siz);
    pots = prod(ots);
    if numel(ots) > numel(siz)
        cpr = [1, cumprod(ots)]';
    else
        cpr = [1, cumprod(siz)]';
    end
    lsz = numel(idx);
    lsx = lsz + 1;
    if length(cpr) > lsx
        cpr = cpr(1:lsx);
    end
catch ne_eo;
    error( ...
        'transio:BadSubsRef', ...
        'Invalid S subsref struct given: ''%s''.', ...
        ne_eo.message ...
    );
end

% initialize output to class
varargout{1} = eval([csn(2:end) '(0)']);

% return empty if nothing is to be read
if pots < 1
    varargout{1}(1) = [];
    varargout{1} = reshape(varargout{1}, ots);
    return;
end

% open file and check
try
    if htio.LittleND
        fid = fopen(htio.FileName, 'rb', 'ieee-le');
    else
        fid = fopen(htio.FileName, 'rb', 'ieee-be');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    fid = 0;
end
if fid < 1
    error( ...
        'transio:FileNotOpen', ...
        'Error opening file: ''%s''.', ...
        htio.FileName ...
    );
end

% get offset within file
ofs = htio.IOOffset;

% take specific care of one length arguments
if lsz == 1

    % treat singular-index input
    idx = idx{1};

    % character index
    if ischar(idx)

        % go to position
        try
            fseek(fid, ofs, -1);
            varargout{1} = fread(fid, [pots, 1], csn);
        catch ne_eo;
            fclose(fid);
            error( ...
                'transio:FileReadError', ...
                'Error reading array(:) from file: %s.', ...
                ne_eo.message ...
            );
        end
        fclose(fid);
        return;
    end

    % if contiguous
    if length(cidx{1}) < 3

        % try read
        try
            fseek(fid, ofs + csz * (idx(1) - 1), -1);
            varargout{1} = reshape(fread(fid, [pots, 1], csn), ots);
        catch ne_eo;
            fclose(fid);
            error( ...
                'transio:FileReadError', ...
                'Error reading array part from file: %s.', ...
                ne_eo.message ...
            );
        end

    % otherwise
    else

        % initialize output
        varargout{1}(pots) = varargout{1}(1);
        varargout{1} = reshape(varargout{1}, ots);

        % get cidx for fast access
        cidx = cidx{1};

        % then ...
        try

            % read as long as breaks are found
            dsz = cidx(1);
            for bc = 1:(length(cidx) - 1)

                % get position of break
                dsz2 = cidx(bc + 1);

                % go to position in file
                fseek(fid, ofs + csz * (idx(dsz) - 1), -1);

                % read as much as needed
                varargout{1}(dsz:(dsz2 - 1)) = fread(fid, [dsz2 - dsz, 1], csn);

                % and copy old break to new start
                dsz = dsz2;
            end

        catch ne_eo;
            error( ...
                'transio:FileReadError', ...
                'Error reading array part from file: %s.', ...
                ne_eo.message ...
            );
        end
    end

    % close file
    fclose(fid);

    % reorder indices if needed
    if ~ischar(ridx{1})
        varargout{1} = reshape(varargout{1}(ridx{1}), ots);
    end

    % return now
    return;

% all are :
elseif lastcol == lsz

    % go to first index and read entire array
    try
        fseek(fid, ofs, -1);
        varargout{1} = reshape(fread(fid, [pots, 1], csn), ots);
    catch ne_eo;
        fclose(fid);
        error( ...
            'transio:FileReadError', ...
            'Error reading array(:,...,:) from file: %s.', ...
            ne_eo.message ...
        );
    end

    % close file and keep track of it!
    fclose(fid);
    fid = 0;

% all but last are colon and last is contiguous
elseif lastcol == (lsz - 1) && ...
    length(cidx{end}) < 3

    % go to first index and read entire array
    try
        fseek(fid, ofs + csz * (idx{end}(1) - 1) * cpr(end - 1), -1);
        varargout{1} = reshape(fread(fid, [pots, 1], csn), ots);
    catch ne_eo;
        fclose(fid);
        error( ...
            'transio:FileReadError', ...
            'Error reading array(:,...,I) from file: %s.', ...
            ne_eo.message ...
        );
    end

    % close file and keep track of it!
    fclose(fid);
    fid = 0;
end

% reordering indices
rS.type = '()';
rS.subs = cell(1, lsz);
userS = false;
for bc = 1:lsz
    rS.subs{bc} = ridx{bc};
    if ~ischar(ridx{bc})
        userS = true;
    end
end

% something has been read already
if fid < 1

    % apply reordering
    if userS
        varargout{1} = reshape(subsref(varargout{1}, rS), ots);
    end

    % return now
    return;
end

% we have true multiple indexing, so initialize output now at the latest
varargout{1}(pots) = varargout{1}(1);
varargout{1} = reshape(varargout{1}, ots);

% rebuilt S
rds = ots;
ofc = ones(1, lsx);
otc = ones(1, lsx);
otf = ones(1, lsx);
S = struct;
S.type = '()';
S.subs = {':'};
rcol = lastcol + 1;
for cc = 1:lsz
    if cc < rcol
        S.subs{cc} = ':';
    else
        S.subs{cc} = 1;
        ofc(cc) = idx{cc}(1);

        % set final block counter
        if length(cidx{cc}) > 1
            otf(cc) = length(cidx{cc}) - 1;
        else
            otf(cc) = length(idx{cc});
        end
        rds(cc) = 1;
    end
end

% read specific parts
cidx = cidx{rcol};
cpsiz = cpr(rcol);
while otc(end) < 2

    % find position
    rpos = ofs + csz * (ofc - 1) * cpr;

    % which indices to read next
    rstr = cidx(otc(rcol)):(cidx(otc(rcol)+1)-1);
    S.subs{rcol} = rstr;

    % number of indices
    rsiz = length(rstr);

    % set into reshape array
    rds(rcol) = rsiz;

    % calculate number of elements
    rdp = cpsiz * rsiz;

    % try to read and directly store in varargout{1}
    try
        fseek(fid, rpos, -1);
        varargout{1} = subsasgn(varargout{1}, S, reshape(fread(fid, [rdp, 1], csn), rds));
    catch ne_eo;
        fclose(fid);
        error( ...
            'transio:FileReadError', ...
            'Error reading array part from file: %s.', ...
            ne_eo.message ...
        );
    end

    % start with rcol
    cc = rcol;

    % and go on until last column reaches
    while cc <= lsx

        % increase block count
        otc(cc) = otc(cc) + 1;

        % more blocks to come for this dim?
        if otc(cc) <= otf(cc)

            % increase non-contiguous counter only for later dims
            if cc > rcol
                S.subs{cc} = S.subs{cc} + 1;
            end

            % reset ofc
            if cc == rcol
                ofc(cc) = idx{cc}(cidx(otc(cc)));
            else
                ofc(cc) = idx{cc}(otc(cc));
            end

            % go on ...
            break;

        % last column (termination criterion)
        elseif cc > lsz

            % means: leave reading loop!
            break;
        end

        % set block counter for cc to 1
        otc(cc) = 1;

        % and also set subs structure
        S.subs{cc} = 1;

        % and ofc
        ofc(cc) = idx{cc}(1);

        % then increase column counter
        cc = cc + 1;
    end
end

% close file !!
fclose(fid);

% reorder indices for rcol ?
if userS
    varargout{1} = reshape(subsref(varargout{1}, rS), ots);
end

function htio = subsasgn(htio, S, V)
% transio::subsasgn  - overloaded method
%
% transio(I) = V       - write values at position(s) I
% transio(XI, YI) = V  - write values at position(s) XI, YI
% ...

% Version:  v0.9c
% Build:    13111411
% Date:     May-2 2011, 5:37 PM EST
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

% class check
if nargin > 2 && ...
    ~isa(htio, 'transio')
    try
        htio = builtin('subsasgn', htio, S, V);
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% not allowed for read-only files
if htio.ReadOnly
    error( ...
        'transio:NoWriteOnReadOnly', ...
        'Write access to read-only files denied.' ...
    );

% not allowed for multiple files
elseif iscell(htio.FileName)
    error( ...
        'transio:NoWriteOnMultiFiles', ...
        'Write access to multi-files transio not implemented.' ...
    );
end

% argument check
if nargin < 3 || ...
   ~isstruct(S) || ...
   (~strcmpi(htio.DataType, class(V)) && ...
    ~isa(V, 'double'))
    error( ...
        'transio:BadSubsAsgn', ...
        'No or bad S struct, or no or bad V given.' ...
    );
end

% get and check sizes
csn = htio.DataType;
csz = htio.TypeSize;
siz = htio.DataDims;
Vsz = size(V);
pvs = prod(Vsz);

% parse S
try

    % use external function to efficiently get indices
    [its, idx, ridx, cidx, lastcol] = ioidx(S, siz);
    pits = prod(its);
    if length(its) > length(siz)
        cpr = [1, cumprod(its)]';
    else
        cpr = [1, cumprod(siz)]';
    end
    lsz = length(idx);
    lsx = lsz + 1;
    if length(cpr) > lsx
        cpr = cpr(1:lsx);
    end
    if numel(its) > numel(Vsz) && ...
        all(its(numel(Vsz)+1:end) == 1)
        Vsz = [Vsz, ones(1, numel(its) - numel(Vsz))];
    end
catch ne_eo;
    error( ...
        'transio:BadSubsAsgn', ...
        'Invalid S subsref struct given: ''%s''.', ...
        ne_eo.message ...
    );
end

% check size(V) with its
if pvs ~= pits && ...
    pvs ~= 1
    error( ...
        'transio:BadSubsAsgn', ...
        'Invalid number of elements, must match S or be 1.' ...
    );

% only one value given but more needed
elseif pvs == 1 && ...
    pits > 1

    % make sure not to waste more space as needed
    if isa(V, 'double') && ...
       ~strcmpi(csn, 'double')
        try
            eval(['V=' csn '(V);']);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            error( ...
                'transio:ConversionFailed', ...
                'Error converting double to %s.', ...
                csn ...
            );
        end
    end

    % blow up to larger array for write operation
    V = repmat(V, its);

elseif length(Vsz) ~= length(its) || ...
    any(Vsz ~= its)
    warning( ...
        'transio:BadArgumentSize', ...
        'Possibly bad argument. Reshaping to good size...' ...
    );
    V = reshape(V, its);
end

% check if sizes match now!
if numel(V) ~= pits
    error( ...
        'transio:SizeMismatch', ...
        'Input value V and sizes in S mismatch.' ...
    );
end

% open file and check
try
    if htio.LittleND
        fid = fopen(htio.FileName, 'r+b', 'ieee-le');
    else
        fid = fopen(htio.FileName, 'r+b', 'ieee-be');
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    fid = 0;
end

% not opened (for writing)
if fid < 1

    % repeat up to 5 times with increased pauses
    rpause = 0.01 * round(100 * exp(-4.385:0.6932:0));
    for crp = 1:numel(rpause)
        drawnow;
        pause(rpause(crp));
        try
            if htio.LittleND
                fid = fopen(htio.FileName, 'r+b', 'ieee-le');
            else
                fid = fopen(htio.FileName, 'r+b', 'ieee-be');
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
        if fid > 0
            break;
        end;
    end

    % give up
    if fid < 1
        error( ...
            'transio:FileNotOpen', ...
            'Error opening file ''%s'' for writing.', ...
            htio.FileName ...
        );
    end
end

% get offset into file
ofs = htio.IOOffset;

% check for reordering indices
rS.type = '()';
rS.subs = cell(1, lsz);
userS = false;
for bc = 1:lsz
    if ~ischar(ridx{bc})
        [sidx{1:2}] = sort(ridx{bc});
        userS = true;
    else
        sidx = {'', ':'};
    end
    rS.subs{bc} = sidx{2};
end

% do resorting
if userS
    V = reshape(subsref(V, rS), Vsz);
end

% take specific care of one length arguments
if lsz == 1

    % treat singular-index input
    idx = idx{1};

    % character indexing (all elements)
    if ischar(idx)

        % go to position
        try
            fseek(fid, ofs, -1);
            fwrite(fid, V, csn);
        catch ne_eo;
            fclose(fid);
            error( ...
                'transio:FileWriteError', ...
                'Error writing array to file: %s.', ...
                ne_eo.message ...
            );
        end
        fclose(fid);
        return;
    end

    % if contiguous
    if length(cidx{1}) < 3

        % go to first index and write array
        try
            fseek(fid, ofs + csz * (idx(1) - 1), -1);
            fwrite(fid, V, csn);
        catch ne_eo;
            fclose(fid);
            error( ...
                'transio:FileWriteError', ...
                'Error writing array part to file: %s.', ...
                ne_eo.message ...
            );
        end

    % otherwise
    else
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

                % write as much as needed
                fwrite(fid, V(dsz:(dsz2 - 1)), csn);

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

    % close file and return
    fclose(fid);
    return;

% all are :
elseif lastcol == lsz

    % go to first index and read entire array
    try
        fseek(fid, ofs, -1);
        fwrite(fid, V, csn);
    catch ne_eo;
        fclose(fid);
        error( ...
            'transio:FileWriteError', ...
            'Error writing array to file: %s.', ...
            ne_eo.message ...
        );
    end

    % close file and return
    fclose(fid);
    return;

% all but last are colon and last is contiguous
elseif lastcol == (lsz - 1) && ...
    length(cidx{end}) < 3

    % go to first index and read entire array
    try
        fseek(fid, ofs + csz * (idx{end}(1) - 1) * cpr(end - 1), -1);
        fwrite(fid, V, csn);
    catch ne_eo;
        fclose(fid);
        error( ...
            'transio:FileReadError', ...
            'Error reading array(:,...,I) from file: %s.', ...
            ne_eo.message ...
        );
    end

    % close file and return
    fclose(fid);
    return
end

% rebuilt S
rds = its;
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

% write specific parts
cidx = cidx{rcol};
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

    % try to write array elements
    try
        fseek(fid, rpos, -1);
        fwrite(fid, subsref(V, S), csn);
    catch ne_eo;
        fclose(fid);
        error( ...
            'transio:FileReadError', ...
            'Error writing array part to file: %s.', ...
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

% close file
fclose(fid);

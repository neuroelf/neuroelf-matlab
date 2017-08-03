function q = eq(A,B)
% transio::eq  - checks for equality
%
% FORMAT:       q = eq(A, B)
%
% Input fields:
%
%       A, B        transio object and variable or another transio object
%
% Output fields:
%
%       q           result of equality test

% Version:  v0.9c
% Build:    11052601
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

% not enough inputs
if nargin < 2
    error( ...
        'transio:BadArgument', ...
        'Too few arguments for EQ.' ...
    );
end

% try block
try

    % only one object
    if isa(A, 'transio') && ...
       ~isa(B, 'transio')
        q = tio_eqso(A, B);
    elseif isa(B, 'transio') && ...
       ~isa(A, 'transio')
        q = tio_eqso(B, A);

    % two objects
    elseif isa(A, 'transio') && ...
        isa(B, 'transio')

        % get structs
        sA = struct(A);
        sB = struct(B);

        % single objects -> compare
        if numel(sA) == 1 && ...
            numel(sB) == 1
            q = tio_eqto(A, B);

        % invalid comparison
        elseif numel(sA) ~= 1 && ...
            numel(sB) ~= 1 && ...
            numel(sA) ~= numel(sB)
            error( ...
                'transio:BadArgument', ...
                'EQ requires two arguments of the same size or one singleton.' ...
            );

        % compare actual objects
        else

            % if either is empty
            if isempty(sA) || ...
                isempty(sB)

                % return empty matrix
                q = logical([]);

            % otherwise
            else

                % create array
                q = false(max(numel(sA), numel(sB)));

                % then compare structs
                if numel(sA) == 1
                    for sc = 1:numel(q)
                        q(sc) = isequal(sA, sB(sc));
                    end
                elseif numel(sB) == 1
                    for sc = 1:numel(q)
                        q(sc) = isequal(sA(sc), sB);
                    end
                else
                    for sc = 1:numel(q)
                        q(sc) = isequal(sA(sc), sB(sc));
                    end
                end
            end
        end

    % no objects?
    else
        error( ...
            'transio:BadArgument', ...
            'One argument to @transio/eq must be of type transio.' ...
        );
    end
catch ne_eo;
    rethrow(ne_eo);
end


% sub-functions (A = transio, B = variable)
function q = tio_eqso(A, B)

% only valid for single transio
sA = struct(A);
if numel(sA) ~= 1
    error( ...
        'transio:BadArgument', ...
        'Call not valid for multiple objects array.' ...
    );
end

% compare sizes
szA = sA.DataDims;
nA = prod(szA);
szB = size(B);
nB = numel(B);
if ~isequal(szA, szB) && ...
    nA ~= 1 && ...
    nB ~= 1
    error( ...
        'transio:BadArgument', ...
        'Matrix dimensions must agree.' ...
    );
end

% get datatype
dt = ['*' sA.DataType];

% single value in object (unlikely, but still)
if nA == 1
    q = (B == double(resolve(A)));

% single value in B
elseif nB == 1

    % single file
    if ischar(sA.FileName)
        fname = {sA.FileName};
        fsize = 1;

    % multiple files
    else
        fname = sA.FileName;
        fsize = szA(end);
        szA(end) = [];
    end

    % size array
    q = false(prod(szA), fsize);

    % inner try block (to close file on error)
    try

        % open files
        for fc = 1:fsize
            nA = prod(szA);
            fid = 0;
            if sA.LittleND
                fid = fopen(fname{fc}, 'rb', 'ieee-le');
            else
                fid = fopen(fname{fc}, 'rb', 'ieee-be');
            end

            % check file ID
            if fid < 1
                error( ...
                    'transio:FileNotOpen', ...
                    'Error opening file: ''%s''.', ...
                    fname{fc} ...
                );
            end

            % seek to correct offset
            fseek(fid, sA.IOOffset(fc), -1);

            % while still a large chunk to work on
            iit = 1;
            iitt = 1e6;
            while nA > 1e6

                % read one million entries and work
                q(iit:iitt, fc) = eq(fread(fid, [1e6, 1], dt), B);
                iit = iit + 1e6;
                iitt = iitt + 1e6;
                nA = nA - 1e6;
            end

            % read rest of buffer
            q(iit:end, fc) = eq(fread(fid, [nA, 1], dt), B);

            % close file
            fclose(fid);
        end

    % error handler
    catch ne_eo;
        if fid > 0
            fclose(fid);
        end
        rethrow(ne_eo);
    end

% same size values
else

    % single file
    if ischar(sA.FileName)
        fname = {sA.FileName};
        fsize = 1;

    % multiple files
    else
        fname = sA.FileName;
        fsize = szA(end);
        szA(end) = [];
    end

    % size array
    q = false(prod(szA), fsize);
    B = reshape(B, size(q));

    % inner try block (to close file on error)
    try

        % open files
        for fc = 1:fsize
            nA = prod(szA);
            fid = 0;
            if sA.LittleND
                fid = fopen(fname{fc}, 'rb', 'ieee-le');
            else
                fid = fopen(fname{fc}, 'rb', 'ieee-be');
            end

            % check file ID
            if fid < 1
                error( ...
                    'transio:FileNotOpen', ...
                    'Error opening file: ''%s''.', ...
                    fname{fc} ...
                );
            end

            % seek to correct offset
            fseek(fid, sA.IOOffset(fc), -1);

            % while still a large chunk to work on
            iit = 1;
            iitt = 1e6;
            while nA > 1e6

                % read one million entries and work
                q(iit:iitt, fc) = eq(fread(fid, [1e6, 1], dt), B(iit:iitt, fc));
                iit = iit + 1e6;
                iitt = iitt + 1e6;
                nA = nA - 1e6;
            end

            % read rest of buffer
            q(iit:end, fc) = eq(fread(fid, [nA, 1], dt), B(iit:end, fc));

            % close file
            fclose(fid);
        end

    % error handler
    catch ne_eo;
        if fid > 0
            fclose(fid);
        end
        rethrow(ne_eo);
    end
end

% reshape correctly
q = reshape(q, sA.DataDims);

% sub-functions (A = transio, B = transio)
function q = tio_eqto(A, B)

% only valid for single transio
sA = struct(A);
sB = struct(B);
if numel(sA) ~= 1 || ...
    numel(sB) ~= 1
    error( ...
        'transio:BadArgument', ...
        'Call not valid for multiple objects array.' ...
    );
end

% compare sizes
szA = sA.DataDims;
nA = prod(szA);
szB = sB.DataDims;
nB = prod(szB);
if ~isequal(szA, szB) && ...
    nA ~= 1 && ...
    nB ~= 1
    error( ...
        'transio:BadArgument', ...
        'Matrix dimensions must agree.' ...
    );
end

% single value in object A and B (unlikely, but still)
if nA == 1 && ...
    nB == 1
    q = (double(resolve(A)) == double(resolve(B)));

% single value in object A (again unlikely, but what the heck)
elseif nA == 1

    % re-use tio_eqso
    try
        q = tio_eqso(B, resolve(A));
    catch ne_eo;
        rethrow(ne_eo);
    end

% single value in object B
elseif nB == 1

    % re-use tio_eqso
    try
        q = tio_eqso(A, resolve(B));
    catch ne_eo;
        rethrow(ne_eo);
    end

% same size values
else

    % quick fix -> same object?
    if isequal(sA, sB)

        % return with true over size
        q = true(szA);
        return;
    end

    % currently only supported with same number of files
    if (ischar(sA.FileName) && ...
        ~ischar(sB.FileName)) || ...
       (~ischar(sA.FileName) && ...
        ischar(sB.FileName))
        error( ...
            'transio:Unsupported', ...
            'Comparing multi-file transio requires same number of files.' ...
        );
    end

    % get datatypes
    dtA = ['*' sA.DataType];
    dtB = ['*' sB.DataType];

    % single file
    if ischar(sA.FileName)
        fnamesA = {sA.FileName};
        fnamesB = {sB.FileName};
        fsize = 1;

    % multiple files
    else
        fnamesA = sA.FileName;
        fnamesB = sB.FileName;
        fsize = szA(end);
        szA(end) = [];
    end

    % size array
    q = false(prod(szA), fsize);

    % inner try block (to close file on error)
    try

        % open files and check fid's
        for fc = 1:fsize
            nA = prod(szA);
            fidA = 0;
            fidB = 0;
            if sA.LittleND
                fidA = fopen(fnamesA{fc}, 'rb', 'ieee-le');
            else
                fidA = fopen(fnamesA{fc}, 'rb', 'ieee-be');
            end
            if fidA < 1
                error( ...
                    'transio:FileNotOpen', ...
                    'Error opening file: ''%s''.', ...
                    fnamesA{fc} ...
                );
            end
            if sB.LittleND
                fidB = fopen(fnamesB{fc}, 'rb', 'ieee-le');
            else
                fidB = fopen(fnamesB{fc}, 'rb', 'ieee-be');
            end
            if fidB < 1
                error( ...
                    'transio:FileNotOpen', ...
                    'Error opening file: ''%s''.', ...
                    fnamesB{fc} ...
                );
            end
            if fidB == fidA
                fidB = 0;
                error( ...
                    'transio:OSError', ...
                    'Operating system doesn''t allow two filepointers on same file.' ...
                );
            end

            % seek to correct offset
            fseek(fidA, sA.IOOffset(fc), -1);
            fseek(fidB, sB.IOOffset(fc), -1);

            % while still a large chunk to work on
            iit = 1;
            iitt = 1e6;
            while nA > 1e6

                % read one million entries and work
                q(iit:iitt, fc) = eq(fread(fidA, [1e6, 1], dtA), fread(fidB, [1e6, 1], dtB));
                iit = iit + 1e6;
                iitt = iitt + 1e6;
                nA = nA - 1e6;
            end

            % read rest of buffer
            q(iit:end, fc) = eq(fread(fidA, [nA, 1], dtA), fread(fidB, [nA, 1], dtB));

            % close file
            fclose(fidA);
            fidA = 0;
            fclose(fidB);
        end

    % error handler
    catch ne_eo;
        if fidA > 0
            fclose(fidA);
        end
        if fidB > 0
            fclose(fidB);
        end
        rethrow(ne_eo);
    end
end

% reshape
q = reshape(q, sA.DataDims);

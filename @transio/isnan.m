function in = isnan(htio)
% transio::isnan  - returns false for other than single and double
%
% FORMAT:       IN = isnan(obj)
%
% Input fields:
%
%       obj         transio object
%
% Output fields:
%
%       IN          isnan

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

% evaluated only for single/double
if any(strcmp(htio.DataType, {'single', 'double'}))

    % if data probably fits into memory
    pd = prod(htio.DataDims);
    if pd <= 1e7

        % use resolve and then the return value
        in = isnan(resolve(htio));

    % or stepwise (one file)
    elseif ischar(htio.FileName)

        % produce array large enough
        in = false(htio.DataDims);

        % get datatype
        dt = ['*' htio.DataType];

        % modified from subsref
        try

            % open file
            fid = 0;
            if htio.LittleND
                fid = fopen(htio.FileName, 'rb', 'ieee-le');
            else
                fid = fopen(htio.FileName, 'rb', 'ieee-be');
            end

            % check file ID
            if fid < 1
                error( ...
                    'transio:FileNotOpen', ...
                    'Error opening file: ''%s''.', ...
                    htio.FileName ...
                );
            end

            % seek to correct offset
            fseek(fid, htio.IOOffset, -1);

            % while still a large chunk to work on
            iit = 1;
            iitt = 1e6;
            while pd > 1e6

                % read one million entries and work
                in(iit:iitt) = isnan(fread(fid, [1e6, 1], dt));
                iit = iit + 1e6;
                iitt = iitt + 1e6;
                pd = pd - 1e6;
            end

            % read rest of buffer
            in(iit:end) = isnan(fread(fid, [pd, 1], dt));

            % close file
            fclose(fid);

        % error handler
        catch ne_eo;
            if fid > 0
                fclose(fid);
            end
            rethrow(ne_eo);
        end

    % multiple files
    else

        % produce array large enough
        ii = false(prod(htio.DataDims(1:end-1)), htio.DataDims(end));

        % get datatype
        dt = ['*' htio.DataType];

        % modified from subsref
        try

            % open file
            for fc = 1:size(ii, 2)
                fid = 0;
                if htio.LittleND
                    fid = fopen(htio.FileName{fc}, 'rb', 'ieee-le');
                else
                    fid = fopen(htio.FileName{fc}, 'rb', 'ieee-be');
                end

                % check file ID
                if fid < 1
                    error( ...
                        'transio:FileNotOpen', ...
                        'Error opening file: ''%s''.', ...
                        htio.FileName{fc} ...
                    );
                end

                % seek to correct offset
                fseek(fid, htio.IOOffset(fc), -1);

                % while still a large chunk to work on
                iit = 1;
                iitt = 1e6;
                while pd > 1e6

                    % read one million entries and work
                    ii(iit:iitt, fc) = isnan(fread(fid, [1e6, 1], dt));
                    iit = iit + 1e6;
                    iitt = iitt + 1e6;
                    pd = pd - 1e6;
                end

                % read rest of buffer
                ii(iit:end, fc) = isnan(fread(fid, [pd, 1], dt));

                % close file
                fclose(fid);
            end

            % reshape
            ii = reshape(ii, htio.DataDims);

        % error handler
        catch ne_eo;
            if fid > 0
                fclose(fid);
            end
            rethrow(ne_eo);
        end
    end

% otherwise
else

    % return false with correct size
    in = false(htio.DataDims);
end

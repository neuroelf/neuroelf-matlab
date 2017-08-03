function y = resolve(htio)
% transio::resolve  - given entire contents in good dims
%
% FORMAT:       y = resolve(tio);
%
% Input fields:
%
%       tio         transio object
%
% Output fields:
%
%       y           data retrieved from transio object

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

% argument check
if nargin < 1 || ...
    numel(struct(htio)) ~= 1
    error( ...
        'transio:BadSubsRef', ...
        'No S subsref struct given.' ...
    );
end

% try retrieval
if ischar(htio.FileName)
    try
        fid = 0;
        if htio.LittleND
            fid = fopen(htio.FileName, 'rb', 'ieee-le');
        else
            fid = fopen(htio.FileName, 'rb', 'ieee-be');
        end
        if fid < 1
            error('FILE_NOT_OPEN');
        end
        fseek(fid, htio.IOOffset, -1);
        if ftell(fid) ~= htio.IOOffset
            error('FILE_SEEK_ERROR');
        end
        y = reshape(fread(fid, ...
            [prod(htio.DataDims), 1], ['*' htio.DataType]), htio.DataDims);
        fclose(fid);
    catch ne_eo;
        if fid > 0
            fclose(fid);
        end
        error( ...
            'transio:FileReadError', ...
            'Error opening file / reading contents: %s.', ...
            ne_eo.message ...
        );
    end
else
    try
        y = eval([htio.DataType '([])']);
        fname = htio.FileName;
        fsize = htio.DataDims(end);
        szA = prod(htio.DataDims(1:end-1));
        y(szA, fsize) = 0;
        for fc = 1:fsize
            fid = 0;
            if htio.LittleND
                fid = fopen(fname{fc}, 'rb', 'ieee-le');
            else
                fid = fopen(fname{fc}, 'rb', 'ieee-be');
            end
            if fid < 1
                error('FILE_NOT_OPEN');
            end
            fseek(fid, htio.IOOffset(fc), -1);
            if ftell(fid) ~= htio.IOOffset(fc)
                error('FILE_SEEK_ERROR');
            end
            y(:, fc) = reshape(fread(fid, [szA, 1], ['*' htio.DataType]), szA, 1);
            fclose(fid);
        end
        y = reshape(y, htio.DataDims);
    catch ne_eo;
        if fid > 0
            fclose(fid);
        end
        error( ...
            'transio:FileReadError', ...
            'Error opening file / reading contents: %s.', ...
            ne_eo.message ...
        );
    end
end

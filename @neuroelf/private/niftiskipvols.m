function niftiskipvols(f, k)
% niftiskipvols  - skip (additional) k volumes from start in NIFTI file f
%
% FORMAT:       niftiskipvols(f, k)
%
% Input fields:
%
%       f           filename of NIFTI file
%       k           number of volumes to skip from beginning (in addition)

% Version:  v1.0
% Build:    15110321
% Date:     Nov-03 2015, 9:56 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, Jochen Weber
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

% check inputs
if nargin < 2 || ...
   ~ischar(f) || ...
    isempty(f) || ...
   ~isa(k, 'double') || ...
    numel(k) ~= 1 || ...
    isinf(k) || ...
    isnan(k) || ...
    k < 1
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% ensure K is integer
k = floor(k);

% ensure file exists
if exist(f, 'file') ~= 2
    error( ...
        'neuroelf:FileNotFound', ...
        'NIFTI file in f doesn''t exist.' ...
    );
end

% open file
fid = fopen(f, 'r+', 'ieee-le');
if fid < 1
    error( ...
        'neuroelf:FileOpenError', ...
        'Error opening NIFTI file.' ...
    );
end

% get filesize
fseek(fid, 0, 1);
filesize = ftell(fid);
if filesize < 400
    fclose(fid);
    error( ...
        'neuroelf:BadFileContent', ...
        'NIFTI file too short.' ...
    );
end
fseek(fid, 0, -1);

% confirm it's a regular NIFTI file
byte_check = fread(fid, [1, 4], 'uint8=>double');
if ~isequal(byte_check, [92, 1, 0, 0])

    % close file
    fclose(fid);

    % could be big-endian, if so, re-open
    if isequal(byte_check, [0, 0, 1, 92])
        fid = fopen(f, 'r+', 'ieee-be');

    % or not...
    else
        error( ...
            'neuroelf:BadFileContent', ...
            'File is not a NIFTI file.' ...
        );
    end
end

% seek forward to position 40, read dimensions
fseek(fid, 40, -1);
d = fread(fid, [1, 8], 'uint16=>double');

% check D(1)
if d(1) ~= 4
    fclose(fid);
    error( ...
        'neuroelf:BadFileContent', ...
        'NIFTI file is not a 4D volume series.' ...
    );
end

% check d(5) (4-th dim) to be > k
if d(5) <= k
    fclose(fid);
    error( ...
        'neuroelf:BadFileContent', ...
        'NIFTI file has too few volumes to skip additional k volumes.' ...
    );
end

% seek forward to position 72, read bits per pixel header info
fseek(fid, 72, -1);
bitsperpixel = fread(fid, [1, 1], 'uint16=>double');
bytesperpixel = ceil(bitsperpixel / 8);

% seek forward to position 108, read current voxel offset
fseek(fid, 108, -1);
voxoffset = fread(fid, [1, 1], 'single=>double');

% check filesize
compsize = voxoffset + bytesperpixel * prod(d(2:5));
if compsize > filesize
    fclose(fid);
    error( ...
        'neuroelf:BadFileSize', ...
        'NIFTI file is too short.' ...
    );
end

% compute new offset
newoffset = voxoffset + bytesperpixel * prod(d(2:4)) * k;

% ensure that new offset can be represented as a single
if newoffset ~= double(single(newoffset))
    fclose(fid);
    error( ...
        'neuroelf:PrecisionError', ...
        'NIFTI header cannot be patched due to precision error in offset.' ...
    );
end

% re-seek to position 48, write new number of volumes
fseek(fid, 48, -1);
fwrite(fid, uint16(d(5) - k), 'uint16');

% re-seek to position 108, write new offset
fseek(fid, 108, -1);
fwrite(fid, single(newoffset), 'single');

% close file
fclose(fid);

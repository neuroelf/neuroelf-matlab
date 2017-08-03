function binwrite(filename, content, cp)
% binwrite  - writes a binary stream from a char/uint8 array to file
%
% FORMAT:       binwrite(filename, content)
%
% Input Fields:
%       filename    name to a file, preferably absolute path
%       content     char/uint8 array to write
%
% See also asciiwrite

% Version:  v0.9c
% Build:    12012316
% Date:     Jan-23 2012, 4:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2012, Jochen Weber
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
  ~ischar(filename) || ...
    isempty(filename) || ...
    (~ischar(content) && ...
     ~isa(content, 'uint8'))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing arguments in call.' ...
    );
end

% filename mangling check
if ispc
    filename = strrep(filename(:)', '/', filesep);
else
    filename = strrep(filename(:)', '\', filesep);
end

% open file and check for handle
ofp = fopen(filename, 'w');
if ofp < 1

    % try for parent dir existance if cp is true
    if nargin > 2 && ...
        islogical(cp) && ...
        numel(cp) == 1 && ...
        cp
        filepath = fileparts(filename);
        if ~isempty(filepath) && ...
            exist(filepath, 'dir') ~= 7
            try
                mkadir(filepath, '-p');
            catch ne_eo;
                rethrow(ne_eo);
            end
            ofp = fopen(filename, 'w');
        else
            error( ...
                'neuroelf:FileNotWritable', ...
                'Couldn''t write to file: %s.', ...
                filename ...
            );
        end
    else
        error( ...
            'neuroelf:FileNotWritable', ...
            'Couldn''t write to file: %s.', ...
            filename ...
        );
    end
end
frewind(ofp);

% rewind file, write content, and close file
if ischar(content)
    fwrite(ofp, content, 'uchar');
else
    fwrite(ofp, content, 'uint8');
end
fclose(ofp);

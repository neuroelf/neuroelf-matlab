function bincont = binread(filename)
% binread  - reads a binary into one uint8 array
%
% FORMAT:       bincont = binread(filename)
%
% Input fields:
%
%       filename    name to a file, preferably absolute path
%
% Output fields:
%
%       bincont     1xN uint8 array
%
% See also asciiread, binwrite

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% enough arguments ?
if nargin < 1 || ...
   ~ischar(filename) || ...
    isempty(filename)
    error( ...
        'neuroelf:BadArgument',...
        'Bad or missing argument.' ...
    );
end

% file exist check
filename = filename(:)';
if exist(filename, 'file') ~= 2
    error( ...
        'neuroelf:FileNotFound',...
        'File not found: %s',...
        filename(1:min(length(filename), 80)) ...
    );
end

% open file check
fp = fopen(filename, 'r');
if fp < 1
    error( ...
        'neuroelf:FileNotReadable',...
        'File not readable: %s',...
        filename ...
    );
end

% read file contents
bincont = fread(fp, [1, Inf], '*uint8');
fclose(fp);

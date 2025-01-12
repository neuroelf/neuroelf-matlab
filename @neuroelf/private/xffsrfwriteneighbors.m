function xffsrfwriteneighbors(fid, varargin)
% xffsrfwriteneighbors  - write neighbors to SRF file
%
% FORMAT:       xffsrfparseneighbors(fid, neighbors)
%
% Input fields:
%
%       fid         input file fid (fopen)
%       neighbors   Nx2 cell array with neighbor numbers and lists
%
% See also xff, xffsrfparseneighbors
%
% Note: this function uses a compiled function, xffsrfwriteneighborsc.

% Version:  v1.1b
% Build:    25011123
% Date:     Jan-11 2025, 11:45 PM EST
% Author:   Jochen Weber, NeuroElf
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2025, Jochen Weber
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
   ~isa(fid, 'double') || ...
    numel(fid) ~= 1 || ...
   ~isreal(fid) || ...
   ~iscell(varargin{1}) || ...
    isempty(varargin{1}) || ...
    length(size(varargin{1})) ~= 2 || ...
    size(varargin{1}, 2) ~= 2
    error( ...
        'xff:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% get neighbors stream
try
    neighstream = xffsrfwriteneighborsc(varargin{1});
catch ne_eo;
    error( ...
        'xff:MEXError', ...
        'Error compiling neighbors stream: %s.', ...
        ne_eo.message ...
    );
end

% write neighbors
try
    fwrite(fid, neighstream, 'uint32');
catch ne_eo;
    error( ...
        'xff:BadFileOrContent', ...
        'Writing of file failed with: ''%s''.', ...
        ne_eo.message ...
    );
end

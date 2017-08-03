function neighbors = xffsrfparseneighbors(fid, numvtx)
% xffsrfparseneighbors  - parse neighbors from SRF file
%
% FORMAT:       neighbors = xffsrfparseneighbors(fid, numvtx)
%
% Input fields:
%
%       fid         input file fid (fopen)
%       numvtx      number of vertices
%
% Output fields
%
%       neighbors   Nx2 cell array with content
%        {N, 1}     1x1 double, number of neighbors for vertex(N)
%        {N, 2}     1xN double, list of neighbors
%
% See also xff
%
% Note: this function uses a compiled function, xffsrfparseneighborsc.

% Version:  v0.9d
% Build:    14072111
% Date:     Jul-21 2014, 11:35 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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
    isempty(fid) || ...
   ~isreal(fid) || ...
   ~any(fopen('all') == fid(1)) || ...
   ~isa(numvtx, 'double') || ...
    isempty(numvtx) || ...
   ~isreal(numvtx) || ...
    isnan(numvtx(1)) || ...
    isinf(numvtx(1)) || ...
    numvtx(1) < 3
    error( ...
        'xff:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% try
try

    % get current position
    cpos = ftell(fid);

    % get remainder of contents
    vcont = fread(fid, [1, Inf], '*uint32');
    vlen  = length(vcont);

    % sanity check: in a mesh, every vertex MUST have neighbors defined
    if vlen < numvtx
        error('Too few DWORDs left.');
    end

    % get neighbors
    [neighbors, skip] = xffsrfparseneighborsc(vcont, numvtx);

    % seek to good position
    fseek(fid, cpos + 4 * skip, -1);

catch ne_eo;
    error( ...
        'xff:BadFileOrContent', ...
        'Reading of file failed with: ''%s''.', ...
        ne_eo.message ...
    );
end

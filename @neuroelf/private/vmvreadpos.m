function [crd, nrm, col] = vmvreadpos(fid, nrc, usecol)
% vmvreadpos  - read one vertex movie position in a VMV file
%
% FORMAT:       [crd, nrm, col] = vmvreadpos(fid, nrc, usecol)
%
% Input fields:
%
%       fid         1x1 file id (fopen, fread, fclose)
%       nrc         number of coordinates
%       usecol      use colors in VMV
%
% Output fields:
%
%       crd         Cx3 coordinates
%       nrm         Cx3 normals
%       col         Cx4 SRF-like color coding (or empty if unused)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
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

% check arguments
if nargin ~= 3 || ...
   ~isa(fid, 'double') || ...
    numel(fid) ~= 1 || ...
    isinf(fid) || ...
    isnan(fid) || ...
    fid < 1 || ...
    fid ~= fix(fid) || ...
   ~isa(nrc, 'double') || ...
    numel(nrc) ~= 1 || ...
    isinf(nrc) || ...
    isnan(nrc) || ...
    nrc < 1 || ...
    nrc ~= fix(nrc) || ...
    numel(usecol) ~= 1
    error( ...
        'neuroelf:BadArgument', ...
        'This call requires exactly 3 arguments.' ...
    );
end
try
    fpos = 0;
    if usecol
        rwwidth = 7;
    else
        rwwidth = 6;
    end
    fpos = ftell(fid);
    tdat = fread(fid, [rwwidth, nrc], '*single')';
    crd = tdat(:, 1:3);
    nrm = tdat(:, 4:6);
    if usecol
        fseek(fid, fpos, -1);
        tdat = fread(fid, [rwwidth, nrc], '*uint32')';
        col = uint322colcode(tdat(:, 7));
    else
        col = [];
    end
catch ne_eo;
    error( ...
        'neuroelf:InternalError', ...
        'Error reading movie position in file id %d at position %d (%s).', ...
        fid, fpos, ne_eo.message ...
    );
end

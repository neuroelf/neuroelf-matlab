function writeok = vmvwritepos(fid, vpos, usecol)
% vmvwritepos  - write one vertex movie position in a VMV file
%
% FORMAT:       writeok = vmvwritepos(fid, vpos, usecol)
%
% Input fields:
%
%       fid         1x1 file id (fopen, fread, fclose)
%       vpos        struct with at least the fields
%        .Coordinates  Cx3 coordinates
%        .Normals      Cx3 normals
%        .Colors       Cx4 SRF-like color coding (must be present if empty)
%       usecol      use colors in VMV
%
% Output fields:
%
%       writeok     boolean flag whether write succeeded

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
   ~isstruct(vpos) || ...
    numel(vpos) ~= 1 || ...
   ~isfield(vpos, 'Coordinates') || ...
   ~isfield(vpos, 'Normals') || ...
   ~isequal(size(vpos.Coordinates), size(vpos.Normals)) || ...
   ~isfield(vpos, 'Colors') || ...
    numel(usecol) ~= 1 || ...
   (usecol && ...
    (size(vpos.Colors, 1) ~= size(vpos.Coordinates, 1) || ...
     size(vpos.Colors, 2) ~= 4 || ...
     ndims(vpos.Colors) ~= 2))
    error( ...
        'neuroelf:BadArgument', ...
        'This call requires exactly 3 arguments.' ...
    );
end

% with error handling
try
    fpos = 0;
    if usecol
        rwwidth = 7;
    else
        rwwidth = 6;
    end
    nrc = size(vpos.Coordinates, 1);
    fpos = ftell(fid);
    tdat = single(zeros(nrc, rwwidth));
    tdat(:, 1:3) = vpos.Coordinates;
    tdat(:, 4:6) = vpos.Normals;
    writeok = (fwrite(fid, tdat', 'single') == numel(tdat));
    if usecol
        fseek(fid, fpos, -1);
        tdat = fread(fid, [rwwidth, nrc], '*uint32')';
        fseek(fid, fpos, -1);
        tdat(:, 7) = colcode2uint32(vpos.Colors);
        writeok = (fwrite(fid, tdat', 'uint32') == numel(tdat));
    end
catch ne_eo;
    error( ...
        'neuroelf:InternalError', ...
        'Error reading movie position in file id %d at position %d (%s).', ...
        fid, fpos, ne_eo.message ...
    );
end

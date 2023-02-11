function xo = srf_SaveAsOBJ(xo, objfile, talorder)
% SRF::SaveAsOBJ  - create an OBJ file from a surface
%
% FORMAT:       srf.SaveAsOBJ(objfile [, talorder]);
%
% Input fields:
%
%       objfile     filename of OBJ file to write
%       talorder    reorder dimensions to be in TAL space (default: false)

% Version:  v1.1
% Build:    23011721
% Date:     Jan-17 2023, 9:22 PM EST
% Author:   Jochen Weber, NeuroElf.net
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2023, Jochen Weber
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

% check input arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
   ~ischar(objfile) || isempty(objfile)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
[fp, fn, fe] = fileparts(xo.F);
objfile = objfile(:)';
if nargin < 3 || ~islogical(talorder) || numel(talorder) ~= 1
    talorder = false;
end
v = bc.VertexCoordinate;
vn = bc.VertexNormal;
t = bc.TriangleVertex;
if talorder
    v = ones(size(v, 1), 1) * bc.MeshCenter - v(:, [3, 1, 2]);
    vn = -vn(:, [3, 1, 2]);
end

% try to open file for output (binary write, big-endian)
fid = fopen(objfile, 'wb', 'ieee-be');
if fid < 1
    error('neuroelf:xff:fileNotWritable', 'Given filename is not writable.');
end

% write header
nl = char(10);
fwrite(fid, ['# OBJ DataFile', nl], 'char');
fwrite(fid, ['g ', fn, nl], 'char');
fwrite(fid, nl, 'char');

% write vertices
fwrite(fid, sprintf('v %f %f %f\n', v'), 'char');
fwrite(fid, nl, 'char');

% write triangles
fwrite(fid, sprintf('f %d %d %d\n', t'), 'char');
fwrite(fid, nl, 'char');

% close file
fclose(fid);

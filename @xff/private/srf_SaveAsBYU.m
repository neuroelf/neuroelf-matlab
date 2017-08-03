function xo = srf_SaveAsBYU(xo, byufile)
% SRF::SaveAsBYU  - create a BYU file from a surface
%
% FORMAT:       srf.SaveAsBYU(byufile);
%
% Input fields:
%
%       byufile     filename of BYU file to write

% Version:  v1.1
% Build:    16021121
% Date:     Feb-11 2016, 9:20 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
   ~ischar(byufile) || isempty(byufile)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
byufile = byufile(:)';

% try to open file for output (binary write, big-endian)
fid = fopen(byufile, 'w');
if fid < 1
    error('neuroelf:xff:fileNotWritable', 'Given filename is not writable.');
end

% write header
fwrite(fid, [sprintf(' %6d', [1, bc.NrOfVertices, bc.NrOfTriangles, 3 * bc.NrOfTriangles, 0]) char(10)], 'char');
fwrite(fid, sprintf(' %6d %6d%c', [1, bc.NrOfTriangles, 10]));

% write coordinates
crd = bc.VertexCoordinate';
crd(4, :) = 10;
fwrite(fid, sprintf('%12.5E%12.5E%12.5E%c', crd(:)'));

% write triangle vertices
tri = bc.TriangleVertex';
tri(3, :) = -tri(3, :);
tri(4, :) = 10;
fwrite(fid, sprintf(' %6d %6d %6d%c', tri(:)'));

% close file
fclose(fid);

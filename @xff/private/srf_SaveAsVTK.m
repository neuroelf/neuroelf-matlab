function xo = srf_SaveAsVTK(xo, vtkfile)
% SRF::SaveAsVTK  - create a VTK file from a surface
%
% FORMAT:       srf.SaveAsVTK(vtkfile);
%
% Input fields:
%
%       vtkfile     filename of VTK file to write

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:53 PM EST
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
   ~ischar(vtkfile) || isempty(vtkfile)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
vtkfile = vtkfile(:)';

% try to open file for output (binary write, big-endian)
fid = fopen(vtkfile, 'wb', 'ieee-be');
if fid < 1
    error('neuroelf:xff:fileNotWritable', 'Given filename is not writable.');
end

% write header
fwrite(fid, ['# vtk DataFile Version 3.0' char(10)], 'char');
fwrite(fid, ['vtk output' char(10)], 'char');
fwrite(fid, ['BINARY' char(10)], 'char');
fwrite(fid, ['DATASET POLYDATA' char(10)], 'char');

% write vertices
numvert = bc.NrOfVertices;
fwrite(fid, [sprintf('POINTS %d float', numvert) char(10)], 'char');
fwrite(fid, (bc.VertexCoordinate' - bc.MeshCenter' * ones(1, numvert)), 'single');
fwrite(fid, char(10), 'char');

% write triangles
numtri = bc.NrOfTriangles;
fwrite(fid, [sprintf('POLYGONS %d %d', numtri, 4 * numtri) char(10)], 'char');
fwrite(fid, [3 * ones(1, numtri); bc.TriangleVertex' - 1], 'int32');
fwrite(fid, char(10), 'char');

% write footer
fwrite(fid, [sprintf('CELL_DATA %d', numtri) char(10)], 'char');
fwrite(fid, [sprintf('POINT_DATA %d', numvert) char(10)], 'char');

% close file
fclose(fid);

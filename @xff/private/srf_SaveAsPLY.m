function xo = srf_SaveAsPLY(xo, plyfile, opts)
% SRF::SaveAsPLY  - create a PLY file from a surface
%
% FORMAT:       srf.SaveAsVTK(plyfile [, opts]);
%
% Input fields:
%
%       plyfile     filename of PLY file to write
%       opts        optional settings
%        .binary    write binary PLY file (default: false)
%
% No output fields.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, 2016, Jochen Weber
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
   ~ischar(plyfile) || isempty(plyfile)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
sf = xo.F;
if isempty(sf)
    sf = 'unsaved.srf';
end
plyfile = plyfile(:)';
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'binary') || ~islogical(opts.binary) || numel(opts.binary) ~= 1
    opts.binary = false;
end

% binary
if opts.binary
    disp('Not yet implemented.');
    return;
end

% try to open file for output (binary write, big-endian)
fid = fopen(plyfile, 'w');
if fid < 1
    error('neuroelf:xff:fileNotWritable', 'Given filename is not writable.');
end

% write header
lf = char(10);
fwrite(fid, ['ply' lf]);
fwrite(fid, ['format ascii 1.0' lf]);
fwrite(fid, ['comment sourcefile ' sf lf]);
fprintf(fid, 'element vertex %d\n', size(bc.VertexCoordinate, 1));
fprintf(fid, 'property float x\nproperty float y\nproperty float z\n');
fprintf(fid, 'element face %d\n', size(bc.TriangleVertex, 1));
fwrite(fid, ['property list uchar int vertex_index' lf]);
fwrite(fid, ['end_header' lf]);

% write data
fprintf(fid, '%g %g %g\n', 128 - bc.VertexCoordinate(:, [3, 1, 2])');
fprintf(fid, '%d %d %d %d\n', ([3 .* ones(size(bc.TriangleVertex, 1), 1), bc.TriangleVertex - 1])');

% close file
fclose(fid);

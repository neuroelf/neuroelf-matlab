function xo = srf_SaveAsJVX(xo, jvxfile)
% SRF::SaveAsJVX  - save SRF to a JavaView JVX file
%
% FORMAT:       srf.SaveAsJVX(jvxfile);
%
% Input fields:
%
%       jvxfile     filename of JVX file to write
%
% No output fields.

% Version:  v1.1
% Build:    16021120
% Date:     Feb-11 2016, 8:57 PM EST
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
   ~ischar(jvxfile) || isempty(jvxfile)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
jvxfile = jvxfile(:)';

% try to open output for writing
jf = fopen(jvxfile, 'w');
if jf < 1
    error('neuroelf:xff:fileNotWritable', 'Destination file not writable.');
end

% print header
fwrite(jf, ['<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>' char(10)]);
fwrite(jf, ['<!DOCTYPE jvx-model SYSTEM "http://www.javaview.de/rsrc/jvx.dtd">' char(10)]);
fwrite(jf, ['<jvx-model>' char([10,9]) '<meta generator="bvsrf2jvx"/>' char([10,9])]);
fwrite(jf, ['<title>' xo.F '</title>' char([10,9]) '<geometries>' char([10,9,9])]);
fwrite(jf, ['<geometry name="' xo.F '">' char(10)]);

% write points/normals
numvert = bc.NrOfVertices;
fwrite(jf, sprintf('\t\t\t<pointSet point="hide" color="hide" dim=3>\n'));
fwrite(jf, sprintf('<points num=%d>\n', numvert));
fprintf(jf, '\t<p>%8.5f %8.5f %8.5f</p>\n', bc.VertexCoordinate' - bc.MeshCenter' * ones(1, numvert));
fwrite(jf, sprintf('</points>\n'));
fwrite(jf, sprintf('<colors num=%d>\n', numvert));
oc0 = fix(255 * bc.ConvexRGBA(1:3));
oc1 = fix(255 * bc.ConcaveRGBA(1:3));
vcol = bc.VertexColor;
col1 = find(~isnan(vcol(:, 1)) & vcol(:, 1) == 0);
col2 = find(~isnan(vcol(:, 1)) & vcol(:, 1) == 1);
vcol(col1, 2:4) = ones(numel(col1), 1) * oc0;
vcol(col2, 2:4) = ones(numel(col2), 1) * oc1;
vcol(vcol < 0) = 0;
vcol(vcol > 255) = 255;
vcol = vcol(:, 2:4);
fprintf(jf, '\t<c>%d %d %d</c>\n', vcol');
fwrite(jf, sprintf('</colors>\n'));
sn = bc.VertexNormal;
fwrite(jf, sprintf('<normals num=%d>\n', numvert));
fprintf(jf,'\t<p>%9.6f %9.6f %9.6f</p>\n', sn');
fwrite(jf, sprintf('\t<length>1</length>\n'));
fwrite(jf, sprintf('</normals>\n'));
fwrite(jf, sprintf('\t\t\t</pointSet>\n'));

% write faces
fwrite(jf, sprintf('\t\t\t<faceSet color="show" colorFromPoints="show" edge="hide" face="show" colorSmooth="show">\n'));
sf = bc.TriangleVertex - 1;
fwrite(jf, sprintf('<faces num=%d>\n', size(sf, 1)));
fprintf(jf, '\t<f>%d %d %d</f>\n', sf');
fwrite(jf, sprintf('</faces>\n'));
fwrite(jf, sprintf('\t\t\t</faceSet>\n'));

% write footer
fwrite(jf, [char([9,9]) '</geometry>' char([10,9])]);
fwrite(jf, ['</geometries>' char(10) '</jvx-model>' char(10)]);

% close file
fclose(jf);

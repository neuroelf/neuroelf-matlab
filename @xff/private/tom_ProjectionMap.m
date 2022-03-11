function smp = tom_ProjectionMap(xo, opts)
% TOM::ProjectionMap  - create SMP with projection maps
%
% FORMAT:       smp = tom.ProjectionMap([opts]);
%
% Input fields:
%
%       opts        optional settings
%
% Output fields:
%
%       smp         SMP object with projection maps

% Version:  v1.1
% Build:    21120216
% Date:     Dec-02 2021, 4:32 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2021, Jochen Weber
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
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'tom')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end

% access TOM properties
tom = xo.C;
ctm = tom.CornerTexVtxMap;
tam = tom.TexVertAMap;
tv = tom.TriangleVertex;

% generate map content
nummaps = max(tam);
m = zeros(tom.NrOfVertices, nummaps);

% iterate over vertices
for c = size(ctm,1):-1:1
    ct = ctm(c, 1);
    cti = ctm(c, 2);
    tnum = tam(ctm(c, 3));
    m(tv(ct, cti), tnum) = m(tv(ct, cti), tnum)+1;
end

% scale results
m = log(1+m);

% generate new SMP object
smp = xff('new:smp');
smpc = smp.C;

% store results
smpc.Map = repmat(smpc.Map(1), 1, nummaps);
smpc.NrOfMaps = nummaps;
smpc.NrOfVertices = tom.NrOfVertices;
for c = 1:nummaps
    smpc.Map(c).SMPData = single(m(:, c));
    smpc.Map(c).LowerThreshold = 1;
    smpc.Map(c).UpperThreshold = 4;
    smpc.Map(c).Name = sprintf('Texture %d', c);
end
smp.C = smpc;

function [y, c] = samplefmrspace(x, c, fmr, trf, imeth, conly)
% samplefmrspace  - samples a data slab at a given set of coordinates
%
% FORMAT:       [y, c] = samplefmrspace(x, c, fmr [, trf [, meth [, conly]])
%
% Input fields:
%
%       x       data in FMR resolution (first matching dims, so STC works)
%       c       Cx3 list of coordinates to sample at, or a cell array
%               suitable to pass as args into ndgrid (in BVsystem
%               coordinates, axes order as TAL, so that BVsys = 128 - Tal)
%       fmr     FMR (or DMR) project object
%       trf     optional 1xT cell array of transformation files,
%               e.g. {ia, fa, acpc, tal}
%       meth    optional, 'cubic', 'lanczos3', {'linear'}, 'nearest'
%       conly   if given, only return coordinates of sampling
%
% Output fields:
%
%       y       Cx1 data values at the C coordinates
%       c       Cx3 coordinates for flexinterpn sampling

% Version:  v0.9b
% Build:    10062313
% Date:     Jun-23 2010, 7:51 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
if nargin < 3 || ...
    isempty(x) || ...
   (~isnumeric(x) && ~istransio(x)) || ...
    ndims(x) < 3 || ...
    numel(fmr) ~= 1 || ...
   (~isxff(fmr, 'fmr') && ...
    ~isxff(fmr, 'dmr')) || ...
   ((~iscell(c) || ...
      numel(c) ~= 3 || ...
      ~isa(c{1}, 'double') || ...
      ~isa(c{2}, 'double') || ...
      ~isa(c{3}, 'double') || ...
      any(isinf([c{1}(:); c{2}(:); c{3}(:)]) | isinf([c{1}(:); c{2}(:); c{3}(:)]))) && ...
    (~isa(c, 'double') || ...
      ndims(c) ~= 2 || ...
      size(c, 2) ~= 3))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
szknown = false;
if iscell(c)
    szknown = [numel(c{1}), numel(c{2}), numel(c{3})];
    [cx{1:3}] = ndgrid(c{1}, c{2}, c{3});
    c = [cx{1}(:), cx{2}(:), cx{3}(:)];
end
if nargin < 4
    trf = {};
end
if any(isxff(trf))
    trfc = cell(1, numel(trf));
    for ci = 1:numel(trf)
        trfc{ci} = trf(ci);
    end
    trf = trfc;
end
if ~iscell(trf)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad trf argument.' ...
    );
end
if nargin < 5 || ...
   ~ischar(imeth) || ...
    isempty(imeth)
    imeth = 'linear';
else
    imeth = lower(imeth(:)');
end
if nargin < 6 || ...
   ~islogical(conly) || ...
    numel(conly) ~= 1
    conly = false;
end

% get dmr/fmr CoordinateFrame
fmrc = fmr.CoordinateFrame();

% get dimensions to work with
sd = size(x);
td = [fmrc.DimX, fmrc.DimY, fmrc.DimZ];
si = 1;
% ti = [1, 2, 3];
while si < numel(sd) && ...
    (sd(si) ~= td(1))
    si = si + 1;
end
if si > numel(sd)
    error( ...
        'neuroelf:BadArgument', ...
        'Dimensions of dataslab and FMR mismatch.' ...
    );
end
% ti(1) = si;
si = si + 1;
while si < numel(sd) && ...
    (sd(si) ~= td(2))
    si = si + 1;
end
if si > numel(sd)
    error( ...
        'neuroelf:BadArgument', ...
        'Dimensions of dataslab and FMR mismatch.' ...
    );
end
% ti(2) = si;
si = si + 1;
while si < numel(sd) && ...
    (sd(si) ~= td(3))
    si = si + 1;
end
if si > numel(sd)
    error( ...
        'neuroelf:BadArgument', ...
        'Dimensions of dataslab and FMR mismatch.' ...
    );
end
% ti(3) = si;

% transform target TAL to source ACPC if necessary
for ci = numel(trf):-1:1
    if isxff(trf{ci}, 'tal')
        c = acpc2tal(c, trf{ci}, true);
        trf(ci) = [];
    end
end

% subtract center coordinate
if fmr.CoordinateSystem ~= 1
    c = 255 - c;
end
c = c - 127.5;

% add forth column ...
c(:, 4) = 1;

% inspect list of transformations with, supposedly, no tal transform
for ci = numel(trf):-1:1
    if ~isxff(trf{ci}, 'trf')
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid trf argument (particle).' ...
        );
    end
    if trf{ci}.TransformationType > 1
        c(:, [2, 3, 1, 4]) = c(:, [2, 3, 1, 4]) * trf{ci}.TFMatrix';
    elseif trf{ci}.AlignmentStep == 2
        c(:, [2, 3, 1, 4]) = c(:, [2, 3, 1, 4]) * inv(trf{ci}.TFMatrix)';
    else
        c = c(:, [2, 3, 1, 4]) * inv(trf{ci}.TFMatrix)';
    end
end

% prepare sampling
c = [c(:, 1) ./ fmrc.ResX + (fmrc.DimX + 1) / 2, ...
     c(:, 2) ./ fmrc.ResY + (fmrc.DimY + 1) / 2, ...
    -c(:, 3) ./ fmrc.ResZ + (fmrc.DimZ + 1) / 2];

% flip slice order (experimental !!)
if ~isempty(trf) && ...
    trf{1}.CreateFMR3DMethod == 2
    c(:, 3) = (fmrc.DimZ + 1) - c(:, 3);
end

% only cooordinates
if conly
    y = zeros(size(c, 1), 1);
    return;
end

% sample x at c
if numel(sd) == 3
    y = flexinterpn_method(x, c, 0, imeth);
else
end

% remove NaNs
y(isnan(y(:))) = 0;

% reshape ?
if numel(szknown) == 3
    y = reshape(y, szknown);
end

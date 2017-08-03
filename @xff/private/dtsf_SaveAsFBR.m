function xo2 = dtsf_SaveAsFBR(xo, fbrfilename, bogus)
% DTSF::SaveAsFBR  - convert a DTI Studio FiberDat into FBR file
%
% FORMAT:       [fbr] = dtsf.SaveAsFBR(fbrfilename [, bogus]);
%
% Input fields:
%
%       fbrfilename FBR filename
%       bogus       if given and true, add bogus group
%
% Output fields:
%
%       fbr         FBR object

% Version:  v1.1
% Build:    16020215
% Date:     Feb-02 2016, 3:20 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'dtsf') || ...
   ~ischar(fbrfilename) || isempty(fbrfilename)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% create empty FBR file in memory
xo2 = xff('new:fbr');

% check file writability
try
    aft_SaveAs(xo2, fbrfilename);
    nbc = xo2.C;
catch xfferror
    neuroelf_lasterr(xfferror);
    delete(xo2);
    error('neuroelf:xff:fileNotWritable', 'File not writable: ''%s''.', fbrfilename);
end

% get source space
bc = xo.C;
pix1 = bc.DimX;
pix2 = bc.DimY;
pix3 = bc.DimZ;
res1 = bc.ResX;
res2 = bc.ResY;
res3 = bc.ResZ;
sori = bc.SliceOrientation;
if bc.SliceSequencing
    res3 = -res3;
end
cnt1 = (res1 * pix1) / 2;
cnt2 = (res2 * pix2) / 2;
cnt3 = (res3 * pix3) / 2;

% create coordinate conversion matrix, slice orientation
switch (sori)

    % coronal slicing (along BV's X axis)
    case 0
        tmat = [ ...
               0,    0, res3, -cnt3 ; ...
               0, res2,    0, -cnt2 ; ...
            res1,    0,    0, -cnt1 ; ...
               0,    0,    0,     1];

    % axial slicing (along BV's Y axis)
    case 1
        tmat = [ ...
               0, res2,    0, -cnt2 ; ...
               0,    0, res3, -cnt3 ; ...
            res1,    0,    0, -cnt1 ; ...
               0,    0,    0,     1];

    % sagittal slicing (along BV's Z axis)
    case 2
        tmat = [ ...
            res1,    0,    0, -cnt1 ; ...
               0, res2,    0, -cnt2 ; ...
               0,    0, res3, -cnt3 ; ...
               0,    0,    0,     1];

    otherwise
        delete(xo2);
        error('neuroelf:xff:invalidField', 'Invalid SliceOrientation field content.');
end

% create new fiber struct
nfb.NrOfPoints = 1;
nfb.FiberPoints = [0 0 0];

% put new fibers into FBR
numf = bc.NrOfFibers;
nbc.Group(1).NrOfFibers = numf;
nbc.Group(1).Fiber = nfb(ones(numf, 1));

% add bogus group ?
if nargin > 2 && (islogical(bogus) || isnumeric(bogus)) && ~isempty(bogus) && bogus(1)
    nbc.NrOfGroups = 2;
    nbc.Group(2).Name = 'BogusGroup';
    nbc.Group(2).Visible = 1;
    nbc.Group(2).Animate = 0;
    nbc.Group(2).Thickness = 0.3;
    nbc.Group(2).Color = [0 0 0];
    nbc.Group(2).NrOfFibers = 1;
    nbc.Group(2).Fiber = nfb;
    nbc.Group(2).Fiber.FiberPoints = [-128 -128 -128];
end

% loop over fibers in FiberDat
for fc = 1:numf
    fiber = bc.Fibers(fc);
    fcoor = fiber.Coord;
    nfb.NrOfPoints = size(fcoor, 1);
    fcoor(:, 4) = 1;
    fcoort = (tmat * fcoor')';
    nfb.FiberPoints = fcoort(:, 1:3) + 128;
    nbc.Group(1).Fiber(fc) = nfb;
end

% put back into content
xo2.C = nbc;

% save to file
aft_SaveAs(xo2, fbrfilename);

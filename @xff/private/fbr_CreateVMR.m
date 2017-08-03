function vmr = fbr_CreateVMR(xo, fgroups, addc, vres, colcode, expand)
% FBR::CreateVMR  - create a VMR from a fibre file
%
% FORMAT:       vmr = fbr.CreateVMR([fgroups, addcrd, vmrres, ccode, ex]);
%
% Input fields:
%
%       fgroups     list of fiber groups to write to VMR (all)
%       addcrd      1x3 double, added to coordinate of fibers ([0, 0, 0])
%       vmrres      VMR resolution (1)
%       ccode       color code (240)
%       expand      number of voxels to expand fibers with (0)
%
% Output fields:
%
%       vmr         VMR object with fibers as voxels

% Version:  v1.1
% Build:    16020309
% Date:     Feb-03 2016, 9:59 AM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fbr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 6 || ~isnumeric(expand) || numel(expand) ~= 1 || ...
    isinf(expand) || isnan(expand) || fix(expand) ~= expand || expand < 1
    expand = 0;
else
    expand = min(3, expand);
end
if nargin < 5 || ~isnumeric(colcode) || numel(colcode) ~= 1 || isinf(colcode) || ...
    isnan(colcode) || fix(colcode) ~= colcode || colcode < 226 || colcode > 255
    colcode = uint8(240);
else
    colcode = uint8(colcode);
end
if nargin < 4 || ~isa(vres, 'double') || numel(vres) ~= 1 || ...
    (vres ~= 1 && vres ~= 0.5 && vres ~= 0.25 && vres ~= 0.2)
    vres = 0.5;
end
if nargin < 3 || ~isa(addc, 'double') || numel(addc) ~= 3 || ...
    any(isinf(addc) | isnan(addc) | addc < -256 | addc > 512)
    addc = [0, 0, 0];
else
    addc = addc(:)';
end
if nargin < 2 || ~isa(fgroups, 'double') || isempty(fgroups) || numel(fgroups) ~= length(fgroups) || ...
    any(isinf(fgroups) | isnan(fgroups) | fgroups < 1)
    fgroups = [];
else
    fgroups = fgroups(:)';
end

% get content
bc = xo.C;

% build expand list
exps = (1 + 2*expand) ^ 3;
expm = zeros(exps, 3);
exps = 0;
for ecx = 0:(2*expand)
    for ecy = 0:(2*expand)
        for ecz = 0:(2*expand)
            exps = exps + 1;
            expm(exps, :) = [ecx, ecy, ecz];
        end
    end
end

% default rvalue
vmr = xff('new:vmr');
vmrc = vmr.C;

% get group list right
fgroups(fgroups > bc.NrOfGroups) = [];
if isempty(fgroups)
    fgroups = 1:bc.NrOfGroups;
end

% get geometry values
ffac = 1 / vres;
fcub = 256 * ffac;
vmrc.SliceThickness = vres;
vmrc.VoxResX = vres;
vmrc.VoxResY = vres;
vmrc.VoxResZ = vres;
vmrc.VoxResInTalairach = 1;
vmrc.FramingCube = ffac * 256;
vmrc.NRows = ffac * 256;
vmrc.NCols = ffac * 256;

% clear VMR data
vmrc.VMRData = uint8([]);

% iterate over groups
pcount = 0;
tpts = zeros(0, 3);
for gc = fgroups(:)'

    % get group
    fg = bc.Group(gc);

    % iterate over fibers
    for fc = 1:fg.NrOfFibers

        % get fiber and points
        fib = fg.Fiber(fc);
        fps = fib.FiberPoints;

        % iterate over points
        pnum = size(fps, 1);
        for pc = 1:pnum

            % get system coordinate
            fpt = round(ffac * (fps(pc,:) + addc));

            % only set point if coordinate is valid
            if all(fpt >= 1 & fpt <= fcub)
                pcount = pcount + 1;
                if size(tpts, 1) < pcount
                    tpts(pcount:pcount+8191,:) = 0;
                end
                tpts(pcount,:) = fpt;
            end
        end
    end
end

% truncate target points
if size(tpts, 1) > pcount
    tpts(pcount+1:end, :) = [];
end

% get framing cube and create VMR data array
fcmin = min(tpts) - 1;
fcmax = max(tpts);
fcsiz = fcmax - fcmin + 2 * expand;
vmrd = uint8(0);
vmrd(fcsiz(1), fcsiz(2), fcsiz(3)) = uint8(0);

% fill VMR data
for pc = 1:pcount
    fpt = tpts(pc,:) - fcmin;

    % expanding
    for ec = 1:exps
        ept = fpt + expm(ec,:);
        vmrd(ept(1), ept(2), ept(3)) = colcode;
    end
end

% update VMR object
vmrc.DimX = fcsiz(1);
vmrc.DimY = fcsiz(2);
vmrc.DimZ = fcsiz(3);
vmrc.OffsetX = fcmin(1);
vmrc.OffsetY = fcmin(2);
vmrc.OffsetZ = fcmin(3);

% put data back into object
vmrc.VMRData = vmrd;
vmr.C = vmrc;

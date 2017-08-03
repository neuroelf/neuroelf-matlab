function sfh = fif_WriteSFH(xo, sfhname)
% FIF::WriteSFH  - write SFH file
%
% FORMAT:       sfh = fif.CreateCTC(sfhname)
%
% Input fields:
%
%       sfhname     name of SFH output file
%
% Output fields:
%
%       sfh         the created object

% Version:  v1.1
% Build:    16021216
% Date:     Feb-12 2016, 4:24 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fif') || ...
   ~ischar(sfhname) || isempty(sfhname)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
sfhname = sfhname(:)';

% get object
bc = xo.C;

% create CTC in memory
sfh = xff('new:sfh');
sfc = sfh.C;

% get fif shortcut
fif = bc.FIFStructure;

% any dig points ?
iel = any(fif.Lookup == 213);
if ~iel
    delete(sfh);
    error('neuroelf:xff:badInputFile', 'FIF file does not contain DigPoint tags.');
end

% get all DigPoints
dps = fif_Value(xo, 'DigPoint');

% parse dig points
cfid = zeros(1, 0);
chpi = zeros(1, 0);
cext = zeros(1, 0);
for pc = 1:numel(dps)
    switch lower(dps(pc).Kind)
        case 'headfiducial'
            if ~any(dps(pc).Ident == [1, 2, 3])
                delete(sfh);
                error('neuroelf:xff:unexpectedToken', 'Only three major fiducials supported.');
            end
            cfid(dps(pc).Ident) = pc;
        case 'hpipoint'
            chpi(dps(pc).Ident) = pc;
        case 'extrapoint'
            cext(dps(pc).Ident) = pc;
        otherwise
            warning('neuroelf:xff:unexpectedToken', ...
                'Unknown kind of DigPoint %d for point %d.', dps(pc).Kind, pc);
    end
end
if numel(cfid) ~= 3 || any([cfid, chpi, cext] == 0)
	delete(sfh);
    error('neuroelf:xff:invalidFile', 'Invalid organization of fiducial points.');
end

% build fiducial list
fstr = struct;
fstr.Fid_T9  = [1000 * dps(cfid(1)).Coord, 3, 255, 128, 255];
fstr.Fid_Nz  = [1000 * dps(cfid(2)).Coord, 3, 255, 128, 255];
fstr.Fid_T10 = [1000 * dps(cfid(3)).Coord, 3, 255, 128, 255];
for fc = 1:numel(chpi)
    fstr.(sprintf('HPI_%d', fc + 3)) = [1000 * dps(chpi(fc)).Coord, 2, 0, 255, 0];
end
for fc = 1:numel(cext)
    fstr.(sprintf('Extra_%d', fc + 3 + numel(chpi))) = ...
        [1000 * dps(cext(fc)).Coord, 2, 0, 255, 255];
end

% set number of points and save
sfc.Fiducials = fstr;
sfc.NrOfPoints = numel(fieldnames(fstr));
sfh.C = sfc;
try
    aft_SaveAs(sfh, sfhname);
catch xfferror
    warning('neuroelf:xff:saveFailed', 'Saving SFH file failed: %s.', xfferror.message);
end

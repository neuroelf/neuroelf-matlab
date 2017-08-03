function ia = initialalignment(vmr, fmr)
% initialalignment  - perform the automatic IA
%
% FORMAT:       ia = initialalignment(vmr, fmr)
%
% Input fields:
%
%       vmr         VMR object
%       fmr         FMR object
%
% Output fields:
%
%       ia          TRF object with IA matrix

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:22 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
if nargin ~= 2 || ...
    numel(vmr) ~= 1 || ...
   ~isxff(vmr, 'vmr') || ...
    numel(fmr) ~= 1 || ...
   ~isxff(fmr, {'fmr', 'vmr'})
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% generate ia
ia = xff('new:trf');

% get coordinate frames
fmrc = fmr.CoordinateFrame;
vmrc = vmr.CoordinateFrame;

% coordinate systems
fmrs = fmr.CoordinateSystem;
vmrs = vmr.CoordinateSystem;
if fmrs ~= 1
    % TODO
end
if vmrs ~= 1
    % TODO
end

% convention
fmrrad = strcmpi(fmr.Convention(1), 'r');
if fmrrad == (vmr.Convention ~= 0)
    % TODO
end

% compute IA matrix
tfm = inv(vmrc.Trans) * fmrc.Trans;
tfm(4, :) = [0, 0, 0, 1];

% other interpretation of systems
tfm(:, 3) = -tfm(:, 3);

% some other settings
ia.CoordinateSystem = vmrc;
if isxff(fmr, 'fmr')
    ia.NSlicesFMRVMR = fmr.NrOfSlices;
else
    slcdiff = diff([ ...
        fmr.Slice1CenterX, fmr.Slice1CenterY, fmr.Slice1CenterZ; ...
        fmr.SliceNCenterX, fmr.SliceNCenterY, fmr.SliceNCenterZ]);
    ia.NSlicesFMRVMR = round(sqrt(sum(slcdiff .* slcdiff)) ./ ...
        (fmr.SliceThickness + fmr.GapThickness));
end
ia.SlThickFMRVMR = fmr.SliceThickness;
ia.SlGapFMRVMR = fmr.GapThickness;
ia.CreateFMR3DMethod = 3;
ia.AlignmentStep = 1;

% vmr transformations applied
if vmr.NrOfPastSpatialTransformations > 0
else
    ia.ExtraVMRTransf = 0;
    ia.ExtraVMRTransfValues = eye(4);
end
ia.SourceFile = fmr.FilenameOnDisk;
ia.TargetFile = vmr.FilenameOnDisk;

% store transformation matrix
ia.TFMatrix = tfm;

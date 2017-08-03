function nii = analyzeto4dnii(anafiles, opts)
% analyzeto4dnii  - read in several 3D analyze files and create a 4D NII
%
% FORMAT:       nii = analyzeto4dnii(anafiles [, opts])
%
% Input fields:
%
%       anafiles    list of analyze (3D) files
%       opts        optional settings
%        .outtype   override datatype (Matlab types)
%        .resample  flag, allow orientations to be different (false)
%        .resmeth   resampling method, ('cubic', see flexinterpn_method)
%        .trans     transform values in 4th dim, {'none'}, 'psc', 'z'
%
% Output fields:
%
%       nii         4D NIftI object with data

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:13 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
if nargin < 1 || ...
   ~iscell(anafiles) || ...
    numel(anafiles) < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
anafiles = anafiles(:)';
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'outtype') || ...
   ~ischar(opts.outtype) || ...
    isempty(opts.outtype) || ...
   ~any(strcmpi(opts.outtype(:)', ...
        {'double',  'int8',  'int16',  'int32', ...
         'single', 'uint8', 'uint16', 'uint32'}))
    opts.outtype = '';
else
    opts.outtype = lower(opts.outtype(:)');
end
if ~isfield(opts, 'resample') || ...
   ~islogical(opts.resample) || ...
    numel(opts.resample) ~= 1
    opts.resample = false;
end
if ~isfield(opts, 'resmeth') || ...
   ~ischar(opts.resmeth) || ...
    isempty(opts.resmeth) || ...
    numel(opts.resmeth) > 16
    resmeth = 'cubic';
else
    resmeth = lower(opts.resmeth(:)');
end
if ~isfield(opts, 'trans') || ...
   ~ischar(opts.trans) || ...
    isempty(opts.trans) || ...
   ~any(lower(opts.trans(1)) == 'npz')
    opts.trans = 'n';
else
    opts.trans = lower(opts.trans(1));
end

% check files
for fc = 1:numel(anafiles)
    try
        if ischar(anafiles{fc}) && ...
            numel(anafiles{fc}) > 4 && ...
            strcmpi(anafiles{fc}(end-3:end), '.img')
            anafiles{fc} = regexprep(anafiles{fc}(:)', ...
                '\.img$', '.hdr', 'preservecase');
        end
        anafiles{fc} = xff(anafiles{fc});
        if ~isxff(anafiles{fc}, 'hdr')
            error('Not an Analyze file.');
        end
    catch ne_eo;
        clearxffobjects(anafiles);
        error( ...
            'neuroelf:FileReadError', ...
            'Error reading file %d: %s.', ...
            fc, ne_eo.message ...
        );
    end
end

% check resolution, etc.
aif = anafiles{1}.Info;
l = aif.Layout;
f = false(numel(anafiles), 1);
fs = zeros(numel(anafiles), 1);
s = zeros(numel(anafiles), 2);
l(numel(anafiles), end) = 0;
for fc = 1:numel(anafiles)
    ai = anafiles{fc}.Info;
    l(fc, :) = ai.Layout;
    f(fc) = ai.FourD;
    fs(fc) = ai.FourDSize;
    s(fc, :) = ai.Scaling;
end

% only one type of input image!
if ~all(f == f(1))
    clearxffobjects(anafiles);
    error( ...
        'neuroelf:BadArgument', ...
        'Input files must be all 3D or all 4D.' ...
    );
end

% no difference in scaling supported
if any(any(diff(s)))
    clearxffobjects(anafiles);
    error( ...
        'neuroelf:BadArgument', ...
        'Input files must not be scaled differently.' ...
    );
end

% difference only allowed with resampling
if any(any(diff(l))) && ...
   ~opts.resample
    clearxffobjects(anafiles);
    error( ...
        'neuroelf:OrientationMismatch', ...
        'Analyze files must match in spatial layout and orientation.' ...
    );
end

% copy first object
nii = anafiles{1}.CopyObject;

% create voxel data to hold result
if isempty(opts.outtype)
    nii.VoxelData = eval([aif.DataType '(0)']);
else
    nii.VoxelData = eval([opts.outtype '(0)']);
    switch (opts.outtype)
        case {'double'}
            nii.ImgDim.DataType = 64;
            nii.ImgDim.BitsPerPixel = 64;
        case {'int16'}
            nii.ImgDim.DataType = 4;
            nii.ImgDim.BitsPerPixel = 16;
        case {'int32'}
            nii.ImgDim.DataType = 8;
            nii.ImgDim.BitsPerPixel = 32;
        case {'int8'}
            nii.ImgDim.DataType = 256;
            nii.ImgDim.BitsPerPixel = 8;
        case {'single'}
            nii.ImgDim.DataType = 16;
            nii.ImgDim.BitsPerPixel = 32;
        case {'uint16'}
            nii.ImgDim.DataType = 512;
            nii.ImgDim.BitsPerPixel = 16;
        case {'uint32'}
            nii.ImgDim.DataType = 768;
            nii.ImgDim.BitsPerPixel = 32;
        case {'uint8'}
            nii.ImgDim.DataType = 2;
            nii.ImgDim.BitsPerPixel = 8;
    end
end
nii.VoxelData(l(1, 1), l(1, 2), l(1, 3), sum(fs)) = 0;
nii.ImgDim.Dim([1, 5]) = [4, sum(fs)];

% patch settings
nii.FileMagic = 'n+1';
nii.NIIFileType = 2;

% correct NIftI frame
trf = aif.CoordinateFrame.Trf;
trf(1:3, 4) = trf(1:3, 4) + trf(1:3, 1:3) * [1;1;1];
nii.DataHist.NIftI.QFormCode = 3;
nii.DataHist.NIftI.SFormCode = 3;
nii.DataHist.NIftI.QuatOffsetX = trf(1, 4);
nii.DataHist.NIftI.QuatOffsetY = trf(2, 4);
nii.DataHist.NIftI.QuatOffsetZ = trf(3, 4);
nii.DataHist.NIftI.AffineTransX = trf(1, :);
nii.DataHist.NIftI.AffineTransY = trf(2, :);
nii.DataHist.NIftI.AffineTransZ = trf(3, :);

% without resampling
if ~opts.resample

    % copy data into object
    tc = 1;
    for fc = 1:numel(anafiles)
        nii.VoxelData(:, :, :, tc:tc+fs(fc)-1) = ...
            anafiles{fc}.VoxelData(:, :, :, :);
        tc = tc + fs(fc);
    end

% with resample
else

    % get coordinate frame correctly again!
    trf = aif.CoordinateFrame.Trf;
    box = [Inf, Inf, Inf; 1, 1, 1; 1, 1, 1; l(1, 1:3)];

    % sample data
    tc = 1;
    for fc = 1:numel(anafiles)
        for vc = 1:fs(fc)
            nii.VoxelData(:, :, :, tc) = anafiles{fc}.SampleData3D(box, ...
                struct('mapvol', vc, 'method', resmeth, 'trans', trf));
            tc = tc + 1;
        end
    end
end

% clear objectc
clearxffobjects(anafiles);

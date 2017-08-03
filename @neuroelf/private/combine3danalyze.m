function h4 = combine3danalyze(files, res, sng)
% combine3danalyze  - combine 3D Analyze/NIftI files into 4D file
%
% FORMAT:       combine3danalyze(files [, res [, sng]])
%
% Input fields:
%
%       files       cell array with list of filenames
%       res         resampling method, if not given, not allowed!
%       sng
%
% Output fields:
%
%       h4          4D NIftI file with data from single volumes

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:14 AM EST
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
   ~iscell(files) || ...
    numel(files) < 2 || ...
   ~ischar(files{1}) || ...
    isempty(files{1}) || ...
    numel(files{1}) ~= size(files{1}, 2)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
files = files(:)';
for fc = 2:numel(files)
    if ~ischar(files{fc}) || ...
        isempty(files{fc}) || ...
        numel(files{fc}) ~= size(files{fc}, 2)
        error( ...
            'neuroelf:BadArgument', ...
            'Bad argument.' ...
        );
    end
end
if nargin < 2 || ...
   ~ischar(res) || ...
   ~any(strcmpi(res(:)', {'cubic', 'linear', 'nearest', 'lanczos3'}))
    res = '';
else
    res = lower(res(:)');
end
if nargin < 3 || ...
   ~islogical(sng) || ...
    numel(sng) ~= 1
    sng = false;
end

% try to load first image
try
    files{1} = xff(files{1});
    if ~isxff(files{1}, 'hdr') || ...
        size(files{1}.VoxelData, 4) ~= 1
        error( ...
            'neuroelf:BadArgument', ...
            'The first file must be a 3D Analyze/NIftI file.' ...
        );
    end
catch ne_eo;
    clearxffobjects(files);
    rethrow(ne_eo);
end

% get data layout, etc.
h4 = files{1}.CopyObject;
h4.VoxelData = h4.VoxelData(:, :, :);
sz = size(h4.VoxelData);
si = h4.ImgDim.ScalingIntercept;
if numel(si) ~= 1 || ...
    isinf(si) || ...
    isnan(si)
    si = 0;
end
ss = h4.ImgDim.ScalingSlope;
if numel(ss) ~= 1 || ...
    isinf(ss) || ...
    isnan(ss) || ...
    ss == 0
    ss = 1;
end
ss = 1 / ss;
if sng
    h4.VoxelData = single(h4.VoxelData);
    if si ~= 0 || ...
        ss ~= 1
        h4.VoxelData = single(double(ss) .* double(h4.VoxelData) + double(si));
        si = 0;
        ss = 1;
    end
    h4.ImgDim.DataType = 16;
    h4.ImgDim.BitsPerPixel = 32;
    h4.ImgDim.ScalingSlope = 1;
    h4.ImgDim.ScalingIntercept = 0;
end
trf = h4.CoordinateFrame.Trf;

% increase storage accordingly
h4.VoxelData(1, 1, 1, numel(files)) = 0;
h4.RunTimeVars.Map = h4.RunTimeVars.Map(ones(1, numel(files)));

% force to 4D compatible file
h4.FileMagic = 'n+1';
h4.NIIFileType = 2;

% check NIftI1 compliance
ntrf = [h4.DataHist.NIftI1.AffineTransX; ...
    h4.DataHist.NIftI1.AffineTransY; ...
    h4.DataHist.NIftI1.AffineTransZ];
if ~isequal(trf(1:3, 1:3), ntrf(:, 1:3))
    ntrf = trf;
    ntrf(4, :) = [0, 0, 0, 1];
    ntrf(1:3, 4) = 0;
    xtrf = ntrf * ones(4, 1);
    ntrf(1:3, 4) = trf(1:3, 4) + xtrf(1:3);
    h4.DataHist.NIftI1.QFormCode = 2;
    h4.DataHist.NIftI1.SFormCode = 2;
    h4.DataHist.NIftI1.QuaternionB = 0;
    h4.DataHist.NIftI1.QuaternionC = 1;
    h4.DataHist.NIftI1.QuaternionD = 0;
    h4.DataHist.NIftI1.QuatOffsetX = ntrf(1, 4);
    h4.DataHist.NIftI1.QuatOffsetY = ntrf(2, 4);
    h4.DataHist.NIftI1.QuatOffsetZ = ntrf(3, 4);
    h4.DataHist.NIftI1.AffineTransX = ntrf(1, :);
    h4.DataHist.NIftI1.AffineTransY = ntrf(2, :);
    h4.DataHist.NIftI1.AffineTransZ = ntrf(3, :);
end

% iterate over additional volumes
for fc = 2:numel(files)

    % try
    try

        % to load object
        fname = files{fc};
        files{fc} = xff(fname);
        hn = files{fc};
        if ~isxff(hn, 'hdr') || ...
            size(hn.VoxelData, 4) ~= 1
            error( ...
                'neuroelf:BadArgument', ...
                'File not an Analyze/NIftI file: %s.', ...
                fname ...
            );
        end

        % resampling required
        if ~isequal(size(hn.VoxelData), sz) || ...
           ~isequal(hn.CoordinateFrame.Trf(:)', trf(:)')

            % but not allowed
            if isempty(res)
                error( ...
                    'neuroelf:BadArgument', ...
                    'No resampling method selected, but images not in same space.' ...
                );
            end

            % resample
            nvd = hn.SampleTalBox(struct('BBox', [1, 1, 1; (sz + 0.25)], ...
                'ResXYZ', [1, 1, 1]), 1, res, trf);

        % no resampling required
        else

            % simply get data
            nvd = hn.VoxelData(:, :, :);
        end

        % scaling parameters
        nsi = hn.ImgDim.ScalingIntercept;
        if numel(nsi) ~= 1 || ...
            isinf(nsi) || ...
            isnan(nsi)
            nsi = 0;
        end
        nss = hn.ImgDim.ScalingSlope;
        if numel(nss) ~= 1 || ...
            isinf(nss) || ...
            isnan(nss) || ...
            nss == 0
            nss = 1;
        end

        % need to scale
        if nss ~= 1 || ...
            nsi ~= 0

            % then do it!
            nvd = nsi + nss .* double(nvd);
        end

        % "un-scale" required
        if ss ~= 1 || ...
            si ~= 0

            % un-scale
            nvd = ss .* (nvd - si);
        end

        % add to 4D object
        h4.VoxelData(:, :, :, fc) = nvd;

        % and copy Map (with description, etc.)
        h4.RunTimeVars.Map(fc) = hn.RunTimeVars.Map(1);

    % error handling
    catch ne_eo;

        % clear objects
        clearxffobjects(files);
        h4.ClearObject;

        % and report error
        rethrow(ne_eo);
    end
end

% clear original objects
clearxffobjects(files);

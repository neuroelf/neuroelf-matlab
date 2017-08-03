function h = newnii(size, opts)
% newnii  - create a new NII object with specific options
%
% FORMAT:       nii = newnii(size [, opts])
%
% Input fields:
%
%       size        1x3 or 1x4 size argument (3D or 4D)
%       opts        optional settings including
%        .dtype     either Matlab (char) or Analyze/NIftI datatype ('int8')
%        .mat       4x4 matrix (default, [0 0 0] in the middle, 1mm)
%        .type      Analyze/NIftI type (default: 'n+1')

if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'mat')
    opts.mat = eye(4);
    opts.mat(1:3, 4) = -0.5 .* (size(1:3)' - 1);
end
opts.dtype = 'int8';
opts.type = 'n+1';

% generate object
h = xff('new:hdr');

% make settings
h.FileMagic = opts.type;
if opts.type(2) == '+'
    h.NIIFileType = 2;
else
    h.NIIFileType = 1;
end
if numel(size) == 3
    h.ImgDim.Dim = [4, size, 1, 0, 0, 0];
else
    h.ImgDim.Dim = [4, size, 0, 0, 0];
end
if ischar(opts.dtype)
    switch (opts.dtype)
        case {'int16'}
            dvalue = int16(0);
            h.ImgDim.DataType = 4;
            h.ImgDim.BitsPerPixel = 16;
        case {'int8'}
            dvalue = int8(0);
            h.ImgDim.DataType = 256;
            h.ImgDim.BitsPerPixel = 8;
        case {'float32', 'single'}
            dvalue = single(0);
            h.ImgDim.DataType = 16;
            h.ImgDim.BitsPerPixel = 32;
    end
else
    h.ImgDim.DataType = opts.dtype;
    switch opts.dtype
        case {4}
            dvalue = int16(0);
            h.ImgDim.BitsPerPixel = 16;
        case {16}
            dvalue = single(0);
            h.ImgDim.BitsPerPixel = 32;
        case {256}
            dvalue = int8(0);
            h.ImgDim.BitsPerPixel = 8;
    end
end
%h.

% set data
h.VoxelData = repmat(dvalue, size);

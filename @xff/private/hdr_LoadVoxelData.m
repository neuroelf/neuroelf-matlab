function xo = hdr_LoadVoxelData(xo)
% HDR::LoadVoxelData  - load Analyze image voxel data
%
% FORMAT:       hdr.LoadVoxelData;
%
% No input/output fields.
%
% Using: analyzetype.

% Version:  v1.1
% Build:    16020516
% Date:     Feb-05 2016, 4:13 PM EST
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

% neuroelf library
global ne_methods;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get HDR file name -> IMG file
hdrc = xo.C;

% already transio
if istransio(hdrc.VoxelData)

    % resolve, store, and return
    xo.C.VoxelData = resolve(hdrc.VoxelData);
    return;
end

% if not empty, return
if ~isempty(hdrc.VoxelData)
    return;
end


% look in same folder
fflp = xo.F;
[ffpn, fflp, fflx] = fileparts(fflp);
ifname = [ffpn filesep fflp '.img'];
ifnamc = [ffpn filesep fflp '.IMG'];
if ~exist(ifname, 'file') == 2 && exist(ifnamc, 'file') == 2
    ifname = ifnamc;
end
if ~exist(ifname, 'file') == 2;
    error('neuroelf:xff:fileNotFound', ...
        'Related %s.img file for %s%s not found.', fflp, fflp, fflx);
end

% get data size and type
try
    dsize = hdrc.ImgDim.Dim(2:hdrc.ImgDim.Dim(1)+1);
    tsize = prod(dsize);
catch xfferror
    neuroelf_lasterr(xfferror);
    error('neuroelf:xff:badFileContent', 'Analyze ImgDim.Dim field error.');
end
dtype = hdrc.ImgDim.DataType;
if isempty(dtype)
    dtype = 0;
end
endian = hdrs.C.Endian;
if dtype > 255
    dtype = fix(dtype / 256);
    switch lower(endian)
        case 'ieee-le'
            endian = 'ieee-be';
        case 'ieee-be'
            endian = 'ieee-le';
        otherwise
            error('neuroelf:xff:internalError', ...
                'Bad machine datatype/encoding syntax combination.');
    end
end
[tmat, stype] = ne_methods.analyzetype(dtype);

% occupy mem
tmat(tsize) = tmat(1);

% open image file
try
    fid = fopen(ifname, 'rb', endian);
    fseek(fid, floor(abs(hdrc.ImgDim.VoxOffset)), 'bof');
    tmat(1:tsize) = fread(fid, [1, tsize], stype);
    fclose(fid);
catch xfferror
    neuroelf_lasterr(xfferror);
    error('neuroelf:xff:fileNotReadable', ...
        'Error opening/reading image file data from %s.', ifname);
end

% reshape voxel data
tmat = reshape(tmat, dsize);

% put data into object
xo.C.VoxelData = tmat;

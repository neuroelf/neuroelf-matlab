function rvalue = vmp_WriteAnalyzeVol(xo, mapno, filename, s2c)
% VMP::WriteAnalyzeVol  - write an Analyze image from a VMP map
%
% FORMAT:       [success = ] vmp.WriteAnalyzeVol(mapno, filename [, s2c])
%
% Input fields:
%
%       mapno       map number
%       filename    volume filename
%       s2c         SPM2 compatibility-flag (default: false)
%
% Output fields:
%
%       success     true if successful
%
% Using: bvcoordconv.

% Version:  v1.1
% Build:    16021314
% Date:     Feb-13 2016, 2:59 PM EST
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
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp') || ...
   ~isa(mapno, 'double') || numel(mapno) ~= 1 || isinf(mapno) || isnan(mapno) || ...
    mapno < 1 || mapno ~= fix(mapno) || ~ischar(filename) || isempty(filename) || ...
    numel(filename) ~= length(filename) || length(filename) < 5 || ...
   (~strcmpi(filename(end-3:end), '.hdr') && ~strcmpi(filename(end-3:end), '.img'))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 4 || ~islogical(s2c) || numel(s2c) ~= 1
    s2c = false;
end
bc = xo.C;
if isempty(bc.Map) || mapno > numel(bc.Map)
    error('neuroelf:xff:badArgument', 'Empty VMP or map number too high.');
end
bcm = bc.Map(mapno);
vr = bc.Resolution;

% default rvalue
rvalue = false;

% filename
filename = filename(:)';

% for SPM2
if s2c

    % ensure filename is valid!
    filename = [filename(1:end-3) 'img'];
    hdrfname = [filename(1:end-3) 'hdr'];
    matfname = [filename(1:end-3) 'mat'];

    % dimension and datatype
    vmpdat = single(permute(bcm.VMPData(end:-1:1, end:-1:1, end:-1:1), [3, 1, 2]));
    vmpdsc = bcm.Name;
    if numel(vmpdsc) > 80
        vmpdsc = [vmpdsc(1:77) '...'];
    end
    isiz = size(vmpdat);
    offx = bc.ZEnd - 128;
    offy = bc.XEnd - 128;
    offz = bc.YEnd - 128;
    orgx = 1 + round(offx / vr);
    orgy = 1 + round(offy / vr);
    orgz = 1 + round(offz / vr);

    % build mat
    tmat = vr .* eye(4);
    tmat(end) = 1;
    tmat(1:3, 4) = -[offx; offy; offz];

    % try volume creation
    hdrc = cell(1, 1);
    try
        hdr = xff('new:hdr');
        hdrc{1} = hdr;
        hdrbc = hdr.C;
    catch xfferror
        clearxffobjects(hdrc);
        error('neuroelf:xff:internalError', ...
            'Error creating Analyze header object: %s.', xfferror.message);
    end

    % set dims and data
    hdrbc.ImgDim.Dim(1:4) = [3, isiz];
    hdrbc.ImgDim.PixSpacing(2:4) = [vr, vr, vr];
    hdrbc.ImgDim.DataType = 16;
    hdrbc.ImgDim.BitsPerPixel = 32;
    hdrbc.ImgDim.CalMaxDisplay = 32767;
    hdrbc.ImgDim.CalMinDisplay = -32768;
    hdrbc.DataHist.Description = vmpdsc;
    hdrbc.DataHist.OriginSPM = [orgx, orgy, orgz, 0, 0];
    hdrbc.VoxelData = vmpdat;
    hdr.C = hdrbc;

    % save hdr/img
    try
        aft_SaveAs(hdr, hdrfname);
        rvalue = hdr_SaveVoxelData(hdr);
    catch xfferror
        warning('neuroelf:xff:internalError', ...
            'Error writing Analyze header/image %s (%s).', hdrfname, xfferror.message);
        return;
    end
    if ~rvalue
        return;
    end

    % save mat
    try
        eval('M=tmat;mat=tmat;save(matfname,''M'',''mat'',''-v6'');')
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
    if exist(matfname, 'file') ~= 2
        rvalue = false;
    end

% for newer versions
else

    % big TRY/CATCH
    nio = {[]};
    try

        % create temporary NII object
        nio{1} = xff('new:nii');
        nii = nio{1}.C;

        % output name extension NII
        if ~isempty(regexpi(filename, '\.nii$'))

            % set magic
            nii.NIIFileType = 2;
            nii.FileMagic = 'n+1';

        % not NII
        else

            % set magic
            nii.NIIFileType = 1;
            nii.FileMagic = 'ni1';
        end

        % set VoxelData
        nii.VoxelData = permute(bcm.VMPData(end:-1:1, end:-1:1, end:-1:1), [3, 1, 2]);
        vs = size(bcm.VMPData);

        % set size, type, etc.
        nii.ImgDim.Dim(1:5) = [3, size(nii.VoxelData), 1];
        nii.ImgDim.DataType = 16;
        nii.ImgDim.BitsPerPixel = 32;
        nii.ImgDim.PixSpacing(2:4) = vr;

        % generate transformation matrix rows
        ti = vr .* (vs - ne_methods.bvcoordconv([0, 0, 0], 'tal2bvc', aft_BoundingBox(xo)));
        nii.DataHist.NIftI1.QFormCode = 2;
        nii.DataHist.NIftI1.SFormCode = 2;
        nii.DataHist.NIftI1.QuaternionB = 0;
        nii.DataHist.NIftI1.QuaternionC = 1;
        nii.DataHist.NIftI1.QuaternionD = 0;
        nii.DataHist.NIftI1.QuatOffsetX = -ti(3);
        nii.DataHist.NIftI1.QuatOffsetY = -ti(1);
        nii.DataHist.NIftI1.QuatOffsetZ = -ti(2);
        nii.DataHist.NIftI1.AffineTransX = [vr,  0,  0, -ti(3)];
        nii.DataHist.NIftI1.AffineTransY = [ 0, vr,  0, -ti(1)];
        nii.DataHist.NIftI1.AffineTransZ = [ 0,  0, vr, -ti(2)];

        % set in object
        nio{1}.C = nii;

        % save file
        aft_SaveAs(nio{1}, filename);
        if nii.NIIFileType == 1
            rvalue = hdr_SaveVoxelData(nio{1});
        end
    catch xfferror
        clearxffobjects(nio);
        rethrow(xfferror);
    end

    % clear object
    clearxffobjects(nio);
end

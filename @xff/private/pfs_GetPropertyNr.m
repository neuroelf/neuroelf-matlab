function p = pfs_GetPropertyNr(xo, idx)
% PFS::GetPropertyNr  - read POIFS property from POIFS file
%
% FORMAT:       p = pfs.GetPropertyNr(idx)
%
% Input fields:
%
%       idx         the index of the POIFS file (1-based)
%
% Output fields:
%
%       p           the obtained property

% Version:  v1.1
% Build:    16021018
% Date:     Feb-10 2016, 6:12 PM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'pfs') || ...
    numel(idx) ~= 1 || isinf(idx) || isnan(idx) || idx < 1
    error('neuroelf:xff:badArgument', 'Invalid argument in call.');
end
bc = xo.C;

% try opening file
fid = fopen(xo.F, 'r', 'ieee-le');
if fid < 1
    error('neuroelf:xff:fileOpenError', 'Error opening POIFS file.');
end
pblidx = bc.FirstPropertyBlock;
oidx = idx;
while (idx > bc.PropertiesPerBlock)
    pblidx = bc.BAT(pblidx);
    if pblidx < 1
        fclose(fid);
        error('neuroelf:xff:POIFSError', 'BAT link error in POIFS file.');
    end
    idx = idx - bc.PropertiesPerBlock;
end
ppos = bc.LargeBlockSize * double(pblidx) + (idx - 1) * 128;
fseek(fid, ppos, -1);
if ftell(fid) ~= ppos
    fclose(fid);
    error('neuroelf:xff:POIFSError', 'Error seeking to property.');
end
p = cell2struct(cell(1, 1, 15), {'Name', 'NameSize', 'Type', 'Color', ...
     'Previous', 'Next', 'FirstChild', 'CLSID', 'UserFlags', 'Created', 'Modified', ...
     'FirstBlock', 'Size', 'BlockAccess', 'Reserved7c'}, 3);
p.Name = fread(fid, [1, 32], 'uint16=>uint16');
p.NameSize = floor(fread(fid, [1, 1], 'uint16=>double') / 2) - 1;
p.Name = char(p.Name(1:p.NameSize));
p.Type = fread(fid, [1, 1], 'uint8=>double');
p.Color = fread(fid, [1, 1], 'uint8=>double');
p.Previous = fread(fid, [1, 1], 'int32=>double') + 1;
p.Next = fread(fid, [1, 1], 'int32=>double') + 1;
p.FirstChild = fread(fid, [1, 1], 'int32=>double') + 1;
p.CLSID = fread(fid, [1, 16], 'uint8=>uint8');
p.UserFlags = fread(fid, [1, 1], 'int32=>int32');
p.Created = fread(fid, [1, 1], 'int64=>int64');
p.Modified = fread(fid, [1, 1], 'int64=>int64');
if p.Created(1) ~= 0
    p.Created = datestr(584389 + double(p.Created) / 8.64e11);
end
if p.Modified ~= 0
    p.Modified = datestr(584389 + double(p.Modified) / 8.64e11);
end
p.FirstBlock = fread(fid, [1, 1], 'int32=>int32') + 1;
p.Size = fread(fid, [1, 1], 'int32=>double');
if p.Size >= bc.SBATThresholdSize
    p.BlockAccess = 'BAT';
else
    p.BlockAccess = 'SBAT';
end
p.Reserved7c = fread(fid, [1, 1], 'uint32=>double');
fclose(fid);

% for Root for Access to BAT
if oidx == 1
    p.BlockAccess = 'BAT';
end

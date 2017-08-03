function [tmat, stype, bsize] = analyzetype(anatype)
% analyzetype  - return properties of analyze image datatype
%
% FORMAT:       [anadata, anareadtype, anareadsize] = analyzetype(anatype)
%
% Input fields:
%
%       anatype     numeric datatype given in Analyze header
%
% Output fields:
%
%       anadata     1x1 data field with correct type and zero content
%       anareadtype type string for fread operations
%       anareadsize number of bytes per pixel/voxel

% Version:  v0.9a
% Build:    10070718
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
   ~isa(anatype, 'double') || ...
    numel(anatype) ~= 1 || ...
    isnan(anatype) || ...
    isinf(anatype) || ...
    anatype < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument for call to %s.', ...
        mfilename ...
    );
end

% double check
while anatype >= 256 && ...
    ~any(anatype == [256, 512, 768])
    anatype = fix(anatype / 256);
end
switch anatype
    case {2},   bsize = 1; tmat = uint8(0);  stype = 'uint8=>uint8';
    case {4},   bsize = 2; tmat = int16(0);  stype = 'int16=>int16';
    case {8},   bsize = 4; tmat = int32(0);  stype = 'int32=>int32';
    case {16},  bsize = 4; tmat = single(0); stype = 'single=>single';
    case {64},  bsize = 8; tmat = 0;         stype = 'double=>double';
    case {130, 256}, bsize = 1; tmat = int8(0);   stype = 'int8=>int8';
    case {132, 512}, bsize = 2; tmat = uint16(0); stype = 'uint16=>uint16';
    case {136, 768}, bsize = 4; tmat = uint32(0); stype = 'uint32=>uint32';
    otherwise
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid anatype value: %d.', ...
            anatype ...
        );
end

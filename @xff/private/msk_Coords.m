function vox = msk_Coords(xo, csc)
% MSK::Coords  - return a list of coordinates in a mask
%
% FORMAT:       vox = msk.Coords([convention])
%
% Input fields:
%
%       convention  either 'bvint', {'bvsys'}, or 'tal'
%
% Output fields:
%
%       vox         Vx3 list of coordinates

% Version:  v1.1
% Build:    16020917
% Date:     Feb-09 2016, 5:06 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'msk')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~ischar(csc) || ~any(strcmpi(csc(:)', {'bvint', 'bvsys', 'int', 'sys', 'tal'}))
    csc = 'bvsys';
else
    csc = lower(csc(:)');
end

% find coordinates
[c1, c2, c3] = ind2sub(size(bc.Mask), find(bc.Mask ~= 0));

% multiply correctly and add ?Start Value
res = bc.Resolution;
c1 = (c1 - 1) .* res + bc.XStart;
c2 = (c2 - 1) .* res + bc.YStart;
c3 = (c3 - 1) .* res + bc.ZStart;

% what convention
switch (csc)

    % BV internal
    case {'bvint', 'int'}
        vox = [c1(:), c2(:), c3(:)];

    % BV-System
    case {'bvsys', 'sys'}
        vox = [c3(:), c1(:), c2(:)];

    % TAL
    case 'tal'
        vox = 128 - [c3(:), c1(:), c2(:)];
end

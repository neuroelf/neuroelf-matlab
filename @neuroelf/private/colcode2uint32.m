function outarray = colcode2uint32(inarray)
% colcode2uint32  - convert from colcode to uint32 array
%
% FORMAT:       outarray = colcode2uint32(inarray)
%
% Input fields:
%
%       inarray     Nx4 input array of type colcolde
%                   (n, 1) = color index or NaN
%                   (n, 2) = red component (0...255)
%                   (n, 3) = green component (0...255)
%                   (n, 4) = blue component (0...255)
%
% Output fields:
%
%       outarray    Nx1 uint32 array
%
% See also uint322colcode

% Version:  v0.9a
% Build:    10051716
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
if nargin ~= 1 || ...
   ~isa(inarray, 'double') || ...
   ~isreal(inarray) || ...
    issparse(inarray) || ...
    ndims(inarray) ~= 2 || ...
    size(inarray, 2) ~= 4
    error( ...
        'neuroelf:BadArgument', ...
        'Bad/missing input argument supplied.' ...
    );
end

% produce output array
outarray = uint32(zeros(size(inarray, 1), 1));

% leave if empty
if isempty(inarray)
    return;
end

% force RGB in RGB(:, [2, 3, 4])
for cc = 2:4
    inarray(isinf(inarray(:, cc)), cc) = 0;
    inarray(isnan(inarray(:, cc)), cc) = 0;
    inarray(inarray(:, cc) < 0   , cc) = 0;
    inarray(inarray(:, cc) > 255 , cc) = 255;
end
inarray = fix(inarray);

% rgbcolormod/min, rfactor, gfactor
rgbcolortot = 256 * 16777216;
rgbcolormin = 63 * 16777216;
rfactor = 65536;
gfactor = 256;

% get colors with indices
rgbcolor = isnan(inarray(:, 1));

% set index color
outarray(~rgbcolor) = uint32(mod(inarray(~rgbcolor, 1), rgbcolortot));

% fill those correctly
outarray(rgbcolor) = uint32(rgbcolormin + ...
    rfactor .* inarray(rgbcolor, 2) + ...
    gfactor .* inarray(rgbcolor, 3) + ...
    inarray(rgbcolor, 4));

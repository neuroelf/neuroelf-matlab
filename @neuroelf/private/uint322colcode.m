function outarray = uint322colcode(inarray)
% uint322colcode  - convert from uint32 to colcode array
%
% FORMAT:       outarray = uint322colcode(inarray)
%
% Input fields:
%
%       inarray     Nx1 (or 1xN) input array of type uint32
%
% Output fields:
%
%       outarray    Nx4 double array
%                   (n, 1) = color index or NaN
%                   (n, 2) = red component (0...255)
%                   (n, 3) = green component (0...255)
%                   (n, 4) = blue component (0...255)
%
% See also colcode2uint32

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
   ~isa(inarray, 'uint32')
    error( ...
        'neuroelf:BadArgument', ...
        'Bad/missing input argument supplied.' ...
    );
end

% linearize array
inarray = double(inarray(:));

% produce output array
outarray = zeros(size(inarray, 1), 4);
outarray(:, 1) = NaN;

% index color range and moduli
idxcolormin = 63 * 16777216;
idxcolormax = 64 * 16777216;
gmodulus = 65536;
bmodulus = 256;

% get colors with indices
rgbcolor = ((inarray >= idxcolormin) & (inarray <  idxcolormax));

% set index
outarray(~rgbcolor, 1) = inarray(~rgbcolor);

% fill RGB correctly
inarray = inarray(rgbcolor) - idxcolormin;
outarray(rgbcolor, 2:4) = floor([inarray / 65536, ...
    mod(inarray, gmodulus) / 256, ...
    mod(inarray, bmodulus)]);

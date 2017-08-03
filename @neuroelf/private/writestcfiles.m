function writestcfiles(fourd, filetmp, sdim)
% writestcfiles  - create STC files from a 4D array
%
% FORMAT:       writestcfiles(fourd, filetmp [, sdim])
%
% Input fields:
%
%       fourd       4-D data array (if not within 0..32767, will be norm.)
%       filetmp     filename template with one %d sequence
%       sdim        dimension of slices (default: 1)

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

% check arguments
if nargin < 2 || ...
   ~isnumeric(fourd) || ...
   ~any([3, 4] == ndims(fourd)) || ...
   ~ischar(filetmp) || ...
    isempty(strfind(filetmp(:)', '%d'))
    error( ...
        'neuroelf:BadArguments', ...
        'Invalid or missing argument given.' ...
    );
end
filetmp = filetmp(:)';
if nargin < 3 || ...
   ~isa(sdim, 'double') || ...
    numel(sdim) ~= 1 || ...
    isinf(sdim) || ...
    isnan(sdim) || ...
    sdim < 1 || ...
    sdim > 4
    sdim = 1;
else
    sdim = real(sdim);
end
if sdim ~= 1
    fourd = permute(fourd, [sdim, setdiff(1:ndims(fourd), sdim)]);
end

% get dimensions
ad = size(fourd);

% check datatype
mna = min(fourd(:));
if mna < 0
    fourd = fourd - mna;
end
mxa = max(fourd(:));
if mxa >= 32767.5
    fourd = fourd / (mxa / 32768);
end
fourd = round(fourd);

% get a slice
slice = xff('new:stc');
slice.NrOfRows = ad(2);
slice.NrOfCols = ad(3);
slice.NrOfVolumes = ad(4);

% iterate over slices
for sc = 1:ad(1)
    slice.STCData = squeeze(fourd(sc, :, :, :));
    slice.SaveAs(sprintf(filetmp, sc));
end

% clear object
slice.ClearObject;

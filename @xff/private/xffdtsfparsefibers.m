function fibers = xffdtsfparsefibers(fid, numfib)
% xffdtsfparsefibers  - parse fibers from DTSF file
%
% FORMAT:       fibers = xffdtsfparsefibers(fid, numfib)
%
% Input fields:
%
%       fid         input file fid (fopen)
%       numfib      number of fibers
%
% Output fields
%
%       fibers      Nx1 struct array with fields
%        .NrOfPoints   number of points for that fiber
%        .Selected     uint8, either 0 or 1
%        .RGB          1x3 uint8 array
%        .FromToPoint  1x2 double array, [0, NrOfPoints]
%        .Coord        Px3 coordinates of fiber points
%
% See also xff

% Version:  v1.1
% Build:    16021216
% Date:     Feb-12 2016, 4:28 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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
if nargin < 2 || ~isa(fid, 'double') || isempty(fid) || ~isreal(fid) || ...
   ~any(fopen('all') == fid(1)) || ~isa(numfib, 'double') || numel(numfib) ~= 1 || ...
   ~isreal(numfib) || isnan(numfib) || isinf(numfib) || numfib < 0
    error('neuroelf:xff:badArgument', 'Bad or missing argument.');
end

% build structure first, so we don't need growing it later!
fibers = struct;
fibers.NrOfPoints  = 1;
fibers.Selected    = 0;
fibers.RGB         = [0, 0, 0];
fibers.FromToPoint = [0, 1];
fibers.Coord       = [0, 0, 0];
if numfib < 1
    fibers(:) = [];
    return;
end
fibers(2:numfib) = fibers;

% try reading
try

    % loop over fibers
    for fc = 1:numfib

        % get number of points
        numpts = fread(fid, [1, 1], 'uint32=>double');
        fibers(fc).NrOfPoints  = numpts;
        fibers(fc).Selected    = fread(fid, [1, 1], '*uint8');
        fibers(fc).RGB         = fread(fid, [1, 3], '*uint8');
        fibers(fc).FromToPoint = fread(fid, [1, 2], 'uint32=>double');
        fibers(fc).Coord       = fread(fid, [3, numpts], 'single=>double')';
    end
catch xfferror
    rethrow(xfferror);
end

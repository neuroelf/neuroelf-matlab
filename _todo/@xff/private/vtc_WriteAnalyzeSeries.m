function rvalue = vtc_WriteAnalyzeSeries(xo, pattern, offset, range)
% VTC::WriteAnalyzeSeries  - writes analyze vols for the (entire) TC
%
% FORMAT:       vtc.WriteAnalyzeSeries(pattern, offset, range)
%
% Input fields:
%
%       pattern     output filename pattern
%       offset      volume offset (for filename)
%       range       range of volumes

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:31 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get number of volumes
bc = xo.C;
numvol = bc.NrOfVolumes;

% check further arguments
if nargin < 4 || ~isa(range, 'double') || isempty(range) || ...
    any(isinf(range(:)') | isnan(range(:)')) || ...
    any(fix(range(:)') ~= range(:)' | range(:)' < 1 | range(:)' > numvol)
    range = 1:numvol;
end
range = range(:)';
if length(range) == 2
    range = range(1):range(2);
end
if nargin < 3 || ~isa(offset, 'double') || numel(offset) ~= 1 || ...
    isinf(offset) || isnan(offset) || fix(offset) ~= offset || offset < 1 || offset > 9999
    offset = 1;
end
if nargin < 2 || ~ischar(pattern) || ~any(pattern(:)' == '%') || ...
    (isempty(strfind(lower(pattern(:)'), '.img')) && isempty(strfind(lower(pattern(:)'), '.hdr')))
    pattern = sprintf('vtcvol_%%0%dd.img', length(num2str(numvol + offset - 1)));
end
pattern = pattern(:)';

% loop over images
for vc = range
    tname = sprintf(pattern, vc + offset - 1);
    succ = vtc_WriteAnalyzeVol(xo, vc, tname);
    if ~succ
        error('neuroelf:xff:internalError', 'Error writing volume %d to file %s.', vc, tname);
    end
end
rvalue = succ;

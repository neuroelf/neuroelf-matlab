function rvalue = fmr_WriteAnalyzeSeries(xo, pattern, offset, range, flip)
% FMR::WriteAnalyzeSeries  - writes analyze vols for the (entire) TC
%
% FORMAT:       [OK] = fmr.WriteAnalyzeSeries([pattern, offset, range, flip]);
%
% Input fields:
%
%       pattern     analyze filename pattern (default: fmrvol_%03d.vol)
%       offset      first number for output file (default: 1)
%       range       which FMR volumes to write (default: all)
%       flip        optional flipping of images (e.g. 'xy', default: '')
%
% Output fields:
%
%       OK          true if write was successful

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:10 AM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neuroelf:xff:BadArgument', 'Invalid call to %s.', mfilename);
end

% get content
bc = xo.C;

% if slice data not loaded at all?
if isempty(bc.Slice) || ~isstruct(bc.Slice) || ~isfield(bc.Slice, 'STCData')

    % make sure data is loaded
    try
        fmr_LoadSTC(xo);
        bc = xo.C;
    catch xfferror
        rethrow(xfferror);
    end
end
try
    if istransio(bc.Slice.STCData)
        bc.Slice.STCData = bc.Slice.STCData(:, :, :, :);
        xo.C = bc;
    end
catch xfferror
    error('neuroelf:xff:InternalError', 'Error resolving transio: %s.', xfferror.message);
end
numvol = bc.NrOfVolumes;

% check further arguments
if nargin < 5 || ~ischar(flip) || numel(flip) > 3
    flip = '';
else
    flip = lower(flip(:)');
end
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
    isnan(offset) || fix(offset) ~= offset || offset < 1 || offset > 9999
    offset = 1;
end
if nargin < 2 || ~ischar(pattern) || ~any(pattern(:)' == '%') || ...
    (isempty(strfind(lower(pattern(:)'), '.img')) && isempty(strfind(lower(pattern(:)'), '.hdr')))
    pattern = sprintf('fmrvol_%%0%dd.img', length(num2str(numvol + offset - 1)));
end
pattern = pattern(:)';

% try to create progress bar
try
    pbar = xprogress;
    xprogress(pbar, 'settitle', 'Writing FMR Analyze series...');
    xprogress(pbar, 0, 'Writing volumes...', 'visible', 1, numel(range));
catch xfferror
    neuroelf_lasterr(xfferror);
    pbar = [];
end

% loop over images
for vc = 1:numel(range)
    tname = sprintf(pattern, range(vc) + offset - 1);
    succ = fmr_WriteAnalyzeVol(xo, range(vc), tname, flip);
    if ~succ
        if ~isempty(pbar)
            closebar(pbar);
        end
        error('neuroelf:xff:internalError', ...
            'Error writing volume %d to file %s.', range(vc), tname);
    end
    if ~isempty(pbar)
        xprogress(pbar, vc, sprintf('Writing volume %d...', vc));
    end
end
if ~isempty(pbar)
    closebar(pbar);
end
rvalue = succ;

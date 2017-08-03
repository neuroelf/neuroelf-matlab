function apmp = averagepmps(pmps, sd)
% averagepmps  - average PMP objects and create new object
%
% FORMAT:       apmp = averagepmps(pmps [, sd])
%
% Input fields:
%
%       pmps        cell array with either objects or filenames
%       sd          if given and true also create std map(s)
%
% Output fields:
%
%       apmp        averaged PMP

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:13 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% check arguments
if nargin < 1 || ...
   ~iscell(pmps) || ...
    numel(pmps) < 2
    error( ...
        'neuroelf:BadArguments', ...
        'Invalid call to %s.', ...
        mfilename ...
    );
end
pnum = numel(pmps);
pldd = false(pnum, 1);
for pc = pnum:-1:1
    if ischar(pmps{pc}) && ...
        exist(pmps{pc}(:)', 'file') == 2
        try
            pldd(pc) = true;
            pmps{pc} = xff(pmps{pc});
            if ~isxff(pmps{pc}, 'pmp')
                error('BADPMP');
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            clearxffobjects(pmps(pldd));
            error( ...
                'neuroelf:BadFilename', ...
                'Invalid PMP filename.' ...
            );
        end
    elseif ~isxff(pmps{pc})
        clearxffobjects(pmps(pldd));
        error( ...
            'neuroelf:BadFilename', ...
            'Invalid PMP filename.' ...
        );
    end
end
r1 = pmps{1}.ThetaResolution;
r2 = pmps{1}.PhiResolution;
nm = pmps{1}.NrOfMaps;
for pc = pnum:-1:2
    if pmps{pc}.ThetaResolution ~= r1 || ...
        pmps{pc}.PhiResolution ~= r2 || ...
        pmps{pc}.NrOfMaps ~= nm
        clearxffobjects(pmps(pldd));
        error( ...
            'neuroelf:ObjectsMismatch', ...
            'All PMPs must have the same number of maps and resolution.' ...
        );
    end
end

% standard deviation maps ?
if nargin > 1 && ...
    numel(sd) == 1 && ...
   (islogical(sd) || ...
    isa(sd, 'double')) && ...
    sd
    sd = true;
else
    sd = false;
end

% create required map field in memory
maps = zeros([r1, r2, pnum]);

% create output map
apmp = xff('new:pmp');
if r1 ~= 360 || ...
    r2 ~= 180
    apmp.FileVersion = apmp.FileVersion + 256;
end
if sd
    apmp.NrOfMaps = 2 * nm;
else
    apmp.NrOfMaps = nm;
end
apmp.ThetaResolution = r1;
apmp.PhiResolution = r2;
apmp.Map(apmp.NrOfMaps).PMPData = [];

% iterate over number of maps, then objects
for mc = 1:nm
    for oc = 1:pnum
        maps(:, :, oc) = pmps{oc}.Map(mc).PMPData(:, :);
    end

    % build average
    apmp.Map(mc).PMPData = mean(maps, 3);

    % also build std map ?
    if sd
        apmp.Map(mc + nm).PMPData = std(maps, [], 3);
    end
end

% unload objects as needed
clearxffobjects(pmps(pldd));

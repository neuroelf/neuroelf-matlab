function bcc = bcc2area(bvfolder)
% bcc2area  - lookup TalairachBrain.vmr color codes
%
% FORMAT:       bcc = bcc2area(BVfolder)
%
% Input fields:
%
%       BVfolder    path name where TalairachBrain.vmr is stored
%
% Output fields:
%
%       bcc         255x1 cell array with text descriptions
%
% See also xff

% Version:  v0.9d
% Build:    14061710
% Date:     Jun-17 2014, 10:00 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% cd to correct folder and load TalairachBrain.vmr
try
    cd(bvfolder);
    tvo = xff('TalairachBrain.vmr');
    tv = tvo.VMRData(:, :, :);
    tvbb = tvo.BoundingBox;
    tvo.ClearObject;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:BadArgument', ...
        'No or bad BrainVoyager installation folder given.' ...
    );
end

% loop over colorcodes
mmm = minmaxmean(tv(tv > 0));
bcc = cell(255, 1);
for cc = mmm(1):mmm(2)

    % say for which color code...
    ccout = {sprintf('Results for color code %d:', cc)};
    disp(ccout{1});

    % find voxels
    tal = bvcoordconv(find(tv(:) == uint8(cc)), 'bvx2tal', tvbb);
    cend = size(tal, 1);

    % prepare arrays
    ttd  = cell(cend, 1);
    ttdr = cell(cend, 5);

    % lookup each coordinate
    for c = 1:cend

        % use tdclient with NGM and 5 mm radius search
        ttd{c} = tdlocal2(5, tal(c, 1), tal(c, 2), tal(c, 3), [5, 1]);

        % keep user updated
        if ~mod(c, 1000)
            disp(sprintf('%d from %d coords looked up...', c, cend));
            drawnow;
        end
    end

    % get unique names
    for c = 1:cend

        % remove non-used part
        ttdr(c, :) = splittocell(regexprep(regexprep( ...
            ttd{c}, '^.*\:\s+', ''), '\s+(\(coo.*)?$', ''), ',');

        % glue rest together again
        ttd{c} = gluetostring(ttdr(c, :), ',');
    end
    u = unique(ttd);

    % print number of voxels and
    ccout{1 + numel(u)} = '';
    for c = 1:numel(u)
        ccout{c + 1} = sprintf('%5d vx: %s', ...
            length(find(strcmp(ttd, u{c}))), u{c});
    end

    % put ccout in correct cell
    bcc{cc} = ccout;
end

% delete empty cells
bcc(cellfun('isempty', bcc)) = [];

function avmr = averagevmrs(v, cv)
% averagevmrs  - average VMRs
%
% FORMAT:       avmr = averagevmrs(v [, centervmr])
%
% Input fields:
%
%       v           cell array with list of VMR (objects or filenames)
%       centervmr   center around middle coordinate (default false)
%
% Output fields:
%
%       avmr        average VMR with maximum bounding box

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

% argument check
if nargin < 1 || ...
   ~iscell(v) || ...
    numel(v) < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~islogical(cv) || ...
    numel(cv) ~= 1
    cv = false;
end
nvmr = numel(v);

% test xprogress
try
    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 640, 36]);
    xprogress(pbar, 'settitle', 'Averaging VMRs...');
    xprogress(pbar, 0, 'Checking VMRs...', 'visible', 0, nvmr + 3);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    pbar = [];
end

% check each cell
xffroot = xff();
v16auto = xffroot.Config('vmr', 'v16auto', false);
vmrtios = xffroot.TransIOSize('vmr', 65536);
v = v(:);
res = zeros(nvmr, 3);
vm8 = true(1, nvmr);
clearloaded = false(1, nvmr);
for cc = 1:nvmr
    if ischar(v{cc}) && ...
        numel(v{cc}) > 4 && ...
        exist(v{cc}(:)', 'file') == 2 && ...
        any(strcmpi(v{cc}(end-3:end), {'.vmr', '.v16'}))
        try
            [v{cc}, clearloaded(cc)] = xff(v{cc}(:)');
        catch ne_eo;
            clearxffobjects(v(clearloaded));
            if ~isempty(pbar)
                closebar(pbar);
            end
            xffroot.Config('vmr', 'v16auto', v16auto);
            xffroot.TransIOSize('vmr', vmrtios);
            error( ...
                'neuroelf:xffError', ...
                'Error loading VMR from file ''%s'': %s.', ...
                v{cc}(:)', ne_eo.message ...
            );
        end
    elseif numel(v{cc}) ~= 1 || ...
       ~isxff(v{cc}, 'vmr')
        clearxffobjects(v(clearloaded));
        if ~isempty(pbar)
            closebar(pbar);
        end
        xffroot.Config('vmr', 'v16auto', v16auto);
        xffroot.TransIOSize('vmr', vmrtios);
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid cell (no VMR filename or object).' ...
        );
    end
    res(cc, :) = v{cc}.BoundingBox.ResXYZ;
    vm8(cc) = v{cc}.VMR8bit;
end
xffroot.Config('vmr', 'v16auto', v16auto);
xffroot.TransIOSize('vmr', vmrtios);

% check resolutions
if any(any(diff(res) ~= 0)) || ...
    any(vm8 ~= vm8(1))
    clearxffobjects(v(clearloaded));
    if ~isempty(pbar)
        closebar(pbar);
    end
    error( ...
        'neuroelf:BadArgument', ...
        'VMRs must match in resolution and bytesize.' ...
    );
end
vm8 = vm8(1);

% get max bounding box
if ~isempty(pbar)
    xprogress(pbar, 1, 'Getting bounding box...');
end
bbox = minmaxbbox(v);
bbox = bbox(3:4, :);

% create array with matching size
vmrsz = diff(bbox) + 1;
if vm8
    avmrd = uint16([]);
    avmrd(vmrsz(1), vmrsz(2), vmrsz(3)) = uint16(0);
else
    avmrd = uint32([]);
    avmrd(vmrsz(1), vmrsz(2), vmrsz(3)) = uint32(0);
end

% iterate over VMRs
if ~isempty(pbar)
    xprogress(pbar, 2, 'Averaging...');
end
for cc = 1:nvmr

    % reframe VMR and get data
    av = v{cc}.CopyObject;
    if istransio(av.VMRData)
        av.VMRData = av.VMRData(:, :, :);
    end
    av.Reframe(bbox);

    % add data to average
    if vm8
        avmrd = avmrd + min(uint16(av.VMRData), uint16(225));
    else
        avmrd = avmrd + av.VMRData;
    end
    if ~isempty(pbar)
        xprogress(pbar, 2 + cc);
    end
    av.ClearObject;
end

% copy first VMR and then clear loaded VMRs
avmr = v{1}.CopyObject;
clearxffobjects(v(clearloaded));

% build average VMR
if max(avmrd(:)) > (225 * nvmr)
    avmrd = uint16(single(avmrd) ./ single(nvmr));
    vm8 = false;
else
    avmrd = uint8(single(avmrd) ./ single(nvmr));
    vm8 = true;
end
if ~isempty(pbar)
    xprogress(pbar, 3 + nvmr);
end

% copy data into avmr
avmr.VMRData = avmrd;
avmr.VMR8bit = vm8;

% reframe around center ?
if cv

    % which resolution
    fc = avmr.FramingCube / 2;
    bb = avmr.BoundingBox.BBox;
    bb(2, :) = bb(2, :) + 1;
    ns = max(abs(fc - bb));
    avmr.Reframe([fc - ns; (fc - 1) + ns]);
end

% clear no longer needed settings
if avmr.FileVersion == 3
    avmr.Slice1CenterX = -127.5;
    avmr.Slice1CenterY = 0;
    avmr.Slice1CenterZ = 0;
    avmr.SliceNCenterX = 127.5;
    avmr.SliceNCenterY = 0;
    avmr.SliceNCenterZ = 0;
    avmr.RowDirX = 0;
    avmr.RowDirY = 1;
    avmr.RowDirZ = 0;
    avmr.ColDirX = 0;
    avmr.ColDirY = 0;
    avmr.ColDirZ = -1;
    avmr.NrOfPastSpatialTransformations = 0;
    avmr.Trf(:) = [];
end

% clear progress bar
if ~isempty(pbar)
    closebar(pbar);
end

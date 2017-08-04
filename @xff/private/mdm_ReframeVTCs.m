function xo = mdm_ReframeVTCs(xo, bbox, newext)
% MDM::ReframeVTCs  - reframe all VTCs referenced in an MDM to the same box
%
% FORMAT:       [mdm = ] mdm.ReframeVTCs(bbox [, newext])
%
% Input fields:
%
%       bbox        output of VTC::BoundingBox from VTC with the desired box
%       newext      string with replacement for '.vtc' (e.g. '_rfbox.vtc')
%
% Output fields:
%
%       mdm         input object
%
% Note: this method will only work on ANY data if ALL data can be framed
%       into the desired box, otherwise an error will be thrown
%
%       data that is outside the requested box will be DISCARDED!! voxels
%       in the requested box not filled by the original data will be set
%       to 0 (for all time points)!
%
%       if newext is not given, VTCs will be overwritten!

% Version:  v1.1
% Build:    17080413
% Date:     Aug-04 2017, 1:53 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2017, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm') || ...
   ~isstruct(bbox) || numel(bbox) ~= 1 || ~isfield(bbox, 'BBox') || ...
   ~isfield(bbox, 'ResXYZ') || ~isequal(size(bbox.BBox), [2, 3]) || ...
    numel(bbox.ResXYZ) ~= 3 || any(mod(bbox.BBox(:), 1) ~= 0) || ...
    any(bbox.ResXYZ ~= bbox.ResXYZ(1)) || ~any(bbox.ResXYZ(1) == (1:12))
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
bc = xo.C;
res = bbox.ResXYZ(1);
bbox = bbox.BBox;
bbox(2, :) = bbox(2, :) + 1;
sbox = diff(bbox) ./ res;
if any(mod(sbox, 1) ~= 0)
    error('neuroelf:xff:badArgument', 'Invalid bounding box requested.');
end
if nargin < 3 || ~ischar(newext) || isempty(newext) || strcmpi(newext(:)', '.vtc')
    newext = '';
end

% all files referenced must be VTC files
vtcfiles = bc.XTC_RTC(:, 1);
if any(cellfun('isempty', regexp(vtcfiles, '\.vtc$')))
    error('neuroelf:xff:badObject', 'All time-course files must be VTCs.');
end

% try to load current bounding box of all linked VTCs
vtco = {[]};
try
    bboxes = zeros(numel(vtcfiles), 7);
    for fc = 1:numel(vtcfiles)
        vtco{1} = xff(vtcfiles{fc}, 't');
        if ~xffisobject(vtco{1}, true, 'vtc')
            error('neuroelf:xff:badObject', 'Not a VTC file: %s.', vtcfiles{fc});
        end
        vtcc = vtco{1}.C;
        bboxes(fc, :) = [vtcc.Resolution, vtcc.XStart, vtcc.XEnd, ...
            vtcc.YStart, vtcc.YEnd, vtcc.ZStart, vtcc.ZEnd];
        delete(vtco{1});
    end
catch xfferror
    if xffisobject(vtco{1}, true)
        delete(vtco{1});
    end
    rethrow(xfferror);
end

% test whether change CAN be made
if any(bboxes(:, 1) ~= res)
    error('neuroelf:xff:badObject', 'Not all VTCs share resolution with bbox.');
end
if any(mod(bboxes(:, 2) - bbox(1, 1), res) ~= 0) || ...
    any(mod(bboxes(:, 4) - bbox(1, 2), res) ~= 0) || ...
    any(mod(bboxes(:, 6) - bbox(1, 3), res) ~= 0)
    error('neuroelf:xff:badObject', 'At least one VTC would require resampling.');
end

% test whether no change is necessary
bboxeq = (bboxes(:, 2:7) == repmat(bbox(:)', numel(vtcfiles), 1));
if all(bboxeq(:))
    return;
end

% which VTCs need handling
workfiles = find(~all(bboxeq, 2));

% process those VTCs
vtco = {[], []};
for fc = 1:numel(workfiles)
    
    % handle errors
    try

        % load VTC (including data!)
        vtco{1} = xff(vtcfiles{workfiles(fc)});

        % copy object
        vtco{2} = xffcopyobject(vtco{1});

        % get data
        vtcd = vtco{1}.C.VTCData;
        if istransio(vtcd)
            vtcd = resolve(transio);
        end
        vtcsz = size(vtcd);
        vv = vtcsz(1);
        vx = vtcsz(2);
        vy = vtcsz(3);
        vz = vtcsz(4);

        % delete source object
        delete(vtco{1});
        vtco{1} = [];

        % create an 1x1 value of the same type with the required size = 0
        newd = vtcd(1);
        newd(:) = 0;

        % figure out how many voxels to add/delete at each dimensions 1/end
        srcbox = reshape(bboxes(workfiles(fc), 2:7), 2, 3);
        difbox = round((srcbox - bbox) / res);
        difbox(2, :) = -difbox(2, :);
        
        % add at XStart
        if difbox(1, 1) > 0
            vtcd = cat(2, repmat(newd, [vv, difbox(1, 1), vtcsz(3:4)]), vtcd);

        % remove at XStart
        elseif difbox(1, 1) < 0
            vtcd(:, 1:-difbox(1, 1), :, :) = [];
        end
        
        % add at XEnd
        if difbox(2, 1) > 0
            vtcd = cat(2, vtcd, repmat(newd, [vv, difbox(2, 1), vtcsz(3:4)]));

        % remove at XEnd
        elseif difbox(2, 1) < 0
            vtcd(:, end+1+difbox(2, 1):end, :, :) = [];
        end
        vx = size(vtcd, 2);
        
        % work on YStart/YEnd
        if difbox(1, 2) > 0
            vtcd = cat(3, repmat(newd, [vv, vx, difbox(1, 2), vz]), vtcd);
        elseif difbox(1, 1) < 0
            vtcd(:, :, 1:-difbox(1, 2), :) = [];
        end
        if difbox(2, 2) > 0
            vtcd = cat(3, vtcd, repmat(newd, [vv, vx, difbox(2, 2), vz]));
        elseif difbox(2, 1) < 0
            vtcd(:, :, end+1+difbox(2, 2):end, :) = [];
        end
        vy = size(vtcd, 3);
        
        % work on ZStart/ZEnd
        if difbox(1, 3) > 0
            vtcd = cat(4, repmat(newd, [vv, vx, vy, difbox(1, 3)]), vtcd);
        elseif difbox(1, 3) < 0
            vtcd(:, :, :, 1:-difbox(1, 3)) = [];
        end
        if difbox(2, 3) > 0
            vtcd = cat(4, vtcd, repmat(newd, [vv, vx, vy, difbox(2, 3)]));
        elseif difbox(2, 3) < 0
            vtcd(:, :, :, end+1+difbox(2, 3):end) = [];
        end
        
        % set data
        vtco{2}.C.XStart = bbox(1, 1);
        vtco{2}.C.XEnd = bbox(2, 1);
        vtco{2}.C.YStart = bbox(1, 2);
        vtco{2}.C.YEnd = bbox(2, 2);
        vtco{2}.C.ZStart = bbox(1, 3);
        vtco{2}.C.ZEnd = bbox(2, 3);
        vtco{2}.C.VTCData = vtcd;

        % add information to RunTimeVars
        rtv = vtco{2}.C.RunTimeVars;
        if isfield(rtv, 'DVARS') && iscell(rtv.DVARS) && numel(rtv.DVARS) == 2 && ...
           ~isempty(rtv.DVARS{1})
            vtco{2}.C.RunTimeVars.DVARSBBox = reshape(bboxes(workfiles(fc), 2:7), 2, 3);
        end

        % save with old filename
        if isempty(newext)
            aft_SaveAs(vtco{2}, vtcfiles{workfiles(fc)});
            aft_SaveRunTimeVars(vtco{2});

        % save with new filename
        else
            aft_SaveAs(vtco{2}, [vtcfiles{workfiles(fc)}(1:end-4) newext]);
            aft_SaveRunTimeVars(vtco{2});
        end
        
    % throw error
    catch xfferror
        if xffisobject(vtco{1}, true)
            delete(vtco{1});
        end
        if xffisobject(vtco{2}, true)
            delete(vtco{2});
        end
        rethrow(xfferror);
    end
end
function varargout = ne_draw(varargin)
% ne_draw  - draw into dataset
%
% FORMAT:       ne_draw(SRC, EVT [, cpos [, ddir]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       cpos        current position (if different from ne_gcfg.fcfg.cpos)
%       ddir        1x2 drawing plane (only for 2D drawing: 1, 2, or 3)
%
% No output fields.
%
% Examples:
%
%       ne_draw(0, 0, [24, 0, -24]);
%
% Notes:
%
%   this function is automatically called upon by the GUI (ne_setslicepos)
%   when drawing is enabled
%
%   for positive codes, the values are written directly into the data
%   for negative values, subtract from data (make darker), with special
%   ranges [-2 .. 0[ : multiply current data with (positive) value;
%   [-4 .. -2[ : multiply original data with (positive) value;

% Version:  v1.1
% Build:    16012316
% Date:     Jan-23 2016, 4:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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

% global variable
global ne_gcfg;

% get handle shortcuts
cc = ne_gcfg.fcfg;
cp = cc.paint;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% currently configured position
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 3 || ...
    any(isinf(varargin{3}) | isnan(varargin{3}) | varargin{3} < -128 | varargin{3} > 128)
    cpos = cc.cpos;
else
    cpos = varargin{3}(:)';
end
if nargin < 4 || ...
   ~isa(varargin{4}, 'double') || ...
    numel(varargin{4}) ~= 2 || ...
    any(isinf(varargin{4}) | isnan(varargin{4})) || ...
   ~any(varargin{4}(1) == [1, 2, 3]) || ...
   ~any(varargin{4}(2) == [1, 2, 3]) || ...
    varargin{4}(1) == varargin{4}(2)
    ddir = cc.ddir;
else
    ddir = varargin{4};
end

% no drawing after all?
svar = cc.SliceVar;
if ~isxff(svar, {'hdr', 'vmr'}) || abs(cp.mode) <= 1 || isempty(ddir)
    return;
end
svartyp = lower(svar.Filetype);

% and transformation values (if not == eye(4))
if ~isequal(cc.strans, eye(4))
    trans = cc.strans;
else
    trans = [];
end

% use Trf to get voxel position
rtv = svar.RunTimeVars;
if isfield(rtv, 'Trf')
    ttrf = rtv.Trf;
else
    ttrf = eye(4);
end
if isfield(rtv, 'TrfPlus')
    tplus = inv(rtv.TrfPlus)';
else
    tplus = eye(4);
end
if ~isempty(trans)
    bvpos = [cpos, 1] * tplus * trans' * ttrf;
else
    bvpos = [cpos, 1] * tplus * ttrf;
end
bvpos(end) = [];

% painting codes
codes = double(cp.code(:)');

% no weights (default)
sw = [];

% undofield name
ufn = 'UndoBuffer';

% for VMRs
if strcmp(svartyp, 'vmr')

    % fieldname
    fn = 'VMRData';

    % V16?
    if isequal(size(svar.(fn)), size(svar.VMRData16)) && ...
        rtv.ShowV16
        fn = [fn, '16'];
        ufn = [ufn, '16'];
    end

    % no transio
    if istransio(svar.(fn))
        svar.(fn) = resolve(svar.(fn));
    end

    % for 2D-painting, use direction from above
    if abs(cp.mode) == 2

        ss = cp.shap2;
        dpos = round(bvpos(ones(size(ss, 1), 1), :));
        dpos(:, ddir) = dpos(:, ddir) + ss;
        if cp.smooth > 0
            sw = cp.shap2w;
        end

    % 3D-painting
    elseif abs(cp.mode) == 3
        ss = cp.shap3;
        dpos = round(bvpos(ones(size(ss, 1), 1), :)) + ss;
        if cp.smooth > 0
            sw = cp.shap3w;
        end
    end

% for Analyze files
elseif strcmp(svartyp, 'hdr') && ...
   (ndims(svar.VoxelData) == 3 || ...
    ndims(svar.VoxelDataRGBA) == 5)

    % resolve if necessary (not drawing to disk!)
    fn = 'VoxelData';
    if isempty(svar.(fn))
        fn = 'VoxelDataRGBA';
    end
    if istransio(svar.(fn))
        svar.(fn) = resolve(svar.(fn));
    end

    % get drawing direction (slice orientation) right
    ddir = 1 + mod(ddir, 3);

    % buffer Trf
    cframe = svar.CoordinateFrame;
    if ~isfield(rtv, 'Trf')
        svar.RunTimeVars.Trf = inv(cframe.Trf(1))';
        rtv = svar.RunTimeVars;
    end

    % apply transformations
    if ~isempty(trans)
        dpos = [cpos, 1] * tplus * trans';
    else
        dpos = [cpos, 1] * tplus;
    end
    dpos(end) = [];

    % for 2D-painting, use direction from above
    if abs(cp.mode) == 2
        ss = cp.shap2 - 0.25;
        ss = [ss; [ss(:, 1), ss(:, 2) + 0.5]; ...
             [ss(:, 1) + 0.5, ss(:, 2)]; ss + 0.5];
        if cp.smooth > 0
            sw = repmat(cp.shap2w, 12, 1);
        end
        dpos = dpos(ones(size(ss, 1), 1), :);
        ddpos = (0.3 * mean(cframe.Resolution)) .* ...
            [-ones(size(dpos)); zeros(size(dpos)); ones(size(dpos))];
        ddpos(:, ddir) = repmat(ss, 3, 1);
        ddpos = ddpos * tplus(1:3, 1:3);
        dpos = repmat(dpos, 3, 1) + ddpos;

    % 3D-painting
    elseif abs(cp.mode) == 3
        ss = cp.shap3 - 0.25;
        ss = [ss; ...
              [ss(:, [1, 2]), ss(:, 3) + 0.5]; ...
              [ss(:, 1), ss(:, 2) + 0.5, ss(:, 3)]; ...
              [ss(:, 1), ss(:, [2, 3]) + 0.5]; ...
              [ss(:, 1) + 0.5, ss(:, [2, 3])]; ...
              [ss(:, 1) + 0.5, ss(:, 2), ss(:, 3) + 0.5]; ...
              [ss(:, [1, 2]) + 0.5, ss(:, 3)]; ...
              ss + 0.5];
        if cp.smooth > 0
            sw = repmat(cp.shap3w, 8, 1);
        end
        dpos = dpos(ones(size(ss, 1), 1), :) + ss;
    end

    % get voxel positions
    dpos(:, 4) = 1;
    dpos = round(dpos * rtv.Trf);

% otherwise no drawing
else
    return;
end

% make sure we have an UndoBuffer
if ~isfield(rtv, ufn) || ...
    isempty(rtv.(ufn))
    svar.RunTimeVars.(ufn) = svar.(fn);
    rtv = svar.RunTimeVars;
end

% remove invalid voxels
sz = size(svar.(fn));
if numel(sz) < 5
    sz(end+1:5) = 1;
end
sz3 = prod(sz(1:3));
if numel(sw) <= 1
    dpos(dpos(:, 1) < 1 | dpos(:, 1) > sz(1) | ...
         dpos(:, 2) < 1 | dpos(:, 2) > sz(2) | ...
         dpos(:, 3) < 1 | dpos(:, 3) > sz(3), :) = [];
    dpos = sub2ind(sz(1:3), dpos(:, 1), dpos(:, 2), dpos(:, 3));
    if isempty(dpos)
        return;
    end
else
    ddpos = (dpos(:, 1) < 1 | dpos(:, 1) > sz(1) | ...
         dpos(:, 2) < 1 | dpos(:, 2) > sz(2) | ...
         dpos(:, 3) < 1 | dpos(:, 3) > sz(3));
    dpos(ddpos, :) = [];
    if isempty(dpos)
        return;
    end
    sw(ddpos) = [];
    dpos = sub2ind(sz(1:3), dpos(:, 1), dpos(:, 2), dpos(:, 3));
    if strcmp(svartyp, 'hdr')
         [udpos, idpos1] = unique(dpos(end:-1:1));
         [dpos, idpos2] = unique(dpos);
         sw = 0.5 .* (sw((numel(sw) + 1) - idpos1) + sw(idpos2));
    end
end

% extend codes
if numel(codes) < sz(5)
    codes = repmat(codes, 1, ceil(sz(5) / numel(codes)));
end

% get current values
cval = svar.(fn)(dpos);

% determine datatype issues
cl = lower(class(cval));
ii = ~isempty(strfind(cl, 'int'));
switch (cl)
    case {'int16'};
        minv = -32768;
        maxv = 32767;
    case {'int8'};
        minv = -128;
        maxv = 127;
    case {'uint16'}
        minv = 0;
        maxv = 65535;
    case {'uint8'}
        minv = 0;
        maxv = 255;
    otherwise
        minv = -Inf;
        maxv = Inf;
end

% paint mode, disable updated
upds = ne_gcfg.c.xff.UpdateState(svartyp, false);

% regular paint mode (> 0)
if cp.mode > 0

    % check overwrite range
    ovr = cp.over;
    if ovr(1) > minv || ...
        ovr(2) < maxv

        % for a single volume
        if sz(5) == 1

            % one check
            cvalr = (cval < ovr(1) | cval > ovr(2));

        % for multiple volumes
        else

            % convert to double
            cval = double(cval);

            % and tally
            for dc = 2:sz(5)
                cval = cval + double(svar.(fn)(dpos + (dc - 1) * sz3));
            end

            % then check
            cvalr = (cval < (sz(5) * ovr(1)) | cval > (sz(5) * ovr(2)));

            % and get values back
            cval = svar.(fn)(dpos);
        end

        % remove from selection
        dpos(cvalr) = [];

        % return if nothing left to do
        if isempty(dpos)
            return;
        end

        % otherwise remove from values and weights
        cval(cvalr) = [];
        if numel(sw) > 1
            sw(cvalr) = [];
        end
    end

    % get negative weights
    if numel(sw) > 1
        nsw = 1 - sw;
    end

    % regular painting (setting code)
    if codes(1) >= 0

        % with weights
        if numel(sw) > 1

            % iterate over 5th dimension
            for dc = 1:sz(5)

                % set in buffer
                if ii
                    svar.(fn)(dpos) = min(maxv, max(minv, round( ...
                        codes(dc) .* sw + nsw .* double(cval))));
                else
                    svar.(fn)(dpos) = ...
                        codes(dc) .* sw + nsw .* double(cval);
                end

                % next volume
                if sz(5) > dc
                    dpos = dpos + sz3;
                    cval = svar.(fn)(dpos);
                end
            end

        % without weights
        else

            % iterate over 5th dimension
            for dc = 1:sz(5)
                svar.(fn)(dpos) = codes(dc);
                if sz(5) > dc
                    dpos = dpos + sz3;
                end
            end
        end

    % subtracting (making darker)
    elseif codes(1) < -4

        % with weights
        if numel(sw) > 1

            % volumes
            for dc = 1:sz(5)
                if ii
                    svar.(fn)(dpos) = max(minv, round( ...
                        sw .* (codes(dc) + double(cval)) + nsw .* double(cval)));
                else
                    svar.(fn)(dpos) = ...
                        sw .* (codes(dc) + double(cval)) + nsw .* double(cval);
                end
                if sz(5) > dc
                    dpos = dpos + sz3;
                    cval = svar.(fn)(dpos);
                end
            end

        % without weights
        else
            for dc = 1:sz(5)
                if ii
                    svar.(fn)(dpos) = max(minv, codes(dc) + cval);
                else
                    svar.(fn)(dpos) = codes(dc) + cval;
                end
                if sz(5) > dc
                    dpos = dpos + sz3;
                    cval = svar.(fn)(dpos);
                end
            end
        end

    % multiplication of buffer (with -code)
    elseif codes(1) < -2

        % adjust codes
        codes = min(2, max(0, -2 - codes));

        % with weights
        if numel(sw) > 1
            for dc = 1:sz(5)
                uval = rtv.(ufn)(dpos);
                if ii
                    svar.(fn)(dpos) = max(minv, min(maxv, round( ...
                        sw .* (codes(dc) .* double(uval)) + nsw .* double(cval))));
                else
                    svar.(fn)(dpos) = ...
                        sw .* (codes(dc) .* double(uval)) + nsw .* double(cval);
                end
                if sz(5) > dc
                    dpos = dpos + sz3;
                    cval = svar.(fn)(dpos);
                end
            end
        else
            for dc = 1:sz(5)
                uval = rtv.(ufn)(dpos);
                if ii
                    svar.(fn)(dpos) = max(minv, min(maxv, round( ...
                        codes(dc) .* double(uval))));
                else
                    svar.(fn)(dpos) = ...
                        codes(dc) .* double(uval);
                end
                if sz(5) > dc
                    dpos = dpos + sz3;
                end
            end
        end

    % multiplication (with -code)
    else

        % adjust codes
        codes = min(2, max(0, -codes));

        % with weights
        if numel(sw) > 1
            for dc = 1:sz(5)
                if ii
                    svar.(fn)(dpos) = max(minv, min(maxv, round( ...
                        sw .* (codes(dc) .* double(cval)) + nsw .* double(cval))));
                else
                    svar.(fn)(dpos) = ...
                        sw .* (codes(dc) .* double(cval)) + nsw .* double(cval);
                end
                if sz(5) > dc
                    dpos = dpos + sz3;
                    cval = svar.(fn)(dpos);
                end
            end
        else
            for dc = 1:sz(5)
                if ii
                    svar.(fn)(dpos) = max(minv, min(maxv, round( ...
                        codes(dc) .* double(cval))));
                else
                    svar.(fn)(dpos) = ...
                        codes(dc) .* double(cval);
                end
                if sz(5) > dc
                    dpos = dpos + sz3;
                    cval = svar.(fn)(dpos);
                end
            end
        end
    end

% undo painting
else
    if numel(sw) > 1
        nsw = 1 - sw;
        for dc = 1:sz(5)
            uval = rtv.(ufn)(dpos);
            if ii
                svar.(fn)(dpos) = max(minv, min(maxv, round( ...
                    sw .* double(uval) + nsw .* double(cval))));
            else
                svar.(fn)(dpos) = ...
                    sw .* double(uval) + nsw .* double(cval);
            end
            if sz(5) > dc
                dpos = dpos + sz3;
                cval = svar.(fn)(dpos);
            end
        end
    else
        for dc = 1:sz(5)
            svar.(fn)(dpos) = rtv.(ufn)(dpos);
            if sz(5) > dc
                dpos = dpos + sz3;
            end
        end
    end
end

% reinstantiate update state
ne_gcfg.c.xff.UpdateState(svartyp, upds);
